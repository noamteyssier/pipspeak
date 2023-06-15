use anyhow::Result;
use disambiseq::Disambibyte;
use hashbrown::{HashMap, HashSet};
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

type BarcodeID = usize;
type EndPos = usize;

#[derive(Debug)]
pub struct Barcodes {
    map: HashMap<Vec<u8>, usize>,
    index: HashMap<usize, Vec<u8>>,
    len: usize,
}
impl Barcodes {
    pub fn from_file(path: &str, exact: bool) -> Result<Self> {
        let reader = File::open(path).map(BufReader::new)?;
        Self::from_buffer(reader, exact)
    }
    pub fn from_file_with_spacer(path: &str, spacer: &Spacer, exact: bool) -> Result<Self> {
        let reader = File::open(path).map(BufReader::new)?;
        Self::from_buffer_with_spacer(reader, spacer, exact)
    }

    pub fn from_buffer<R: BufRead>(reader: R, exact: bool) -> Result<Self> {
        Self::parse_buffer(reader, None, exact)
    }

    pub fn from_buffer_with_spacer<R: BufRead>(
        reader: R,
        spacer: &Spacer,
        exact: bool,
    ) -> Result<Self> {
        Self::parse_buffer(reader, Some(spacer), exact)
    }

    /// Parses a buffer and returns a Barcodes object
    /// If a spacer is given, it is appended to each barcode
    pub fn parse_buffer<R: BufRead>(
        reader: R,
        spacer: Option<&Spacer>,
        exact: bool,
    ) -> Result<Self> {
        let mut map = HashMap::new();
        let mut index = HashMap::new();
        let mut sizes = HashSet::new();

        for (idx, line) in reader.lines().enumerate() {
            let barcode = line.map(|l| Self::read_sequence(&l, spacer))?;
            sizes.insert(barcode.len());
            map.entry(barcode.clone()).or_insert(idx);
            index.entry(idx).or_insert(barcode);
        }

        if !exact {
            let parent_barcodes = map.keys().cloned().collect::<Vec<_>>();
            let dsb = Disambibyte::from_slice(&parent_barcodes);
            dsb.unambiguous().iter().for_each(|(child, parent)| {
                map.insert(
                    child.sequence().to_owned(),
                    *map.get(parent.sequence()).unwrap(),
                );
            });
        }

        let len = if sizes.len() == 1 {
            sizes.into_iter().next().unwrap()
        } else {
            anyhow::bail!("Barcodes have different lengths");
        };

        Ok(Self { map, index, len })
    }

    /// Reads a sequence from a line and appends a spacer if given
    /// Returns the sequence as a vector of integer nucleotides
    fn read_sequence(line: &str, spacer: Option<&Spacer>) -> Vec<u8> {
        let barcode = line.trim().as_bytes().to_vec();
        if let Some(spacer) = spacer {
            let mut barcode_with_spacer = barcode.clone();
            barcode_with_spacer.extend_from_slice(spacer.seq());
            return barcode_with_spacer;
        } else {
            return barcode;
        }
    }

    /// Checks if a sequence contains a barcode as a substring
    /// and returns the position of the first nucleotide after the barcode
    /// as well as the barcode index
    pub fn match_sequence(&self, sequence: &[u8]) -> Option<(EndPos, BarcodeID)> {
        if sequence.len() < self.len {
            return None;
        }
        sequence
            .windows(self.len)
            .position(|window| self.map.contains_key(window))
            .map(|pos| {
                (
                    pos + self.len,
                    *self.map.get(&sequence[pos..pos + self.len]).unwrap(),
                )
            })
    }

    /// Matches a subsequence of a sequence
    /// and returns the position of the first nucleotide after the barcode
    /// as well as the barcode index
    pub fn match_subsequence(
        &self,
        sequence: &[u8],
        start: usize,
        end: usize,
    ) -> Option<(EndPos, BarcodeID)> {
        if start > sequence.len() || end > sequence.len() || start > end {
            return None;
        }
        self.match_sequence(&sequence[start..end])
    }

    /// Returns the barcode sequence for a given index
    pub fn get_barcode(&self, idx: usize) -> Option<&[u8]> {
        self.index.get(&idx).map(|bc| &bc[..])
    }

    /// Returns the barcode index for a given sequence
    #[allow(dead_code)]
    pub fn get_id(&self, barcode: &[u8]) -> Option<usize> {
        self.map.get(barcode).map(|id| *id)
    }

    /// Returns the length of each barcode
    pub fn len(&self) -> usize {
        self.len
    }
}

pub struct Spacer {
    seq: Vec<u8>,
}
impl Spacer {
    pub fn from_str(seq: &str) -> Self {
        Self {
            seq: seq.as_bytes().to_vec(),
        }
    }
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
}

#[cfg(test)]
mod testing {
    use super::*;

    const TEST_FILE: &str = "data/barcodes_v3/fb_v3_bc1.tsv";
    const TEST_BUFFER: &[u8] = b"AGAAACCA\nGATTTCCC\nAAGTCCAA\nGAGAAACC";
    const MALFORMED_BUFFER: &[u8] = b"AGAAACCA\nGATTTCCC\nAAGTCCAA\nGAGAAACCC";
    const TEST_SPACER: &str = "ATG";
    const NOMATCH_SEQ: &[u8] = b"SHOULDNOTMATCHANYTHING";
    const ENDMATCH_SEQ: &[u8] = b"OFFSETXAGAAACCA";
    const ENDMATCH_SEQ_1D: &[u8] = b"OFFSETXTGAAACCA";
    const STARTMATCH_SEQ: &[u8] = b"AGAAACCAANDSOMETHINGELSE";
    const STARTMATCH_SEQ_1D: &[u8] = b"TGAAACCAANDSOMETHINGELSE";
    const OFFSETMATCH_SEQ: &[u8] = b"123AGAAACCASOMETHINGELSE";
    const OFFSETMATCH_SEQ_1D: &[u8] = b"123TGAAACCASOMETHINGELSE";

    #[test]
    fn from_file() {
        let barcodes = Barcodes::from_file(TEST_FILE, false).unwrap();
        assert_eq!(barcodes.len(), 8);
        assert_eq!(barcodes.map.len(), 2360);
        assert_eq!(barcodes.index.len(), 96);
    }

    #[test]
    fn from_file_exact() {
        let barcodes = Barcodes::from_file(TEST_FILE, true).unwrap();
        assert_eq!(barcodes.len(), 8);
        assert_eq!(barcodes.map.len(), 96);
        assert_eq!(barcodes.index.len(), 96);
    }

    #[test]
    fn from_buffer_exact() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, true).unwrap();
        assert_eq!(barcodes.len(), 8);
        assert_eq!(barcodes.map.len(), 4);
        assert_eq!(barcodes.index.len(), 4);

        assert_eq!(barcodes.get_barcode(0).unwrap(), b"AGAAACCA");
        assert_eq!(barcodes.get_barcode(1).unwrap(), b"GATTTCCC");
        assert_eq!(barcodes.get_barcode(2).unwrap(), b"AAGTCCAA");
        assert_eq!(barcodes.get_barcode(3).unwrap(), b"GAGAAACC");

        assert_eq!(barcodes.get_id(b"AGAAACCA").unwrap(), 0);
        assert_eq!(barcodes.get_id(b"GATTTCCC").unwrap(), 1);
        assert_eq!(barcodes.get_id(b"AAGTCCAA").unwrap(), 2);
        assert_eq!(barcodes.get_id(b"GAGAAACC").unwrap(), 3);
    }

    #[test]
    fn from_buffer() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, false).unwrap();
        assert_eq!(barcodes.len(), 8);
        assert_eq!(barcodes.map.len(), 100);
        assert_eq!(barcodes.index.len(), 4);

        assert_eq!(barcodes.get_barcode(0).unwrap(), b"AGAAACCA");
        assert_eq!(barcodes.get_barcode(1).unwrap(), b"GATTTCCC");
        assert_eq!(barcodes.get_barcode(2).unwrap(), b"AAGTCCAA");
        assert_eq!(barcodes.get_barcode(3).unwrap(), b"GAGAAACC");

        // no mismatch
        assert_eq!(barcodes.get_id(b"AGAAACCA").unwrap(), 0);
        assert_eq!(barcodes.get_id(b"GATTTCCC").unwrap(), 1);
        assert_eq!(barcodes.get_id(b"AAGTCCAA").unwrap(), 2);
        assert_eq!(barcodes.get_id(b"GAGAAACC").unwrap(), 3);

        // with mismatch
        assert_eq!(barcodes.get_id(b"TGAAACCA").unwrap(), 0);
        assert_eq!(barcodes.get_id(b"CATTTCCC").unwrap(), 1);
        assert_eq!(barcodes.get_id(b"TAGTCCAA").unwrap(), 2);
        assert_eq!(barcodes.get_id(b"CAGAAACC").unwrap(), 3);

        // hamming distance of 2 should fail
        assert_eq!(barcodes.get_id(b"TCAAACCA"), None);
        assert_eq!(barcodes.get_id(b"CTTTTCCC"), None);
        assert_eq!(barcodes.get_id(b"TTGTCCAA"), None);
        assert_eq!(barcodes.get_id(b"CCGAAACC"), None);
    }

    #[test]
    fn from_file_with_spacer() {
        let spacer = Spacer::from_str(TEST_SPACER);
        let barcodes = Barcodes::from_file_with_spacer(TEST_FILE, &spacer, false).unwrap();
        assert_eq!(barcodes.len(), 11);
        assert_eq!(barcodes.map.len(), 3224);
        assert_eq!(barcodes.index.len(), 96);
    }

    #[test]
    fn from_file_with_spacer_exact() {
        let spacer = Spacer::from_str(TEST_SPACER);
        let barcodes = Barcodes::from_file_with_spacer(TEST_FILE, &spacer, true).unwrap();
        assert_eq!(barcodes.len(), 11);
        assert_eq!(barcodes.map.len(), 96);
        assert_eq!(barcodes.index.len(), 96);
    }

    #[test]
    fn from_buffer_with_spacer() {
        let spacer = Spacer::from_str(TEST_SPACER);
        let barcodes = Barcodes::from_buffer_with_spacer(TEST_BUFFER, &spacer, false).unwrap();
        assert_eq!(barcodes.len(), 11);
        assert_eq!(barcodes.map.len(), 136);
        assert_eq!(barcodes.index.len(), 4);

        assert_eq!(barcodes.get_barcode(0).unwrap(), b"AGAAACCAATG");
        assert_eq!(barcodes.get_barcode(1).unwrap(), b"GATTTCCCATG");
        assert_eq!(barcodes.get_barcode(2).unwrap(), b"AAGTCCAAATG");
        assert_eq!(barcodes.get_barcode(3).unwrap(), b"GAGAAACCATG");

        // no mismatch
        assert_eq!(barcodes.get_id(b"AGAAACCAATG").unwrap(), 0);
        assert_eq!(barcodes.get_id(b"GATTTCCCATG").unwrap(), 1);
        assert_eq!(barcodes.get_id(b"AAGTCCAAATG").unwrap(), 2);
        assert_eq!(barcodes.get_id(b"GAGAAACCATG").unwrap(), 3);

        // with mismatch
        assert_eq!(barcodes.get_id(b"TGAAACCAATG").unwrap(), 0);
        assert_eq!(barcodes.get_id(b"TATTTCCCATG").unwrap(), 1);
        assert_eq!(barcodes.get_id(b"TAGTCCAAATG").unwrap(), 2);
        assert_eq!(barcodes.get_id(b"TAGAAACCATG").unwrap(), 3);

        // with mismatch of 2 should fail
        assert_eq!(barcodes.get_id(b"TTAAACCAATG"), None);
        assert_eq!(barcodes.get_id(b"TTTTTCCCATG"), None);
        assert_eq!(barcodes.get_id(b"TTGTCCAAATG"), None);
        assert_eq!(barcodes.get_id(b"TTGAAACCATG"), None);
    }

    #[test]
    fn from_buffer_with_spacer_exact() {
        let spacer = Spacer::from_str(TEST_SPACER);
        let barcodes = Barcodes::from_buffer_with_spacer(TEST_BUFFER, &spacer, true).unwrap();
        assert_eq!(barcodes.len(), 11);
        assert_eq!(barcodes.map.len(), 4);
        assert_eq!(barcodes.index.len(), 4);

        assert_eq!(barcodes.get_barcode(0).unwrap(), b"AGAAACCAATG");
        assert_eq!(barcodes.get_barcode(1).unwrap(), b"GATTTCCCATG");
        assert_eq!(barcodes.get_barcode(2).unwrap(), b"AAGTCCAAATG");
        assert_eq!(barcodes.get_barcode(3).unwrap(), b"GAGAAACCATG");

        assert_eq!(barcodes.get_id(b"AGAAACCAATG").unwrap(), 0);
        assert_eq!(barcodes.get_id(b"GATTTCCCATG").unwrap(), 1);
        assert_eq!(barcodes.get_id(b"AAGTCCAAATG").unwrap(), 2);
        assert_eq!(barcodes.get_id(b"GAGAAACCATG").unwrap(), 3);
    }

    #[test]
    fn size_variance() {
        let barcodes = Barcodes::from_buffer(MALFORMED_BUFFER, false);
        assert!(barcodes.is_err());
    }

    #[test]
    fn size_variance_exact() {
        let barcodes = Barcodes::from_buffer(MALFORMED_BUFFER, true);
        assert!(barcodes.is_err());
    }

    #[test]
    fn size_variance_with_spacer() {
        let spacer = Spacer::from_str(TEST_SPACER);
        let barcodes = Barcodes::from_buffer_with_spacer(MALFORMED_BUFFER, &spacer, false);
        assert!(barcodes.is_err());
    }

    #[test]
    fn size_variance_with_spacer_exact() {
        let spacer = Spacer::from_str(TEST_SPACER);
        let barcodes = Barcodes::from_buffer_with_spacer(MALFORMED_BUFFER, &spacer, true);
        assert!(barcodes.is_err());
    }

    #[test]
    fn match_sequence() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, false).unwrap();

        // no mismatch
        assert_eq!(barcodes.match_sequence(NOMATCH_SEQ), None);
        assert_eq!(
            barcodes.match_sequence(ENDMATCH_SEQ),
            Some((7 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_sequence(STARTMATCH_SEQ),
            Some((0 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_sequence(OFFSETMATCH_SEQ),
            Some((3 + barcodes.len(), 0))
        );

        // with mismatch
        assert_eq!(
            barcodes.match_sequence(ENDMATCH_SEQ_1D),
            Some((7 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_sequence(STARTMATCH_SEQ_1D),
            Some((0 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_sequence(OFFSETMATCH_SEQ_1D),
            Some((3 + barcodes.len(), 0))
        );
    }

    #[test]
    fn match_sequence_exact() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, true).unwrap();

        // no mismatch
        assert_eq!(barcodes.match_sequence(NOMATCH_SEQ), None);
        assert_eq!(
            barcodes.match_sequence(ENDMATCH_SEQ),
            Some((7 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_sequence(STARTMATCH_SEQ),
            Some((0 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_sequence(OFFSETMATCH_SEQ),
            Some((3 + barcodes.len(), 0))
        );

        // with mismatch
        assert_eq!(barcodes.match_sequence(ENDMATCH_SEQ_1D), None);
        assert_eq!(barcodes.match_sequence(STARTMATCH_SEQ_1D), None);
        assert_eq!(barcodes.match_sequence(OFFSETMATCH_SEQ_1D), None);
    }

    #[test]
    fn match_subsequence() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, false).unwrap();
        let start_pos = 7;
        let end_pos = start_pos + barcodes.len();

        // no mismatch
        assert_eq!(
            barcodes.match_subsequence(NOMATCH_SEQ, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(ENDMATCH_SEQ, start_pos, end_pos),
            Some((0 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_subsequence(STARTMATCH_SEQ, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(OFFSETMATCH_SEQ, start_pos, end_pos),
            None
        );

        // with mismatch
        assert_eq!(
            barcodes.match_subsequence(ENDMATCH_SEQ_1D, start_pos, end_pos),
            Some((0 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_subsequence(STARTMATCH_SEQ_1D, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(OFFSETMATCH_SEQ_1D, start_pos, end_pos),
            None
        );
    }

    #[test]
    fn match_subsequence_exact() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, true).unwrap();
        let start_pos = 7;
        let end_pos = start_pos + barcodes.len();

        // no mismatch
        assert_eq!(
            barcodes.match_subsequence(NOMATCH_SEQ, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(ENDMATCH_SEQ, start_pos, end_pos),
            Some((0 + barcodes.len(), 0))
        );
        assert_eq!(
            barcodes.match_subsequence(STARTMATCH_SEQ, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(OFFSETMATCH_SEQ, start_pos, end_pos),
            None
        );

        // with mismatch
        assert_eq!(
            barcodes.match_subsequence(ENDMATCH_SEQ_1D, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(STARTMATCH_SEQ_1D, start_pos, end_pos),
            None
        );
        assert_eq!(
            barcodes.match_subsequence(OFFSETMATCH_SEQ_1D, start_pos, end_pos),
            None
        );
    }

    #[test]
    fn match_empty() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, false).unwrap();
        assert_eq!(barcodes.match_sequence(b""), None);
        assert_eq!(barcodes.match_subsequence(b"", 0, barcodes.len()), None);
    }

    #[test]
    fn match_empty_exact() {
        let barcodes = Barcodes::from_buffer(TEST_BUFFER, true).unwrap();
        assert_eq!(barcodes.match_sequence(b""), None);
        assert_eq!(barcodes.match_subsequence(b"", 0, barcodes.len()), None);
    }
}
