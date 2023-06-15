use anyhow::Result;
use hashbrown::{HashMap, HashSet};
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

#[derive(Debug)]
pub struct Barcodes {
    map: HashMap<Vec<u8>, usize>,
    index: HashMap<usize, Vec<u8>>,
    len: usize,
}
impl Barcodes {
    pub fn from_file(path: &str) -> Result<Self> {
        let reader = File::open(path).map(BufReader::new)?;
        Self::from_buffer(reader)
    }
    pub fn from_file_with_spacer(path: &str, spacer: &Spacer) -> Result<Self> {
        let reader = File::open(path).map(BufReader::new)?;
        Self::from_buffer_with_spacer(reader, spacer)
    }

    pub fn from_buffer<R: BufRead>(reader: R) -> Result<Self> {
        Self::parse_buffer(reader, None)
    }

    pub fn from_buffer_with_spacer<R: BufRead>(reader: R, spacer: &Spacer) -> Result<Self> {
        Self::parse_buffer(reader, Some(spacer))
    }

    pub fn parse_buffer<R: BufRead>(reader: R, spacer: Option<&Spacer>) -> Result<Self> {
        let mut map = HashMap::new();
        let mut index = HashMap::new();
        let mut sizes = HashSet::new();

        for (idx, line) in reader.lines().enumerate() {
            let barcode = line.map(|l| Self::read_sequence(&l, spacer))?;
            sizes.insert(barcode.len());
            map.entry(barcode.clone()).or_insert(idx);
            index.entry(idx).or_insert(barcode);
        }

        let len = if sizes.len() == 1 {
            sizes.into_iter().next().unwrap()
        } else {
            anyhow::bail!("Barcodes have different lengths");
        };

        Ok(Self { map, index, len })
    }

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
    pub fn match_sequence(&self, sequence: &[u8]) -> Option<(usize, usize)> {
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
    ) -> Option<(usize, usize)> {
        self.match_sequence(&sequence[start..end])
    }

    /// Returns the barcode sequence for a given index
    pub fn get_barcode(&self, idx: usize) -> Option<&[u8]> {
        self.index.get(&idx).map(|bc| &bc[..])
    }

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
