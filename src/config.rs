use crate::barcodes::{Barcodes, Spacer};
use anyhow::Result;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct ConfigYaml {
    barcodes: ConfigBarcodes,
    spacers: ConfigSpacers,
}

#[derive(Debug, Deserialize)]
pub struct ConfigBarcodes {
    bc1: String,
    bc2: String,
    bc3: String,
    bc4: String,
}

#[derive(Debug, Deserialize)]
pub struct ConfigSpacers {
    s1: String,
    s2: String,
    s3: String,
}

pub struct Config {
    bc1: Barcodes,
    bc2: Barcodes,
    bc3: Barcodes,
    bc4: Barcodes,
    linkers: bool,
}
impl Config {
    pub fn from_file(path: &str, exact: bool, linkers: bool) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        let yaml = serde_yaml::from_str::<ConfigYaml>(&contents)?;
        Self::from_yaml(yaml, exact, linkers)
    }

    pub fn from_yaml(yaml: ConfigYaml, exact: bool, linkers: bool) -> Result<Self> {
        let spacer1 = Spacer::from_str(&yaml.spacers.s1);
        let spacer2 = Spacer::from_str(&yaml.spacers.s2);
        let spacer3 = Spacer::from_str(&yaml.spacers.s3);
        let bc1 = Self::load_barcode(&yaml.barcodes.bc1, Some(&spacer1), exact)?;
        let bc2 = Self::load_barcode(&yaml.barcodes.bc2, Some(&spacer2), exact)?;
        let bc3 = Self::load_barcode(&yaml.barcodes.bc3, Some(&spacer3), exact)?;
        let bc4 = Self::load_barcode(&yaml.barcodes.bc4, None, exact)?;
        Ok(Self {
            bc1,
            bc2,
            bc3,
            bc4,
            linkers,
        })
    }

    fn load_barcode(path: &str, spacer: Option<&Spacer>, exact: bool) -> Result<Barcodes> {
        if let Some(spacer) = spacer {
            Barcodes::from_file_with_spacer(path, spacer, exact)
        } else {
            Barcodes::from_file(path, exact)
        }
    }

    /// Matches a subsequence starting from `pos` against one of the barcode sets.
    /// Returns the end nucleotide position of the match and the within-set barcode index
    pub fn match_subsequence(
        &self,
        seq: &[u8],
        set_idx: usize,
        pos: usize,
        offset: Option<usize>,
    ) -> Option<(usize, usize)> {
        let bc = match set_idx {
            0 => &self.bc1,
            1 => &self.bc2,
            2 => &self.bc3,
            3 => &self.bc4,
            _ => panic!("Invalid set index: {}", set_idx),
        };
        if let Some(off) = offset {
            bc.match_subsequence(seq, pos, pos + bc.len() + off)
        } else {
            bc.match_subsequence(seq, pos, pos + bc.len())
        }
    }

    /// Builds a full barcode from the 4 barcode indices
    pub fn build_barcode(
        &self,
        b1_idx: usize,
        b2_idx: usize,
        b3_idx: usize,
        b4_idx: usize,
    ) -> Vec<u8> {
        let mut bc =
            Vec::with_capacity(self.bc1.len() + self.bc2.len() + self.bc3.len() + self.bc4.len());
        bc.extend_from_slice(
            self.bc1
                .get_barcode(b1_idx, self.linkers)
                .expect("Invalid barcode index in bc1"),
        );
        bc.extend_from_slice(
            self.bc2
                .get_barcode(b2_idx, self.linkers)
                .expect("Invalid barcode index in bc2"),
        );
        bc.extend_from_slice(
            self.bc3
                .get_barcode(b3_idx, self.linkers)
                .expect("Invalid barcode index in bc3"),
        );
        bc.extend_from_slice(
            self.bc4
                .get_barcode(b4_idx, self.linkers)
                .expect("Invalid barcode index in bc4"),
        );
        bc
    }
}

#[cfg(test)]
mod testing {

    use super::*;

    const TEST_PATH: &str = "data/config_v3.yaml";

    #[test]
    fn load_yaml() {
        let config = Config::from_file(TEST_PATH, false, false);
        assert!(config.is_ok());
    }

    #[test]
    fn load_yaml_exact() {
        let config = Config::from_file(TEST_PATH, true, false);
        assert!(config.is_ok());
    }

    #[test]
    fn barcode_lengths() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        assert_eq!(config.bc1.len(), 8 + 3);
        assert_eq!(config.bc2.len(), 6 + 3);
        assert_eq!(config.bc3.len(), 6 + 5);
        assert_eq!(config.bc4.len(), 8);
    }

    #[test]
    fn barcode_lengths_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();
        assert_eq!(config.bc1.len(), 8 + 3);
        assert_eq!(config.bc2.len(), 6 + 3);
        assert_eq!(config.bc3.len(), 6 + 5);
        assert_eq!(config.bc4.len(), 8);
    }

    #[test]
    fn barcode_sequences() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();

        assert_eq!(config.bc1.get_barcode(0, true).unwrap(), b"AGAAACCAATG");
        assert_eq!(config.bc1.get_barcode(95, true).unwrap(), b"TCTTTGACATG");
        assert_eq!(config.bc1.get_barcode(96, true), None);

        assert_eq!(config.bc1.get_barcode(0, false).unwrap(), b"AGAAACCA");
        assert_eq!(config.bc1.get_barcode(95, false).unwrap(), b"TCTTTGAC");
        assert_eq!(config.bc1.get_barcode(96, false), None);

        assert_eq!(config.bc2.get_barcode(0, true).unwrap(), b"TCTGTGGAG");
        assert_eq!(config.bc2.get_barcode(95, true).unwrap(), b"GTAATCGAG");
        assert_eq!(config.bc2.get_barcode(96, true), None);

        assert_eq!(config.bc2.get_barcode(0, false).unwrap(), b"TCTGTG");
        assert_eq!(config.bc2.get_barcode(95, false).unwrap(), b"GTAATC");
        assert_eq!(config.bc2.get_barcode(96, false), None);

        assert_eq!(config.bc3.get_barcode(0, true).unwrap(), b"AAAGTGTCGAG");
        assert_eq!(config.bc3.get_barcode(95, true).unwrap(), b"CTGAAGTCGAG");
        assert_eq!(config.bc3.get_barcode(96, false), None);

        assert_eq!(config.bc3.get_barcode(0, false).unwrap(), b"AAAGTG");
        assert_eq!(config.bc3.get_barcode(95, false).unwrap(), b"CTGAAG");
        assert_eq!(config.bc3.get_barcode(96, false), None);

        assert_eq!(config.bc4.get_barcode(0, true).unwrap(), b"CTGGGTAT");
        assert_eq!(config.bc4.get_barcode(95, true).unwrap(), b"AAACTACA");
        assert_eq!(config.bc4.get_barcode(96, true), None);

        assert_eq!(config.bc4.get_barcode(0, false).unwrap(), b"CTGGGTAT");
        assert_eq!(config.bc4.get_barcode(95, false).unwrap(), b"AAACTACA");
        assert_eq!(config.bc4.get_barcode(96, false), None);
    }

    #[test]
    fn barcode_sequences_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();

        assert_eq!(config.bc1.get_barcode(0, true).unwrap(), b"AGAAACCAATG");
        assert_eq!(config.bc1.get_barcode(95, true).unwrap(), b"TCTTTGACATG");
        assert_eq!(config.bc1.get_barcode(96, true), None);

        assert_eq!(config.bc2.get_barcode(0, true).unwrap(), b"TCTGTGGAG");
        assert_eq!(config.bc2.get_barcode(95, true).unwrap(), b"GTAATCGAG");
        assert_eq!(config.bc2.get_barcode(96, true), None);

        assert_eq!(config.bc3.get_barcode(0, true).unwrap(), b"AAAGTGTCGAG");
        assert_eq!(config.bc3.get_barcode(95, true).unwrap(), b"CTGAAGTCGAG");
        assert_eq!(config.bc3.get_barcode(96, true), None);

        assert_eq!(config.bc4.get_barcode(0, true).unwrap(), b"CTGGGTAT");
        assert_eq!(config.bc4.get_barcode(95, true).unwrap(), b"AAACTACA");
        assert_eq!(config.bc4.get_barcode(96, true), None);
    }

    #[test]
    fn construct_building_a() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        let bc = config.build_barcode(0, 0, 0, 0);
        let exp = [
            "AGAAACCA".as_bytes(),
            "TCTGTG".as_bytes(),
            "AAAGTG".as_bytes(),
            "CTGGGTAT".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }

    #[test]
    fn construct_building_b() {
        let config = Config::from_file(TEST_PATH, false, false).unwrap();
        let bc = config.build_barcode(0, 95, 0, 95);
        let exp = [
            "AGAAACCA".as_bytes(),
            "GTAATC".as_bytes(),
            "AAAGTG".as_bytes(),
            "AAACTACA".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }

    #[test]
    fn construct_building_a_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();
        let bc = config.build_barcode(0, 0, 0, 0);
        let exp = [
            "AGAAACCA".as_bytes(),
            "TCTGTG".as_bytes(),
            "AAAGTG".as_bytes(),
            "CTGGGTAT".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }

    #[test]
    fn construct_building_b_exact() {
        let config = Config::from_file(TEST_PATH, true, false).unwrap();
        let bc = config.build_barcode(0, 95, 0, 95);
        let exp = [
            "AGAAACCA".as_bytes(),
            "GTAATC".as_bytes(),
            "AAAGTG".as_bytes(),
            "AAACTACA".as_bytes(),
        ]
        .concat();
        assert_eq!(bc, exp);
    }
}
