use anyhow::{anyhow, Result};
use yaml_rust::{Yaml, YamlLoader};

use crate::barcodes::{Barcodes, Spacer};

pub struct Config {
    bc1: Barcodes,
    bc2: Barcodes,
    bc3: Barcodes,
    bc4: Barcodes,
}
impl Config {
    pub fn from_file(path: &str) -> Result<Self> {
        let contents = std::fs::read_to_string(path)?;
        let yaml = YamlLoader::load_from_str(&contents)?;
        Self::from_yaml(yaml)
    }

    pub fn from_yaml(yaml: Vec<Yaml>) -> Result<Self> {
        let spacer1 = Self::load_spacer(&yaml[0], "s1")?;
        let spacer2 = Self::load_spacer(&yaml[0], "s2")?;
        let spacer3 = Self::load_spacer(&yaml[0], "s3")?;
        let bc1 = Self::load_barcode(&yaml[0], "bc1", Some(&spacer1))?;
        let bc2 = Self::load_barcode(&yaml[0], "bc2", Some(&spacer2))?;
        let bc3 = Self::load_barcode(&yaml[0], "bc3", Some(&spacer3))?;
        let bc4 = Self::load_barcode(&yaml[0], "bc4", None)?;
        Ok(Self { bc1, bc2, bc3, bc4 })
    }

    fn load_spacer(yaml: &Yaml, spacer: &str) -> Result<Spacer> {
        yaml["spacers"][spacer]
            .as_str()
            .map(Spacer::from_str)
            .ok_or(anyhow!("Spacer {} not found", spacer))
    }

    fn load_barcode(yaml: &Yaml, barcode: &str, spacer: Option<&Spacer>) -> Result<Barcodes> {
        if let Some(bc) = yaml["barcodes"][barcode].as_str() {
            if let Some(spacer) = spacer {
                Barcodes::from_file_with_spacer(bc, spacer)
            } else {
                Barcodes::from_file(bc)
            }
        } else {
            Err(anyhow!("Barcode {} not found", barcode))
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
                .get_barcode(b1_idx)
                .expect("Invalid barcode index in bc1"),
        );
        bc.extend_from_slice(
            self.bc2
                .get_barcode(b2_idx)
                .expect("Invalid barcode index in bc2"),
        );
        bc.extend_from_slice(
            self.bc3
                .get_barcode(b3_idx)
                .expect("Invalid barcode index in bc3"),
        );
        bc.extend_from_slice(
            self.bc4
                .get_barcode(b4_idx)
                .expect("Invalid barcode index in bc4"),
        );
        bc
    }
}
