use std::{
    fs::File,
    io::{BufWriter, Write},
};

use anyhow::Result;
use hashbrown::HashSet;
use serde::Serialize;

#[derive(Debug, Default, Serialize, Clone)]
pub struct Statistics {
    pub total_reads: usize,
    pub passing_reads: usize,
    pub fraction_passing: f64,
    pub whitelist_size: usize,
    pub num_filtered_1: usize,
    pub num_filtered_2: usize,
    pub num_filtered_3: usize,
    pub num_filtered_4: usize,
    #[serde(skip)]
    pub whitelist: HashSet<Vec<u8>>,
}
impl Statistics {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn calculate_metrics(&mut self) {
        self.fraction_passing = self.passing_reads as f64 / self.total_reads as f64;
        self.whitelist_size = self.whitelist.len();
    }
    pub fn whitelist_to_file(&self, file: &str) -> Result<()> {
        let mut writer = File::create(file).map(BufWriter::new)?;
        for seq in &self.whitelist {
            writer.write(seq)?;
            writer.write(b"\n")?;
        }
        Ok(())
    }
}

#[derive(Debug, Serialize)]
pub struct Timing {
    pub timestamp: String,
    pub elapsed_time: f64,
}

#[derive(Debug, Serialize)]
pub struct FileIO {
    pub readpath_r1: String,
    pub readpath_r2: String,
    pub writepath_r1: String,
    pub writepath_r2: String,
    pub whitelist_path: String,
}

#[derive(Debug, Serialize)]
pub struct Parameters {
    pub offset: usize,
    pub umi_len: usize,
    pub exact_matching: bool,
    pub write_linkers: bool,
    pub pipspeak_version: String,
}

#[derive(Debug, Serialize)]
/// A struct to hold the information about the run
pub struct Log {
    pub parameters: Parameters,
    pub file_io: FileIO,
    pub statistics: Statistics,
    pub timing: Timing,
}
impl Log {
    pub fn stderr(&self) -> Result<()> {
        let yaml = serde_yaml::to_string(&self)?;
        eprint!("{}", yaml);
        Ok(())
    }

    pub fn to_file(&self, path: &str) -> Result<()> {
        let yaml = serde_yaml::to_string(&self)?;
        std::fs::write(path, yaml)?;
        Ok(())
    }
}
