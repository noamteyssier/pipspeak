use anyhow::Result;
use serde::Serialize;

#[derive(Debug, Default, Serialize, Copy, Clone)]
pub struct Statistics {
    pub total_reads: usize,
    pub passing_reads: usize,
    pub fraction_passing: f64,
    pub num_filtered_1: usize,
    pub num_filtered_2: usize,
    pub num_filtered_3: usize,
    pub num_filtered_4: usize,
}
impl Statistics {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn calculate_fraction_passing(&mut self) {
        self.fraction_passing = self.passing_reads as f64 / self.total_reads as f64;
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
}

#[derive(Debug, Serialize)]
pub struct Parameters {
    pub offset: usize,
    pub umi_len: usize,
    pub exact_matching: bool,
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
