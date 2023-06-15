use anyhow::Result;
use serde::Serialize;

#[derive(Debug, Serialize)]
/// A struct to hold the information about the run
pub struct Log {
    timestamp: String,
    elapsed_time: f64,
    total_reads: usize,
    passing_reads: usize,
    fraction_passing: f64,
    num_filtered_1: usize,
    num_filtered_2: usize,
    num_filtered_3: usize,
    num_filtered_4: usize,
    readpath_r1: String,
    readpath_r2: String,
    writepath_r1: String,
    writepath_r2: String,
}
impl Log {
    pub fn new(
        total_reads: usize,
        passing_reads: usize,
        num_filtered_1: usize,
        num_filtered_2: usize,
        num_filtered_3: usize,
        num_filtered_4: usize,
        readpath_r1: String,
        readpath_r2: String,
        writepath_r1: String,
        writepath_r2: String,
        timestamp: String,
        elapsed_time: f64,
    ) -> Self {
        let fraction_passing = passing_reads as f64 / total_reads as f64;
        Self {
            total_reads,
            passing_reads,
            fraction_passing,
            num_filtered_1,
            num_filtered_2,
            num_filtered_3,
            num_filtered_4,
            readpath_r1,
            readpath_r2,
            writepath_r1,
            writepath_r2,
            timestamp,
            elapsed_time,
        }
    }

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
