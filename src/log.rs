use anyhow::Result;
use clap::crate_version;
use serde::Serialize;

#[derive(Debug, Serialize)]
struct Statistics {
    total_reads: usize,
    passing_reads: usize,
    fraction_passing: f64,
    num_filtered_1: usize,
    num_filtered_2: usize,
    num_filtered_3: usize,
    num_filtered_4: usize,
}

#[derive(Debug, Serialize)]
struct Timing {
    timestamp: String,
    elapsed_time: f64,
}

#[derive(Debug, Serialize)]
struct FileIO {
    readpath_r1: String,
    readpath_r2: String,
    writepath_r1: String,
    writepath_r2: String,
}

#[derive(Debug, Serialize)]
struct Parameters {
    offset: usize,
    umi_len: usize,
    exact_matching: bool,
    pipspeak_version: String,
}

#[derive(Debug, Serialize)]
/// A struct to hold the information about the run
pub struct Log {
    parameters: Parameters,
    file_io: FileIO,
    statistics: Statistics,
    timing: Timing,
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
        offset: usize,
        umi_len: usize,
        exact_matching: bool,
    ) -> Self {
        let pipspeak_version = crate_version!().to_owned();
        let fraction_passing = passing_reads as f64 / total_reads as f64;
        let statistics = Statistics {
            total_reads,
            passing_reads,
            fraction_passing,
            num_filtered_1,
            num_filtered_2,
            num_filtered_3,
            num_filtered_4,
        };
        let file_io = FileIO {
            readpath_r1,
            readpath_r2,
            writepath_r1,
            writepath_r2,
        };
        let parameters = Parameters {
            offset,
            umi_len,
            exact_matching,
            pipspeak_version,
        };
        let timing = Timing {
            timestamp,
            elapsed_time,
        };
        Self {
            parameters,
            file_io,
            statistics,
            timing,
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
