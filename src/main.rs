mod barcodes;
mod cli;
mod config;
mod log;

use anyhow::Result;
use chrono::Local;
use clap::Parser;
use cli::Cli;
use config::Config;
use flate2::{write::GzEncoder, Compression};
use fxread::{initialize_reader, FastxRead, Record};
use indicatif::ProgressBar;
use log::{FileIO, Log, Parameters, Statistics, Timing};
use std::{
    fs::File,
    io::Write,
    time::{Duration, Instant},
};

/// Writes a record to a gzip fastq file
fn write_to_fastq<W: Write>(writer: &mut W, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn parse_records(
    r1: Box<dyn FastxRead<Item = Record>>,
    r2: Box<dyn FastxRead<Item = Record>>,
    r1_out: &mut GzEncoder<File>,
    r2_out: &mut GzEncoder<File>,
    config: &Config,
    offset: usize,
    umi_len: usize,
) -> Result<Statistics> {
    let mut statistics = Statistics::new();
    let pb = ProgressBar::new_spinner();
    pb.enable_steady_tick(Duration::from_millis(100));
    let record_iter = r1
        .zip(r2)
        .inspect(|_| statistics.total_reads += 1)
        .enumerate()
        .map(|(idx, pair)| {
            if idx % 125 == 0 {
                pb.set_message(format!("Processed {} reads", idx));
            }
            pair
        })
        .filter_map(|(rec1, rec2)| {
            if let Some((pos, b1_idx)) = config.match_subsequence(rec1.seq(), 0, 0, Some(offset)) {
                Some((rec1, rec2, pos, b1_idx))
            } else {
                statistics.num_filtered_1 += 1;
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx)| {
            if let Some((new_pos, b2_idx)) = config.match_subsequence(rec1.seq(), 1, pos, None) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx))
            } else {
                statistics.num_filtered_2 += 1;
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx, b2_idx)| {
            if let Some((new_pos, b3_idx)) = config.match_subsequence(&rec1.seq(), 2, pos, None) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx, b3_idx))
            } else {
                statistics.num_filtered_3 += 1;
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx, b2_idx, b3_idx)| {
            if let Some((new_pos, b4_idx)) = config.match_subsequence(&rec1.seq(), 3, pos, None) {
                statistics.passing_reads += 1;
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx, b3_idx, b4_idx))
            } else {
                statistics.num_filtered_4 += 1;
                None
            }
        })
        .map(|(rec1, rec2, pos, b1_idx, b2_idx, b3_idx, b4_idx)| {
            let umi = &rec1.seq()[pos..pos + umi_len];
            (
                b1_idx,
                b2_idx,
                b3_idx,
                b4_idx,
                umi.to_vec(),
                pos + umi_len,
                rec1,
                rec2,
            )
        })
        .map(|(b1_idx, b2_idx, b3_idx, b4_idx, umi, pos, rec1, rec2)| {
            let mut construct_seq = config.build_barcode(b1_idx, b2_idx, b3_idx, b4_idx);
            construct_seq.extend_from_slice(&umi);
            let construct_qual = rec1.qual().unwrap()[pos - construct_seq.len()..pos].to_vec();
            (construct_seq, construct_qual, rec1, rec2)
        });

    for (c_seq, c_qual, rec1, rec2) in record_iter {
        statistics.whitelist.insert(c_seq.clone());
        write_to_fastq(r1_out, rec1.id(), &c_seq, &c_qual)?;
        write_to_fastq(r2_out, rec2.id(), rec2.seq(), rec2.qual().unwrap())?;
    }
    statistics.calculate_metrics();
    pb.finish_with_message(format!(
        "Processed {} reads, {} passed filters ({:.4}%)",
        statistics.total_reads,
        statistics.passing_reads,
        statistics.fraction_passing * 100.0
    ));
    Ok(statistics)
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let config = Config::from_file(&args.config, args.exact, args.linkers)?;
    let r1 = initialize_reader(&args.r1)?;
    let r2 = initialize_reader(&args.r2)?;

    let r1_filename = args.prefix.clone() + "_R1.fq.gz";
    let r2_filename = args.prefix.clone() + "_R2.fq.gz";
    let log_filename = args.prefix.clone() + "_log.yaml";
    let whitelist_filename = args.prefix.clone() + "_whitelist.txt";

    let mut r1_writer = GzEncoder::new(File::create(&r1_filename)?, Compression::default());
    let mut r2_writer = GzEncoder::new(File::create(&r2_filename)?, Compression::default());

    let timestamp = Local::now().to_string();
    let start_time = Instant::now();

    let statistics = parse_records(
        r1,
        r2,
        &mut r1_writer,
        &mut r2_writer,
        &config,
        args.offset,
        args.umi_len,
    )?;
    statistics.whitelist_to_file(&whitelist_filename)?;

    let elapsed_time = start_time.elapsed().as_secs_f64();
    let timing = Timing {
        timestamp,
        elapsed_time,
    };

    let parameters = Parameters {
        offset: args.offset,
        umi_len: args.umi_len,
        exact_matching: args.exact,
        pipspeak_version: env!("CARGO_PKG_VERSION").to_string(),
    };

    let file_io = FileIO {
        readpath_r1: args.r1,
        readpath_r2: args.r2,
        writepath_r1: r1_filename,
        writepath_r2: r2_filename,
        whitelist_path: whitelist_filename,
    };

    let log = Log {
        parameters,
        timing,
        statistics,
        file_io,
    };

    if !args.quiet {
        log.stderr()?;
    }
    log.to_file(&log_filename)?;

    Ok(())
}
