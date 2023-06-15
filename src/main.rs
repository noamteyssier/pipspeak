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
use fxread::initialize_reader;
use log::Log;
use std::{fs::File, io::Write, time::Instant};

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

fn main() -> Result<()> {
    let args = Cli::parse();
    let config = Config::from_file(&args.config, args.exact)?;
    let r1 = initialize_reader(&args.r1)?;
    let r2 = initialize_reader(&args.r2)?;

    let r1_filename = args.prefix.clone() + "_R1.fq.gz";
    let r2_filename = args.prefix.clone() + "_R2.fq.gz";
    let log_filename = args.prefix.clone() + "_log.yaml";

    let mut num_filtered_1 = 0;
    let mut num_filtered_2 = 0;
    let mut num_filtered_3 = 0;
    let mut num_filtered_4 = 0;
    let mut total_reads = 0;
    let mut passing_reads = 0;
    let start_time = Instant::now();
    let timestamp = Local::now().to_string();

    let record_iter = r1
        .zip(r2)
        .inspect(|_| total_reads += 1)
        .filter_map(|(rec1, rec2)| {
            if let Some((pos, b1_idx)) =
                config.match_subsequence(rec1.seq(), 0, 0, Some(args.offset))
            {
                Some((rec1, rec2, pos, b1_idx))
            } else {
                num_filtered_1 += 1;
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx)| {
            if let Some((new_pos, b2_idx)) = config.match_subsequence(rec1.seq(), 1, pos, None) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx))
            } else {
                num_filtered_2 += 1;
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx, b2_idx)| {
            if let Some((new_pos, b3_idx)) = config.match_subsequence(&rec1.seq(), 2, pos, None) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx, b3_idx))
            } else {
                num_filtered_3 += 1;
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx, b2_idx, b3_idx)| {
            if let Some((new_pos, b4_idx)) = config.match_subsequence(&rec1.seq(), 3, pos, None) {
                passing_reads += 1;
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx, b3_idx, b4_idx))
            } else {
                num_filtered_4 += 1;
                None
            }
        })
        .map(|(rec1, rec2, pos, b1_idx, b2_idx, b3_idx, b4_idx)| {
            let umi = &rec1.seq()[pos..pos + args.umi_len];
            (
                b1_idx,
                b2_idx,
                b3_idx,
                b4_idx,
                umi.to_vec(),
                pos + args.umi_len,
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

    let mut r1_writer = GzEncoder::new(File::create(&r1_filename)?, Compression::default());
    let mut r2_writer = GzEncoder::new(File::create(&r2_filename)?, Compression::default());

    for (construct_seq, construct_qual, rec1, rec2) in record_iter {
        write_to_fastq(&mut r1_writer, rec1.id(), &construct_seq, &construct_qual)?;
        write_to_fastq(&mut r2_writer, rec2.id(), rec2.seq(), rec2.qual().unwrap())?;
    }

    let elapsed = start_time.elapsed().as_secs_f64();

    let log = Log::new(
        total_reads,
        passing_reads,
        num_filtered_1,
        num_filtered_2,
        num_filtered_3,
        num_filtered_4,
        args.r1.clone(),
        args.r2.clone(),
        r1_filename.clone(),
        r2_filename.clone(),
        timestamp,
        elapsed,
    );

    if !args.quiet {
        log.stderr()?;
    }
    log.to_file(&log_filename)?;

    Ok(())
}
