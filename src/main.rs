use anyhow::Result;
use fxread::initialize_reader;
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

fn read_barcodes(filename: &str) -> Result<Vec<Vec<u8>>> {
    let reader = File::open(filename).map(BufReader::new)?;
    let mut barcodes = Vec::new();
    for line in reader.lines() {
        let barcode = line.map(|l| l.trim().as_bytes().to_vec())?;
        barcodes.push(barcode);
    }
    Ok(barcodes)
}

/// Check if a sequence contains a barcode as a substring
fn contains_barcode(seq: &[u8], barcode: &[u8]) -> Option<usize> {
    seq.windows(barcode.len())
        .position(|window| window == barcode)
        .map(|pos| pos + barcode.len())
}

fn find_barcode(
    seq: &[u8],
    barcodes: &[Vec<u8>],
    start_pos: usize,
    end_pos: usize,
) -> Option<(usize, usize)> {
    for (bc_idx, barcode) in barcodes.iter().enumerate() {
        if let Some(pos) = contains_barcode(&seq[start_pos..end_pos], barcode) {
            return Some((pos, bc_idx));
        }
    }
    None
}

fn main() -> Result<()> {
    let b1_fn = "data/barcodes_v3/fb_v3_bc1.tsv";
    let b2_fn = "data/barcodes_v3/fb_v3_bc2.tsv";
    let b3_fn = "data/barcodes_v3/fb_v3_bc3.tsv";
    let b4_fn = "data/barcodes_v3/fb_v3_bc4.tsv";

    let b1 = read_barcodes(b1_fn)?;
    let b2 = read_barcodes(b2_fn)?;
    let b3 = read_barcodes(b3_fn)?;
    let b4 = read_barcodes(b4_fn)?;
    let spacers: Vec<&[u8]> = vec![b"ATG", b"GAG", b"TCGAG"];
    let b1 = b1
        .iter()
        .map(|bc| [bc, spacers[0]].concat())
        .collect::<Vec<Vec<u8>>>();
    let b2 = b2
        .iter()
        .map(|bc| [bc, spacers[1]].concat())
        .collect::<Vec<Vec<u8>>>();
    let b3 = b3
        .iter()
        .map(|bc| [bc, spacers[2]].concat())
        .collect::<Vec<Vec<u8>>>();
    let offset = 3;

    let b1_size = b1[0].len();
    let b2_size = b2[0].len();
    let b3_size = b3[0].len();
    let b4_size = b4[0].len();

    let r1_fn = "data/subset_R1.fastq.gz";
    let r2_fn = "data/subset_R2.fastq.gz";
    let r1 = initialize_reader(r1_fn)?;
    let r2 = initialize_reader(r2_fn)?;

    let record_iter = r1
        .zip(r2)
        .filter_map(|(rec1, rec2)| {
            if let Some((pos, b1_idx)) = find_barcode(&rec1.seq(), &b1, 0, b1_size + offset) {
                Some((rec1, rec2, pos, b1_idx))
            } else {
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx)| {
            if let Some((new_pos, b2_idx)) = find_barcode(&rec1.seq(), &b2, pos, pos + b2_size) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx))
            } else {
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx, b2_idx)| {
            if let Some((new_pos, b3_idx)) = find_barcode(&rec1.seq(), &b3, pos, pos + b3_size) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx, b3_idx))
            } else {
                None
            }
        })
        .filter_map(|(rec1, rec2, pos, b1_idx, b2_idx, b3_idx)| {
            if let Some((new_pos, b4_idx)) = find_barcode(&rec1.seq(), &b4, pos, pos + b4_size) {
                Some((rec1, rec2, pos + new_pos, b1_idx, b2_idx, b3_idx, b4_idx))
            } else {
                None
            }
        })
        .map(|(rec1, rec2, pos, b1_idx, b2_idx, b3_idx, b4_idx)| {
            let umi = &rec1.seq()[pos..pos + 12];
            (b1_idx, b2_idx, b3_idx, b4_idx, umi.to_vec(), rec2)
        })
        .map(|(b1_idx, b2_idx, b3_idx, b4_idx, umi, rec2)| {
            let construct = [
                b1[b1_idx].clone(),
                b2[b2_idx].clone(),
                b3[b3_idx].clone(),
                b4[b4_idx].clone(),
                umi,
            ]
            .concat();
            (construct, rec2)
        });

    for (construct, rec2) in record_iter {
        println!(
            "{}\t{}",
            String::from_utf8(construct).unwrap(),
            rec2.seq().len()
        );
    }

    Ok(())
}
