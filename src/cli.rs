use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about)]
pub struct Cli {
    /// Input file for R1
    #[clap(short = 'i', long, value_parser)]
    pub r1: String,

    /// Input file for R2
    #[clap(short = 'I', long, value_parser)]
    pub r2: String,

    /// Output file prefix (output files will be named <prefix>_R[12].fq.gz)
    #[clap(short = 'p', long, value_parser)]
    pub prefix: String,

    /// The amount of nucleotides away from the start of R1 to accept a barcode
    #[clap(short = 's', long, default_value = "5")]
    pub offset: usize,

    /// The yaml config file describing the file paths of the 4 barcodes and the spacers
    #[clap(short = 'c', long, value_parser)]
    pub config: String,

    /// The length of the UMI
    #[clap(short = 'u', long, default_value = "12")]
    pub umi_len: usize,
}
