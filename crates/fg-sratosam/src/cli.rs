//! Command-line argument definitions for fg-sratosam.

use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

/// High-performance SRA-to-SAM/BAM converter.
///
/// Converts NCBI SRA archives to SAM or BAM format, replacing `sam-dump`
/// with multi-threaded processing for significantly higher throughput.
#[derive(Debug, Parser)]
#[command(name = "fg-sratosam", version, about)]
pub struct Cli {
    /// SRA accession(s) or file path(s) to convert.
    #[arg(required = true)]
    pub accessions: Vec<String>,

    // ── Core options ──────────────────────────────────────────────────

    /// Output unaligned reads along with aligned reads.
    #[arg(short = 'u', long = "unaligned")]
    pub unaligned: bool,

    /// Output only primary alignments.
    #[arg(short = '1', long = "primary")]
    pub primary: bool,

    /// Filter by genomic region (repeatable). Format: name[:from-to]
    #[arg(long = "aligned-region")]
    pub aligned_region: Vec<String>,

    /// Minimum MAPQ to output.
    #[arg(long = "min-mapq")]
    pub min_mapq: Option<u32>,

    /// Suppress SAM header in output.
    #[arg(short = 'n', long = "no-header")]
    pub no_header: bool,

    /// Reconstruct header from metadata.
    #[arg(short = 'r', long = "header")]
    pub header: bool,

    /// Use external header file.
    #[arg(long = "header-file")]
    pub header_file: Option<PathBuf>,

    /// Add @CO comment line(s) to header (repeatable).
    #[arg(long = "header-comment")]
    pub header_comment: Vec<String>,

    /// Use SEQ_ID instead of NAME for RNAME.
    #[arg(short = 's', long = "seqid")]
    pub seqid: bool,

    /// Output only unaligned spots (spots with no alignments).
    #[arg(long = "unaligned-spots-only")]
    pub unaligned_spots_only: bool,

    // ── Output options ────────────────────────────────────────────────

    /// Write to file instead of stdout.
    #[arg(long = "output-file")]
    pub output_file: Option<PathBuf>,

    /// Output format.
    #[arg(long = "output-format", default_value = "sam")]
    pub output_format: OutputFormat,

    /// Compress SAM output with gzip.
    #[arg(long = "gzip")]
    pub gzip: bool,

    /// Compress SAM output with bzip2.
    #[arg(long = "bzip2")]
    pub bzip2: bool,

    /// Output in FASTA format.
    #[arg(long = "fasta")]
    pub fasta: bool,

    /// Output in FASTQ format.
    #[arg(long = "fastq")]
    pub fastq: bool,

    /// Omit quality values.
    #[arg(short = 'o', long = "omit-quality")]
    pub omit_quality: bool,

    // ── Formatting options ────────────────────────────────────────────

    /// Use long CIGAR form.
    #[arg(short = 'c', long = "cigar-long")]
    pub cigar_long: bool,

    /// Output `=` for bases matching reference.
    #[arg(long = "hide-identical")]
    pub hide_identical: bool,

    /// Append .SPOT_GROUP to QNAME.
    #[arg(short = 'g', long = "spot-group")]
    pub spot_group: bool,

    /// Prepend prefix to QNAME.
    #[arg(short = 'p', long = "prefix")]
    pub prefix: Option<String>,

    /// Reverse unaligned reads per read type.
    #[arg(long = "reverse")]
    pub reverse: bool,

    /// Compute and output MD tag.
    #[arg(long = "with-md-flag")]
    pub with_md_flag: bool,

    /// Quality score quantization (e.g. "1:10,10:20,20:30,30:40").
    #[arg(short = 'Q', long = "qual-quant")]
    pub qual_quant: Option<String>,

    /// Output alignment ID in XI:i tag.
    #[arg(long = "XI")]
    pub xi_tag: bool,

    /// Detect RNA splicing (replace D with N in CIGAR, add XS:A tag).
    #[arg(long = "rna-splicing")]
    pub rna_splicing: bool,

    /// Mismatch tolerance for splice site detection (0, 1, or 2).
    #[arg(long = "rna-splice-level", default_value = "0")]
    pub rna_splice_level: u8,

    /// Log splice events to file.
    #[arg(long = "rna-splice-log")]
    pub rna_splice_log: Option<PathBuf>,

    // ── Performance options ───────────────────────────────────────────

    /// Number of worker threads (default: available cores).
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,
}

/// Output format for converted records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum OutputFormat {
    Sam,
    Bam,
}

impl Cli {
    /// Run the conversion with the parsed CLI options.
    pub fn execute(&self) -> Result<()> {
        // TODO: implement conversion pipeline
        anyhow::bail!("not yet implemented")
    }
}
