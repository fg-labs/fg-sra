//! Command-line argument definitions for fg-sra.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::manager::VdbManager;

/// High-performance SRA toolkit.
#[derive(Debug, Parser)]
#[command(name = "fg-sra", version, about)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

/// Available subcommands.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Convert SRA archives to SAM or BAM format.
    #[command(name = "tosam")]
    ToSam(ToSam),

    /// Pre-populate the local VDB reference sequence cache.
    #[command(name = "cache-refs")]
    CacheRefs(CacheRefs),
}

/// Convert NCBI SRA archives to SAM or BAM format, replacing `sam-dump`
/// with multi-threaded processing for significantly higher throughput.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Parser)]
pub struct ToSam {
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

    /// Use `SEQ_ID` instead of NAME for RNAME.
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

    /// Append `.SPOT_GROUP` to QNAME.
    #[arg(short = 'g', long = "spot-group")]
    pub spot_group: bool,

    /// Prepend prefix to QNAME.
    #[arg(short = 'p', long = "prefix")]
    pub prefix: Option<String>,

    /// Reverse unaligned reads per read type.
    #[arg(long = "reverse")]
    pub reverse: bool,

    /// Quality score quantization (e.g. "1:10,10:20,20:30,30:40").
    #[arg(short = 'Q', long = "qual-quant")]
    pub qual_quant: Option<String>,

    /// Output alignment ID in XI:i tag.
    #[arg(long = "XI")]
    pub xi_tag: bool,

    // ── Performance options ───────────────────────────────────────────
    /// Number of worker threads (default: available cores).
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,

    /// Explicit VDB cursor pool size (number of cursors).
    /// Default: one cursor per thread.
    #[arg(long = "pool-size")]
    pub pool_size: Option<usize>,
}

/// Output format for converted records.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum OutputFormat {
    Sam,
    Bam,
}

/// Pre-populate the local VDB reference sequence cache for one or more
/// SRA accessions. This resolves and caches reference sequences serially,
/// avoiding SDL resolver failures that occur under heavy concurrent load.
#[derive(Debug, Parser)]
pub struct CacheRefs {
    /// SRA accession(s) or file path(s) to resolve references for.
    #[arg(required = true)]
    pub accessions: Vec<String>,
}

impl Cli {
    /// Dispatch to the appropriate subcommand.
    pub fn execute(&self) -> Result<()> {
        match &self.command {
            Command::ToSam(cmd) => cmd.execute(),
            Command::CacheRefs(cmd) => cmd.execute(),
        }
    }
}

/// Open a VDB database for reading from an accession or file path.
///
/// Creates a VDB manager, disables the pagemap thread (best-effort), and
/// opens the database for reading.
fn open_database(accession: &str) -> Result<VDatabase> {
    let mgr = VdbManager::make_read().context("failed to create VDB manager")?;
    mgr.disable_pagemap_thread().ok();
    mgr.open_db_read(accession).with_context(|| format!("failed to open database: {accession}"))
}

impl CacheRefs {
    /// Resolve and cache reference dependencies for each accession.
    pub fn execute(&self) -> Result<()> {
        let mut num_failed = 0usize;

        for accession in &self.accessions {
            if let Err(e) = self.process_accession(accession) {
                eprintln!("[cache-refs] {accession}: ERROR: {e:#}");
                num_failed += 1;
            }
        }

        let num_ok = self.accessions.len() - num_failed;
        eprintln!("[cache-refs] Done: {num_ok} accession(s) processed, {num_failed} failed");

        if num_failed > 0 {
            anyhow::bail!(
                "{num_failed} of {} accession(s) failed to resolve",
                self.accessions.len()
            );
        }

        Ok(())
    }

    /// Process a single accession: open database, list dependencies, report results.
    #[allow(clippy::unused_self)]
    fn process_accession(&self, accession: &str) -> Result<()> {
        eprintln!("[cache-refs] {accession}: resolving dependencies...");

        let db = open_database(accession)?;

        let deps = db
            .list_dependencies(false)
            .with_context(|| format!("failed to list dependencies: {accession}"))?;

        let infos = deps
            .all_info()
            .with_context(|| format!("failed to read dependency info: {accession}"))?;

        let num_local = infos.iter().filter(|d| d.local).count();
        let num_resolved = infos.len() - num_local;

        eprintln!(
            "[cache-refs] {accession}: {} dependenc{} ({} local, {} resolved)",
            infos.len(),
            if infos.len() == 1 { "y" } else { "ies" },
            num_local,
            num_resolved,
        );

        Ok(())
    }
}

impl ToSam {
    /// Run the conversion with the parsed CLI options.
    pub fn execute(&self) -> Result<()> {
        for accession in &self.accessions {
            self.process_accession(accession)?;
        }
        Ok(())
    }

    /// Process a single SRA accession or file path.
    fn process_accession(&self, accession: &str) -> Result<()> {
        use crate::aligned::{AlignConfig, process_aligned_table};
        use crate::header::generate_header;
        use crate::output::OutputWriter;
        use crate::progress::ProgressLogger;
        use crate::record::FormatOptions;
        use crate::unaligned::process_unaligned_reads;

        const PROGRESS_INTERVAL: u64 = 1_000_000;

        let db = open_database(accession)?;

        let output_mode = if self.output_format == OutputFormat::Bam {
            crate::record::OutputMode::Bam
        } else if self.fasta {
            crate::record::OutputMode::Fasta
        } else if self.fastq {
            crate::record::OutputMode::Fastq
        } else {
            crate::record::OutputMode::Sam
        };

        let mut writer = if output_mode == crate::record::OutputMode::Bam {
            match &self.output_file {
                Some(path) => OutputWriter::bam_from_path(path)?,
                None => OutputWriter::bam_stdout(),
            }
        } else {
            let compression = if self.gzip {
                crate::output::CompressionMode::Gzip
            } else if self.bzip2 {
                crate::output::CompressionMode::Bzip2
            } else {
                crate::output::CompressionMode::None
            };
            match &self.output_file {
                Some(path) => OutputWriter::from_path_with_compression(path, compression)?,
                None => OutputWriter::stdout_with_compression(compression),
            }
        };

        // Generate the header for SAM and BAM modes.
        let header_text = if !self.no_header
            && matches!(
                output_mode,
                crate::record::OutputMode::Sam | crate::record::OutputMode::Bam
            ) {
            let h = generate_header(
                &db,
                self.header,
                self.seqid,
                &self.header_comment,
                self.header_file.as_deref(),
            )?;
            writer.write_header(&h)?;
            Some(h)
        } else {
            None
        };

        // Build ref_name → ref_id map for BAM output.
        let ref_name_to_id = if output_mode == crate::record::OutputMode::Bam {
            header_text.as_deref().map(|h| crate::header::build_ref_id_map(&db, h)).transpose()?
        } else {
            None
        };

        let qual_table =
            self.qual_quant.as_deref().map(crate::quality::parse_qual_quant).transpose()?;
        let opts = FormatOptions {
            prefix: self.prefix.as_deref(),
            spot_group_in_name: self.spot_group,
            xi_tag: self.xi_tag,
            reverse_unaligned: self.reverse,
            omit_quality: self.omit_quality,
            qual_quant: qual_table.as_ref(),
            output_mode,
            ref_name_to_id: ref_name_to_id.as_ref(),
        };

        let align_config = AlignConfig {
            use_seqid: self.seqid,
            use_long_cigar: self.cigar_long,
            primary_only: self.primary,
            min_mapq: self.min_mapq,
            num_threads: self.threads.unwrap_or_else(|| {
                std::thread::available_parallelism().map(std::num::NonZero::get).unwrap_or(1)
            }),
            pool_size_override: self.pool_size,
            opts: &opts,
            regions: &self.aligned_region,
        };

        // Aligned reads (unless --unaligned-spots-only).
        if !self.unaligned_spots_only {
            process_aligned_table(&db, &mut writer, &align_config, PROGRESS_INTERVAL)?;
        }

        // Unaligned reads (if requested).
        if self.unaligned || self.unaligned_spots_only {
            let unaligned_progress = ProgressLogger::new(0, PROGRESS_INTERVAL);
            process_unaligned_reads(
                &db,
                &mut writer,
                &opts,
                self.unaligned_spots_only,
                &unaligned_progress,
            )?;
        }

        writer.finish()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::*;

    /// Parse a `ToSam` subcommand from a slice of arguments (program name auto-prepended).
    fn parse(args: &[&str]) -> ToSam {
        let mut full = vec!["fg-sra", "tosam"];
        full.extend_from_slice(args);
        let cli = Cli::parse_from(full);
        match cli.command {
            Command::ToSam(cmd) => cmd,
            Command::CacheRefs(_) => panic!("expected ToSam command"),
        }
    }

    #[test]
    fn test_minimal_args() {
        let cmd = parse(&["SRR123456"]);
        assert_eq!(cmd.accessions, vec!["SRR123456"]);
        assert!(!cmd.unaligned);
        assert!(!cmd.primary);
        assert_eq!(cmd.output_format, OutputFormat::Sam);
    }

    #[test]
    fn test_multiple_accessions() {
        let cmd = parse(&["SRR111", "SRR222", "SRR333"]);
        assert_eq!(cmd.accessions, vec!["SRR111", "SRR222", "SRR333"]);
    }

    #[test]
    fn test_primary_and_unaligned() {
        let cmd = parse(&["-1", "-u", "SRR123456"]);
        assert!(cmd.primary);
        assert!(cmd.unaligned);
    }

    #[test]
    fn test_bam_output_format() {
        let cmd = parse(&["--output-format", "bam", "SRR123456"]);
        assert_eq!(cmd.output_format, OutputFormat::Bam);
    }

    #[test]
    fn test_output_file() {
        let cmd = parse(&["--output-file", "/tmp/out.sam", "SRR123456"]);
        assert_eq!(cmd.output_file.as_deref(), Some(std::path::Path::new("/tmp/out.sam")));
    }

    #[test]
    fn test_compression_flags() {
        let cmd = parse(&["--gzip", "SRR123456"]);
        assert!(cmd.gzip);
        assert!(!cmd.bzip2);

        let cmd = parse(&["--bzip2", "SRR123456"]);
        assert!(!cmd.gzip);
        assert!(cmd.bzip2);
    }

    #[test]
    fn test_cigar_and_formatting() {
        let cmd = parse(&["--cigar-long", "--spot-group", "--prefix", "PRE", "SRR123456"]);
        assert!(cmd.cigar_long);
        assert!(cmd.spot_group);
        assert_eq!(cmd.prefix.as_deref(), Some("PRE"));
    }

    #[test]
    fn test_header_options() {
        let cmd = parse(&[
            "--no-header",
            "--header-comment",
            "line one",
            "--header-comment",
            "line two",
            "SRR123456",
        ]);
        assert!(cmd.no_header);
        assert_eq!(cmd.header_comment, vec!["line one", "line two"]);
    }

    #[test]
    fn test_aligned_region() {
        let cmd =
            parse(&["--aligned-region", "chr1:1000-2000", "--aligned-region", "chr2", "SRR123456"]);
        assert_eq!(cmd.aligned_region, vec!["chr1:1000-2000", "chr2"]);
    }

    #[test]
    fn test_min_mapq() {
        let cmd = parse(&["--min-mapq", "30", "SRR123456"]);
        assert_eq!(cmd.min_mapq, Some(30));
    }

    #[test]
    fn test_qual_quant() {
        let cmd = parse(&["--qual-quant", "1:10,10:20,20:30,30:40", "SRR123456"]);
        assert_eq!(cmd.qual_quant.as_deref(), Some("1:10,10:20,20:30,30:40"));
    }

    #[test]
    fn test_threads() {
        let cmd = parse(&["-t", "8", "SRR123456"]);
        assert_eq!(cmd.threads, Some(8));
    }

    #[test]
    fn test_fasta_fastq() {
        let cmd = parse(&["--fasta", "SRR123456"]);
        assert!(cmd.fasta);
        assert!(!cmd.fastq);

        let cmd = parse(&["--fastq", "SRR123456"]);
        assert!(!cmd.fasta);
        assert!(cmd.fastq);
    }

    #[test]
    fn test_short_flags() {
        let cmd = parse(&["-n", "-r", "-s", "-c", "-g", "-o", "-t", "4", "SRR123456"]);
        assert!(cmd.no_header);
        assert!(cmd.header);
        assert!(cmd.seqid);
        assert!(cmd.cigar_long);
        assert!(cmd.spot_group);
        assert!(cmd.omit_quality);
        assert_eq!(cmd.threads, Some(4));
    }

    #[test]
    fn test_pool_size_flag() {
        let cmd = parse(&["--pool-size", "4", "SRR123456"]);
        assert_eq!(cmd.pool_size, Some(4));

        let cmd = parse(&["SRR123456"]);
        assert_eq!(cmd.pool_size, None);
    }

    fn parse_cache_refs(args: &[&str]) -> CacheRefs {
        let mut full = vec!["fg-sra", "cache-refs"];
        full.extend_from_slice(args);
        let cli = Cli::parse_from(full);
        match cli.command {
            Command::CacheRefs(cmd) => cmd,
            Command::ToSam(_) => panic!("expected CacheRefs command"),
        }
    }

    #[test]
    fn test_cache_refs_single_accession() {
        let cmd = parse_cache_refs(&["SRR123456"]);
        assert_eq!(cmd.accessions, vec!["SRR123456"]);
    }

    #[test]
    fn test_cache_refs_multiple_accessions() {
        let cmd = parse_cache_refs(&["SRR111", "SRR222", "SRR333"]);
        assert_eq!(cmd.accessions, vec!["SRR111", "SRR222", "SRR333"]);
    }

    #[test]
    fn test_cache_refs_missing_accession_fails() {
        let result = Cli::try_parse_from(["fg-sra", "cache-refs"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_missing_subcommand_fails() {
        let result = Cli::try_parse_from(["fg-sra"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_missing_accession_fails() {
        let result = Cli::try_parse_from(["fg-sra", "tosam"]);
        assert!(result.is_err());
    }
}
