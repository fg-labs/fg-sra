//! Command-line argument definitions for fg-sra.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use fg_sra_vdb::database::VDatabase;
use fg_sra_vdb::manager::{VdbManager, disable_remote_access};

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

    /// Convert an SRA archive to FASTQ, in spot order, with mates paired.
    #[command(name = "fastq")]
    Fastq(crate::fastq::Fastq),

    /// Describe SRA archives: kind, totals, qualities, names and read layout.
    #[command(name = "info")]
    Info(crate::info::Info),
}

/// Convert NCBI SRA archives to SAM or BAM format, replacing `sam-dump`
/// with multi-threaded processing for significantly higher throughput.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Parser)]
pub struct ToSam {
    /// SRA accession(s) or file path(s) to convert. Several accessions are
    /// written, in order, to one output; BAM output takes a single accession,
    /// and SAM output of several accessions requires `--no-header`.
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

    /// Use external header file. As with the stored header, an `@HD SO:coordinate`
    /// in it is written as `SO:unsorted` when the output is not coordinate-sorted.
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
    /// Write to file instead of stdout. The output is written under a temporary
    /// name in the same directory and renamed into place only on success, so an
    /// existing file is kept (and needs room alongside) until the run completes.
    /// Written in place instead: symlinks, devices, pipes, hard-linked files,
    /// files owned by another user or group, and when the temporary file
    /// cannot be created (e.g. a read-only directory).
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

    // ── Access options ────────────────────────────────────────────────
    /// Never use the network: find each accession, and an aligned run's
    /// references, only locally (e.g. beside the archive, as `prefetch` puts
    /// them).
    #[arg(long)]
    pub offline: bool,
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
            Command::Fastq(cmd) => cmd.execute(),
            Command::Info(cmd) => cmd.execute(),
        }
    }
}

/// Write `header`, demoting its coordinate sort order (see
/// [`crate::header::demote_coordinate_sort_order`]) when the planned aligned
/// output (`None` when no aligned reads are output) is not in coordinate order
/// with respect to its `@SQ` lines. `ref_name_to_id` is the header's reference
/// id map when already built (BAM output); otherwise it is built only if needed.
fn write_header_for_output(
    db: &VDatabase,
    writer: &mut crate::output::OutputWriter,
    header: &str,
    aligned_plan: Option<&crate::aligned::AlignedPlan>,
    ref_name_to_id: Option<&HashMap<String, i32>>,
) -> Result<()> {
    let demoted = header_for_output(header, || {
        let Some(plan) = aligned_plan else {
            return Ok(true); // No aligned records: only unaligned ones, at the end.
        };
        Ok(match ref_name_to_id {
            Some(ids) => plan.is_coordinate_sorted(ids),
            None => plan.is_coordinate_sorted(&crate::header::build_ref_id_map(db, header)?),
        })
    })?;
    match demoted {
        Some(demoted) => {
            eprintln!(
                "[header] output is not coordinate-sorted; writing @HD SO:unsorted instead of \
                 SO:coordinate"
            );
            writer.write_header(&demoted)
        }
        None => writer.write_header(header),
    }
}

/// The header to write in place of `header`, or `None` to write it unchanged:
/// demoted when it claims coordinate order and `is_sorted` (only evaluated
/// then) says the output is not coordinate-sorted.
fn header_for_output(
    header: &str,
    is_sorted: impl FnOnce() -> Result<bool>,
) -> Result<Option<String>> {
    match crate::header::demote_coordinate_sort_order(header) {
        Some(demoted) if !is_sorted()? => Ok(Some(demoted)),
        _ => Ok(None),
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
    ///
    /// All accessions are written, in order, to one output (`--output-file` or
    /// stdout). A header cannot follow records, so several accessions need output
    /// without one: BAM output takes a single accession, and SAM output requires
    /// `--no-header` (FASTA/FASTQ have no header).
    pub fn execute(&self) -> Result<()> {
        self.validate()?;
        if self.offline {
            disable_remote_access().context("failed to turn off remote access")?;
        }
        if let Some(output) = &self.output_file {
            for accession in &self.accessions {
                crate::pending_file::refuse_overwriting_input(output, Path::new(accession))?;
            }
        }
        let Some((first, rest)) = self.accessions.split_first() else {
            anyhow::bail!("no accessions given");
        };
        // Open the first accession before creating the output, so that an
        // accession that cannot be opened does not truncate an existing file.
        let first_db = open_database(first)?;
        let mut writer = self.create_writer()?;
        self.process_database(&first_db, &mut writer)?;
        drop(first_db);
        for accession in rest {
            self.process_database(&open_database(accession)?, &mut writer)?;
        }
        writer.finish()
    }

    /// Reject option combinations that cannot produce valid output, before any
    /// work is done.
    fn validate(&self) -> Result<()> {
        let num_accessions = self.accessions.len();
        if num_accessions == 1 {
            return Ok(());
        }
        match self.output_mode() {
            crate::record::OutputMode::Bam => anyhow::bail!(
                "BAM output takes a single accession ({num_accessions} given): concatenated BAM \
                 streams are not a valid BAM; convert each accession separately"
            ),
            crate::record::OutputMode::Sam if !self.no_header => anyhow::bail!(
                "SAM output of several accessions ({num_accessions} given) would put each \
                 accession's header after the previous accession's records; pass --no-header or \
                 convert each accession separately"
            ),
            _ => Ok(()),
        }
    }

    /// The record format selected by the output options.
    fn output_mode(&self) -> crate::record::OutputMode {
        if self.output_format == OutputFormat::Bam {
            crate::record::OutputMode::Bam
        } else if self.fasta {
            crate::record::OutputMode::Fasta
        } else if self.fastq {
            crate::record::OutputMode::Fastq
        } else {
            crate::record::OutputMode::Sam
        }
    }

    /// Open the output (`--output-file` or stdout) for the selected format.
    fn create_writer(&self) -> Result<crate::output::OutputWriter> {
        use crate::output::{CompressionMode, OutputWriter};
        if self.output_mode() == crate::record::OutputMode::Bam {
            return Ok(match &self.output_file {
                Some(path) => OutputWriter::bam_from_path(path)?,
                None => OutputWriter::bam_stdout(),
            });
        }
        let compression = if self.gzip {
            CompressionMode::Gzip
        } else if self.bzip2 {
            CompressionMode::Bzip2
        } else {
            CompressionMode::None
        };
        Ok(match &self.output_file {
            Some(path) => OutputWriter::from_path_with_compression(path, compression)?,
            None => OutputWriter::stdout_with_compression(compression),
        })
    }

    /// Convert a single opened SRA database, writing to `writer`.
    fn process_database(
        &self,
        db: &VDatabase,
        writer: &mut crate::output::OutputWriter,
    ) -> Result<()> {
        use crate::aligned::{AlignConfig, plan_aligned_tables, process_aligned_plan};
        use crate::header::generate_header;
        use crate::progress::ProgressLogger;
        use crate::record::FormatOptions;
        use crate::unaligned::process_unaligned_reads;

        const PROGRESS_INTERVAL: u64 = 1_000_000;

        let output_mode = self.output_mode();

        // Generate the header for SAM and BAM modes.
        let header_text = if !self.no_header
            && matches!(
                output_mode,
                crate::record::OutputMode::Sam | crate::record::OutputMode::Bam
            ) {
            // Written once the aligned output order is known (below).
            Some(generate_header(
                db,
                self.header,
                self.seqid,
                &self.header_comment,
                self.header_file.as_deref(),
            )?)
        } else {
            None
        };

        // Build ref_name → ref_id map for BAM output.
        let ref_name_to_id = if output_mode == crate::record::OutputMode::Bam {
            header_text.as_deref().map(|h| crate::header::build_ref_id_map(db, h)).transpose()?
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

        // Plan the aligned reads (unless --unaligned-spots-only) before writing the
        // header, whose @HD SO must not claim coordinate order the output lacks.
        // Unaligned records (RNAME `*`) follow all aligned ones, as coordinate
        // order requires, so they never affect it.
        let aligned_plan = if self.unaligned_spots_only {
            None
        } else {
            Some(plan_aligned_tables(db, &align_config)?)
        };
        if let Some(header) = &header_text {
            write_header_for_output(
                db,
                writer,
                header,
                aligned_plan.as_ref(),
                ref_name_to_id.as_ref(),
            )?;
        }

        if let Some(plan) = &aligned_plan {
            process_aligned_plan(db, plan, writer, &align_config, PROGRESS_INTERVAL)?;
        }

        // Unaligned reads (if requested).
        if self.unaligned || self.unaligned_spots_only {
            let unaligned_progress = ProgressLogger::new(0, PROGRESS_INTERVAL);
            process_unaligned_reads(
                db,
                writer,
                &opts,
                self.unaligned_spots_only,
                &unaligned_progress,
            )?;
        }

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
            Command::CacheRefs(_) | Command::Fastq(_) | Command::Info(_) => {
                panic!("expected ToSam command")
            }
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

    #[test]
    fn test_offline_flag() {
        assert!(parse(&["--offline", "SRR123456"]).offline);
        assert!(!parse(&["SRR123456"]).offline);
    }

    fn parse_cache_refs(args: &[&str]) -> CacheRefs {
        let mut full = vec!["fg-sra", "cache-refs"];
        full.extend_from_slice(args);
        let cli = Cli::parse_from(full);
        match cli.command {
            Command::CacheRefs(cmd) => cmd,
            Command::ToSam(_) | Command::Fastq(_) | Command::Info(_) => {
                panic!("expected CacheRefs command")
            }
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
    fn test_header_for_output_demotes_only_unsorted_coordinate_claims() {
        let coordinate = "@HD\tVN:1.4\tSO:coordinate\n";
        let demoted = header_for_output(coordinate, || Ok(false)).unwrap();
        assert_eq!(demoted.as_deref(), Some("@HD\tVN:1.4\tSO:unsorted\n"));
        assert_eq!(header_for_output(coordinate, || Ok(true)).unwrap(), None);
        // Errors deciding sortedness propagate.
        assert!(header_for_output(coordinate, || anyhow::bail!("no REFERENCE")).is_err());
    }

    #[test]
    fn test_header_for_output_checks_order_only_for_coordinate_claims() {
        // Without a coordinate claim, sortedness (which may need the database) is
        // never evaluated, so a failing check cannot fail the run.
        for header in ["@HD\tVN:1.4\tSO:unsorted\n", "@HD\tVN:1.3\n", ""] {
            let result = header_for_output(header, || panic!("must not be evaluated"));
            assert_eq!(result.unwrap(), None);
        }
    }

    #[test]
    fn test_missing_accession_fails() {
        let result = Cli::try_parse_from(["fg-sra", "tosam"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_bam_with_several_accessions_is_rejected() {
        for args in [
            &["--output-format", "bam", "SRR1", "SRR2"][..],
            &["--output-format", "bam", "--output-file", "out.bam", "SRR1", "SRR2"][..],
            &["--output-format", "bam", "--gzip", "SRR1", "SRR2"][..],
            &["--output-format", "bam", "--no-header", "SRR1", "SRR2"][..],
        ] {
            let err = parse(args).validate().unwrap_err();
            assert!(err.to_string().contains("BAM output takes a single accession (2 given)"));
        }
    }

    #[test]
    fn test_sam_with_several_accessions_needs_no_header() {
        for args in [&["SRR1", "SRR2"][..], &["--output-file", "out.sam", "SRR1", "SRR2"][..]] {
            let err = parse(args).validate().unwrap_err();
            assert!(err.to_string().contains("pass --no-header"), "{err}");
        }
        parse(&["--no-header", "--output-file", "out.sam", "SRR1", "SRR2"]).validate().unwrap();
    }

    #[test]
    fn test_single_or_headerless_accessions_are_accepted() {
        parse(&["--output-format", "bam", "--output-file", "out.bam", "SRR1"]).validate().unwrap();
        parse(&["--output-file", "out.sam", "SRR1"]).validate().unwrap();
        parse(&["--fastq", "--output-file", "out.fq", "SRR1", "SRR2"]).validate().unwrap();
        parse(&["--fasta", "SRR1", "SRR2", "SRR3"]).validate().unwrap();
    }

    #[test]
    fn test_unopenable_accession_leaves_existing_output_file_untouched() {
        let dir = std::env::temp_dir().join(format!("fg-sra-cli-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let output = dir.join("existing.sam");
        std::fs::write(&output, "keep me\n").unwrap();
        let missing = dir.join("missing.sra");
        let result = parse(&["--output-file", output.to_str().unwrap(), missing.to_str().unwrap()])
            .execute();
        let contents = std::fs::read_to_string(&output).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(result.is_err());
        assert_eq!(contents, "keep me\n");
    }

    #[test]
    fn test_output_that_is_an_input_archive_is_refused() {
        let dir = std::env::temp_dir().join(format!("fg-sra-cli-input-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let archive = dir.join("in.sra");
        std::fs::write(&archive, b"archive").unwrap();
        let other = dir.join("other.sra");
        let result = parse(&[
            "--no-header",
            "--output-file",
            archive.to_str().unwrap(),
            other.to_str().unwrap(),
            archive.to_str().unwrap(),
        ])
        .execute();
        let contents = std::fs::read(&archive).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        let err = result.unwrap_err();
        assert!(err.to_string().contains("is the input archive"), "{err:#}");
        assert_eq!(contents, b"archive");
    }

    #[test]
    fn test_output_mode() {
        use crate::record::OutputMode;
        assert_eq!(parse(&["SRR1"]).output_mode(), OutputMode::Sam);
        assert_eq!(parse(&["--output-format", "bam", "SRR1"]).output_mode(), OutputMode::Bam);
        assert_eq!(parse(&["--fasta", "SRR1"]).output_mode(), OutputMode::Fasta);
        assert_eq!(parse(&["--fastq", "SRR1"]).output_mode(), OutputMode::Fastq);
    }
}
