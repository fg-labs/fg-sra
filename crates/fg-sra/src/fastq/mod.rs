//! `fg-sra fastq`: SRA archives to FASTQ, in spot order, with mates paired.
//!
//! - [`reader`]: reading spots from an archive's reads table ([`crate::archive`]).
//! - [`prefetch`]: rebuilding, ahead of conversion, aligned reads stored far from their batch's.
//! - [`spot`]: a spot's reads, and which output each spot goes to.
//! - [`defline`] and [`format`]: read names and records.
//! - [`pipeline`]: converting batches on worker threads and writing them in order.
//! - [`counts`]: tallies of what was read and written, and the metrics file.

pub mod counts;
pub mod defline;
pub mod format;
pub mod pipeline;
pub mod prefetch;
pub mod reader;
pub mod spot;

use std::collections::HashSet;
use std::ops::RangeInclusive;
use std::os::unix::fs::MetadataExt;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, anyhow, bail};
use bgzf::CompressionLevel;
use clap::{ArgGroup, Parser, ValueEnum};
use fg_sra_vdb::database::{VDatabase, VTable};
use fg_sra_vdb::manager::disable_remote_access;
use fg_sra_vdb::reference::{ReferenceList, reflist_options};

use crate::archive::{Archive, CONSENSUS_TABLE, QualitySource, default_accession};
use crate::progress::{ProgressLogger, format_count};
use crate::refstore::{ReferenceRows, ReferenceStore, preload_references, split_by_length};
use counts::FastqMetrics;
use defline::Defline;
use format::RecordFormatter;
use pipeline::{
    BATCH_SPOTS, Encoding, Output, OutputLayout, OutputTarget, PipelineConfig,
    VDB_THREAD_STACK_BYTES,
};
use prefetch::PrefetchedReads;
use reader::{ReadColumns, References, SpotSource, VdbSpotReader};
use spot::{ReadFilter, SpotRouter};

/// `--table` value that lets [`crate::archive::default_table`] choose.
const AUTO_TABLE: &str = "auto";

/// Placeholder for the technical read number in a `--technical` path template.
const TECHNICAL_PLACEHOLDER: &str = "{i}";

/// Spots between progress lines on stderr.
const PROGRESS_INTERVAL: u64 = 10_000_000;

/// At most this many threads load references, each with its own VDB manager.
const MAX_REFERENCE_LOADERS: usize = 8;

/// Convert an SRA archive to FASTQ (or FASTA), in spot order, with mates paired.
///
/// Each spot's non-empty biological reads decide where it goes: two go to the paired outputs
/// (`--r1`/`--r2`, or `--interleaved`), one goes to `--unpaired`, and spots with none or
/// more than two are dropped and counted. A spot whose output wasn't given is dropped and
/// counted too; it is an error if nothing at all is written. Filters test biological reads
/// and fail the whole spot. `--technical` writes the technical reads of the spots written,
/// one file per technical read, in step with the biological output.
///
/// Outputs ending `.gz` or `.bgz` are BGZF-compressed; any output may be `-` for stdout.
/// Output bytes are the same at any thread count. Aligned (cSRA) archives are supported: their
/// references are loaded into memory first (about a byte per reference base), and aligned
/// reads are rebuilt from their alignments.
#[derive(Debug, Parser)]
#[command(group(
    ArgGroup::new("biological_outputs")
        .required(true)
        .multiple(true)
        .args(["r1", "interleaved", "unpaired"])
))]
pub struct Fastq {
    /// SRA archive: a local `.sra`/`.sralite` file (recommended), an accession, or an https URL.
    pub input: String,

    /// Table to read within a database. `auto` reads CONSENSUS where an unaligned database has
    /// one (PacBio and Oxford Nanopore native loads, whose SEQUENCE holds each molecule's
    /// subreads or strands), as fasterq-dump does, and SEQUENCE otherwise.
    #[arg(long, default_value = AUTO_TABLE, help_heading = "Input", value_name = "NAME")]
    pub table: String,

    /// First spot to convert (1-based).
    #[arg(long, help_heading = "Input", value_name = "N")]
    pub min_spot_id: Option<i64>,

    /// Last spot to convert (inclusive).
    #[arg(long, help_heading = "Input", value_name = "N")]
    pub max_spot_id: Option<i64>,

    /// First mate of each pair.
    #[arg(
        short = '1',
        long = "r1",
        requires = "r2",
        help_heading = "Outputs",
        value_name = "PATH"
    )]
    pub r1: Option<PathBuf>,

    /// Second mate of each pair.
    #[arg(
        short = '2',
        long = "r2",
        requires = "r1",
        help_heading = "Outputs",
        value_name = "PATH"
    )]
    pub r2: Option<PathBuf>,

    /// Pairs, with mates adjacent (instead of --r1/--r2).
    #[arg(short = 'p', long, conflicts_with = "r1", help_heading = "Outputs", value_name = "PATH")]
    pub interleaved: Option<PathBuf>,

    /// Spots with one biological read: single-end runs, and orphans of empty mates.
    #[arg(short = 'u', long, help_heading = "Outputs", value_name = "PATH")]
    pub unpaired: Option<PathBuf>,

    /// Technical reads of the spots written. Either one path containing `{i}`, replaced by
    /// the technical read number, or the option repeated with one path per technical read.
    /// Needs exactly one kind of biological output: pairs, or --unpaired.
    #[arg(long, help_heading = "Outputs", value_name = "PATH")]
    pub technical: Vec<String>,

    /// Minimum biological read length; a spot with a shorter biological read is dropped.
    #[arg(long, help_heading = "Filters", value_name = "N")]
    pub min_read_len: Option<u32>,

    /// `READ_FILTER` values to keep (comma-separated); a spot with a biological read of any
    /// other value is dropped. Default: all. Runs loaded from BAM mark duplicates `criteria`, so
    /// `pass` also drops them.
    #[arg(long, value_enum, value_delimiter = ',', help_heading = "Filters", value_name = "VALUE")]
    pub read_filter: Vec<ReadFilter>,

    /// Read-name template: `$ac` accession, `$si` spot id, `$sn` original name (else the
    /// spot id), `$sg` spot group, `$ri` read number within its type, `$rl` read length.
    #[arg(long, default_value = "$ac.$si", help_heading = "Format", value_name = "TEMPLATE")]
    pub defline: Defline,

    /// Value of `$ac` [default: the input's name, without .sra/.sralite].
    #[arg(long, help_heading = "Format", value_name = "ACC")]
    pub accession: Option<String>,

    /// Write FASTA instead of FASTQ; qualities aren't read.
    #[arg(long, help_heading = "Format")]
    pub fasta: bool,

    /// Output compression: `auto` writes BGZF to paths ending `.gz`/`.bgz` and plain text
    /// otherwise (including stdout); `bgzf` and `none` apply to every output.
    #[arg(long, value_enum, default_value_t = OutputCompression::Auto, help_heading = "Format")]
    pub output_compression: OutputCompression,

    /// BGZF compression level (1-12).
    #[arg(
        short = 'c',
        long,
        default_value_t = 1,
        value_parser = clap::value_parser!(u8).range(1..=12),
        help_heading = "Format"
    )]
    pub compression_level: u8,

    /// Threads for reading, formatting and compressing; writing uses one more per output.
    #[arg(
        short = 't',
        long,
        default_value_t = 4,
        value_parser = clap::value_parser!(u16).range(1..),
        help_heading = "Run"
    )]
    pub threads: u16,

    /// Metrics TSV; a summary always goes to stderr.
    #[arg(short = 'm', long, help_heading = "Run", value_name = "PATH")]
    pub metrics: Option<PathBuf>,

    /// Never use the network: find the archive, and an aligned archive's references, only
    /// locally (e.g. beside the archive, as `prefetch` puts them).
    #[arg(long, help_heading = "Run")]
    pub offline: bool,
}

impl Fastq {
    /// Run the conversion.
    pub fn execute(&self) -> Result<()> {
        self.check_outputs()?;
        if self.offline {
            disable_remote_access().context("failed to turn off remote access")?;
        }
        let requested_table = (self.table != AUTO_TABLE).then_some(self.table.as_str());
        let archive = Archive::inspect(&self.input, requested_table)?;
        if requested_table.is_none() && archive.table.as_deref() == Some(CONSENSUS_TABLE) {
            eprintln!(
                "[fastq] {}: reading consensus reads (CONSENSUS); --table SEQUENCE reads each \
                 molecule's subreads or strands instead",
                self.input
            );
        }
        // Without alignments, READ is safe to read even when it is virtual, as it is for
        // colour-space runs decoded from CSREAD. With alignments, a virtual READ reaches into
        // the reference cache, which isn't thread-safe, so those reads are restored instead.
        let restore_aligned_reads = match (archive.aligned, archive.stores_cmp_read) {
            (false, _) => false,
            (true, true) => true,
            (true, false) if archive.stores_bases => false,
            (true, false) => bail!(
                "{}: the reads table stores neither READ nor CMP_READ; this layout is not \
                 supported",
                self.input
            ),
        };
        let accession = self.accession.clone().unwrap_or_else(|| default_accession(&self.input));
        self.warn_about_qualities(&archive, &accession)?;
        let spots = self.spot_range(&archive)?;
        let columns = ReadColumns {
            qualities: !self.fasta,
            names: archive.has_names && self.defline.uses_spot_name(),
            spot_groups: archive.has_spot_groups && self.defline.uses_spot_group(),
        };
        // One open table serves every worker, each with its own cursor; it stays open until
        // the pipeline's threads have finished.
        let table = archive.open_table()?;
        let technical_paths = self.technical_paths(&table.table, *spots.start(), columns)?;
        let (outputs, layout) = self.outputs(&technical_paths);
        refuse_overwriting_input(&outputs, Path::new(&archive.location))?;
        let removable = removable_outputs(&outputs);

        let config = PipelineConfig {
            router: self.router(layout.technical.len()),
            formatter: RecordFormatter::new(self.defline.clone(), &accession, self.fasta),
            layout,
            with_qualities: !self.fasta,
            compression_level: CompressionLevel::new(self.compression_level)?,
            batch_spots: BATCH_SPOTS,
        };
        let aligned_parts = match table.database() {
            Some(database) if restore_aligned_reads => {
                let alignments = database
                    .open_table_read("PRIMARY_ALIGNMENT")
                    .context("failed to open the PRIMARY_ALIGNMENT table")?;
                let (store, rows) =
                    load_references(&archive, database, &accession, usize::from(self.threads))?;
                Some((alignments, store, rows))
            }
            _ => None,
        };
        let prefetched = match &aligned_parts {
            Some((alignments, store, rows)) => Some(prefetch_far_reads(
                &table.table,
                alignments,
                References { store, rows },
                spots.clone(),
                usize::from(self.threads),
                &accession,
            )?),
            None => None,
        };
        let sources = (0..self.threads)
            .map(|_| match &aligned_parts {
                Some((alignments, store, rows)) => VdbSpotReader::new_aligned(
                    &table.table,
                    alignments,
                    columns,
                    References { store, rows },
                    prefetched.as_ref(),
                ),
                None => VdbSpotReader::new(&table.table, columns),
            })
            .collect::<Result<Vec<_>>>()?;
        let progress = ProgressLogger::new(0, PROGRESS_INTERVAL);
        let result = pipeline::run(sources, spots.clone(), &outputs, &config, &progress)
            .and_then(|summary| self.finish(&archive, &accession, &spots, &summary));
        if result.is_err() {
            remove_outputs(&removable);
        }
        result
    }

    /// Check the outputs requested make sense together, before opening the archive.
    fn check_outputs(&self) -> Result<()> {
        let pairs = self.r1.is_some() || self.interleaved.is_some();
        if !self.technical.is_empty() && pairs == self.unpaired.is_some() {
            bail!(
                "--technical needs exactly one kind of biological output, so its files line up \
                 record for record: pairs (--r1/--r2 or --interleaved), or --unpaired"
            );
        }
        let templates = self.technical.iter().filter(|p| p.contains(TECHNICAL_PLACEHOLDER));
        if templates.count() > 0 && self.technical.len() > 1 {
            bail!(
                "--technical takes one path containing {TECHNICAL_PLACEHOLDER}, or one path per \
                 technical read, not both"
            );
        }
        if let (Some(min), Some(max)) = (self.min_spot_id, self.max_spot_id)
            && min > max
        {
            bail!("--min-spot-id {min} is greater than --max-spot-id {max}");
        }

        let paths: Vec<&Path> = [&self.r1, &self.r2, &self.interleaved, &self.unpaired]
            .into_iter()
            .flatten()
            .map(PathBuf::as_path)
            .chain(self.technical.iter().map(Path::new))
            .collect();
        if paths.iter().filter(|&&p| p == Path::new("-")).count() > 1 {
            bail!("at most one output may be - (stdout)");
        }
        let mut seen = HashSet::new();
        if let Some(duplicate) = paths.iter().find(|&&p| p != Path::new("-") && !seen.insert(p)) {
            bail!("{} is given as more than one output", duplicate.display());
        }
        Ok(())
    }

    /// Warn about synthesised or constant qualities; fail if FASTQ needs qualities there are
    /// none of.
    fn warn_about_qualities(&self, archive: &Archive, accession: &str) -> Result<()> {
        if self.fasta {
            return Ok(());
        }
        if !archive.has_qualities {
            bail!("{} has no qualities to write as FASTQ; use --fasta", self.input);
        }
        match archive.quality_source {
            QualitySource::Stored => {}
            QualitySource::Lite => eprintln!(
                "[fastq] {accession}: warning: SRA Lite archive; qualities are synthesised \
                 (Q30, or Q3 for reads that fail the filter)"
            ),
            QualitySource::None => eprintln!(
                "[fastq] {accession}: warning: the archive stores no qualities; those written \
                 are synthesised"
            ),
        }
        if let Some(phred) = archive.constant_quality {
            eprintln!("[fastq] {accession}: warning: every base has quality {phred}");
        }
        Ok(())
    }

    /// The spot ids to convert: the archive's, narrowed by --min-spot-id/--max-spot-id.
    fn spot_range(&self, archive: &Archive) -> Result<RangeInclusive<i64>> {
        let (first, last) = (archive.first_spot, archive.last_spot());
        let start = self.min_spot_id.unwrap_or(first).max(first);
        let end = self.max_spot_id.unwrap_or(last).min(last);
        if archive.spot_count == 0 {
            bail!("{} has no spots", self.input);
        }
        if start > end {
            bail!("the requested spots are outside {}'s spots {first}-{last}", self.input);
        }
        Ok(start..=end)
    }

    /// The technical output paths: the template expanded for each technical read, or the
    /// paths as given. Their number must match the technical reads of the first spot
    /// converted.
    fn technical_paths(
        &self,
        table: &VTable,
        first_spot: i64,
        columns: ReadColumns,
    ) -> Result<Vec<PathBuf>> {
        if self.technical.is_empty() {
            return Ok(Vec::new());
        }
        let mut reader = VdbSpotReader::new(table, columns)?;
        let spot = reader.read(first_spot)?;
        spot.validate(columns.qualities)?;
        let technical_reads = spot.technical_reads().count();
        if technical_reads == 0 {
            bail!("spot {first_spot} of {} has no technical reads to write", self.input);
        }
        match self.technical.as_slice() {
            [template] if template.contains(TECHNICAL_PLACEHOLDER) => Ok((1..=technical_reads)
                .map(|number| {
                    PathBuf::from(template.replace(TECHNICAL_PLACEHOLDER, &number.to_string()))
                })
                .collect()),
            paths if paths.len() == technical_reads => {
                Ok(paths.iter().map(PathBuf::from).collect())
            }
            paths => bail!(
                "--technical gives {} path(s), but spot {first_spot} of {} has {technical_reads} \
                 technical read(s)",
                paths.len(),
                self.input
            ),
        }
    }

    /// The outputs, in a fixed order, and which record kinds go to each.
    fn outputs(&self, technical_paths: &[PathBuf]) -> (Vec<Output>, OutputLayout) {
        let mut outputs = Vec::new();
        let mut add = |path: &Path| {
            outputs.push(Output {
                target: OutputTarget::from_arg(path),
                encoding: self.output_compression.encoding(path),
            });
            outputs.len() - 1
        };
        let layout = OutputLayout {
            first_mates: self.r1.as_deref().map(&mut add),
            second_mates: self.r2.as_deref().map(&mut add),
            interleaved: self.interleaved.as_deref().map(&mut add),
            unpaired: self.unpaired.as_deref().map(&mut add),
            technical: technical_paths.iter().map(|path| add(path)).collect(),
        };
        (outputs, layout)
    }

    /// The router for the filters and outputs requested, with `technical_outputs` technical
    /// read outputs.
    fn router(&self, technical_outputs: usize) -> SpotRouter {
        let kept_filters = if self.read_filter.is_empty() {
            [true; 4]
        } else {
            ReadFilter::ALL.map(|filter| self.read_filter.contains(&filter))
        };
        SpotRouter {
            min_read_len: self.min_read_len.unwrap_or(0),
            kept_filters,
            write_pairs: self.r1.is_some() || self.interleaved.is_some(),
            write_unpaired: self.unpaired.is_some(),
            technical_reads: (technical_outputs > 0).then_some(technical_outputs),
        }
    }

    /// Check the finished conversion, report it, and write the metrics.
    ///
    /// Fails if the run was complete but its totals don't match the archive's stored ones,
    /// or if nothing was written.
    fn finish(
        &self,
        archive: &Archive,
        accession: &str,
        spots: &RangeInclusive<i64>,
        summary: &pipeline::PipelineSummary,
    ) -> Result<()> {
        let counts = &summary.counts;
        let whole_run = !summary.stopped_early
            && *spots.start() == archive.first_spot
            && *spots.end() == archive.last_spot();
        let metrics = FastqMetrics::new(
            accession,
            archive,
            counts,
            summary.technical_reads_written,
            whole_run,
        );
        metrics.report();
        if let Some(path) = &self.metrics {
            fgoxide::io::DelimFile::default()
                .write_tsv(path, std::iter::once(&metrics))
                .with_context(|| format!("failed to write metrics to {}", path.display()))?;
        }

        if summary.stopped_early {
            eprintln!(
                "[fastq] {accession}: an output's reader went away, so conversion stopped early"
            );
            return Ok(());
        }
        if let Some(mismatch) = metrics.integrity_mismatch() {
            bail!("{accession}: {mismatch}; the conversion is incomplete or the archive damaged");
        }
        if counts.spots_written() == 0 {
            bail!("{accession}: no spots were written. {}", counts.nothing_written_hint());
        }
        Ok(())
    }
}

/// How outputs are compressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OutputCompression {
    /// BGZF for paths ending `.gz` or `.bgz`, plain text otherwise.
    Auto,
    /// BGZF for every output.
    Bgzf,
    /// Plain text for every output.
    None,
}

impl OutputCompression {
    /// How an output at `path` is encoded.
    fn encoding(self, path: &Path) -> Encoding {
        let gzip_extension = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("bgz"));
        match self {
            Self::Auto if gzip_extension => Encoding::Bgzf,
            Self::Bgzf => Encoding::Bgzf,
            Self::Auto | Self::None => Encoding::Plain,
        }
    }
}

/// Load every reference of `database`, `archive`'s, into memory on up to `threads` threads,
/// and which `REFERENCE` rows each covers, before any worker reads alignments.
///
/// A reference list and the references it opens aren't safe to share between threads, so each
/// loader opens the archive again, with its own manager, and loads its share of the references,
/// balanced by length, through its own list.
fn load_references(
    archive: &Archive,
    database: &VDatabase,
    accession: &str,
    threads: usize,
) -> Result<(ReferenceStore, ReferenceRows)> {
    let started = Instant::now();
    let reflist = ReferenceList::make_database(database, reflist_options::READ_4NA, 0)
        .context("failed to list the archive's references")?;
    let rows = ReferenceRows::new(&reflist, reference_row_length(database)?)?;
    let lengths = (0..reflist.count()?)
        .map(|idx| Ok((idx, u64::from(reflist.get(idx)?.seq_length()?))))
        .collect::<Result<Vec<_>>>()?;
    let groups = split_by_length(&lengths, threads.min(MAX_REFERENCE_LOADERS));
    let store = if groups.len() <= 1 {
        preload_references(&reflist, &groups.concat())
    } else {
        let opened = groups.iter().map(|_| archive.open_table()).collect::<Result<Vec<_>>>()?;
        std::thread::scope(|scope| {
            let loaders = opened
                .into_iter()
                .zip(&groups)
                .map(|(open, indices)| {
                    std::thread::Builder::new().stack_size(VDB_THREAD_STACK_BYTES).spawn_scoped(
                        scope,
                        move || {
                            let load = || -> Result<ReferenceStore> {
                                let database =
                                    open.database().context("the archive has no database")?;
                                let reflist = ReferenceList::make_database(
                                    database,
                                    reflist_options::READ_4NA,
                                    0,
                                )?;
                                preload_references(&reflist, indices)
                            };
                            (load(), open)
                        },
                    )
                })
                .collect::<std::io::Result<Vec<_>>>()
                .context("failed to start a reference loader")?;
            // Managers are released only once every loader has finished: in ncbi-vdb before
            // VDB-6280 (the vendored copy is patched), releasing one frees a process-wide string
            // that the other loaders' resolvers read.
            let (stores, _managers): (Vec<_>, Vec<_>) = loaders
                .into_iter()
                .map(|loader| match loader.join() {
                    Ok((store, open)) => (store, Some(open)),
                    Err(_) => (Err(anyhow!("a reference loader panicked")), None),
                })
                .unzip();
            stores.into_iter().collect::<Result<Vec<_>>>()
        })
        .map(ReferenceStore::merge)
    };
    let store = store.context(
        "failed to load the archive's references; external ones must be available locally \
         (e.g. beside the archive, as `prefetch` puts them) or, without --offline, over the \
         network",
    )?;
    eprintln!(
        "[fastq] {accession}: loaded {} references ({} bases) in {:.1}s",
        store.num_references(),
        format_count(store.total_bases()),
        started.elapsed().as_secs_f64()
    );
    Ok((store, rows))
}

/// Rebuild, ahead of conversion, the aligned reads of `spots` whose alignments lie far from the
/// rest of their batch's (see [`prefetch`]).
fn prefetch_far_reads(
    table: &VTable,
    alignments: &VTable,
    references: References<'_>,
    spots: RangeInclusive<i64>,
    threads: usize,
    accession: &str,
) -> Result<PrefetchedReads> {
    let started = Instant::now();
    let reads =
        PrefetchedReads::far_alignments(table, alignments, references, spots, BATCH_SPOTS, threads)
            .context("failed to prefetch aligned reads")?;
    eprintln!(
        "[fastq] {accession}: prefetched {} aligned reads far from their batches' others ({} \
         bases) in {:.1}s",
        format_count(reads.len() as u64),
        format_count(reads.total_bases() as u64),
        started.elapsed().as_secs_f64()
    );
    Ok(reads)
}

/// Bases per row of the `REFERENCE` table (`MAX_SEQ_LEN`, read from its first row, as the
/// schema's `REF_POS` reads it), or 0 when the table has no such column.
fn reference_row_length(database: &VDatabase) -> Result<u32> {
    let table = database.open_table_read("REFERENCE").context("failed to open REFERENCE")?;
    let cursor = table.create_cursor_read()?;
    let column = match cursor.add_column("(U32)MAX_SEQ_LEN") {
        Ok(column) => column,
        Err(e) if e.is_not_found() => return Ok(0),
        Err(e) => return Err(e).context("failed to read REFERENCE.MAX_SEQ_LEN"),
    };
    cursor.open()?;
    cursor.read_u32(1, column).context("failed to read REFERENCE.MAX_SEQ_LEN")
}

/// Fail if an output is the input archive itself, however it is named (another spelling, a
/// symbolic or a hard link): writing it would truncate the archive, and a failed run removes it.
fn refuse_overwriting_input(outputs: &[Output], input: &Path) -> Result<()> {
    let Ok(input_meta) = std::fs::metadata(input) else { return Ok(()) };
    for output in outputs {
        if let OutputTarget::Path(path) = &output.target
            && let Ok(meta) = std::fs::metadata(path)
            && (meta.dev(), meta.ino()) == (input_meta.dev(), input_meta.ino())
        {
            bail!("{} is the input archive, which would be overwritten", path.display());
        }
    }
    Ok(())
}

/// The outputs that are regular files this run creates or truncates, which a failed run
/// removes; FIFOs, devices (e.g. `/dev/null`) and stdout are left alone.
fn removable_outputs(outputs: &[Output]) -> Vec<PathBuf> {
    outputs
        .iter()
        .filter_map(|output| match &output.target {
            OutputTarget::Path(path) => Some(path),
            OutputTarget::Stdout => None,
        })
        .filter(|path| std::fs::metadata(path).map_or(true, |meta| meta.is_file()))
        .cloned()
        .collect()
}

/// Remove outputs of a failed run, so a partial file can't be mistaken for a complete one.
fn remove_outputs(paths: &[PathBuf]) {
    for path in paths {
        // Best effort: the file may never have been created.
        let _ = std::fs::remove_file(path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Fastq, clap::Error> {
        Fastq::try_parse_from(std::iter::once("fastq").chain(args.iter().copied()))
    }

    #[test]
    fn pairs_with_defaults_parse() {
        let cmd = parse(&["SRR1.sra", "-1", "a.fq.gz", "-2", "b.fq.gz"]).unwrap();
        assert_eq!(cmd.r1.as_deref(), Some(Path::new("a.fq.gz")));
        assert_eq!(cmd.threads, 4);
        assert_eq!(cmd.compression_level, 1);
        assert_eq!(cmd.defline, Defline::parse("$ac.$si").unwrap());
        assert_eq!(cmd.table, "auto");
        assert!(cmd.read_filter.is_empty());
    }

    #[test]
    fn some_biological_output_is_required() {
        assert!(parse(&["SRR1.sra"]).is_err());
    }

    #[test]
    fn r1_requires_r2() {
        assert!(parse(&["SRR1.sra", "-1", "a.fq"]).is_err());
    }

    #[test]
    fn r2_requires_r1() {
        assert!(parse(&["SRR1.sra", "-2", "b.fq"]).is_err());
    }

    #[test]
    fn interleaved_conflicts_with_r1() {
        assert!(parse(&["SRR1.sra", "-p", "-", "-1", "a.fq", "-2", "b.fq"]).is_err());
    }

    #[test]
    fn unpaired_alone_is_enough() {
        assert!(parse(&["SRR1.sra", "-u", "single.fq"]).is_ok());
    }

    #[test]
    fn technical_can_be_repeated() {
        let cmd = parse(&["SRR1.sra", "-u", "r.fq", "--technical", "a.fq", "--technical", "b.fq"])
            .unwrap();
        assert_eq!(cmd.technical, vec!["a.fq", "b.fq"]);
    }

    #[test]
    fn read_filter_takes_a_comma_separated_list() {
        let cmd = parse(&["SRR1.sra", "-u", "r.fq", "--read-filter", "pass,criteria"]).unwrap();
        assert_eq!(cmd.read_filter, vec![ReadFilter::Pass, ReadFilter::Criteria]);
    }

    #[test]
    fn invalid_defline_is_a_parse_error() {
        assert!(parse(&["SRR1.sra", "-u", "r.fq", "--defline", "$ac.$zz"]).is_err());
    }

    #[test]
    fn compression_level_above_12_is_a_parse_error() {
        assert!(parse(&["SRR1.sra", "-u", "r.fq", "-c", "13"]).is_err());
    }

    #[test]
    fn zero_threads_is_a_parse_error() {
        assert!(parse(&["SRR1.sra", "-u", "r.fq", "-t", "0"]).is_err());
    }

    #[test]
    fn technical_with_both_pairs_and_unpaired_is_refused() {
        let cmd =
            parse(&["S.sra", "-1", "a", "-2", "b", "-u", "c", "--technical", "t{i}"]).unwrap();
        assert!(cmd.check_outputs().is_err());
    }

    #[test]
    fn technical_with_only_pairs_is_accepted() {
        let cmd = parse(&["S.sra", "-1", "a", "-2", "b", "--technical", "t{i}"]).unwrap();
        assert!(cmd.check_outputs().is_ok());
    }

    #[test]
    fn technical_template_among_several_paths_is_refused() {
        let cmd = parse(&["S.sra", "-u", "c", "--technical", "t{i}", "--technical", "x"]).unwrap();
        assert!(cmd.check_outputs().is_err());
    }

    #[test]
    fn two_outputs_to_stdout_are_refused() {
        let cmd = parse(&["S.sra", "-p", "-", "-u", "-"]).unwrap();
        assert!(cmd.check_outputs().is_err());
    }

    #[test]
    fn the_same_path_twice_is_refused() {
        let cmd = parse(&["S.sra", "-1", "x.fq", "-2", "x.fq"]).unwrap();
        assert!(cmd.check_outputs().is_err());
    }

    #[test]
    fn min_spot_above_max_spot_is_refused() {
        let cmd =
            parse(&["S.sra", "-u", "c", "--min-spot-id", "10", "--max-spot-id", "5"]).unwrap();
        assert!(cmd.check_outputs().is_err());
    }

    /// A fresh directory holding `archive.sra`.
    fn dir_with_archive(name: &str) -> (PathBuf, PathBuf) {
        let dir = std::env::temp_dir().join(format!("fg_sra_fastq_{}_{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let archive = dir.join("archive.sra");
        std::fs::write(&archive, b"archive").unwrap();
        (dir, archive)
    }

    fn file_output(path: &Path) -> Output {
        Output { target: OutputTarget::from_arg(path), encoding: Encoding::Plain }
    }

    #[test]
    fn output_spelling_the_input_differently_is_refused() {
        let (dir, archive) = dir_with_archive("spelling");
        let outputs =
            [file_output(&dir.join("r1.fq")), file_output(&dir.join(".").join("archive.sra"))];
        assert!(refuse_overwriting_input(&outputs, &archive).is_err());
        assert_eq!(std::fs::read(&archive).unwrap(), b"archive");
    }

    #[test]
    fn output_hard_linked_to_the_input_is_refused() {
        let (dir, archive) = dir_with_archive("hard_link");
        let link = dir.join("link.fq");
        std::fs::hard_link(&archive, &link).unwrap();
        assert!(refuse_overwriting_input(&[file_output(&link)], &archive).is_err());
    }

    #[test]
    fn outputs_other_than_the_input_are_accepted() {
        let (dir, archive) = dir_with_archive("others");
        std::fs::write(dir.join("old.fq"), b"old").unwrap();
        let outputs = [
            file_output(&dir.join("old.fq")),
            file_output(&dir.join("new.fq")),
            file_output(Path::new("-")),
        ];
        assert!(refuse_overwriting_input(&outputs, &archive).is_ok());
    }

    #[test]
    fn an_input_that_is_not_a_local_file_has_no_outputs_to_refuse() {
        let (dir, _) = dir_with_archive("remote");
        let outputs = [file_output(&dir.join("SRR1"))];
        assert!(refuse_overwriting_input(&outputs, Path::new("SRR1")).is_ok());
    }

    #[test]
    fn gz_and_bgz_paths_are_bgzf_under_auto() {
        assert_eq!(OutputCompression::Auto.encoding(Path::new("r1.fq.gz")), Encoding::Bgzf);
        assert_eq!(OutputCompression::Auto.encoding(Path::new("r1.fq.BGZ")), Encoding::Bgzf);
    }

    #[test]
    fn other_paths_and_stdout_are_plain_under_auto() {
        assert_eq!(OutputCompression::Auto.encoding(Path::new("r1.fq")), Encoding::Plain);
        assert_eq!(OutputCompression::Auto.encoding(Path::new("-")), Encoding::Plain);
    }

    #[test]
    fn explicit_compression_applies_whatever_the_extension() {
        assert_eq!(OutputCompression::Bgzf.encoding(Path::new("-")), Encoding::Bgzf);
        assert_eq!(OutputCompression::None.encoding(Path::new("r1.fq.gz")), Encoding::Plain);
    }

    #[test]
    fn router_keeps_every_filter_by_default() {
        let cmd = parse(&["S.sra", "-u", "c"]).unwrap();
        assert_eq!(cmd.router(0).kept_filters, [true; 4]);
    }

    #[test]
    fn router_keeps_only_the_filters_listed() {
        let cmd = parse(&["S.sra", "-u", "c", "--read-filter", "pass"]).unwrap();
        assert_eq!(cmd.router(0).kept_filters, [true, false, false, false]);
    }
}
