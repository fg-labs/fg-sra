//! What an input SRA archive is, and opening its reads table.

use std::path::Path;

use anyhow::{Context, Result, bail};
use fg_sra_vdb::database::{KMetadata, VDatabase, VTable};
use fg_sra_vdb::manager::{PathType, VdbManager};

/// Name of the reads table in an SRA database.
const SEQUENCE_TABLE: &str = "SEQUENCE";

/// Name of the table of consensus reads in PacBio and Oxford Nanopore native databases.
pub const CONSENSUS_TABLE: &str = "CONSENSUS";

/// Name of the alignments table of an aligned (cSRA) database.
const ALIGNMENT_TABLE: &str = "PRIMARY_ALIGNMENT";

/// An input archive, inspected once before it is read.
// Each flag is an independent property of the archive, not an encoded state.
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone)]
pub struct Archive {
    /// What VDB is given to open: an absolute path for a local file, else the input as given.
    pub location: String,
    /// Reads table within a database; `None` for a flat table.
    pub table: Option<String>,
    /// Whether the archive has alignments (a `PRIMARY_ALIGNMENT` table).
    pub aligned: bool,
    /// Whether the reads table stores bases directly (a physical `READ` column).
    pub stores_bases: bool,
    /// Whether the reads table stores only its unaligned reads' bases (a physical
    /// `CMP_READ` column), as aligned archives do, the rest being rebuilt from alignments.
    pub stores_cmp_read: bool,
    /// Sequencing platform of the first spot, e.g. `ILLUMINA`.
    pub platform: Option<&'static str>,
    /// First spot (row) id.
    pub first_spot: i64,
    /// Number of spots (rows).
    pub spot_count: u64,
    pub quality_source: QualitySource,
    /// Whether `QUALITY` can be read at all.
    pub has_qualities: bool,
    /// The single phred value of every base, when the quality histogram has only one bin.
    pub constant_quality: Option<u8>,
    /// Whether the archive kept original spot names (rather than none, or serial numbers).
    pub has_names: bool,
    pub has_spot_groups: bool,
    pub totals: StoredTotals,
}

impl Archive {
    /// Inspect `input` (a path, accession or URL), choosing `table` within a database, or by
    /// default the one [`default_table`] picks.
    ///
    /// Fails if the input isn't an SRA archive, or if `table` is given for a flat table.
    pub fn inspect(input: &str, table: Option<&str>) -> Result<Self> {
        let location = vdb_location(input)?;
        let manager = new_manager()?;
        let (database, table_name) = match manager.path_type(&location)? {
            PathType::Table => {
                if let Some(table) = table {
                    bail!("--table {table} applies only to databases, and {input} is a flat table");
                }
                (None, None)
            }
            PathType::Database => {
                let database = manager
                    .open_db_read(&location)
                    .with_context(|| format!("failed to open database {input}"))?;
                let name = match table {
                    Some(table) => table.to_owned(),
                    None => default_table(&database.list_tables()?).to_owned(),
                };
                (Some(database), Some(name))
            }
            PathType::NotFound => bail!("no SRA archive found at {input}"),
            PathType::Other(_) => bail!("{input} is not an SRA archive (a VDB table or database)"),
        };
        let opened = open_table(manager, database, &location, table_name.as_deref(), input)?;
        let table = &opened.table;

        let physical = table.list_physical_columns()?;
        let readable = table.list_readable_columns()?;
        let has = |columns: &[String], name: &str| columns.iter().any(|c| c == name);
        let tables = match &opened.database {
            Some(database) => database.list_tables()?,
            None => Vec::new(),
        };

        // The Lite marker sits in the top-level metadata: the database's, or the flat table's.
        let top_metadata = match &opened.database {
            Some(database) => database.open_metadata_read()?,
            None => table.open_metadata_read()?,
        };
        let table_metadata = table.open_metadata_read()?;
        let has_qualities = has(&readable, "QUALITY");
        let quality_source =
            if has_qualities && (has(&physical, "QUALITY") || has(&physical, "ORIGINAL_QUALITY")) {
                QualitySource::Stored
            } else if is_sra_lite(&top_metadata) {
                QualitySource::Lite
            } else {
                QualitySource::None
            };

        let (first_spot, spot_count) = id_range(table)?;
        let has_names =
            spot_count > 0 && has(&readable, "NAME") && names_are_original(table, first_spot)?;
        let platform = if spot_count > 0 && has(&readable, "PLATFORM") {
            platform(table, first_spot)?
        } else {
            None
        };
        Ok(Self {
            location,
            table: table_name,
            aligned: tables.iter().any(|t| t == ALIGNMENT_TABLE),
            stores_bases: has(&physical, "READ"),
            stores_cmp_read: has(&physical, "CMP_READ"),
            platform,
            first_spot,
            spot_count,
            quality_source,
            has_qualities,
            constant_quality: constant_quality(&table_metadata),
            has_names,
            has_spot_groups: has(&readable, "SPOT_GROUP"),
            totals: StoredTotals::read(&table_metadata),
        })
    }

    /// Last spot (row) id.
    pub fn last_spot(&self) -> i64 {
        self.first_spot + self.spot_count as i64 - 1
    }

    /// Open the reads table with a fresh manager.
    pub fn open_table(&self) -> Result<OpenTable> {
        let manager = new_manager()?;
        let database = match &self.table {
            Some(_) => Some(manager.open_db_read(&self.location)?),
            None => None,
        };
        open_table(manager, database, &self.location, self.table.as_deref(), &self.location)
    }
}

/// A reads table opened with its own manager. Cursors on it may be used on other threads
/// while it stays open.
///
/// Fields drop in declaration order, releasing the table before what it was opened from.
pub struct OpenTable {
    pub table: VTable,
    database: Option<VDatabase>,
    _manager: VdbManager,
}

impl OpenTable {
    /// The database the table is in; `None` for a flat table.
    pub fn database(&self) -> Option<&VDatabase> {
        self.database.as_ref()
    }
}

/// How an archive's qualities came to be, classified as `sra-info` does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum QualitySource {
    /// Qualities as submitted.
    Stored,
    /// Removed to make an SRA Lite object; readers see synthesised Q30, or Q3 for reads that
    /// fail the filter.
    Lite,
    /// Never stored (or removed without SRA Lite's marker); any qualities read are synthesised.
    None,
}

/// Totals the archive stores in `STATS/TABLE`, when present (old tables lack them).
#[derive(Debug, Clone, Copy, Default)]
pub struct StoredTotals {
    pub spots: Option<u64>,
    pub bases: Option<u64>,
    pub biological_bases: Option<u64>,
}

impl StoredTotals {
    /// The totals in `metadata`, the reads table's.
    fn read(metadata: &KMetadata) -> Self {
        let read = |name: &str| {
            metadata
                .open_node_read(&format!("STATS/TABLE/{name}"))
                .and_then(|node| node.read_u64())
                .ok()
        };
        Self {
            spots: read("SPOT_COUNT"),
            bases: read("BASE_COUNT"),
            biological_bases: read("BIO_BASE_COUNT"),
        }
    }
}

/// The accession `$ac` defaults to: the input's last path component, without `.sra`,
/// `.sralite` or `.lite.sra` (so `/data/SRR1.sra` and `SRR1` both give `SRR1`), and without
/// a URL's query. A run accession also loses the version that objects fetched from NCBI or
/// the cloud carry (`SRR1.1`, `SRR1.sralite.1` and `SRR1.lite.1` all give `SRR1`).
pub fn default_accession(input: &str) -> String {
    let path =
        if input.contains("://") { input.split(['?', '#']).next().unwrap_or(input) } else { input };
    let name = path.trim_end_matches('/').rsplit('/').next().unwrap_or(path);
    let name = without_archive_extension(name);
    if let Some((unversioned, version)) = name.rsplit_once('.')
        && !version.is_empty()
        && version.bytes().all(|b| b.is_ascii_digit())
    {
        let unversioned = without_archive_extension(unversioned);
        if is_run_accession(unversioned) {
            return unversioned.to_owned();
        }
    }
    name.to_owned()
}

/// `name` without a trailing `.sra`, `.sralite` or `.lite.sra`, or a bare `.lite`.
fn without_archive_extension(name: &str) -> &str {
    let name = name.strip_suffix(".sralite").or_else(|| name.strip_suffix(".sra")).unwrap_or(name);
    name.strip_suffix(".lite").unwrap_or(name)
}

/// Whether `name` is an SRA run accession: `SRR`, `ERR` or `DRR` and digits.
fn is_run_accession(name: &str) -> bool {
    let digits = name
        .strip_prefix("SRR")
        .or_else(|| name.strip_prefix("ERR"))
        .or_else(|| name.strip_prefix("DRR"));
    digits.is_some_and(|digits| !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()))
}

/// The reads table of a database with `tables`: `CONSENSUS` in an unaligned database that has
/// one, as fasterq-dump reads (PacBio and Oxford Nanopore native loads, whose `SEQUENCE` holds
/// each molecule's subreads or strands as separate reads), and `SEQUENCE` otherwise.
pub fn default_table(tables: &[String]) -> &'static str {
    let has = |name: &str| tables.iter().any(|t| t == name);
    if has(CONSENSUS_TABLE) && !has(ALIGNMENT_TABLE) { CONSENSUS_TABLE } else { SEQUENCE_TABLE }
}

/// A local path made absolute, else the input as given. VDB tries a bare name without `/`
/// as an accession first, asking the network, so local files must not be passed bare.
fn vdb_location(input: &str) -> Result<String> {
    let path = Path::new(input);
    if !path.exists() {
        return Ok(input.to_owned());
    }
    let absolute =
        std::fs::canonicalize(path).with_context(|| format!("failed to resolve {input}"))?;
    absolute
        .to_str()
        .map(str::to_owned)
        .with_context(|| format!("{} is not valid UTF-8", absolute.display()))
}

/// A read manager, without the pagemap thread.
fn new_manager() -> Result<VdbManager> {
    let manager = VdbManager::make_read().context("failed to create VDB manager")?;
    // Best effort: the pagemap thread only contends with our own workers.
    manager.disable_pagemap_thread().ok();
    Ok(manager)
}

/// Open the reads table: `table_name` within `database`, or the flat table at `location`.
fn open_table(
    manager: VdbManager,
    database: Option<VDatabase>,
    location: &str,
    table_name: Option<&str>,
    input: &str,
) -> Result<OpenTable> {
    let table = match (&database, table_name) {
        (Some(database), Some(name)) => database.open_table_read(name).with_context(|| {
            let tables = database.list_tables().unwrap_or_default().join(", ");
            format!("{input} has no table {name} (tables: {tables})")
        })?,
        _ => manager
            .open_table_read(location)
            .with_context(|| format!("failed to open table {input}"))?,
    };
    Ok(OpenTable { table, database, _manager: manager })
}

/// First row id and row count of `table`.
fn id_range(table: &VTable) -> Result<(i64, u64)> {
    let cursor = table.create_cursor_read()?;
    let column = cursor.add_column("READ_LEN")?;
    cursor.open()?;
    Ok(cursor.id_range(column)?)
}

/// Whether `NAME` holds original names. Some loads that dropped names still serve `NAME`,
/// as the serial spot number; that shows in the first spot's name being its id.
fn names_are_original(table: &VTable, first_spot: i64) -> Result<bool> {
    let cursor = table.create_cursor_read()?;
    let column = cursor.add_column("(ascii)NAME")?;
    cursor.open()?;
    let name = cursor.read_str(first_spot, column)?;
    Ok(!name.is_empty() && name != first_spot.to_string())
}

/// Name of the sequencing platform of spot `id`, from `PLATFORM` (`INSDC:SRA:platform_id`).
fn platform(table: &VTable, id: i64) -> Result<Option<&'static str>> {
    let cursor = table.create_cursor_read()?;
    let column = cursor.add_column("(INSDC:SRA:platform_id)PLATFORM")?;
    cursor.open()?;
    Ok(platform_name(cursor.read_u8(id, column)?))
}

/// The `SRA_PLATFORM_*` name of a platform id (`insdc/sra.vschema`), without its prefix;
/// `None` for `UNDEFINED` or an id newer than this list.
fn platform_name(id: u8) -> Option<&'static str> {
    const NAMES: [&str; 21] = [
        "UNDEFINED",
        "454",
        "ILLUMINA",
        "ABSOLID",
        "COMPLETE_GENOMICS",
        "HELICOS",
        "PACBIO_SMRT",
        "ION_TORRENT",
        "CAPILLARY",
        "OXFORD_NANOPORE",
        "ELEMENT_BIO",
        "TAPESTRI",
        "VELA_DIAG",
        "GENAPSYS",
        "ULTIMA",
        "GENEMIND",
        "BGISEQ",
        "DNBSEQ",
        "SINGULAR_GENOMICS",
        "GENEUS_TECH",
        "SALUS",
    ];
    NAMES.get(usize::from(id)).copied().filter(|&name| name != "UNDEFINED")
}

/// Whether the metadata carries SRA Lite's marker, `SOFTWARE/delite@name = "delite"`.
fn is_sra_lite(metadata: &KMetadata) -> bool {
    metadata
        .open_node_read("SOFTWARE/delite")
        .and_then(|node| node.read_attr("name"))
        .is_ok_and(|name| name == "delite")
}

/// The phred value of every base when `STATS/QUALITY` has a single `PHRED_<n>` bin.
fn constant_quality(metadata: &KMetadata) -> Option<u8> {
    let bins = metadata.open_node_read("STATS/QUALITY").ok()?.list_children().ok()?;
    match bins.as_slice() {
        [only] => only.strip_prefix("PHRED_")?.parse().ok(),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accession_is_the_file_name_without_sra() {
        assert_eq!(default_accession("/data/runs/SRR123.sra"), "SRR123");
    }

    #[test]
    fn accession_drops_sralite() {
        assert_eq!(default_accession("SRR123.sralite"), "SRR123");
    }

    #[test]
    fn accession_drops_lite_sra() {
        assert_eq!(default_accession("ERR015558.lite.sra"), "ERR015558");
    }

    #[test]
    fn bare_accession_is_kept() {
        assert_eq!(default_accession("SRR123"), "SRR123");
    }

    #[test]
    fn url_accession_is_the_last_path_segment() {
        let url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR123/SRR123";
        assert_eq!(default_accession(url), "SRR123");
    }

    #[test]
    fn trailing_slash_is_ignored() {
        assert_eq!(default_accession("/data/SRR123/"), "SRR123");
    }

    #[test]
    fn accession_drops_the_version_of_a_downloaded_object() {
        assert_eq!(default_accession("SRR10984982.1"), "SRR10984982");
    }

    #[test]
    fn accession_drops_the_version_of_a_downloaded_lite_object() {
        assert_eq!(default_accession("/data/SRR6819246.sralite.1"), "SRR6819246");
        assert_eq!(default_accession("SRR18323606.lite.1"), "SRR18323606");
    }

    #[test]
    fn url_query_is_dropped() {
        let url = "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-11/SRR341/SRR341578.1?ncbi_phid=1";
        assert_eq!(default_accession(url), "SRR341578");
    }

    #[test]
    fn numeric_suffix_of_a_name_that_is_not_an_accession_is_kept() {
        assert_eq!(default_accession("sample.2"), "sample.2");
        assert_eq!(default_accession("SRRx.1"), "SRRx.1");
    }

    #[test]
    fn other_extensions_are_kept() {
        assert_eq!(default_accession("run.vdb"), "run.vdb");
    }

    fn tables(names: &[&str]) -> Vec<String> {
        names.iter().map(|&name| name.to_owned()).collect()
    }

    #[test]
    fn default_table_is_consensus_in_a_native_long_read_database() {
        let pacbio = tables(&["CONSENSUS", "PASSES", "SEQUENCE", "ZMW_METRICS"]);
        assert_eq!(default_table(&pacbio), "CONSENSUS");
        assert_eq!(default_table(&tables(&["CONSENSUS", "SEQUENCE"])), "CONSENSUS");
    }

    #[test]
    fn default_table_is_sequence_without_a_consensus_table() {
        assert_eq!(default_table(&tables(&["SEQUENCE"])), "SEQUENCE");
        assert_eq!(
            default_table(&tables(&["PRIMARY_ALIGNMENT", "REFERENCE", "SEQUENCE"])),
            "SEQUENCE"
        );
    }

    #[test]
    fn default_table_is_sequence_in_an_aligned_database() {
        let aligned = tables(&["CONSENSUS", "PRIMARY_ALIGNMENT", "REFERENCE", "SEQUENCE"]);
        assert_eq!(default_table(&aligned), "SEQUENCE");
    }

    #[test]
    fn platform_names_follow_the_schema() {
        assert_eq!(platform_name(2), Some("ILLUMINA"));
        assert_eq!(platform_name(9), Some("OXFORD_NANOPORE"));
        assert_eq!(platform_name(20), Some("SALUS"));
    }

    #[test]
    fn undefined_and_unknown_platforms_have_no_name() {
        assert_eq!(platform_name(0), None);
        assert_eq!(platform_name(200), None);
    }

    #[test]
    fn missing_local_file_is_passed_through_as_given() {
        assert_eq!(vdb_location("SRR0000000001").unwrap(), "SRR0000000001");
    }

    #[test]
    fn existing_local_file_is_made_absolute() {
        let dir = std::env::temp_dir().join(format!("fg_sra_fastq_loc_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("SRR1.sra");
        std::fs::write(&file, b"").unwrap();
        let location = vdb_location(file.to_str().unwrap()).unwrap();
        assert!(Path::new(&location).is_absolute());
        assert!(location.ends_with("SRR1.sra"));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
