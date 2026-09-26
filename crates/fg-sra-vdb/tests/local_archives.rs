//! Tests of path types, flat-table access and table metadata that need no network.
//!
//! Tests that need a real archive are opt-in, so no test data is committed. Set:
//! - `FG_SRA_TEST_TABLE_SRA` to a flat-table (unaligned) SRA file, e.g. SRR2584863;
//! - `FG_SRA_TEST_ALIGNED_SRA` to an aligned (cSRA) SRA file, e.g. SRR390728.
//!
//! Without them, those tests print a note and pass.

use std::path::PathBuf;

use fg_sra_vdb::cursor::{BlobColumn, VCursor};
use fg_sra_vdb::database::VTable;
use fg_sra_vdb::error::VdbError;
use fg_sra_vdb::manager::{PathType, VdbManager};

/// Returns the archive named by env var `key`, or `None` (after a note) if it is unset.
fn archive_from_env(key: &str) -> Option<String> {
    let path = std::env::var(key).ok();
    if path.is_none() {
        eprintln!("skipping: set {key} to an SRA file to run this test");
    }
    path
}

/// Returns a fresh directory under the system temp dir, namespaced by process and `name`.
fn scratch_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fg_sra_vdb_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Returns the number of rows in `table`, from a cursor over its `READ_LEN` column.
fn row_count(table: &VTable) -> u64 {
    let cursor = table.create_cursor_read().unwrap();
    let col = cursor.add_column("READ_LEN").unwrap();
    cursor.open().unwrap();
    cursor.id_range(col).unwrap().1
}

fn manager() -> VdbManager {
    VdbManager::make_read().unwrap()
}

/// Rows read by the blob-read tests: enough to span several blobs of any column.
const BLOB_TEST_ROWS: std::ops::RangeInclusive<i64> = 1..=20_000;

/// Asserts that reading `column` of `table` through its blobs, at each of `rows` in turn, gives
/// the cells that reading it cell by cell does. Each way reads through its own cursor.
fn assert_blob_reads_match_cell_reads<T: PartialEq + std::fmt::Debug>(
    table: &VTable,
    column: &str,
    rows: impl Iterator<Item = i64>,
    cell_read: impl Fn(&VCursor, i64, u32, &mut Vec<T>) -> Result<(), VdbError>,
    blob_read: impl Fn(&mut BlobColumn, &VCursor, i64, &mut Vec<T>) -> Result<(), VdbError>,
) {
    let open = || {
        let cursor = table.create_cursor_read().unwrap();
        let col = cursor.add_column(column).unwrap();
        cursor.open().unwrap();
        (cursor, col)
    };
    let (cell_cursor, cell_col) = open();
    let (blob_cursor, blob_col) = open();
    let mut blobs = BlobColumn::new(blob_col);
    let (mut expected, mut actual) = (Vec::new(), Vec::new());
    for row in rows {
        cell_read(&cell_cursor, row, cell_col, &mut expected).unwrap();
        blob_read(&mut blobs, &blob_cursor, row, &mut actual).unwrap();
        assert_eq!(actual, expected, "{column} row {row}");
    }
}

#[test]
fn path_type_of_a_missing_file_is_not_found() {
    let path = scratch_dir("missing").join("no-such-run.sra");
    assert_eq!(manager().path_type(path.to_str().unwrap()).unwrap(), PathType::NotFound);
}

#[test]
fn path_type_of_a_plain_directory_is_other() {
    let dir = scratch_dir("plain-dir");
    assert!(matches!(manager().path_type(dir.to_str().unwrap()).unwrap(), PathType::Other(_)));
}

#[test]
fn paths_containing_a_percent_sign_are_refused() {
    let path = scratch_dir("percent-open").join("run%s%s%s%d.sra");
    let path = path.to_str().unwrap();
    let mgr = manager();
    assert_eq!(mgr.path_type(path).err(), Some(VdbError::PercentInPath));
    assert_eq!(mgr.open_db_read(path).err(), Some(VdbError::PercentInPath));
    assert_eq!(mgr.open_table_read(path).err(), Some(VdbError::PercentInPath));
}

#[test]
fn flat_table_is_reported_as_a_table() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    assert_eq!(manager().path_type(&sra).unwrap(), PathType::Table);
}

#[test]
fn flat_table_cannot_be_opened_as_a_database() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    assert!(manager().open_db_read(&sra).is_err());
}

#[test]
fn flat_table_spot_count_statistic_matches_its_row_count() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let table = manager().open_table_read(&sra).unwrap();
    let meta = table.open_metadata_read().unwrap();
    let spot_count = meta.open_node_read("STATS/TABLE/SPOT_COUNT").unwrap().read_u64().unwrap();
    assert!(spot_count > 0);
    assert_eq!(spot_count, row_count(&table));
}

#[test]
fn flat_table_physical_columns_include_read_and_quality() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let columns = manager().open_table_read(&sra).unwrap().list_physical_columns().unwrap();
    assert!(columns.iter().any(|c| c == "READ"), "physical columns: {columns:?}");
    assert!(columns.iter().any(|c| c == "QUALITY"), "physical columns: {columns:?}");
}

#[test]
fn flat_table_statistics_node_lists_its_children() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let meta = manager().open_table_read(&sra).unwrap().open_metadata_read().unwrap();
    let children = meta.open_node_read("STATS/TABLE").unwrap().list_children().unwrap();
    assert!(children.iter().any(|c| c == "SPOT_COUNT"), "children: {children:?}");
    assert!(children.iter().any(|c| c == "BASE_COUNT"), "children: {children:?}");
}

#[test]
fn missing_metadata_node_is_not_found() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let meta = manager().open_table_read(&sra).unwrap().open_metadata_read().unwrap();
    let err = meta.open_node_read("NO/SUCH/NODE").err().expect("node should not exist");
    assert!(err.is_not_found(), "{err}");
}

#[test]
fn loader_name_attribute_is_readable() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let meta = manager().open_table_read(&sra).unwrap().open_metadata_read().unwrap();
    let loader = meta.open_node_read("SOFTWARE/loader").unwrap();
    assert!(!loader.read_attr("name").unwrap().is_empty());
}

#[test]
fn missing_attribute_is_not_found() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let meta = manager().open_table_read(&sra).unwrap().open_metadata_read().unwrap();
    let loader = meta.open_node_read("SOFTWARE/loader").unwrap();
    let err = loader.read_attr("no-such-attribute").expect_err("attribute should not exist");
    assert!(err.is_not_found(), "{err}");
}

#[test]
fn aligned_archive_is_reported_as_a_database() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    assert_eq!(manager().path_type(&sra).unwrap(), PathType::Database);
}

#[test]
fn aligned_archive_cannot_be_opened_as_a_table() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    assert!(manager().open_table_read(&sra).is_err());
}

#[test]
fn aligned_sequence_table_stores_cmp_read_not_read() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    let db = manager().open_db_read(&sra).unwrap();
    let columns = db.open_table_read("SEQUENCE").unwrap().list_physical_columns().unwrap();
    assert!(columns.iter().any(|c| c == "CMP_READ"), "physical columns: {columns:?}");
    assert!(!columns.iter().any(|c| c == "READ"), "physical columns: {columns:?}");
}

#[test]
fn aligned_sequence_spot_count_statistic_matches_its_row_count() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    let db = manager().open_db_read(&sra).unwrap();
    let sequence = db.open_table_read("SEQUENCE").unwrap();
    let meta = sequence.open_metadata_read().unwrap();
    let spot_count = meta.open_node_read("STATS/TABLE/SPOT_COUNT").unwrap().read_u64().unwrap();
    assert_eq!(spot_count, row_count(&sequence));
}

#[test]
fn flat_table_blob_reads_match_cell_reads() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let table = manager().open_table_read(&sra).unwrap();
    for column in
        ["(INSDC:4na:bin)READ", "(INSDC:quality:phred)QUALITY", "(INSDC:SRA:xread_type)READ_TYPE"]
    {
        assert_blob_reads_match_cell_reads(
            &table,
            column,
            BLOB_TEST_ROWS,
            VCursor::read_u8_slice_into,
            BlobColumn::read_u8_slice_into,
        );
    }
    assert_blob_reads_match_cell_reads(
        &table,
        "(INSDC:coord:zero)READ_START",
        BLOB_TEST_ROWS,
        VCursor::read_i32_slice_into,
        BlobColumn::read_i32_slice_into,
    );
    assert_blob_reads_match_cell_reads(
        &table,
        "(INSDC:coord:len)READ_LEN",
        BLOB_TEST_ROWS,
        VCursor::read_u32_slice_into,
        BlobColumn::read_u32_slice_into,
    );
}

#[test]
fn blob_reads_in_descending_row_order_match_cell_reads() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let table = manager().open_table_read(&sra).unwrap();
    assert_blob_reads_match_cell_reads(
        &table,
        "(INSDC:quality:phred)QUALITY",
        BLOB_TEST_ROWS.rev(),
        VCursor::read_u8_slice_into,
        BlobColumn::read_u8_slice_into,
    );
}

#[test]
fn aligned_sequence_blob_reads_match_cell_reads() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    let db = manager().open_db_read(&sra).unwrap();
    let sequence = db.open_table_read("SEQUENCE").unwrap();
    assert_blob_reads_match_cell_reads(
        &sequence,
        "(INSDC:4na:bin)CMP_READ",
        BLOB_TEST_ROWS,
        VCursor::read_u8_slice_into,
        BlobColumn::read_u8_slice_into,
    );
    assert_blob_reads_match_cell_reads(
        &sequence,
        "(I64)PRIMARY_ALIGNMENT_ID",
        BLOB_TEST_ROWS,
        VCursor::read_i64_slice_into,
        BlobColumn::read_i64_slice_into,
    );
}

#[test]
fn blob_reads_are_refused_on_a_cached_cursor() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_TABLE_SRA") else { return };
    let table = manager().open_table_read(&sra).unwrap();
    let cursor = table.create_cached_cursor_read(1 << 20).unwrap();
    let mut blobs = BlobColumn::new(cursor.add_column("(INSDC:coord:len)READ_LEN").unwrap());
    cursor.open().unwrap();
    let mut cell = Vec::new();
    assert_eq!(
        blobs.read_u32_slice_into(&cursor, 1, &mut cell).err(),
        Some(VdbError::BlobReadOnCachedCursor)
    );
}
