//! End-to-end test that several accessions are written, in order, to one
//! `--output-file`.
//!
//! Uses a small aligned database shipped with the vendored ncbi-vdb test data
//! (its references are embedded, so no network access is needed). Converting
//! it twice in one run must produce exactly two copies of the single-run
//! output: nothing truncated, dropped, or finished early between accessions.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The vendored aligned database used as input.
fn test_database() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../vendor/ncbi-vdb/test/vdb/db/VDB-3418.sra")
}

/// Run `fg-sra tosam --fastq` on `accessions`, writing to `output`.
fn convert_to_fastq(output: &Path, accessions: &[&Path]) {
    let result = Command::new(env!("CARGO_BIN_EXE_fg-sra"))
        .args(["tosam", "--fastq", "--output-file"])
        .arg(output)
        .args(accessions)
        .output()
        .expect("failed to run fg-sra");
    assert!(
        result.status.success(),
        "fg-sra failed with {}: {}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn multi_accession_output_is_concatenated() {
    let database = test_database();
    let dir = std::env::temp_dir().join(format!("fg-sra-multi-accession-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let single = dir.join("single.fq");
    let combined = dir.join("combined.fq");

    convert_to_fastq(&single, &[&database]);
    convert_to_fastq(&combined, &[&database, &database]);
    let single_records = std::fs::read(&single).unwrap();
    let combined_records = std::fs::read(&combined).unwrap();
    std::fs::remove_dir_all(&dir).ok();

    assert!(!single_records.is_empty(), "the test database produced no records");
    assert!(
        combined_records == [single_records.as_slice(), single_records.as_slice()].concat(),
        "combined output ({} bytes) is not the single-accession output ({} bytes) twice",
        combined_records.len(),
        single_records.len()
    );
}
