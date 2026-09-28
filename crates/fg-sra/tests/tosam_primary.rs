//! End-to-end test that `tosam -u` writes every read exactly one primary record.
//!
//! Uses a small aligned database shipped with the vendored ncbi-vdb test data. Its spot 4
//! has a single read whose only alignment is secondary (`PRIMARY_ALIGNMENT_ID` 0,
//! `ALIGNMENT_COUNT` 1). `tosam -u` writes that read's secondary alignment and, since the
//! read has no primary alignment, an unmapped primary record, so tools that collect a read's
//! primary record (e.g. `samtools fastq`) still see it. sam-dump, which counts secondary
//! alignments when deciding whether a spot is aligned, writes only the secondary record.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const UNMAPPED: u32 = 0x4;
const FIRST: u32 = 0x40;
const LAST: u32 = 0x80;
const SECONDARY: u32 = 0x100;
const SUPPLEMENTARY: u32 = 0x800;

/// The vendored aligned database used as input.
fn test_database() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../vendor/ncbi-vdb/test/vdb/db/blob_val_inv_chsum.sra")
}

/// Run `fg-sra tosam -u --no-header` on the test database, returning each record's QNAME and
/// FLAG.
fn convert() -> Vec<(String, u32)> {
    let result = Command::new(env!("CARGO_BIN_EXE_fg-sra"))
        .args(["tosam", "-u", "--no-header"])
        .arg(test_database())
        .output()
        .expect("failed to run fg-sra");
    assert!(
        result.status.success(),
        "fg-sra failed with {}: {}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );
    String::from_utf8(result.stdout)
        .unwrap()
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            (fields[0].to_string(), fields[1].parse().unwrap())
        })
        .collect()
}

#[test]
fn every_read_has_one_primary_record() {
    let records = convert();
    // Primary records per read, by QNAME and which segment it is; secondary records' reads
    // are counted too, so a read with only secondary records shows up with none.
    let mut primaries: HashMap<(String, u32), usize> = HashMap::new();
    for (qname, flags) in &records {
        let count = primaries.entry((qname.clone(), flags & (FIRST | LAST))).or_default();
        if flags & (SECONDARY | SUPPLEMENTARY) == 0 {
            *count += 1;
        }
    }
    for (read, count) in &primaries {
        assert_eq!(*count, 1, "{read:?}: {count} primary records");
    }

    // Spot 4's read: its secondary alignment, and an unmapped primary record.
    let spot4: Vec<u32> =
        records.iter().filter(|(qname, _)| qname == "4").map(|&(_, flags)| flags).collect();
    assert_eq!(spot4.len(), 2, "spot 4: {spot4:?}");
    assert!(spot4.iter().any(|flags| flags & SECONDARY != 0 && flags & UNMAPPED == 0));
    assert!(spot4.contains(&UNMAPPED), "spot 4 has no unmapped primary record: {spot4:?}");
}
