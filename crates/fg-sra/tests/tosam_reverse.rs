//! End-to-end test of `tosam --reverse` on unaligned reads.
//!
//! Uses a small aligned database shipped with the vendored ncbi-vdb test data (its
//! references are embedded, so no network access is needed). Of its 47 unaligned
//! spots, each a single biological read, 27 have `READ_TYPE` `BIOLOGICAL|REVERSE`
//! and 20 `BIOLOGICAL|FORWARD`. `--reverse` must reverse-complement exactly the
//! REVERSE ones and set their 0x10 flag, as `sam-dump --reverse` does, and leave
//! FASTQ output as stored.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The vendored aligned database used as input.
fn test_database() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../fg-sra-vdb-sys/vendor/ncbi-vdb/test/vdb/db/VDB-3418.sra")
}

/// Run `fg-sra tosam --unaligned-spots-only --no-header` with `args` on the test
/// database, returning each SAM record's fields.
fn convert(args: &[&str]) -> Vec<Vec<String>> {
    let result = Command::new(env!("CARGO_BIN_EXE_fg-sra"))
        .args(["tosam", "--unaligned-spots-only", "--no-header"])
        .args(args)
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
        .map(|line| line.split('\t').map(str::to_string).collect())
        .collect()
}

/// The reverse complement of `bases`.
fn reverse_complement(bases: &str) -> String {
    bases
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other,
        })
        .collect()
}

#[test]
fn reverse_flips_only_reads_typed_reverse() {
    let stored = convert(&[]);
    let reversed = convert(&["--reverse"]);
    assert_eq!(stored.len(), 47);
    assert_eq!(reversed.len(), stored.len());

    let mut num_reversed = 0;
    for (as_stored, record) in stored.iter().zip(&reversed) {
        assert_eq!(as_stored[0], record[0], "records are in a different order");
        assert_eq!(as_stored[1], "4", "{}: flag without --reverse", as_stored[0]);
        match record[1].as_str() {
            "20" => {
                num_reversed += 1;
                assert_eq!(record[9], reverse_complement(&as_stored[9]), "{}: SEQ", record[0]);
                let quality: String = as_stored[10].chars().rev().collect();
                assert_eq!(record[10], quality, "{}: QUAL", record[0]);
            }
            "4" => assert_eq!(record, as_stored),
            flag => panic!("{}: unexpected flag {flag}", record[0]),
        }
    }
    assert_eq!(num_reversed, 27);
}

#[test]
fn reverse_leaves_fastq_as_stored() {
    let stored = convert(&["--fastq"]);
    assert_eq!(stored.len(), 47 * 4);
    assert_eq!(convert(&["--fastq", "--reverse"]), stored);
}
