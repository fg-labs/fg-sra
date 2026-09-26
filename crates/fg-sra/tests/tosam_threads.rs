//! End-to-end test that `tosam` output doesn't depend on `--threads`.
//!
//! Uses a small aligned database with two references, shipped with the vendored
//! ncbi-vdb test data (its references are embedded, so no network access is
//! needed). With several threads the references are preloaded in parallel and the
//! alignments converted on worker threads; the output must match a one-thread run
//! byte for byte.

use std::path::{Path, PathBuf};
use std::process::Command;

/// The vendored aligned database used as input.
fn test_database() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../vendor/ncbi-vdb/test/vdb/db/VDB-3418.sra")
}

/// Run `fg-sra tosam` with `args` on the test database, returning its output file's bytes.
fn convert(dir: &Path, name: &str, args: &[&str]) -> Vec<u8> {
    let output = dir.join(name);
    let result = Command::new(env!("CARGO_BIN_EXE_fg-sra"))
        .arg("tosam")
        .args(args)
        .arg("--output-file")
        .arg(&output)
        .arg(test_database())
        .output()
        .expect("failed to run fg-sra");
    assert!(
        result.status.success(),
        "fg-sra failed with {}: {}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );
    std::fs::read(&output).unwrap()
}

#[test]
fn output_is_identical_at_any_thread_count() {
    let dir = std::env::temp_dir().join(format!("fg-sra-tosam-threads-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let formats: [(&str, &[&str]); 3] = [
        ("sam", &["--unaligned"]),
        ("sam.gz", &["--unaligned", "--gzip"]),
        ("bam", &["--unaligned", "--output-format", "bam"]),
    ];
    for (extension, args) in formats {
        let one = convert(&dir, &format!("t1.{extension}"), &[args, &["-t", "1"]].concat());
        let eight = convert(&dir, &format!("t8.{extension}"), &[args, &["-t", "8"]].concat());
        assert!(!one.is_empty(), "{extension}: no output");
        assert!(one == eight, "{extension}: output differs between 1 and 8 threads");
    }
    // The compressed outputs hold the SAM itself.
    let sam = std::fs::read(dir.join("t1.sam")).unwrap();
    let mut gunzipped = Vec::new();
    std::io::Read::read_to_end(
        &mut flate2::read::MultiGzDecoder::new(std::fs::File::open(dir.join("t1.sam.gz")).unwrap()),
        &mut gunzipped,
    )
    .unwrap();
    assert!(gunzipped == sam, "gzip output is not the SAM output");
    std::fs::remove_dir_all(&dir).ok();
}
