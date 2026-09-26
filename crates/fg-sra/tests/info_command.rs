//! End-to-end tests of `fg-sra info`, driving the built binary.
//!
//! Tests that need a real archive are opt-in, so no test data is committed. Set:
//! - `FG_SRA_TEST_PAIRED_SRA` to an unaligned paired-end SRA file, e.g. SRR2584863;
//! - `FG_SRA_TEST_ALIGNED_SRA` to an aligned (cSRA) SRA file, e.g. SRR390728.
//!
//! Without them, those tests print a note and pass.

use std::process::{Command, Output};

/// Returns the archive named by env var `key`, or `None` (after a note) if it is unset.
fn archive_from_env(key: &str) -> Option<String> {
    let path = std::env::var(key).ok();
    if path.is_none() {
        eprintln!("skipping: set {key} to an SRA file to run this test");
    }
    path
}

fn info(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fg-sra")).arg("info").args(args).output().unwrap()
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

#[test]
fn missing_input_is_reported_and_fails() {
    let dir = std::env::temp_dir().join(format!("fg_sra_info_{}", std::process::id()));
    let output = info(&[dir.join("no-such-run.sra").to_str().unwrap()]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("no SRA archive found"));
}

#[test]
fn paired_run_is_described_as_pairs() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let output = info(&[&sra, "--layout-spots", "100"]);
    assert!(output.status.success());
    let text = stdout(&output);
    assert!(text.contains("read layout       first 100 of "), "{text}");
    assert!(text.contains(" spots only (--layout-spots samples more)"), "{text}");
    assert!(text.contains("biological reads  2 (pairs): 100 spots"), "{text}");
}

#[test]
fn aligned_run_reports_its_alignments() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    let output = info(&[&sra, "--layout-spots", "100"]);
    assert!(output.status.success());
    let text = stdout(&output);
    assert!(text.contains("database, aligned (cSRA)"), "{text}");
    assert!(text.contains("  alignments  "), "{text}");
}

#[test]
fn one_failing_input_does_not_stop_the_others() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let output = info(&["/no/such/run.sra", &sra, "--layout-spots", "10"]);
    assert!(!output.status.success());
    assert!(stdout(&output).contains("biological reads"));
}
