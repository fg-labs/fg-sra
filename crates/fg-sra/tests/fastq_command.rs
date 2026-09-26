//! End-to-end tests of `fg-sra fastq`, driving the built binary.
//!
//! Tests that need a real archive are opt-in, so no test data is committed. Set:
//! - `FG_SRA_TEST_PAIRED_SRA` to an unaligned paired-end SRA file, e.g. SRR2584863;
//! - `FG_SRA_TEST_ALIGNED_SRA` to an aligned (cSRA) SRA file whose references are available
//!   locally, e.g. SRR390728 fetched with `prefetch SRR390728`, which puts them beside it (the
//!   missing-reference test needs one with external references, and skips others);
//! - `FG_SRA_TEST_COLOR_SPACE_SRA` to a colour-space run stored as `CSREAD` with no physical
//!   `READ`, e.g. ERR048905.
//!
//! Without them, those tests print a note and pass. Tests that need only some archive use
//! a small aligned database from the vendored ncbi-vdb tests.

use std::fs::{File, Permissions};
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// A small aligned database with embedded references, from the vendored ncbi-vdb tests.
fn vendored_archive() -> String {
    path(Path::new(env!("CARGO_MANIFEST_DIR")), "../../vendor/ncbi-vdb/test/vdb/db/VDB-3418.sra")
}

/// Returns the archive named by env var `key`, or `None` (after a note) if it is unset.
fn archive_from_env(key: &str) -> Option<String> {
    let path = std::env::var(key).ok();
    if path.is_none() {
        eprintln!("skipping: set {key} to an SRA file to run this test");
    }
    path
}

/// A fresh directory for one test's outputs.
fn scratch_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fg_sra_fastq_cmd_{}_{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn fastq(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fg-sra")).arg("fastq").args(args).output().unwrap()
}

fn path(dir: &Path, name: &str) -> String {
    dir.join(name).to_str().unwrap().to_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn fastq_records(file: &str) -> usize {
    std::fs::read_to_string(file).unwrap().lines().count() / 4
}

/// The value of `column` in a one-row metrics TSV.
fn metric(file: &str, column: &str) -> String {
    let text = std::fs::read_to_string(file).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let values: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let index = header.iter().position(|&h| h == column).unwrap();
    values[index].to_owned()
}

#[test]
fn missing_input_fails_with_a_clear_error() {
    let dir = scratch_dir("missing");
    let output = fastq(&[&path(&dir, "no-such-run.sra"), "-u", &path(&dir, "u.fq")]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("no SRA archive found"), "{}", stderr(&output));
}

#[test]
fn input_path_with_a_percent_sign_fails_cleanly() {
    let dir = scratch_dir("percent%s");
    let output = fastq(&[&path(&dir, "run.sra"), "-u", &path(&dir, "u.fq")]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("'%'"), "{}", stderr(&output));
}

#[test]
fn paired_run_writes_equal_numbers_of_mates_and_passes_integrity() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let dir = scratch_dir("paired");
    let (r1, r2, metrics) = (path(&dir, "r1.fq"), path(&dir, "r2.fq"), path(&dir, "m.tsv"));
    let output = fastq(&[&sra, "-1", &r1, "-2", &r2, "-m", &metrics]);
    assert!(output.status.success(), "{}", stderr(&output));
    let pairs = fastq_records(&r1);
    assert!(pairs > 0);
    assert_eq!(pairs, fastq_records(&r2));
    assert_eq!(metric(&metrics, "pairs_written"), pairs.to_string());
    assert_eq!(metric(&metrics, "integrity"), "pass");
}

#[test]
fn bgzf_output_is_identical_at_one_and_eight_threads() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let convert = |threads: &str| {
        let dir = scratch_dir(&format!("threads{threads}"));
        let (r1, r2) = (path(&dir, "r1.fq.gz"), path(&dir, "r2.fq.gz"));
        let output = fastq(&[&sra, "-1", &r1, "-2", &r2, "-t", threads, "--max-spot-id", "300000"]);
        assert!(output.status.success(), "{}", stderr(&output));
        (std::fs::read(r1).unwrap(), std::fs::read(r2).unwrap())
    };
    assert!(convert("1") == convert("8"), "BGZF output differs between 1 and 8 threads");
}

#[test]
fn spot_range_converts_only_those_spots_and_skips_integrity() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let dir = scratch_dir("range");
    let (r1, r2, metrics) = (path(&dir, "r1.fq"), path(&dir, "r2.fq"), path(&dir, "m.tsv"));
    let range = ["--min-spot-id", "11", "--max-spot-id", "1010"];
    let output =
        fastq(&[&[sra.as_str(), "-1", &r1, "-2", &r2, "-m", &metrics][..], &range].concat());
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(metric(&metrics, "spots_converted"), "1000");
    assert_eq!(metric(&metrics, "integrity"), "n/a");
    let first_name = std::fs::read_to_string(&r1).unwrap().lines().next().unwrap().to_owned();
    assert!(first_name.ends_with(".11"), "{first_name}");
}

#[test]
fn paired_run_sent_only_to_unpaired_fails_and_leaves_no_output() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let dir = scratch_dir("nothing-written");
    let unpaired = path(&dir, "u.fq");
    let output = fastq(&[&sra, "-u", &unpaired, "--max-spot-id", "1000"]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("did you mean --r1/--r2?"), "{}", stderr(&output));
    assert!(!Path::new(&unpaired).exists());
}

#[test]
fn failed_run_keeps_an_existing_output_it_could_not_open() {
    let dir = scratch_dir("read-only-output");
    let (r1, r2, unpaired) = (path(&dir, "r1.fq"), path(&dir, "r2.fq"), path(&dir, "u.fq"));
    std::fs::write(&unpaired, "kept\n").unwrap();
    std::fs::set_permissions(&unpaired, Permissions::from_mode(0o444)).unwrap();
    if File::options().append(true).open(&unpaired).is_ok() {
        eprintln!("skipping: permissions don't stop this user writing {unpaired}");
        return;
    }
    let output = fastq(&[&vendored_archive(), "-1", &r1, "-2", &r2, "-u", &unpaired]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("failed to create"), "{}", stderr(&output));
    assert!(!Path::new(&r1).exists());
    assert!(!Path::new(&r2).exists());
    assert_eq!(std::fs::read_to_string(&unpaired).unwrap(), "kept\n");
}

#[test]
fn defline_template_names_reads() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_PAIRED_SRA") else { return };
    let dir = scratch_dir("defline");
    let (r1, r2) = (path(&dir, "r1.fa"), path(&dir, "r2.fa"));
    let args = ["--fasta", "--accession", "RUN", "--defline", "$ac:$si/$ri", "--max-spot-id", "1"];
    let output = fastq(&[&[sra.as_str(), "-1", &r1, "-2", &r2][..], &args].concat());
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(std::fs::read_to_string(&r1).unwrap().starts_with(">RUN:1/1\n"));
    assert!(std::fs::read_to_string(&r2).unwrap().starts_with(">RUN:1/2\n"));
}

#[test]
fn offline_run_does_not_resolve_an_accession() {
    let dir = scratch_dir("offline");
    let output = fastq(&["SRR0000000001", "-u", &path(&dir, "u.fq"), "--offline"]);
    assert!(!output.status.success());
    assert!(stderr(&output).contains("no SRA archive found"), "{}", stderr(&output));
}

#[test]
fn aligned_archive_is_converted_offline_from_its_local_references() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    let dir = scratch_dir("aligned");
    let (r1, r2) = (path(&dir, "r1.fq"), path(&dir, "r2.fq"));
    let output = fastq(&[&sra, "-1", &r1, "-2", &r2, "--offline", "--max-spot-id", "1000"]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert!(stderr(&output).contains("references"), "{}", stderr(&output));
    assert_eq!(fastq_records(&r1), fastq_records(&r2));
    let text = std::fs::read_to_string(&r1).unwrap();
    let bases_are_dna = text
        .lines()
        .skip(1)
        .step_by(4)
        .all(|line| !line.is_empty() && line.bytes().all(|b| b"ACGTNMRWSYKVHDB".contains(&b)));
    assert!(bases_are_dna, "{}", &text[..text.len().min(400)]);
}

/// The number of references `fg-sra info` reports as external (not embedded) for `sra`.
fn external_references(sra: &str) -> u32 {
    let output = Command::new(env!("CARGO_BIN_EXE_fg-sra"))
        .args(["info", sra, "--layout-spots", "1"])
        .output()
        .unwrap();
    assert!(output.status.success(), "{}", stderr(&output));
    let text = String::from_utf8_lossy(&output.stdout);
    let (_, after) = text.split_once(" references (").expect(&text);
    let (external, _) = after.split_once(" external)").expect(&text);
    external.parse().unwrap()
}

#[test]
fn aligned_archive_without_its_references_names_the_missing_one_offline() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_ALIGNED_SRA") else { return };
    if external_references(&sra) == 0 {
        eprintln!("skipping: {sra} embeds all its references, so none can be missing");
        return;
    }
    let dir = scratch_dir("aligned-no-references");
    let alone = dir.join("alone.sra");
    if std::fs::hard_link(&sra, &alone).is_err() {
        std::fs::copy(&sra, &alone).unwrap();
    }
    // Run where no reference can be found: a home (and so reference cache) and working
    // directory of the scratch directory alone.
    let output = Command::new(env!("CARGO_BIN_EXE_fg-sra"))
        .args(["fastq", alone.to_str().unwrap(), "-u", &path(&dir, "u.fq"), "--offline"])
        .current_dir(&dir)
        .env("HOME", &dir)
        .env_remove("NCBI_SETTINGS")
        .output()
        .unwrap();
    assert!(!output.status.success());
    let message = stderr(&output);
    assert!(message.contains("failed to read reference "), "{message}");
}

#[test]
fn colour_space_run_is_converted_to_base_space() {
    let Some(sra) = archive_from_env("FG_SRA_TEST_COLOR_SPACE_SRA") else { return };
    let dir = scratch_dir("colour-space");
    let (r1, r2, metrics) = (path(&dir, "r1.fq"), path(&dir, "r2.fq"), path(&dir, "m.tsv"));
    let output = fastq(&[&sra, "-1", &r1, "-2", &r2, "-m", &metrics, "-t", "8"]);
    assert!(output.status.success(), "{}", stderr(&output));
    assert_eq!(metric(&metrics, "integrity"), "pass");
    let text = std::fs::read_to_string(&r1).unwrap();
    let bases_are_dna =
        text.lines().skip(1).step_by(4).all(|line| line.bytes().all(|b| b"ACGTN".contains(&b)));
    assert!(bases_are_dna, "{}", &text[..text.len().min(400)]);
}
