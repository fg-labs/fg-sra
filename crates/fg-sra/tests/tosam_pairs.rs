//! End-to-end test that `tosam -u` writes a partly aligned spot's reads as a pair, as
//! sam-dump does: the aligned read keeps its pairing flags (with 0x8, mate unmapped), and
//! the unaligned read is paired with it, its RNEXT and PNEXT giving the aligned read's
//! position and its 0x20 the aligned read's strand (which is its `READ_TYPE`'s REVERSE
//! bit, which 0x20 comes from, in archives loaded by bam-load).
//!
//! Only pairs are checked: a record that is neither first nor last in its template (the
//! middle read of three) is skipped.
//!
//! No vendored archive has partly aligned spots, so this is opt-in: set
//! `FG_SRA_TEST_ALIGNED_SRA` to an aligned archive (or accession) with some unaligned mates.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::process::{Child, Command, Stdio};

const PAIRED: u32 = 0x1;
const UNMAPPED: u32 = 0x4;
const MATE_UNMAPPED: u32 = 0x8;
const REVERSE: u32 = 0x10;
const MATE_REVERSE: u32 = 0x20;
const FIRST: u32 = 0x40;
const LAST: u32 = 0x80;
/// Secondary and supplementary alignments: not the aligned read of a spot's pair.
const NOT_PRIMARY: u32 = 0x100 | 0x800;

/// A child process, killed and reaped if the test panics before waiting for it.
struct Reaped(Child);

impl Drop for Reaped {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

/// An aligned read whose mate is unaligned: its reference, 1-based position and flags.
struct AlignedMate {
    rname: String,
    pos: String,
    flags: u32,
}

#[test]
fn partly_aligned_spots_are_written_as_pairs() {
    let Ok(archive) = std::env::var("FG_SRA_TEST_ALIGNED_SRA") else {
        eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA to an aligned SRA to run this test");
        return;
    };
    let mut child = Reaped(
        Command::new(env!("CARGO_BIN_EXE_fg-sra"))
            .args(["tosam", "-u", "--no-header", &archive])
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("failed to run fg-sra"),
    );

    // Aligned reads with an unaligned mate, by QNAME and which segment they are (the
    // aligned records all come before the unaligned ones).
    let mut aligned: HashMap<(String, u32), AlignedMate> = HashMap::new();
    let mut pairs = 0;
    for line in BufReader::new(child.0.stdout.take().unwrap()).lines() {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let flags: u32 = fields[1].parse().unwrap();
        let segment = flags & (FIRST | LAST);
        if segment != FIRST && segment != LAST {
            continue;
        }
        if flags & UNMAPPED == 0 {
            if flags & MATE_UNMAPPED != 0 && flags & NOT_PRIMARY == 0 {
                assert_ne!(flags & PAIRED, 0, "{}: 0x8 without 0x1", fields[0]);
                let mate =
                    AlignedMate { rname: fields[2].to_string(), pos: fields[3].to_string(), flags };
                aligned.insert((fields[0].to_string(), segment), mate);
            }
            continue;
        }
        if fields[6] == "*" {
            continue; // The mate is unaligned too, or there is none.
        }
        let mate_segment = segment ^ (FIRST | LAST);
        let mate = aligned
            .remove(&(fields[0].to_string(), mate_segment))
            .unwrap_or_else(|| panic!("{}: no aligned mate for an unaligned read", fields[0]));
        assert_ne!(flags & PAIRED, 0, "{}: unaligned read not paired", fields[0]);
        assert_eq!(flags & MATE_UNMAPPED, 0, "{}: aligned mate called unmapped", fields[0]);
        assert_eq!((fields[6], fields[7]), (mate.rname.as_str(), mate.pos.as_str()));
        assert_eq!(flags & MATE_REVERSE != 0, mate.flags & REVERSE != 0, "{}: 0x20", fields[0]);
        pairs += 1;
    }
    assert!(child.0.wait().unwrap().success(), "fg-sra failed");
    assert!(pairs > 0, "the archive has no partly aligned spots");
    assert!(
        aligned.is_empty(),
        "{} aligned reads' unaligned mates were not written",
        aligned.len()
    );
}
