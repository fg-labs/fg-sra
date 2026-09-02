//! Integration tests for the `ReferenceObj` reference-base reads.
//!
//! Opt-in: set `FG_SRA_TEST_ALIGNED_SRA` to an aligned SRA path (one with a
//! `REFERENCE` table). Without it every test skips, so no test data is
//! committed and CI still compiles the file.

use fg_sra_vdb::manager::VdbManager;
use fg_sra_vdb::reference::{ReferenceList, reflist_options};

/// Open the env-provided SRA and build a 4na `ReferenceList`, or `None` to skip.
fn open_reflist() -> Option<ReferenceList> {
    let sra = std::env::var("FG_SRA_TEST_ALIGNED_SRA").ok()?;
    let mgr = VdbManager::make_read().expect("make_read");
    mgr.disable_pagemap_thread().ok();
    let db = mgr.open_db_read(&sra).expect("open_db_read");
    let opts = reflist_options::READ_4NA | reflist_options::USE_PRIMARY_IDS;
    Some(ReferenceList::make_database(&db, opts, 0).expect("make_database"))
}

#[test]
fn reference_read_into_returns_valid_4na() {
    let Some(reflist) = open_reflist() else {
        eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA");
        return;
    };
    let obj = reflist.get(0).expect("get(0)");
    let mut buf = [0u8; 100];
    let written = obj.read_into(1_000_000, &mut buf).expect("read_into") as usize;
    assert!(written > 0 && written <= buf.len());
    assert!(buf[..written].iter().all(|&b| b < 16), "READ_4NA bytes must be 0..=15");
}

#[test]
fn reference_read_spans_row_boundary() {
    // REFERENCE rows are MAX_SEQ_LEN (5000) bases; a read that crosses a row
    // boundary must equal the concatenation of the two halves.
    let Some(reflist) = open_reflist() else {
        eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA");
        return;
    };
    let obj = reflist.get(0).expect("get(0)");
    let mut whole = [0u8; 20];
    let mut first = [0u8; 10];
    let mut second = [0u8; 10];
    let nw = obj.read_into(4990, &mut whole).expect("read whole") as usize;
    let n1 = obj.read_into(4990, &mut first).expect("read first") as usize;
    let n2 = obj.read_into(5000, &mut second).expect("read second") as usize;
    assert_eq!(nw, n1 + n2);
    assert_eq!(&whole[..nw], &[&first[..n1], &second[..n2]].concat()[..]);
}

#[test]
fn reference_read_all_matches_seq_length() {
    let Some(reflist) = open_reflist() else {
        eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA");
        return;
    };
    // chr21 is a small autosome present in hg38-aligned runs.
    let obj = reflist.find("21").expect("find chr21");
    let bases = obj.read_all().expect("read_all");
    assert_eq!(bases.len(), obj.seq_length().expect("seq_length") as usize);
    assert!(bases.iter().all(|&b| b < 16));
}

#[test]
fn reference_flags_external_noncircular() {
    let Some(reflist) = open_reflist() else {
        eprintln!("skipping: set FG_SRA_TEST_ALIGNED_SRA");
        return;
    };
    let obj = reflist.get(0).expect("get(0)");
    assert!(obj.external().expect("external"), "hg38 references are external refseq objects");
    assert!(!obj.circular().expect("circular"), "reference 0 is not circular");
}
