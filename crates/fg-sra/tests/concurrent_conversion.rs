//! Regression test for the aligned-read concurrency crash.
//!
//! Multi-threaded `tosam` on aligned runs used to read the virtual
//! `PRIMARY_ALIGNMENT.READ` column, which makes libncbi-vdb reconstruct the
//! read from the `REFERENCE` table through an internal sub-cursor whose MRU
//! blob cache is unsynchronized. Concurrent workers corrupted that cache and
//! crashed intermittently with a native SIGSEGV/SIGBUS. The fix reconstructs
//! `READ` in fg-sra from the physically stored alignment deltas plus a
//! reference preloaded once on the main thread, so no `REFERENCE` access
//! happens while workers run.
//!
//! The race is timing-dependent, so a single well-behaved run rarely triggers
//! it. This test reproduces it the way it was originally observed: several
//! `fg-sra tosam` processes run at once (oversubscribing the CPU widens the
//! thread-scheduling race windows), repeated for a few rounds. On the buggy
//! build a large fraction of runs abort with SIGSEGV; on the fixed build every
//! run succeeds.
//!
//! Opt-in: set `FG_SRA_TEST_ALIGNED_SRA` to an aligned SRA (one with a
//! `PRIMARY_ALIGNMENT` table). Without it the test is skipped, so no test data is
//! committed. By default it drives the binary built by this crate; set
//! `FG_SRA_TEST_BIN` to point at a different `fg-sra` to stress that instead.
//!
//! Tunables (env): `FG_SRA_TEST_STRESS_ROUNDS` (default 4),
//! `FG_SRA_TEST_STRESS_PARALLEL` (default 3), `FG_SRA_TEST_STRESS_THREADS`
//! (default 8).

use std::path::PathBuf;
use std::process::Command;

fn env_usize(key: &str, default: usize) -> usize {
    // Only accept a positive value: a zero (or unparseable) setting would make
    // the stress test vacuous — zero rounds/parallelism launches no conversion,
    // and zero threads takes the sequential path instead of exercising the race.
    std::env::var(key)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&n| n > 0)
        .unwrap_or(default)
}

#[test]
fn concurrent_bam_conversion_is_crash_free() {
    let Ok(sra) = std::env::var("FG_SRA_TEST_ALIGNED_SRA") else {
        eprintln!(
            "skipping: set FG_SRA_TEST_ALIGNED_SRA to an aligned SRA path to run this stress test"
        );
        return;
    };
    let bin =
        std::env::var("FG_SRA_TEST_BIN").unwrap_or_else(|_| env!("CARGO_BIN_EXE_fg-sra").into());
    let rounds = env_usize("FG_SRA_TEST_STRESS_ROUNDS", 4);
    let parallel = env_usize("FG_SRA_TEST_STRESS_PARALLEL", 3);
    let threads = env_usize("FG_SRA_TEST_STRESS_THREADS", 8);

    // Namespace outputs by pid so concurrent test binaries don't collide.
    let pid = std::process::id();
    let out_path = |round: usize, w: usize| -> PathBuf {
        std::env::temp_dir().join(format!("fg_sra_stress_{pid}_{round}_{w}.bam"))
    };

    for round in 1..=rounds {
        // Launch `parallel` conversions concurrently, each to its own output.
        let children: Vec<_> = (0..parallel)
            .map(|w| {
                let out = out_path(round, w);
                let child = Command::new(&bin)
                    .args([
                        "tosam",
                        "--output-format",
                        "bam",
                        "-t",
                        &threads.to_string(),
                        "--output-file",
                        out.to_str().unwrap(),
                        &sra,
                    ])
                    .spawn()
                    .expect("failed to spawn fg-sra");
                (child, out)
            })
            .collect();

        // Drain every worker before asserting: waiting for and cleaning up each
        // child regardless of outcome avoids leaking processes/output files when
        // one fails. Record the first failure and assert once the round is drained.
        let mut first_failure: Option<String> = None;
        for (w, (mut child, out)) in children.into_iter().enumerate() {
            let status = child.wait().expect("failed to wait for fg-sra");
            // Clean up the output regardless of outcome.
            let _ = std::fs::remove_file(&out);
            if !status.success() && first_failure.is_none() {
                first_failure = Some(format!(
                    "round {round}/{rounds} worker {w} of concurrent `fg-sra tosam -t {threads}` \
                     failed with {status} (a crash here indicates the concurrent-cursor race)"
                ));
            }
        }
        assert!(first_failure.is_none(), "{}", first_failure.unwrap());
    }
}
