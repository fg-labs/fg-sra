//! Regression test for the concurrent-VDB-cursor data race.
//!
//! Multi-threaded `tosam` used to create every worker cursor from a single
//! shared `VDatabase`. Cursors opened on the same table share libncbi-vdb's
//! internal KColumn/blob caches, which are not thread-safe for concurrent
//! reads, so conversions crashed intermittently with a native SIGSEGV. The fix
//! gives each worker its own `VDatabase` (hence its own `VDBManager` and cache
//! graph).
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
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
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

    for round in 1..=rounds {
        // Launch `parallel` conversions concurrently, each to its own output.
        let children: Vec<_> = (0..parallel)
            .map(|w| {
                let out: PathBuf = std::env::temp_dir().join(format!("fg_sra_stress_{w}.bam"));
                Command::new(&bin)
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
                    .expect("failed to spawn fg-sra")
            })
            .collect();

        for (w, mut child) in children.into_iter().enumerate() {
            let status = child.wait().expect("failed to wait for fg-sra");
            assert!(
                status.success(),
                "round {round}/{rounds} worker {w} of concurrent `fg-sra tosam -t {threads}` \
                 failed with {status} (a crash here indicates the concurrent-cursor race)",
            );
        }
    }
}
