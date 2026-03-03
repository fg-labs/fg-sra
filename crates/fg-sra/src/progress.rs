//! Progress logging for long-running SRA conversions.
//!
//! Logs periodic status lines to stderr showing record counts, reference
//! completion, elapsed time, and throughput rate. Thread-safe via atomics.

use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use std::time::{Duration, Instant};

/// Logs progress to stderr at record-count interval boundaries and on
/// reference completion.
pub struct ProgressLogger {
    record_count: AtomicU64,
    refs_completed: AtomicU32,
    total_refs: u32,
    interval: u64,
    start: Instant,
}

impl ProgressLogger {
    /// Create a new progress logger.
    ///
    /// `total_refs` is the number of references to process (0 for unaligned).
    /// `interval` is the number of records between periodic log lines.
    pub fn new(total_refs: u32, interval: u64) -> Self {
        Self {
            record_count: AtomicU64::new(0),
            refs_completed: AtomicU32::new(0),
            total_refs,
            interval,
            start: Instant::now(),
        }
    }

    /// Record `additional` processed records.  Emits a log line for each
    /// interval boundary crossed (e.g. if interval=1_000_000 and we go from
    /// 999_998 to 1_000_002, one line is emitted for the 1_000_000 crossing).
    pub fn record(&self, additional: u64) {
        let prev = self.record_count.fetch_add(additional, Ordering::Relaxed);
        if self.interval == 0 {
            return;
        }
        let new = prev + additional;
        // How many interval boundaries did we cross?
        let prev_bucket = prev / self.interval;
        let new_bucket = new / self.interval;
        for _ in prev_bucket..new_bucket {
            self.log(None);
        }
    }

    /// Mark one reference as complete and unconditionally emit a log line.
    pub fn reference_done(&self) {
        self.refs_completed.fetch_add(1, Ordering::Relaxed);
        self.log(None);
    }

    /// Emit a final log line with "(complete)" suffix, unless the current
    /// record count falls exactly on an interval boundary (already logged).
    pub fn complete(&self) {
        let count = self.record_count.load(Ordering::Relaxed);
        if self.interval == 0 || count % self.interval != 0 {
            self.log(Some("(complete)"));
        }
    }

    /// Format and emit a progress line to stderr.
    fn log(&self, suffix: Option<&str>) {
        let count = self.record_count.load(Ordering::Relaxed);
        let refs = self.refs_completed.load(Ordering::Relaxed);
        let elapsed = self.start.elapsed();

        let ref_part = if self.total_refs > 0 {
            format!(" across {} of {} work items", refs, self.total_refs)
        } else {
            " unaligned".to_owned()
        };

        let rate = format_rate(count, elapsed);
        let suffix_part = match suffix {
            Some(s) => format!(" {s}"),
            None => String::new(),
        };

        eprintln!(
            "[progress] Processed {} records{} in {} ({}){}",
            format_count(count),
            ref_part,
            format_duration(elapsed),
            rate,
            suffix_part,
        );
    }
}

/// Format a `Duration` as a human-readable string: "0s", "45s", "2m 15s", "1h 30m".
fn format_duration(d: Duration) -> String {
    let secs = d.as_secs();
    if secs < 60 {
        format!("{secs}s")
    } else if secs < 3600 {
        let m = secs / 60;
        let s = secs % 60;
        if s == 0 { format!("{m}m") } else { format!("{m}m {s}s") }
    } else {
        let h = secs / 3600;
        let m = (secs % 3600) / 60;
        if m == 0 { format!("{h}h") } else { format!("{h}h {m}m") }
    }
}

/// Format a `u64` with thousands separators: 1234567 -> "1,234,567".
fn format_count(n: u64) -> String {
    if n < 1_000 {
        return n.to_string();
    }
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result
}

/// Format a throughput rate: "22,222 records/s" or "50.0 records/min" for slow rates.
fn format_rate(count: u64, elapsed: Duration) -> String {
    let secs = elapsed.as_secs_f64();
    if secs < f64::EPSILON {
        return "-- records/s".to_owned();
    }
    let per_sec = count as f64 / secs;
    if per_sec >= 1.0 {
        format!("{} records/s", format_count(per_sec as u64))
    } else {
        format!("{:.1} records/min", per_sec * 60.0)
    }
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use super::*;

    // ── format_duration ──────────────────────────────────────────────

    #[test]
    fn test_format_duration_zero() {
        assert_eq!(format_duration(Duration::from_secs(0)), "0s");
    }

    #[test]
    fn test_format_duration_seconds() {
        assert_eq!(format_duration(Duration::from_secs(45)), "45s");
    }

    #[test]
    fn test_format_duration_exact_minute() {
        assert_eq!(format_duration(Duration::from_secs(60)), "1m");
    }

    #[test]
    fn test_format_duration_minutes_and_seconds() {
        assert_eq!(format_duration(Duration::from_secs(135)), "2m 15s");
    }

    #[test]
    fn test_format_duration_exact_hour() {
        assert_eq!(format_duration(Duration::from_secs(3600)), "1h");
    }

    #[test]
    fn test_format_duration_hour_and_minutes() {
        assert_eq!(format_duration(Duration::from_secs(5400)), "1h 30m");
    }

    // ── format_count ─────────────────────────────────────────────────

    #[test]
    fn test_format_count_zero() {
        assert_eq!(format_count(0), "0");
    }

    #[test]
    fn test_format_count_small() {
        assert_eq!(format_count(999), "999");
    }

    #[test]
    fn test_format_count_thousand() {
        assert_eq!(format_count(1000), "1,000");
    }

    #[test]
    fn test_format_count_millions() {
        assert_eq!(format_count(1_234_567), "1,234,567");
    }

    // ── format_rate ──────────────────────────────────────────────────

    #[test]
    fn test_format_rate_fast() {
        assert_eq!(format_rate(100_000, Duration::from_secs(5)), "20,000 records/s");
    }

    #[test]
    fn test_format_rate_slow() {
        // 10 records in 60s = ~0.167/s = 10.0/min
        assert_eq!(format_rate(10, Duration::from_secs(60)), "10.0 records/min");
    }

    #[test]
    fn test_format_rate_zero_duration() {
        assert_eq!(format_rate(100, Duration::ZERO), "-- records/s");
    }

    // ── ProgressLogger::record ───────────────────────────────────────

    #[test]
    fn test_record_no_log_before_interval() {
        // Just verifying it doesn't panic and the count is correct.
        let logger = ProgressLogger::new(10, 1_000_000);
        logger.record(100);
        assert_eq!(logger.record_count.load(Ordering::Relaxed), 100);
    }

    #[test]
    fn test_record_crosses_interval() {
        // Logs are emitted to stderr; we mainly check the count is correct.
        let logger = ProgressLogger::new(10, 100);
        logger.record(50);
        logger.record(60); // crosses 100 boundary
        assert_eq!(logger.record_count.load(Ordering::Relaxed), 110);
    }

    #[test]
    fn test_record_crosses_multiple_intervals() {
        let logger = ProgressLogger::new(10, 100);
        logger.record(250); // crosses 100 and 200
        assert_eq!(logger.record_count.load(Ordering::Relaxed), 250);
    }

    // ── ProgressLogger::reference_done ───────────────────────────────

    #[test]
    fn test_reference_done_increments() {
        let logger = ProgressLogger::new(10, 1_000_000);
        logger.reference_done();
        logger.reference_done();
        assert_eq!(logger.refs_completed.load(Ordering::Relaxed), 2);
    }

    // ── ProgressLogger::complete ─────────────────────────────────────

    #[test]
    fn test_complete_emits_when_not_on_boundary() {
        // Count 50 with interval 100 -> not on boundary -> logs.
        let logger = ProgressLogger::new(10, 100);
        logger.record(50);
        logger.complete(); // should emit (no panic)
    }

    #[test]
    fn test_complete_skips_when_on_boundary() {
        // Count 100 with interval 100 -> on boundary -> skips.
        let logger = ProgressLogger::new(10, 100);
        logger.record(100); // crosses boundary, log emitted
        logger.complete(); // should NOT emit a duplicate
    }

    #[test]
    fn test_zero_interval_record_does_not_panic() {
        let logger = ProgressLogger::new(0, 0);
        logger.record(100);
        assert_eq!(logger.record_count.load(Ordering::Relaxed), 100);
    }

    #[test]
    fn test_zero_interval_complete_emits() {
        let logger = ProgressLogger::new(0, 0);
        logger.record(50);
        logger.complete(); // interval=0, so always emits
    }
}
