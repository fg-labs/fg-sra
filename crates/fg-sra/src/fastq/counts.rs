//! Tallies of spots, reads and bases over a conversion, and the metrics reported from them.

use serde::Serialize;

use super::spot::{DropReason, ReadFilter, Spot, SpotOutcome};
use crate::archive::{Archive, QualitySource};
use crate::progress::format_count;

/// Counts of what a conversion read and wrote. Workers keep one per batch and the totals
/// are summed with [`SpotCounts::merge`].
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SpotCounts {
    /// Spots read, whatever became of them.
    pub spots: u64,
    /// Spots the subsampler skipped, unread.
    pub spots_skipped_by_subsampling: u64,
    pub pairs_written: u64,
    pub unpaired_written: u64,
    pub spots_no_bio_reads: u64,
    pub spots_too_many_bio_reads: u64,
    pub spots_filtered: u64,
    pub spots_technical_mismatch: u64,
    pub pairs_unrouted: u64,
    pub unpaired_unrouted: u64,
    /// Non-empty reads, by `[technical, biological][READ_FILTER]`.
    pub reads: [[u64; 4]; 2],
    /// Bases, by `[technical, biological][READ_FILTER]`.
    pub bases: [[u64; 4]; 2],
}

impl SpotCounts {
    /// Count `spot`, whose outcome was `outcome`.
    ///
    /// Reads and bases are counted for every spot, written or not, so the totals can be
    /// checked against the archive's stored statistics.
    pub fn record(&mut self, spot: &Spot<'_>, outcome: SpotOutcome) {
        self.spots += 1;
        for slot in 0..spot.read_count() {
            let len = spot.read_len(slot);
            if len == 0 {
                continue;
            }
            let read_type = usize::from(spot.is_biological(slot));
            let filter = spot.filter(slot).index();
            self.reads[read_type][filter] += 1;
            self.bases[read_type][filter] += u64::from(len);
        }
        let counter = match outcome {
            SpotOutcome::Pair(..) => &mut self.pairs_written,
            SpotOutcome::Unpaired(_) => &mut self.unpaired_written,
            SpotOutcome::Dropped(DropReason::NoBioReads) => &mut self.spots_no_bio_reads,
            SpotOutcome::Dropped(DropReason::TooManyBioReads) => &mut self.spots_too_many_bio_reads,
            SpotOutcome::Dropped(DropReason::Filtered) => &mut self.spots_filtered,
            SpotOutcome::Dropped(DropReason::PairUnrouted) => &mut self.pairs_unrouted,
            SpotOutcome::Dropped(DropReason::UnpairedUnrouted) => &mut self.unpaired_unrouted,
            SpotOutcome::Dropped(DropReason::TechnicalMismatch) => {
                &mut self.spots_technical_mismatch
            }
        };
        *counter += 1;
    }

    /// Add `other`'s counts to these.
    pub fn merge(&mut self, other: &SpotCounts) {
        self.spots += other.spots;
        self.spots_skipped_by_subsampling += other.spots_skipped_by_subsampling;
        self.pairs_written += other.pairs_written;
        self.unpaired_written += other.unpaired_written;
        self.spots_no_bio_reads += other.spots_no_bio_reads;
        self.spots_too_many_bio_reads += other.spots_too_many_bio_reads;
        self.spots_filtered += other.spots_filtered;
        self.spots_technical_mismatch += other.spots_technical_mismatch;
        self.pairs_unrouted += other.pairs_unrouted;
        self.unpaired_unrouted += other.unpaired_unrouted;
        for read_type in 0..2 {
            for filter in ReadFilter::ALL {
                let filter = filter.index();
                self.reads[read_type][filter] += other.reads[read_type][filter];
                self.bases[read_type][filter] += other.bases[read_type][filter];
            }
        }
    }

    /// Spots written to any output.
    pub fn spots_written(&self) -> u64 {
        self.pairs_written + self.unpaired_written
    }

    /// Bases of all reads, technical ones included; comparable to `STATS/TABLE/BASE_COUNT`.
    pub fn total_bases(&self) -> u64 {
        self.bases.iter().flatten().sum()
    }

    /// Bases of biological reads; comparable to `STATS/TABLE/BIO_BASE_COUNT`.
    pub fn biological_bases(&self) -> u64 {
        self.bases[1].iter().sum()
    }

    /// Why nothing may have been written, from the reasons spots were dropped.
    pub fn nothing_written_hint(&self) -> String {
        let hints = [
            (self.unpaired_unrouted, "had one biological read; did you mean --unpaired?"),
            (self.pairs_unrouted, "had two biological reads; did you mean --r1/--r2?"),
            (self.spots_too_many_bio_reads, "had more than two biological reads"),
            (self.spots_no_bio_reads, "had no biological reads"),
            (self.spots_filtered, "failed the filters"),
            (self.spots_technical_mismatch, "had a different number of technical reads"),
            (self.spots_skipped_by_subsampling, "were skipped by subsampling"),
        ];
        let hints: Vec<String> = hints
            .into_iter()
            .filter(|&(count, _)| count > 0)
            .map(|(count, what)| format!("{} spots {what}", format_count(count)))
            .collect();
        if hints.is_empty() { "No spots were converted.".to_owned() } else { hints.join("; ") }
    }
}

/// One row of the `--metrics` TSV.
#[derive(Debug, Clone, Serialize)]
pub struct FastqMetrics {
    pub accession: String,
    /// Reads table within a database; empty for a flat table.
    pub table: Option<String>,
    pub platform: Option<&'static str>,
    pub quality_source: QualitySource,
    pub original_names: bool,
    /// Spots in the archive.
    pub spots_total: u64,
    pub spots_converted: u64,
    pub spots_skipped_by_subsampling: u64,
    pub pairs_written: u64,
    pub unpaired_written: u64,
    pub technical_reads_written: u64,
    pub spots_no_bio_reads: u64,
    pub spots_too_many_bio_reads: u64,
    pub spots_filtered: u64,
    pub spots_technical_mismatch: u64,
    pub pairs_unrouted: u64,
    pub unpaired_unrouted: u64,
    pub biological_pass_reads: u64,
    pub biological_pass_bases: u64,
    pub biological_reject_reads: u64,
    pub biological_reject_bases: u64,
    pub biological_criteria_reads: u64,
    pub biological_criteria_bases: u64,
    pub biological_redacted_reads: u64,
    pub biological_redacted_bases: u64,
    pub technical_pass_reads: u64,
    pub technical_pass_bases: u64,
    pub technical_reject_reads: u64,
    pub technical_reject_bases: u64,
    pub technical_criteria_reads: u64,
    pub technical_criteria_bases: u64,
    pub technical_redacted_reads: u64,
    pub technical_redacted_bases: u64,
    /// Bases converted, technical ones included.
    pub bases_converted: u64,
    pub biological_bases_converted: u64,
    /// `STATS/TABLE/SPOT_COUNT`, if stored.
    pub stored_spots: Option<u64>,
    /// `STATS/TABLE/BASE_COUNT`, if stored.
    pub stored_bases: Option<u64>,
    /// `STATS/TABLE/BIO_BASE_COUNT`, if stored.
    pub stored_biological_bases: Option<u64>,
    pub integrity: Integrity,
}

impl FastqMetrics {
    /// Metrics for a conversion of `archive` that counted `counts`. Totals are checked
    /// against the archive's only when `whole_run`, i.e. every spot was converted.
    pub fn new(
        accession: &str,
        archive: &Archive,
        counts: &SpotCounts,
        technical_reads_written: u64,
        whole_run: bool,
    ) -> Self {
        let [technical, biological] = counts.reads;
        let [technical_bases, biological_bases] = counts.bases;
        let totals = archive.totals;
        let mut metrics = Self {
            accession: accession.to_owned(),
            table: archive.table.clone(),
            platform: archive.platform,
            quality_source: archive.quality_source,
            original_names: archive.has_names,
            spots_total: archive.spot_count,
            spots_converted: counts.spots,
            spots_skipped_by_subsampling: counts.spots_skipped_by_subsampling,
            pairs_written: counts.pairs_written,
            unpaired_written: counts.unpaired_written,
            technical_reads_written,
            spots_no_bio_reads: counts.spots_no_bio_reads,
            spots_too_many_bio_reads: counts.spots_too_many_bio_reads,
            spots_filtered: counts.spots_filtered,
            spots_technical_mismatch: counts.spots_technical_mismatch,
            pairs_unrouted: counts.pairs_unrouted,
            unpaired_unrouted: counts.unpaired_unrouted,
            biological_pass_reads: biological[0],
            biological_pass_bases: biological_bases[0],
            biological_reject_reads: biological[1],
            biological_reject_bases: biological_bases[1],
            biological_criteria_reads: biological[2],
            biological_criteria_bases: biological_bases[2],
            biological_redacted_reads: biological[3],
            biological_redacted_bases: biological_bases[3],
            technical_pass_reads: technical[0],
            technical_pass_bases: technical_bases[0],
            technical_reject_reads: technical[1],
            technical_reject_bases: technical_bases[1],
            technical_criteria_reads: technical[2],
            technical_criteria_bases: technical_bases[2],
            technical_redacted_reads: technical[3],
            technical_redacted_bases: technical_bases[3],
            bases_converted: counts.total_bases(),
            biological_bases_converted: counts.biological_bases(),
            stored_spots: totals.spots,
            stored_bases: totals.bases,
            stored_biological_bases: totals.biological_bases,
            integrity: Integrity::NotChecked,
        };
        let stores_totals =
            totals.spots.is_some() || totals.bases.is_some() || totals.biological_bases.is_some();
        if whole_run && stores_totals {
            metrics.integrity = if metrics.integrity_mismatches().is_empty() {
                Integrity::Pass
            } else {
                Integrity::Fail
            };
        }
        metrics
    }

    /// Print a summary of a conversion to stderr.
    pub fn report(&self) {
        let accession = &self.accession;
        eprintln!(
            "[fastq] {accession}: converted {} of {} spots; wrote {} pairs, {} unpaired reads and \
             {} technical reads",
            format_count(self.spots_converted),
            format_count(self.spots_total),
            format_count(self.pairs_written),
            format_count(self.unpaired_written),
            format_count(self.technical_reads_written),
        );
        if self.spots_skipped_by_subsampling > 0 {
            eprintln!(
                "[fastq] {accession}: subsampling skipped {} spots",
                format_count(self.spots_skipped_by_subsampling)
            );
        }
        let dropped = [
            (self.spots_no_bio_reads, "with no biological reads"),
            (self.spots_too_many_bio_reads, "with more than two biological reads"),
            (self.spots_filtered, "failing the filters"),
            (self.spots_technical_mismatch, "with a different number of technical reads"),
            (self.pairs_unrouted, "pairs with no paired output"),
            (self.unpaired_unrouted, "single reads with no unpaired output"),
        ];
        let dropped: Vec<String> = dropped
            .into_iter()
            .filter(|&(count, _)| count > 0)
            .map(|(count, what)| format!("{} {what}", format_count(count)))
            .collect();
        if !dropped.is_empty() {
            eprintln!("[fastq] {accession}: dropped {}", dropped.join(", "));
        }
        let integrity = match self.integrity {
            Integrity::Pass => "totals match the archive's",
            Integrity::Fail => "totals DO NOT match the archive's",
            Integrity::NotChecked => "totals not checked",
        };
        eprintln!("[fastq] {accession}: {integrity}");
    }

    /// A description of how the totals failed to match, if they were checked and didn't.
    pub fn integrity_mismatch(&self) -> Option<String> {
        (self.integrity == Integrity::Fail).then(|| self.integrity_mismatches().join("; "))
    }

    /// Each converted total that differs from the archive's.
    fn integrity_mismatches(&self) -> Vec<String> {
        let checks = [
            ("spots", self.spots_converted, Some(self.spots_total), "rows"),
            ("spots", self.spots_converted, self.stored_spots, "SPOT_COUNT"),
            ("bases", self.bases_converted, self.stored_bases, "BASE_COUNT"),
            (
                "biological bases",
                self.biological_bases_converted,
                self.stored_biological_bases,
                "BIO_BASE_COUNT",
            ),
        ];
        checks
            .into_iter()
            .filter_map(|(what, converted, stored, name)| match stored {
                Some(stored) if stored != converted => Some(format!(
                    "converted {} {what} but the archive has {} ({name})",
                    format_count(converted),
                    format_count(stored)
                )),
                _ => None,
            })
            .collect()
    }
}

/// Whether a conversion's totals matched the archive's stored ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Integrity {
    #[serde(rename = "pass")]
    Pass,
    #[serde(rename = "fail")]
    Fail,
    /// Not checked: only part of the archive was converted, or it stores no totals.
    #[serde(rename = "n/a")]
    NotChecked,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::spot::SelectedRead;
    use crate::fastq::spot::tests::TestSpot;
    use crate::record::READ_TYPE_BIOLOGICAL as B;

    const T: u8 = 0;
    const PASS: u8 = 0;
    const REJECT: u8 = 1;
    const REDACTED: u8 = 3;

    fn unpaired() -> SpotOutcome {
        SpotOutcome::Unpaired(SelectedRead { slot: 0, number: 1 })
    }

    fn counts_of(reads: &[(&str, u8, u8)], outcome: SpotOutcome) -> SpotCounts {
        let mut counts = SpotCounts::default();
        counts.record(&TestSpot::new(reads).view(), outcome);
        counts
    }

    #[test]
    fn reads_and_bases_are_counted_by_type_and_filter() {
        let counts =
            counts_of(&[("ACGT", B, PASS), ("GG", T, REJECT), ("CCC", B, REDACTED)], unpaired());
        assert_eq!(counts.reads[1][0], 1);
        assert_eq!(counts.bases[1][0], 4);
        assert_eq!(counts.reads[0][1], 1);
        assert_eq!(counts.bases[0][1], 2);
        assert_eq!(counts.reads[1][3], 1);
        assert_eq!(counts.bases[1][3], 3);
    }

    #[test]
    fn empty_reads_are_not_counted() {
        let counts = counts_of(&[("", B, PASS), ("ACGT", B, PASS)], unpaired());
        assert_eq!(counts.reads[1][0], 1);
    }

    #[test]
    fn reads_of_dropped_spots_are_still_counted() {
        let outcome = SpotOutcome::Dropped(DropReason::Filtered);
        let counts = counts_of(&[("ACGT", B, PASS)], outcome);
        assert_eq!(counts.spots_filtered, 1);
        assert_eq!(counts.total_bases(), 4);
    }

    #[test]
    fn total_bases_include_technical_reads() {
        let counts = counts_of(&[("ACGT", B, PASS), ("GG", T, PASS)], unpaired());
        assert_eq!(counts.total_bases(), 6);
        assert_eq!(counts.biological_bases(), 4);
    }

    #[test]
    fn each_outcome_has_its_own_counter() {
        let mut counts = SpotCounts::default();
        let spot = TestSpot::new(&[("A", B, PASS)]);
        let read = SelectedRead { slot: 0, number: 1 };
        for outcome in [
            SpotOutcome::Pair(read, read),
            SpotOutcome::Unpaired(read),
            SpotOutcome::Dropped(DropReason::NoBioReads),
            SpotOutcome::Dropped(DropReason::TooManyBioReads),
            SpotOutcome::Dropped(DropReason::Filtered),
            SpotOutcome::Dropped(DropReason::PairUnrouted),
            SpotOutcome::Dropped(DropReason::UnpairedUnrouted),
            SpotOutcome::Dropped(DropReason::TechnicalMismatch),
        ] {
            counts.record(&spot.view(), outcome);
        }
        let per_outcome = [
            counts.pairs_written,
            counts.unpaired_written,
            counts.spots_no_bio_reads,
            counts.spots_too_many_bio_reads,
            counts.spots_filtered,
            counts.pairs_unrouted,
            counts.unpaired_unrouted,
            counts.spots_technical_mismatch,
        ];
        assert_eq!(per_outcome, [1; 8]);
        assert_eq!(counts.spots, 8);
        assert_eq!(counts.spots_written(), 2);
    }

    #[test]
    fn merged_counts_are_the_sum_of_both() {
        let mut first = counts_of(&[("ACGT", B, PASS)], unpaired());
        let second = counts_of(&[("GG", T, REJECT)], SpotOutcome::Dropped(DropReason::NoBioReads));
        first.merge(&second);
        assert_eq!(first.spots, 2);
        assert_eq!(first.unpaired_written, 1);
        assert_eq!(first.spots_no_bio_reads, 1);
        assert_eq!(first.total_bases(), 6);
        assert_eq!(first.reads[0][1], 1);
    }

    #[test]
    fn merged_spots_skipped_by_subsampling_are_the_sum_of_both() {
        let mut first = SpotCounts { spots_skipped_by_subsampling: 3, ..SpotCounts::default() };
        first.merge(&SpotCounts { spots_skipped_by_subsampling: 4, ..SpotCounts::default() });
        assert_eq!(first.spots_skipped_by_subsampling, 7);
    }

    #[test]
    fn nothing_written_hint_counts_spots_skipped_by_subsampling() {
        let counts = SpotCounts { spots_skipped_by_subsampling: 5, ..SpotCounts::default() };
        assert_eq!(counts.nothing_written_hint(), "5 spots were skipped by subsampling");
    }
}
