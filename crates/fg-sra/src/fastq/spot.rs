//! One spot's reads, and which FASTQ output each spot goes to.

use anyhow::{Result, bail};

use crate::record::READ_TYPE_BIOLOGICAL;

/// Decides each spot's outcome from the filters and the outputs requested.
#[derive(Debug, Clone)]
pub struct SpotRouter {
    /// Biological reads shorter than this fail the spot; 0 means no minimum.
    pub min_read_len: u32,
    /// Which `READ_FILTER` values a biological read may have, indexed by [`ReadFilter::index`].
    pub kept_filters: [bool; 4],
    /// Whether pairs have an output.
    pub write_pairs: bool,
    /// Whether single reads have an output.
    pub write_unpaired: bool,
    /// When technical reads are output, the number each written spot must have.
    pub technical_reads: Option<usize>,
}

impl SpotRouter {
    /// Decide where `spot` goes. The spot must have passed [`Spot::validate`].
    ///
    /// Checks run in order: the number of non-empty biological reads, then the filters,
    /// then whether an output was requested for the spot, then its technical reads.
    pub fn route(&self, spot: &Spot<'_>) -> SpotOutcome {
        let mut mates = [SelectedRead { slot: 0, number: 0 }; 2];
        let mut num_mates = 0;
        let mut bio_number = 0;
        for slot in 0..spot.read_count() {
            if !spot.is_biological(slot) {
                continue;
            }
            bio_number += 1;
            if spot.read_len(slot) == 0 {
                continue;
            }
            if num_mates == mates.len() {
                return SpotOutcome::Dropped(DropReason::TooManyBioReads);
            }
            mates[num_mates] = SelectedRead { slot, number: bio_number };
            num_mates += 1;
        }
        let mates = &mates[..num_mates];
        if mates.is_empty() {
            return SpotOutcome::Dropped(DropReason::NoBioReads);
        }
        if !mates.iter().all(|&read| self.passes_filters(spot, read)) {
            return SpotOutcome::Dropped(DropReason::Filtered);
        }

        let outcome = if let [first, second] = *mates {
            if !self.write_pairs {
                return SpotOutcome::Dropped(DropReason::PairUnrouted);
            }
            SpotOutcome::Pair(first, second)
        } else {
            if !self.write_unpaired {
                return SpotOutcome::Dropped(DropReason::UnpairedUnrouted);
            }
            SpotOutcome::Unpaired(mates[0])
        };
        match self.technical_reads {
            Some(expected) if spot.technical_reads().count() != expected => {
                SpotOutcome::Dropped(DropReason::TechnicalMismatch)
            }
            _ => outcome,
        }
    }

    /// Whether biological `read` passes `--min-read-len` and `--read-filter`.
    fn passes_filters(&self, spot: &Spot<'_>, read: SelectedRead) -> bool {
        spot.read_len(read.slot) >= self.min_read_len
            && self.kept_filters[spot.filter(read.slot).index()]
    }
}

/// What happens to a spot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpotOutcome {
    /// Written to the paired outputs, as first and second mate.
    Pair(SelectedRead, SelectedRead),
    /// Written to the unpaired output.
    Unpaired(SelectedRead),
    /// Not written.
    Dropped(DropReason),
}

/// Why a spot wasn't written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropReason {
    /// No non-empty biological read.
    NoBioReads,
    /// More than two non-empty biological reads.
    TooManyBioReads,
    /// A biological read failed `--min-read-len` or `--read-filter`.
    Filtered,
    /// A pair, with no paired output requested.
    PairUnrouted,
    /// A single read, with no unpaired output requested.
    UnpairedUnrouted,
    /// Technical reads are output, and this spot has a different number of them.
    TechnicalMismatch,
}

/// One spot's data, borrowed from the buffers it was read into.
///
/// `bases` holds every read of the spot end to end, and `qualities` their phred values, or
/// nothing when qualities weren't read (FASTA output). Read `i` is `read_lens[i]` bases from
/// `read_starts[i]`. Call [`Spot::validate`] before using a spot from an archive.
#[derive(Debug, Clone, Copy)]
pub struct Spot<'a> {
    /// Spot (row) id, 1-based.
    pub id: i64,
    /// Original spot name; empty when the archive kept none.
    pub name: &'a [u8],
    /// Spot group (e.g. read group or barcode); empty when there is none.
    pub group: &'a [u8],
    /// Bases of every read, end to end.
    pub bases: &'a [u8],
    /// Phred value of each base; empty when qualities weren't read.
    pub qualities: &'a [u8],
    /// Offset of each read within `bases`.
    pub read_starts: &'a [i32],
    /// Length of each read.
    pub read_lens: &'a [u32],
    /// `READ_TYPE` bits of each read; bit 0 set means biological.
    pub read_types: &'a [u8],
    /// Raw `READ_FILTER` value of each read.
    pub read_filters: &'a [u8],
}

impl<'a> Spot<'a> {
    /// Check that the spot's columns agree with each other, so later accessors can't go out
    /// of bounds and output can't silently be wrong.
    ///
    /// `with_qualities` says whether qualities were read and so must match the bases.
    pub fn validate(&self, with_qualities: bool) -> Result<()> {
        let num_reads = self.read_lens.len();
        if self.read_starts.len() != num_reads
            || self.read_types.len() != num_reads
            || self.read_filters.len() != num_reads
        {
            bail!(
                "spot {}: READ_LEN has {} values but READ_START has {}, READ_TYPE {} and \
                 READ_FILTER {}",
                self.id,
                num_reads,
                self.read_starts.len(),
                self.read_types.len(),
                self.read_filters.len()
            );
        }
        let total_len: u64 = self.read_lens.iter().map(|&len| u64::from(len)).sum();
        if total_len != self.bases.len() as u64 {
            bail!(
                "spot {}: READ_LEN sums to {total_len} but READ has {} bases",
                self.id,
                self.bases.len()
            );
        }
        if with_qualities && self.qualities.len() != self.bases.len() {
            bail!(
                "spot {}: QUALITY has {} values but READ has {} bases",
                self.id,
                self.qualities.len(),
                self.bases.len()
            );
        }
        for (i, (&start, &len)) in self.read_starts.iter().zip(self.read_lens).enumerate() {
            let end = i64::from(start) + i64::from(len);
            if start < 0 || end > self.bases.len() as i64 {
                bail!(
                    "spot {}: read {} spans {start}..{end}, outside the spot's {} bases",
                    self.id,
                    i + 1,
                    self.bases.len()
                );
            }
        }
        if let Some(&raw) =
            self.read_filters.iter().find(|&&raw| ReadFilter::from_raw(raw).is_none())
        {
            bail!("spot {}: unknown READ_FILTER value {raw}", self.id);
        }
        Ok(())
    }

    /// Number of read slots in the spot, empty ones included.
    pub fn read_count(&self) -> usize {
        self.read_lens.len()
    }

    /// Whether read `slot` is biological rather than technical.
    pub fn is_biological(&self, slot: usize) -> bool {
        self.read_types[slot] & READ_TYPE_BIOLOGICAL != 0
    }

    /// Length of read `slot`.
    pub fn read_len(&self, slot: usize) -> u32 {
        self.read_lens[slot]
    }

    /// `READ_FILTER` of read `slot`.
    pub fn filter(&self, slot: usize) -> ReadFilter {
        ReadFilter::from_raw(self.read_filters[slot]).expect("READ_FILTER checked by validate")
    }

    /// Bases of read `slot`.
    pub fn read_bases(&self, slot: usize) -> &'a [u8] {
        let start = self.read_starts[slot] as usize;
        &self.bases[start..start + self.read_lens[slot] as usize]
    }

    /// Phred qualities of read `slot`, or an empty slice when qualities weren't read.
    pub fn read_qualities(&self, slot: usize) -> &'a [u8] {
        if self.qualities.is_empty() {
            return &[];
        }
        let start = self.read_starts[slot] as usize;
        &self.qualities[start..start + self.read_lens[slot] as usize]
    }

    /// The spot's non-empty technical reads, numbered 1, 2, … in slot order.
    pub fn technical_reads(&self) -> impl Iterator<Item = SelectedRead> + '_ {
        (0..self.read_count())
            .filter(|&slot| !self.is_biological(slot) && self.read_len(slot) > 0)
            .zip(1..)
            .map(|(slot, number)| SelectedRead { slot, number })
    }
}

/// A read chosen for output: its slot in the spot, and its number within its read type.
///
/// Biological reads are numbered by biological slot, empty ones included, so the survivor of
/// a `[0, 150]` pair is biological read 2. Technical reads are numbered among the non-empty
/// technical slots.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SelectedRead {
    /// Index of the read within the spot.
    pub slot: usize,
    /// Number of the read within its type, from 1.
    pub number: u32,
}

/// A read's `READ_FILTER` value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum ReadFilter {
    /// Passed the instrument's or submitter's filters.
    Pass,
    /// Rejected: e.g. failed Illumina chastity filtering, SRA Lite's quality test, or marked
    /// QC-fail in a submitted BAM.
    Reject,
    /// Filtered by other criteria: e.g. marked duplicate in a submitted BAM.
    Criteria,
    /// Masked by SRA (all `N`), e.g. by human-read scrubbing.
    Redacted,
}

impl ReadFilter {
    /// Every value, in `READ_FILTER` order.
    pub const ALL: [ReadFilter; 4] = [Self::Pass, Self::Reject, Self::Criteria, Self::Redacted];

    /// Interpret a raw `READ_FILTER` value; `None` for one outside the schema's range.
    pub fn from_raw(raw: u8) -> Option<Self> {
        Self::ALL.get(usize::from(raw)).copied()
    }

    /// The raw `READ_FILTER` value, usable as an index into per-filter tallies.
    pub fn index(self) -> usize {
        self as usize
    }

    /// The value's name, as `--read-filter` takes it.
    pub fn name(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Reject => "reject",
            Self::Criteria => "criteria",
            Self::Redacted => "redacted",
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    const T: u8 = 0;
    const B: u8 = READ_TYPE_BIOLOGICAL;
    /// Biological, reverse-oriented (`BIOLOGICAL | REVERSE`), as bam-load writes.
    const B_REV: u8 = READ_TYPE_BIOLOGICAL | 4;
    /// Technical, forward-oriented (`TECHNICAL | FORWARD`).
    const T_FWD: u8 = 2;

    const PASS: u8 = 0;
    const REJECT: u8 = 1;

    /// Owned data for a test spot; [`TestSpot::view`] borrows it as a [`Spot`].
    pub(crate) struct TestSpot {
        pub id: i64,
        pub name: Vec<u8>,
        pub spot_group: Vec<u8>,
        pub bases: Vec<u8>,
        pub qualities: Vec<u8>,
        pub read_starts: Vec<i32>,
        pub read_lens: Vec<u32>,
        pub read_types: Vec<u8>,
        pub read_filters: Vec<u8>,
    }

    impl TestSpot {
        /// A spot with one read per `(bases, read type, read filter)`, stored end to end, with
        /// qualities counting up from 10 across the spot.
        pub(crate) fn new(reads: &[(&str, u8, u8)]) -> Self {
            let mut spot = Self {
                id: 1,
                name: Vec::new(),
                spot_group: Vec::new(),
                bases: Vec::new(),
                qualities: Vec::new(),
                read_starts: Vec::new(),
                read_lens: Vec::new(),
                read_types: Vec::new(),
                read_filters: Vec::new(),
            };
            for &(bases, read_type, read_filter) in reads {
                spot.read_starts.push(spot.bases.len() as i32);
                spot.read_lens.push(bases.len() as u32);
                spot.read_types.push(read_type);
                spot.read_filters.push(read_filter);
                spot.bases.extend_from_slice(bases.as_bytes());
            }
            spot.qualities = (0..spot.bases.len()).map(|i| 10 + i as u8).collect();
            spot
        }

        pub(crate) fn view(&self) -> Spot<'_> {
            Spot {
                id: self.id,
                name: &self.name,
                group: &self.spot_group,
                bases: &self.bases,
                qualities: &self.qualities,
                read_starts: &self.read_starts,
                read_lens: &self.read_lens,
                read_types: &self.read_types,
                read_filters: &self.read_filters,
            }
        }
    }

    /// A router writing pairs and unpaired reads, with no filters or technical output.
    fn router() -> SpotRouter {
        SpotRouter {
            min_read_len: 0,
            kept_filters: [true; 4],
            write_pairs: true,
            write_unpaired: true,
            technical_reads: None,
        }
    }

    fn read(slot: usize, number: u32) -> SelectedRead {
        SelectedRead { slot, number }
    }

    fn route(router: &SpotRouter, reads: &[(&str, u8, u8)]) -> SpotOutcome {
        let spot = TestSpot::new(reads);
        spot.view().validate(true).unwrap();
        router.route(&spot.view())
    }

    // ── validate ────────────────────────────────────────────────────────

    #[test]
    fn well_formed_spot_is_valid() {
        let spot = TestSpot::new(&[("ACGT", B, PASS), ("GG", T, PASS)]);
        assert!(spot.view().validate(true).is_ok());
    }

    #[test]
    fn read_lengths_not_summing_to_the_bases_are_invalid() {
        let mut spot = TestSpot::new(&[("ACGT", B, PASS)]);
        spot.read_lens[0] = 5;
        assert!(spot.view().validate(false).is_err());
    }

    #[test]
    fn qualities_shorter_than_the_bases_are_invalid() {
        let mut spot = TestSpot::new(&[("ACGT", B, PASS)]);
        spot.qualities.pop();
        assert!(spot.view().validate(true).is_err());
    }

    #[test]
    fn missing_qualities_are_valid_when_not_read() {
        let mut spot = TestSpot::new(&[("ACGT", B, PASS)]);
        spot.qualities.clear();
        assert!(spot.view().validate(false).is_ok());
    }

    #[test]
    fn read_start_outside_the_spot_is_invalid() {
        let mut spot = TestSpot::new(&[("AC", B, PASS), ("GT", B, PASS)]);
        spot.read_starts[1] = 3;
        assert!(spot.view().validate(true).is_err());
    }

    #[test]
    fn negative_read_start_is_invalid() {
        let mut spot = TestSpot::new(&[("AC", B, PASS)]);
        spot.read_starts[0] = -1;
        assert!(spot.view().validate(true).is_err());
    }

    #[test]
    fn mismatched_per_read_column_lengths_are_invalid() {
        let mut spot = TestSpot::new(&[("AC", B, PASS), ("GT", B, PASS)]);
        spot.read_types.pop();
        assert!(spot.view().validate(true).is_err());
    }

    #[test]
    fn unknown_read_filter_value_is_invalid() {
        let spot = TestSpot::new(&[("AC", B, 7)]);
        let err = spot.view().validate(true).unwrap_err();
        assert!(err.to_string().contains("READ_FILTER value 7"), "{err}");
    }

    // ── read access ─────────────────────────────────────────────────────

    #[test]
    fn read_bases_and_qualities_come_from_the_read_start() {
        let spot = TestSpot::new(&[("AC", T, PASS), ("GTT", B, PASS)]);
        let view = spot.view();
        assert_eq!(view.read_bases(1), b"GTT");
        assert_eq!(view.read_qualities(1), &[12, 13, 14]);
    }

    #[test]
    fn read_qualities_are_empty_when_not_read() {
        let mut spot = TestSpot::new(&[("AC", B, PASS)]);
        spot.qualities.clear();
        assert!(spot.view().read_qualities(0).is_empty());
    }

    #[test]
    fn technical_reads_are_numbered_among_non_empty_technical_slots() {
        let spot =
            TestSpot::new(&[("", T, PASS), ("AC", B, PASS), ("GG", T, PASS), ("T", T_FWD, PASS)]);
        let technical: Vec<_> = spot.view().technical_reads().collect();
        assert_eq!(technical, vec![read(2, 1), read(3, 2)]);
    }

    // ── routing ─────────────────────────────────────────────────────────

    #[test]
    fn two_biological_reads_are_a_pair() {
        let outcome = route(&router(), &[("ACGT", B, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Pair(read(0, 1), read(1, 2)));
    }

    #[test]
    fn one_biological_read_is_unpaired() {
        let outcome = route(&router(), &[("ACGT", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(0, 1)));
    }

    #[test]
    fn orientation_bits_do_not_hide_a_biological_read() {
        let outcome = route(&router(), &[("ACGT", B_REV, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Pair(read(0, 1), read(1, 2)));
    }

    #[test]
    fn forward_technical_read_is_not_biological() {
        let outcome = route(&router(), &[("ACGT", T_FWD, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(1, 1)));
    }

    #[test]
    fn technical_reads_are_skipped_when_numbering_mates() {
        let outcome = route(
            &router(),
            &[("ACGT", T, PASS), ("AAA", B, PASS), ("CC", T, PASS), ("GGG", B, PASS)],
        );
        assert_eq!(outcome, SpotOutcome::Pair(read(1, 1), read(3, 2)));
    }

    #[test]
    fn empty_first_mate_makes_the_second_an_orphan_numbered_two() {
        let outcome = route(&router(), &[("", B, PASS), ("ACGT", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(1, 2)));
    }

    #[test]
    fn empty_second_mate_makes_the_first_an_orphan_numbered_one() {
        let outcome = route(&router(), &[("ACGT", B, PASS), ("", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(0, 1)));
    }

    #[test]
    fn spot_with_only_technical_reads_has_no_bio_reads() {
        let outcome = route(&router(), &[("ACGT", T, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::NoBioReads));
    }

    #[test]
    fn spot_with_only_empty_reads_has_no_bio_reads() {
        let outcome = route(&router(), &[("", B, PASS), ("", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::NoBioReads));
    }

    #[test]
    fn three_biological_reads_are_too_many() {
        let outcome = route(&router(), &[("A", B, PASS), ("C", B, PASS), ("G", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::TooManyBioReads));
    }

    #[test]
    fn empty_third_biological_read_still_leaves_a_pair() {
        let outcome = route(&router(), &[("A", B, PASS), ("C", B, PASS), ("", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Pair(read(0, 1), read(1, 2)));
    }

    // ── filters ─────────────────────────────────────────────────────────

    #[test]
    fn short_mate_fails_the_whole_spot() {
        let router = SpotRouter { min_read_len: 4, ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS), ("TTG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::Filtered));
    }

    #[test]
    fn mates_at_the_minimum_length_pass() {
        let router = SpotRouter { min_read_len: 4, ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Pair(read(0, 1), read(1, 2)));
    }

    #[test]
    fn short_technical_read_does_not_fail_the_spot() {
        let router = SpotRouter { min_read_len: 4, ..router() };
        let outcome = route(&router, &[("AC", T, PASS), ("ACGT", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(1, 1)));
    }

    #[test]
    fn empty_mate_is_absent_rather_than_too_short() {
        let router = SpotRouter { min_read_len: 4, ..router() };
        let outcome = route(&router, &[("", B, PASS), ("ACGT", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(1, 2)));
    }

    #[test]
    fn rejected_mate_fails_the_whole_spot_when_reject_is_not_kept() {
        let router = SpotRouter { kept_filters: [true, false, true, true], ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS), ("TTGG", B, REJECT)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::Filtered));
    }

    #[test]
    fn rejected_mate_is_kept_by_default() {
        let outcome = route(&router(), &[("ACGT", B, PASS), ("TTGG", B, REJECT)]);
        assert_eq!(outcome, SpotOutcome::Pair(read(0, 1), read(1, 2)));
    }

    #[test]
    fn rejected_technical_read_does_not_fail_the_spot() {
        let router = SpotRouter { kept_filters: [true, false, true, true], ..router() };
        let outcome = route(&router, &[("ACGT", T, REJECT), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(1, 1)));
    }

    // ── outputs requested ───────────────────────────────────────────────

    #[test]
    fn pair_without_paired_output_is_unrouted() {
        let router = SpotRouter { write_pairs: false, ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::PairUnrouted));
    }

    #[test]
    fn single_read_without_unpaired_output_is_unrouted() {
        let router = SpotRouter { write_unpaired: false, ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::UnpairedUnrouted));
    }

    #[test]
    fn filtered_spot_is_reported_as_filtered_even_without_its_output() {
        let router = SpotRouter { min_read_len: 10, write_pairs: false, ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::Filtered));
    }

    // ── technical reads ─────────────────────────────────────────────────

    #[test]
    fn spot_with_the_expected_technical_reads_is_written() {
        let router = SpotRouter { technical_reads: Some(2), write_pairs: false, ..router() };
        let outcome = route(&router, &[("ACGT", T, PASS), ("TTGG", B, PASS), ("CC", T, PASS)]);
        assert_eq!(outcome, SpotOutcome::Unpaired(read(1, 1)));
    }

    #[test]
    fn spot_with_fewer_technical_reads_than_expected_is_a_mismatch() {
        let router = SpotRouter { technical_reads: Some(2), write_pairs: false, ..router() };
        let outcome = route(&router, &[("ACGT", T, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::TechnicalMismatch));
    }

    #[test]
    fn empty_technical_read_counts_as_missing() {
        let router = SpotRouter { technical_reads: Some(2), write_pairs: false, ..router() };
        let outcome = route(&router, &[("ACGT", T, PASS), ("TTGG", B, PASS), ("", T, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::TechnicalMismatch));
    }

    #[test]
    fn technical_reads_are_not_checked_on_unrouted_spots() {
        let router = SpotRouter { technical_reads: Some(2), write_pairs: false, ..router() };
        let outcome = route(&router, &[("ACGT", B, PASS), ("TTGG", B, PASS)]);
        assert_eq!(outcome, SpotOutcome::Dropped(DropReason::PairUnrouted));
    }

    // ── ReadFilter ──────────────────────────────────────────────────────

    #[test]
    fn read_filter_values_match_the_schema() {
        assert_eq!(ReadFilter::from_raw(0), Some(ReadFilter::Pass));
        assert_eq!(ReadFilter::from_raw(1), Some(ReadFilter::Reject));
        assert_eq!(ReadFilter::from_raw(2), Some(ReadFilter::Criteria));
        assert_eq!(ReadFilter::from_raw(3), Some(ReadFilter::Redacted));
        assert_eq!(ReadFilter::from_raw(4), None);
    }

    #[test]
    fn read_filter_names_are_the_values_read_filter_takes() {
        use clap::ValueEnum;
        for filter in ReadFilter::ALL {
            assert_eq!(ReadFilter::from_str(filter.name(), false), Ok(filter));
        }
    }
}
