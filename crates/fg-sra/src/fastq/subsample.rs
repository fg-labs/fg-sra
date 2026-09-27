//! Choosing a random subset of spots, the same at any thread count.

/// `SplitMix64`'s increment: 2^64 divided by the golden ratio, rounded to odd.
const GOLDEN_GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;

/// Keeps each spot with probability `fraction`, chosen by a seed.
///
/// A spot's choice depends only on the seed and its spot id, never on which thread converts
/// it or when, so output is the same at any thread count and all of a spot's reads are kept or
/// skipped together. Choices are independent, so any range of spots keeps about `fraction` of
/// its spots, and with one seed a smaller fraction keeps a subset of a larger one's spots.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subsampler {
    fraction: f64,
    /// Where the seed's `SplitMix64` sequence starts: the seed hashed, so that nearby seeds
    /// start far apart.
    start: u64,
}

impl Subsampler {
    /// Keep `fraction` of spots, which [`parse_fraction`] accepts, chosen by `seed`.
    #[must_use]
    pub fn new(fraction: f64, seed: u64) -> Self {
        Self { fraction, start: mix(seed) }
    }

    /// Whether every spot is kept.
    #[must_use]
    pub fn keeps_all(&self) -> bool {
        self.fraction >= 1.0
    }

    /// Whether spot `id` is kept.
    #[must_use]
    pub fn keeps(&self, id: i64) -> bool {
        // Output `id` of a `SplitMix64` generator seeded with `start`, whose top 53 bits make
        // a value uniform in [0, 1): always below a fraction of 1.
        let hash = mix(self.start.wrapping_add((id as u64).wrapping_mul(GOLDEN_GAMMA)));
        let uniform = (hash >> 11) as f64 / (1u64 << 53) as f64;
        uniform < self.fraction
    }
}

/// `SplitMix64`'s output function, which spreads every bit of `z` over every bit of the result.
fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Parse a fraction of spots to keep: more than 0, and at most 1.
pub fn parse_fraction(text: &str) -> Result<f64, String> {
    let fraction: f64 = text.parse().map_err(|_| "not a number".to_owned())?;
    if fraction > 0.0 && fraction <= 1.0 {
        Ok(fraction)
    } else {
        Err("must be more than 0 and at most 1".to_owned())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// How many of `ids` `subsampler` keeps.
    fn kept(subsampler: &Subsampler, ids: std::ops::RangeInclusive<i64>) -> usize {
        ids.filter(|&id| subsampler.keeps(id)).count()
    }

    /// Assert `count` of `trials` is within five standard deviations of `fraction` of them.
    fn assert_about(count: usize, trials: usize, fraction: f64) {
        let expected = trials as f64 * fraction;
        let tolerance = 5.0 * (expected * (1.0 - fraction)).sqrt();
        let count = count as f64;
        assert!((count - expected).abs() <= tolerance, "{count} of {trials}, expected {expected}");
    }

    #[test]
    fn a_fraction_of_one_keeps_every_spot() {
        let subsampler = Subsampler::new(1.0, 42);
        assert!(subsampler.keeps_all());
        assert_eq!(kept(&subsampler, 1..=100_000), 100_000);
    }

    #[test]
    fn a_fraction_below_one_does_not_keep_all() {
        assert!(!Subsampler::new(0.999, 42).keeps_all());
    }

    #[test]
    fn about_the_fraction_of_spots_is_kept() {
        assert_about(kept(&Subsampler::new(0.1, 42), 1..=1_000_000), 1_000_000, 0.1);
    }

    #[test]
    fn about_the_fraction_of_a_small_range_is_kept() {
        assert_about(kept(&Subsampler::new(0.25, 42), 5_001..=6_000), 1_000, 0.25);
    }

    #[test]
    fn about_the_fraction_of_a_range_of_large_spot_ids_is_kept() {
        let first = 4_000_000_000;
        assert_about(kept(&Subsampler::new(0.5, 7), first..=first + 99_999), 100_000, 0.5);
    }

    #[test]
    fn neighbouring_spots_are_kept_independently() {
        let subsampler = Subsampler::new(0.3, 42);
        let both = (1..=200_000).filter(|&id| subsampler.keeps(id) && subsampler.keeps(id + 1));
        assert_about(both.count(), 200_000, 0.3 * 0.3);
    }

    #[test]
    fn the_same_seed_keeps_the_same_spots() {
        let (a, b) = (Subsampler::new(0.2, 42), Subsampler::new(0.2, 42));
        assert!((1..=100_000).all(|id| a.keeps(id) == b.keeps(id)));
    }

    #[test]
    fn neighbouring_seeds_keep_spots_independently() {
        let (a, b) = (Subsampler::new(0.5, 42), Subsampler::new(0.5, 43));
        let both = (1..=100_000).filter(|&id| a.keeps(id) && b.keeps(id));
        assert_about(both.count(), 100_000, 0.25);
    }

    #[test]
    fn a_smaller_fraction_keeps_a_subset_of_a_larger_ones_spots() {
        let (smaller, larger) = (Subsampler::new(0.1, 42), Subsampler::new(0.2, 42));
        assert!((1..=100_000).all(|id| !smaller.keeps(id) || larger.keeps(id)));
    }

    #[test]
    fn fractions_up_to_one_are_accepted() {
        assert_eq!(parse_fraction("1"), Ok(1.0));
        assert_eq!(parse_fraction("0.25"), Ok(0.25));
        assert_eq!(parse_fraction("1e-6"), Ok(1e-6));
    }

    #[test]
    fn fractions_of_zero_or_less_are_refused() {
        assert!(parse_fraction("0").is_err());
        assert!(parse_fraction("-0.5").is_err());
    }

    #[test]
    fn fractions_over_one_are_refused() {
        assert!(parse_fraction("1.5").is_err());
    }

    #[test]
    fn a_fraction_that_is_not_a_number_is_refused() {
        assert!(parse_fraction("NaN").is_err());
        assert!(parse_fraction("half").is_err());
    }
}
