//! Quality score quantization.
//!
//! Supports the `--qual-quant` option for binning quality scores into
//! user-defined ranges (e.g. `"0:10,10:20,20:30,30:40,40:94"`).

use anyhow::{Context, Result, bail, ensure};

/// A 256-byte lookup table for quality score quantization.
/// Index by raw Phred value, returns quantized Phred value.
pub type QuantTable = [u8; 256];

/// The highest Phred quality that Phred+33 (SAM, FASTQ) can encode (`~`).
pub const MAX_PHRED: u8 = 93;

/// Parse a quantization spec like `"0:10,10:20=16,20:94"`.
///
/// Each comma-separated `low:high[=value]` range maps Phred values in `[low, high)` to
/// `value`, or without `=value` to the range's probability-averaged quality (see
/// [`average_phred`]). `high` may be up to 94, so a range can include Q93. Values outside
/// every range are left unchanged. Ranges may be given in any order, but must not be
/// empty or overlap.
pub fn parse_qual_quant(spec: &str) -> Result<QuantTable> {
    let mut ranges = Vec::new();
    for range in spec.split(',') {
        let range = range.trim();
        let (bounds, value) = match range.split_once('=') {
            Some((bounds, value)) => (bounds, Some(value)),
            None => (range, None),
        };
        let (low, high) = bounds
            .split_once(':')
            .with_context(|| format!("invalid qual-quant range {range:?}: expected low:high"))?;
        let low = parse_number(low, "low value", range)?;
        let high = high.trim();
        if high == "-" {
            bail!(
                "invalid qual-quant range {range:?}: ranges are low:high[=value], not \
                 sam-dump's value:limit with a trailing `-`; to include Q93, end a range at 94"
            );
        }
        let high = parse_number(high, "high value", range)?;
        ensure!(low < high, "invalid qual-quant range {range:?}: low must be below high");
        ensure!(
            high <= MAX_PHRED + 1,
            "invalid qual-quant range {range:?}: high must be at most {}",
            MAX_PHRED + 1
        );
        let value = match value {
            Some(value) => {
                let value = parse_number(value, "value", range)?;
                ensure!(
                    value <= MAX_PHRED,
                    "invalid qual-quant range {range:?}: value must be at most {MAX_PHRED}"
                );
                value
            }
            None => average_phred(low, high),
        };
        ranges.push((low, high, value));
    }

    ranges.sort_unstable();
    for pair in ranges.windows(2) {
        let ((low, high, _), (next_low, next_high, _)) = (pair[0], pair[1]);
        ensure!(
            high <= next_low,
            "invalid qual-quant spec: ranges {low}:{high} and {next_low}:{next_high} overlap"
        );
    }

    let mut table: QuantTable = std::array::from_fn(|i| i as u8);
    for (low, high, value) in ranges {
        table[low as usize..high as usize].fill(value);
    }
    Ok(table)
}

/// Parse the number `what` (e.g. "low value") of the qual-quant range `range`.
fn parse_number(text: &str, what: &str, range: &str) -> Result<u8> {
    text.trim().parse().with_context(|| format!("invalid qual-quant {what} in {range:?}"))
}

/// The probability-averaged quality of Phred values in `[low, high)`: the mean of their
/// error probabilities, converted back to Phred and rounded. Every value in the range is
/// weighted equally, whether or not it occurs in the data.
fn average_phred(low: u8, high: u8) -> u8 {
    let mean_error =
        (low..high).map(|q| 10f64.powf(-f64::from(q) / 10.0)).sum::<f64>() / f64::from(high - low);
    // The mean lies between the range's lowest and highest error probabilities, so the
    // result lies in [low, high) and fits in a u8.
    (-10.0 * mean_error.log10()).round() as u8
}

/// Quantize a Phred+33 quality character using the lookup table.
///
/// Strips the +33 offset, applies the table, and re-adds +33.
pub fn quantize_phred33(qual_char: u8, table: &QuantTable) -> u8 {
    let phred = qual_char.saturating_sub(33);
    table[phred as usize] + 33
}

/// Quantize a raw Phred value using the lookup table.
pub fn quantize_phred(phred: u8, table: &QuantTable) -> u8 {
    table[phred as usize]
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[test]
    fn test_parse_identity_outside_ranges() {
        let table = parse_qual_quant("10:20").unwrap();
        assert_eq!(table[0], 0);
        assert_eq!(table[9], 9);
        // `high` is excluded: [10, 20).
        assert_eq!(table[20], 20);
        assert_eq!(table[93], 93);
    }

    // Without `=value`, a range bins to its probability-averaged quality, e.g. Q0-9 to
    // -10 * log10(mean(10^(-q/10))) = 3.59, so Q4.
    #[test]
    fn test_parse_averages_ranges() {
        let table = parse_qual_quant("0:10,10:20,20:30,30:40,40:94").unwrap();
        for (low, high, value) in
            [(0, 10, 4), (10, 20, 14), (20, 30, 24), (30, 40, 34), (40, 94, 50)]
        {
            assert!(table[low..high].iter().all(|&q| q == value), "{low}:{high}");
        }
        assert_eq!(table[94], 94);
    }

    #[test]
    fn test_parse_single_quality_range() {
        let table = parse_qual_quant("30:31").unwrap();
        assert_eq!(table[29], 29);
        assert_eq!(table[30], 30);
        assert_eq!(table[31], 31);
    }

    #[test]
    fn test_parse_explicit_values() {
        let table = parse_qual_quant("1:10=7,10:20=16,40:94=40").unwrap();
        assert_eq!(table[0], 0);
        assert!(table[1..10].iter().all(|&q| q == 7));
        assert!(table[10..20].iter().all(|&q| q == 16));
        assert_eq!(table[20], 20);
        assert!(table[40..94].iter().all(|&q| q == 40));
    }

    #[test]
    fn test_parse_ranges_in_any_order_with_whitespace() {
        let table = parse_qual_quant(" 10:20 = 16 , 0 : 10 ").unwrap();
        assert_eq!(table[5], 4);
        assert_eq!(table[15], 16);
    }

    // The limits themselves are accepted: a value of 93, and a range ending at 94.
    #[test]
    fn test_parse_accepts_limits() {
        let table = parse_qual_quant("0:10=93,93:94").unwrap();
        assert!(table[0..10].iter().all(|&q| q == 93));
        assert_eq!(table[92], 92);
        assert_eq!(table[93], 93);
    }

    // The README's examples: `0:20` bins Q0-19 to Q6 and leaves Q20 and above alone, and
    // `0:10=1,10:20=10,20:30=20,30:94=30` bins Q0-93 as sam-dump's `1:10,10:20,20:30,30:-` does.
    #[test]
    fn test_parse_readme_examples() {
        let table = parse_qual_quant("0:20").unwrap();
        assert!(table[0..20].iter().all(|&q| q == 6));
        assert_eq!(table[20], 20);

        let table = parse_qual_quant("0:10=1,10:20=10,20:30=20,30:94=30").unwrap();
        let sam_dump = |q: usize| match q {
            0..10 => 1,
            10..20 => 10,
            20..30 => 20,
            _ => 30,
        };
        assert!((0..=93).all(|q| table[q] == sam_dump(q)));
    }

    // A range's average is always one of its own qualities, so it never leaves the range.
    #[test]
    fn test_average_phred_within_range() {
        for high in 1..=94 {
            for low in 0..high {
                let value = average_phred(low, high);
                assert!((low..high).contains(&value), "{low}:{high} averaged to {value}");
            }
        }
    }

    #[rstest]
    #[case::empty_spec("", "expected low:high")]
    #[case::empty_range("0:10,,10:20", "expected low:high")]
    #[case::no_colon("not-a-range", "expected low:high")]
    #[case::bad_low("abc:10", "invalid qual-quant low value")]
    #[case::negative_low("-1:10", "invalid qual-quant low value")]
    #[case::bad_high("10:abc", "invalid qual-quant high value")]
    #[case::bad_value("0:10=x", "invalid qual-quant value")]
    #[case::missing_value("0:10=", "invalid qual-quant value")]
    #[case::empty_bounds("10:10", "low must be below high")]
    #[case::reversed_bounds("20:10", "low must be below high")]
    #[case::high_above_94("0:95", "high must be at most 94")]
    #[case::value_above_93("0:10=94", "value must be at most 93")]
    #[case::overlapping("0:10,5:15", "ranges 0:10 and 5:15 overlap")]
    #[case::overlapping_out_of_order("5:15,0:10", "ranges 0:10 and 5:15 overlap")]
    #[case::same_range_twice("0:10,0:10", "ranges 0:10 and 0:10 overlap")]
    fn test_parse_rejects(#[case] spec: &str, #[case] expected: &str) {
        let err = parse_qual_quant(spec).expect_err(spec);
        let msg = format!("{err:#}");
        assert!(msg.contains(expected), "{spec:?}: {msg}");
    }

    // sam-dump's `value:limit,...,value:-` spec is rejected with a pointer to tosam's syntax.
    // Its `v:-` gives a value, not a bound, so the message doesn't suggest a range from it.
    #[test]
    fn test_parse_rejects_sam_dump_spec() {
        let err = parse_qual_quant("2:10,15:-").unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("sam-dump"), "{msg}");
        assert!(msg.contains("end a range at 94"), "{msg}");
        assert!(!msg.contains("15:94"), "{msg}");
    }

    #[test]
    fn test_quantize_phred33() {
        let table = parse_qual_quant("0:10,10:20,20:30").unwrap();
        assert_eq!(quantize_phred33(5 + 33, &table), 4 + 33);
        assert_eq!(quantize_phred33(15 + 33, &table), 14 + 33);
        // Phred 35 is outside every range.
        assert_eq!(quantize_phred33(35 + 33, &table), 35 + 33);
    }

    #[test]
    fn test_quantize_phred() {
        let table = parse_qual_quant("0:10,10:20").unwrap();
        assert_eq!(quantize_phred(5, &table), 4);
        assert_eq!(quantize_phred(15, &table), 14);
        assert_eq!(quantize_phred(25, &table), 25);
    }
}
