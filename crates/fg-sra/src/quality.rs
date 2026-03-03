//! Quality score quantization.
//!
//! Supports the `--qual-quant` option for binning quality scores into
//! user-defined ranges (e.g. `"1:10,10:20,20:30,30:40"`).

use anyhow::{Context, Result};

/// A 256-byte lookup table for quality score quantization.
/// Index by raw Phred value, returns quantized Phred value.
pub type QuantTable = [u8; 256];

/// Parse a quantization spec like `"1:10,10:20,20:30,30:40"`.
///
/// Each `low:high` pair means: Phred values in `[low, high)` are mapped to `high`.
/// Values outside any range are left unchanged.
pub fn parse_qual_quant(spec: &str) -> Result<QuantTable> {
    let mut table: QuantTable = std::array::from_fn(|i| i as u8);
    for pair in spec.split(',') {
        let (low, high) =
            pair.split_once(':').context("invalid qual-quant format: expected low:high")?;
        let low: u8 = low.trim().parse().context("invalid qual-quant low value")?;
        let high: u8 = high.trim().parse().context("invalid qual-quant high value")?;
        for q in low..high {
            table[q as usize] = high;
        }
    }
    Ok(table)
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
    use super::*;

    #[test]
    fn test_parse_identity() {
        // Values outside any range should be unchanged.
        let table = parse_qual_quant("10:20").unwrap();
        assert_eq!(table[0], 0);
        assert_eq!(table[9], 9);
        assert_eq!(table[20], 20);
        assert_eq!(table[30], 30);
    }

    #[test]
    fn test_parse_single_range() {
        let table = parse_qual_quant("10:20").unwrap();
        assert_eq!(table[10], 20);
        assert_eq!(table[15], 20);
        assert_eq!(table[19], 20);
        // Boundary: 20 itself is NOT included in [10, 20).
        assert_eq!(table[20], 20); // stays at 20 (identity)
    }

    #[test]
    fn test_parse_multiple_ranges() {
        let table = parse_qual_quant("1:10,10:20,20:30,30:40").unwrap();
        assert_eq!(table[0], 0); // below first range
        assert_eq!(table[1], 10);
        assert_eq!(table[9], 10);
        assert_eq!(table[10], 20);
        assert_eq!(table[19], 20);
        assert_eq!(table[20], 30);
        assert_eq!(table[29], 30);
        assert_eq!(table[30], 40);
        assert_eq!(table[39], 40);
        assert_eq!(table[40], 40); // identity
        assert_eq!(table[50], 50); // identity
    }

    #[test]
    fn test_parse_invalid_format() {
        assert!(parse_qual_quant("not-a-range").is_err());
        assert!(parse_qual_quant("10:abc").is_err());
        assert!(parse_qual_quant("abc:10").is_err());
    }

    #[test]
    fn test_quantize_phred33() {
        let table = parse_qual_quant("0:10,10:20,20:30").unwrap();
        // Phred 5 → quantized to 10, re-encoded as 10+33 = 43
        assert_eq!(quantize_phred33(5 + 33, &table), 10 + 33);
        // Phred 15 → quantized to 20, re-encoded as 20+33 = 53
        assert_eq!(quantize_phred33(15 + 33, &table), 20 + 33);
        // Phred 35 → identity (no range), re-encoded as 35+33 = 68
        assert_eq!(quantize_phred33(35 + 33, &table), 35 + 33);
    }

    #[test]
    fn test_quantize_phred() {
        let table = parse_qual_quant("0:10,10:20").unwrap();
        assert_eq!(quantize_phred(5, &table), 10);
        assert_eq!(quantize_phred(15, &table), 20);
        assert_eq!(quantize_phred(25, &table), 25);
    }

    #[test]
    fn test_empty_range() {
        // Same low and high means no values in [low, high) → identity.
        let table = parse_qual_quant("10:10").unwrap();
        assert_eq!(table[10], 10);
    }
}
