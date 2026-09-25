//! FASTQ and FASTA record formatting.

use super::defline::{Defline, DeflineFields};
use super::spot::{SelectedRead, Spot};

/// Highest phred value FASTQ's phred+33 encoding can hold (`~`); higher values are clamped.
const MAX_PHRED: u8 = 93;

/// Offset of the phred+33 (Sanger) quality encoding.
const PHRED_OFFSET: u8 = 33;

/// Formats reads as FASTQ or FASTA records.
#[derive(Debug, Clone)]
pub struct RecordFormatter {
    defline: Defline,
    accession: Vec<u8>,
    fasta: bool,
}

impl RecordFormatter {
    /// A formatter naming reads with `defline`, where `$ac` is `accession`. With `fasta`
    /// set, records are unwrapped FASTA; otherwise FASTQ with a bare `+` line.
    pub fn new(defline: Defline, accession: &str, fasta: bool) -> Self {
        Self { defline, accession: accession.as_bytes().to_vec(), fasta }
    }

    /// Append the record for `read` of `spot` to `out`.
    pub fn write(&self, out: &mut Vec<u8>, spot: &Spot<'_>, read: SelectedRead) {
        let bases = spot.read_bases(read.slot);
        let fields = DeflineFields {
            accession: &self.accession,
            spot_id: spot.id,
            spot_name: spot.name,
            spot_group: spot.group,
            read_number: read.number,
            read_length: spot.read_len(read.slot),
        };
        out.push(if self.fasta { b'>' } else { b'@' });
        self.defline.render(out, &fields);
        out.push(b'\n');
        out.extend_from_slice(bases);
        if self.fasta {
            out.push(b'\n');
            return;
        }
        out.extend_from_slice(b"\n+\n");
        out.extend(spot.read_qualities(read.slot).iter().map(|&q| q.min(MAX_PHRED) + PHRED_OFFSET));
        out.push(b'\n');
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::spot::tests::TestSpot;
    use crate::record::READ_TYPE_BIOLOGICAL as B;

    fn formatter(template: &str, fasta: bool) -> RecordFormatter {
        RecordFormatter::new(Defline::parse(template).unwrap(), "SRR1", fasta)
    }

    fn format(formatter: &RecordFormatter, spot: &TestSpot, read: SelectedRead) -> String {
        let mut out = Vec::new();
        formatter.write(&mut out, &spot.view(), read);
        String::from_utf8(out).unwrap()
    }

    fn read(slot: usize, number: u32) -> SelectedRead {
        SelectedRead { slot, number }
    }

    #[test]
    fn fastq_record_has_a_bare_plus_line_and_phred_33_qualities() {
        let mut spot = TestSpot::new(&[("ACGT", B, 0)]);
        spot.qualities = vec![0, 10, 30, 40];
        let record = format(&formatter("$ac.$si", false), &spot, read(0, 1));
        assert_eq!(record, "@SRR1.1\nACGT\n+\n!+?I\n");
    }

    #[test]
    fn qualities_above_93_are_clamped_to_tilde() {
        let mut spot = TestSpot::new(&[("ACG", B, 0)]);
        spot.qualities = vec![93, 94, 255];
        let record = format(&formatter("$ac.$si", false), &spot, read(0, 1));
        assert_eq!(record, "@SRR1.1\nACG\n+\n~~~\n");
    }

    #[test]
    fn record_holds_only_the_selected_read() {
        let spot = TestSpot::new(&[("AAAA", B, 0), ("CCG", B, 0)]);
        let record = format(&formatter("$ac.$si/$ri", false), &spot, read(1, 2));
        assert_eq!(record, "@SRR1.1/2\nCCG\n+\n/01\n");
    }

    #[test]
    fn read_length_variable_is_the_selected_read_length() {
        let spot = TestSpot::new(&[("AAAA", B, 0), ("CCG", B, 0)]);
        let record = format(&formatter("$si length=$rl", true), &spot, read(1, 2));
        assert_eq!(record, ">1 length=3\nCCG\n");
    }

    #[test]
    fn fasta_record_has_no_qualities() {
        let mut spot = TestSpot::new(&[("ACGT", B, 0)]);
        spot.qualities.clear();
        let record = format(&formatter("$ac.$si", true), &spot, read(0, 1));
        assert_eq!(record, ">SRR1.1\nACGT\n");
    }

    #[test]
    fn spot_name_and_group_come_from_the_spot() {
        let mut spot = TestSpot::new(&[("A", B, 0)]);
        spot.id = 7;
        spot.name = b"M0:1:FC:1:1:2:3".to_vec();
        spot.spot_group = b"RG1".to_vec();
        let record = format(&formatter("$sn $sg", true), &spot, read(0, 1));
        assert_eq!(record, ">M0:1:FC:1:1:2:3 RG1\nA\n");
    }

    #[test]
    fn records_append_to_the_buffer() {
        let spot = TestSpot::new(&[("A", B, 0), ("C", B, 0)]);
        let formatter = formatter("$si", true);
        let mut out = b"existing\n".to_vec();
        formatter.write(&mut out, &spot.view(), read(0, 1));
        formatter.write(&mut out, &spot.view(), read(1, 2));
        assert_eq!(String::from_utf8(out).unwrap(), "existing\n>1\nA\n>1\nC\n");
    }
}
