//! Read-name templates (`--defline`), in the variable syntax of sra-tools' `--seq-defline`.

use anyhow::{Result, bail};

/// A parsed header template: the text after a record's `@` (FASTQ) or `>` (FASTA).
///
/// Variables are `$ac` (accession), `$si` (spot id), `$sn` (original spot name, or the spot
/// id when there is none), `$sg` (spot group), `$ri` (read number within its type) and `$rl`
/// (read length); `$$` is a literal `$`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Defline {
    tokens: Vec<Token>,
}

impl Defline {
    /// Parse a template.
    ///
    /// A leading `@` or `>` is dropped, so templates written for sra-tools work unchanged.
    /// Fails on an empty template, a line break, or an unknown `$` variable.
    pub fn parse(template: &str) -> Result<Self> {
        let body = template.strip_prefix(['@', '>']).unwrap_or(template);
        if body.is_empty() {
            bail!("the defline template is empty");
        }
        if body.contains(['\n', '\r']) {
            bail!("the defline template {template:?} contains a line break");
        }

        let mut tokens = Vec::new();
        let mut literal = Vec::new();
        let mut rest = body;
        while let Some(dollar) = rest.find('$') {
            literal.extend_from_slice(&rest.as_bytes()[..dollar]);
            let after = &rest[dollar + 1..];
            if let Some(tail) = after.strip_prefix('$') {
                literal.push(b'$');
                rest = tail;
                continue;
            }
            let name = after.get(..2).unwrap_or(after);
            let token = match name {
                "ac" => Token::Accession,
                "si" => Token::SpotId,
                "sn" => Token::SpotName,
                "sg" => Token::SpotGroup,
                "ri" => Token::ReadNumber,
                "rl" => Token::ReadLength,
                _ => bail!(
                    "unknown variable ${name} in the defline template {template:?}; expected \
                     one of $ac, $si, $sn, $sg, $ri, $rl or $$"
                ),
            };
            if !literal.is_empty() {
                tokens.push(Token::Literal(std::mem::take(&mut literal)));
            }
            tokens.push(token);
            rest = &after[2..];
        }
        literal.extend_from_slice(rest.as_bytes());
        if !literal.is_empty() {
            tokens.push(Token::Literal(literal));
        }
        Ok(Self { tokens })
    }

    /// Append the header for `fields` to `out`, without the leading `@`/`>` or line end.
    pub fn render(&self, out: &mut Vec<u8>, fields: &DeflineFields<'_>) {
        let mut numbers = itoa::Buffer::new();
        for token in &self.tokens {
            match token {
                Token::Literal(text) => out.extend_from_slice(text),
                Token::Accession => out.extend_from_slice(fields.accession),
                Token::SpotId => out.extend_from_slice(numbers.format(fields.spot_id).as_bytes()),
                Token::SpotName if fields.spot_name.is_empty() => {
                    out.extend_from_slice(numbers.format(fields.spot_id).as_bytes());
                }
                Token::SpotName => out.extend_from_slice(fields.spot_name),
                Token::SpotGroup => out.extend_from_slice(fields.spot_group),
                Token::ReadNumber => {
                    out.extend_from_slice(numbers.format(fields.read_number).as_bytes());
                }
                Token::ReadLength => {
                    out.extend_from_slice(numbers.format(fields.read_length).as_bytes());
                }
            }
        }
    }

    /// Whether the template uses `$sn`, so that spot names need reading.
    pub fn uses_spot_name(&self) -> bool {
        self.tokens.contains(&Token::SpotName)
    }

    /// Whether the template uses `$sg`, so that spot groups need reading.
    pub fn uses_spot_group(&self) -> bool {
        self.tokens.contains(&Token::SpotGroup)
    }
}

impl std::str::FromStr for Defline {
    type Err = anyhow::Error;

    fn from_str(template: &str) -> Result<Self> {
        Self::parse(template)
    }
}

/// The values a [`Defline`]'s variables are filled from, for one record.
#[derive(Debug, Clone, Copy)]
pub struct DeflineFields<'a> {
    pub accession: &'a [u8],
    pub spot_id: i64,
    /// Original spot name; empty when the archive kept none.
    pub spot_name: &'a [u8],
    pub spot_group: &'a [u8],
    pub read_number: u32,
    pub read_length: u32,
}

/// One piece of a parsed template: literal text or a variable.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    Literal(Vec<u8>),
    Accession,
    SpotId,
    SpotName,
    SpotGroup,
    ReadNumber,
    ReadLength,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fields() -> DeflineFields<'static> {
        DeflineFields {
            accession: b"SRR123",
            spot_id: 42,
            spot_name: b"M0:1:FC:1:1101:10:20",
            spot_group: b"GRP",
            read_number: 2,
            read_length: 151,
        }
    }

    fn render(template: &str, fields: &DeflineFields<'_>) -> String {
        let mut out = Vec::new();
        Defline::parse(template).unwrap().render(&mut out, fields);
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn default_template_is_accession_dot_spot() {
        assert_eq!(render("$ac.$si", &fields()), "SRR123.42");
    }

    #[test]
    fn every_variable_is_substituted() {
        let rendered = render("$ac $si $sn $sg $ri $rl", &fields());
        assert_eq!(rendered, "SRR123 42 M0:1:FC:1:1101:10:20 GRP 2 151");
    }

    #[test]
    fn fasterq_dump_default_template_renders_as_it_does() {
        let rendered = render("$ac.$si $sn length=$rl", &fields());
        assert_eq!(rendered, "SRR123.42 M0:1:FC:1:1101:10:20 length=151");
    }

    #[test]
    fn spot_name_falls_back_to_spot_id_when_there_is_none() {
        let fields = DeflineFields { spot_name: b"", ..fields() };
        assert_eq!(render("$sn/$ri", &fields), "42/2");
    }

    #[test]
    fn empty_spot_group_renders_as_nothing() {
        let fields = DeflineFields { spot_group: b"", ..fields() };
        assert_eq!(render("$ac.$si[$sg]", &fields), "SRR123.42[]");
    }

    #[test]
    fn double_dollar_is_a_literal_dollar() {
        assert_eq!(render("$$$ac$$", &fields()), "$SRR123$");
    }

    #[test]
    fn variable_directly_followed_by_text_is_recognised() {
        assert_eq!(render("$ac$si", &fields()), "SRR12342");
    }

    #[test]
    fn leading_at_sign_is_dropped() {
        assert_eq!(render("@$ac.$si/$ri", &fields()), "SRR123.42/2");
    }

    #[test]
    fn leading_greater_than_sign_is_dropped() {
        assert_eq!(render(">$ac.$si", &fields()), "SRR123.42");
    }

    #[test]
    fn only_one_leading_at_sign_is_dropped() {
        assert_eq!(render("@@$si", &fields()), "@42");
    }

    #[test]
    fn unknown_variable_is_an_error() {
        let err = Defline::parse("$ac.$ix").unwrap_err();
        assert!(err.to_string().contains("$ix"), "{err}");
    }

    #[test]
    fn dollar_at_the_end_is_an_error() {
        assert!(Defline::parse("$ac.$").is_err());
    }

    #[test]
    fn truncated_variable_is_an_error() {
        assert!(Defline::parse("$ac.$s").is_err());
    }

    #[test]
    fn empty_template_is_an_error() {
        assert!(Defline::parse("").is_err());
    }

    #[test]
    fn template_of_only_an_at_sign_is_an_error() {
        assert!(Defline::parse("@").is_err());
    }

    #[test]
    fn line_break_is_an_error() {
        assert!(Defline::parse("$ac\n+").is_err());
    }

    #[test]
    fn multibyte_text_after_a_dollar_is_an_error_not_a_panic() {
        assert!(Defline::parse("$é").is_err());
    }

    #[test]
    fn spot_name_use_is_detected() {
        assert!(Defline::parse("$ac.$si $sn").unwrap().uses_spot_name());
        assert!(!Defline::parse("$ac.$si").unwrap().uses_spot_name());
    }

    #[test]
    fn spot_group_use_is_detected() {
        assert!(Defline::parse("$sg:$si").unwrap().uses_spot_group());
        assert!(!Defline::parse("$ac.$si").unwrap().uses_spot_group());
    }
}
