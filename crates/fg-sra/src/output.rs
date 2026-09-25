//! Output dispatch for SAM, BAM, FASTA, and FASTQ formats.
//!
//! Manages the output pipeline including optional gzip (as BGZF) or bzip2
//! compression, and BGZF for BAM output, compressed on several threads.

use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use bzip2::write::BzEncoder;

use crate::bgzf_writer::ParallelBgzfWriter;
use crate::pending_file::PendingFile;

/// Default buffer size for the text output writer (256 KB).
const OUTPUT_BUF_SIZE: usize = 256 * 1024;

/// Compression mode for text output.
#[derive(Clone, Copy)]
pub enum CompressionMode {
    None,
    /// BGZF: gzip that any gzip reader reads, compressed on several threads.
    Gzip,
    Bzip2,
}

/// Abstraction over output destinations (stdout or file, with optional compression).
///
/// For SAM/FASTA/FASTQ text output, uses a `BufWriter` with optional gzip (BGZF)
/// or bzip2 compression. For BAM output, uses a BGZF writer. BGZF blocks are
/// compressed on several threads.
pub struct OutputWriter {
    inner: WriterInner,
    /// For file output, the temporary file being written and its final path.
    pending: Option<PendingFile>,
}

/// A text output stream, optionally gzip/bzip2-compressed.
///
/// Kept as concrete encoder types (not a `Box<dyn Write>`) so [`Self::finish`]
/// can end a compressed stream and report errors writing its final block and
/// trailer, which the encoders' `Drop` would silently discard.
enum TextSink {
    Plain(Box<dyn Write + Send>),
    Gzip(ParallelBgzfWriter),
    Bzip2(BzEncoder<Box<dyn Write + Send>>),
}

impl TextSink {
    /// A sink writing to `raw`, compressing gzip on `threads` threads.
    fn new(raw: Box<dyn Write + Send>, mode: CompressionMode, threads: usize) -> Self {
        match mode {
            CompressionMode::None => Self::Plain(raw),
            CompressionMode::Gzip => Self::Gzip(ParallelBgzfWriter::new(raw, threads)),
            CompressionMode::Bzip2 => {
                Self::Bzip2(BzEncoder::new(raw, bzip2::Compression::default()))
            }
        }
    }

    /// End the stream (writing any compressed trailer) and flush it.
    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(w) => w.finish(),
            Self::Bzip2(w) => w.finish()?.flush(),
        }
    }
}

impl Write for TextSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gzip(w) => w.write(buf),
            Self::Bzip2(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
            Self::Bzip2(w) => w.flush(),
        }
    }
}

/// Internal writer variant: text (buffered) or BAM (BGZF-compressed).
enum WriterInner {
    Text(BufWriter<TextSink>),
    Bgzf(ParallelBgzfWriter),
}

impl OutputWriter {
    /// Create a text writer to stdout with the given compression, compressing on
    /// `threads` threads.
    pub fn stdout_with_compression(mode: CompressionMode, threads: usize) -> Self {
        let sink = TextSink::new(Box::new(io::stdout()), mode, threads);
        Self {
            inner: WriterInner::Text(BufWriter::with_capacity(OUTPUT_BUF_SIZE, sink)),
            pending: None,
        }
    }

    /// Create a text writer to a file with the given compression, compressing on
    /// `threads` threads.
    pub fn from_path_with_compression(
        path: &Path,
        mode: CompressionMode,
        threads: usize,
    ) -> Result<Self> {
        let (file, pending) = PendingFile::create(path)?;
        let sink = TextSink::new(Box::new(file), mode, threads);
        Ok(Self {
            inner: WriterInner::Text(BufWriter::with_capacity(OUTPUT_BUF_SIZE, sink)),
            pending,
        })
    }

    /// Create a BAM writer to stdout, compressing on `threads` threads.
    pub fn bam_stdout(threads: usize) -> Self {
        let bgzf = ParallelBgzfWriter::new(Box::new(io::stdout()), threads);
        Self { inner: WriterInner::Bgzf(bgzf), pending: None }
    }

    /// Create a BAM writer to a file, compressing on `threads` threads.
    pub fn bam_from_path(path: &Path, threads: usize) -> Result<Self> {
        let (file, pending) = PendingFile::create(path)?;
        let bgzf = ParallelBgzfWriter::new(Box::new(file), threads);
        Ok(Self { inner: WriterInner::Bgzf(bgzf), pending })
    }

    /// Write a SAM header (text mode) or BAM header (BAM mode).
    ///
    /// For text mode, writes the header string verbatim.
    /// For BAM mode, encodes the BAM header: magic bytes, SAM header text,
    /// and reference sequence dictionary parsed from `@SQ` lines.
    pub fn write_header(&mut self, header: &str) -> Result<()> {
        match &mut self.inner {
            WriterInner::Text(w) => {
                w.write_all(header.as_bytes())?;
            }
            WriterInner::Bgzf(w) => {
                write_bam_header(w, header)?;
            }
        }
        Ok(())
    }

    /// Write pre-formatted bytes to the output.
    pub fn write_bytes(&mut self, data: &[u8]) -> Result<()> {
        match &mut self.inner {
            WriterInner::Text(w) => w.write_all(data)?,
            WriterInner::Bgzf(w) => w.write_all(data)?,
        }
        Ok(())
    }

    /// Flush and finalize the output.
    ///
    /// For BAM mode, writes the BGZF EOF marker. For file output, then renames
    /// the temporary file to the requested path. A writer dropped without
    /// calling `finish` (e.g. on an error) removes its temporary file.
    pub fn finish(self) -> Result<()> {
        let Self { inner, pending } = self;
        match inner {
            WriterInner::Text(w) => {
                // The stream must be complete (including any compressed trailer)
                // before the file is renamed into place.
                let sink = w
                    .into_inner()
                    .map_err(io::IntoInnerError::into_error)
                    .context("failed to flush output")?;
                sink.finish().context("failed to finalize output")?;
            }
            WriterInner::Bgzf(w) => {
                w.finish().context("failed to finalize BAM output")?;
            }
        }
        if let Some(pending) = pending {
            pending.commit()?;
        }
        Ok(())
    }
}

impl io::Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match &mut self.inner {
            WriterInner::Text(w) => w.write(buf),
            WriterInner::Bgzf(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match &mut self.inner {
            WriterInner::Text(w) => w.flush(),
            WriterInner::Bgzf(w) => w.flush(),
        }
    }
}

/// Write the BAM header section through a BGZF writer.
///
/// Encodes: magic (`BAM\1`), SAM header text, and reference sequence
/// dictionary (from `@SQ` lines in the header text).
fn write_bam_header(writer: &mut impl Write, header_text: &str) -> Result<()> {
    // BAM magic.
    writer.write_all(b"BAM\x01")?;

    // SAM header text.
    let text_bytes = header_text.as_bytes();
    writer.write_all(&(text_bytes.len() as i32).to_le_bytes())?;
    writer.write_all(text_bytes)?;

    // Parse @SQ lines for the reference dictionary.
    let refs = parse_sq_lines(header_text);

    writer.write_all(&(refs.len() as i32).to_le_bytes())?;
    for (name, length) in &refs {
        let name_bytes = name.as_bytes();
        writer.write_all(&((name_bytes.len() + 1) as i32).to_le_bytes())?;
        writer.write_all(name_bytes)?;
        writer.write_all(&[0])?; // null terminator
        writer.write_all(&length.to_le_bytes())?;
    }

    Ok(())
}

/// Parse `@SQ` lines from a SAM header, returning `(name, length)` pairs.
fn parse_sq_lines(header: &str) -> Vec<(&str, i32)> {
    let mut refs = Vec::new();
    for line in header.lines() {
        if !line.starts_with("@SQ\t") {
            continue;
        }
        let mut name = "";
        let mut length: i32 = 0;
        for field in line.split('\t').skip(1) {
            if let Some(val) = field.strip_prefix("SN:") {
                name = val;
            } else if let Some(val) = field.strip_prefix("LN:") {
                length = val.parse().unwrap_or(0);
            }
        }
        if !name.is_empty() {
            refs.push((name, length));
        }
    }
    refs
}

/// Build a map from reference sequence name to BAM reference ID (0-based index).
///
/// Parses `@SQ` lines from the SAM header to establish the mapping.
pub fn build_ref_name_to_id(header: &str) -> std::collections::HashMap<String, i32> {
    parse_sq_lines(header)
        .into_iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_owned(), i as i32))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::path::PathBuf;

    use crate::pending_file::signal_cleanup;

    /// A fresh, empty directory under the system temp dir for one test.
    fn test_dir(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("fg-sra-output-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn test_file_output_appears_only_after_finish() {
        let dir = test_dir("finish");
        let path = dir.join("out.sam");
        let tmp = PendingFile::tmp_path_for(&path);
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"record\n").unwrap();
        assert!(!path.exists(), "final path must not exist before finish");
        assert!(tmp.exists());

        writer.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"record\n");
        assert!(!tmp.exists());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_unfinished_output_leaves_no_file() {
        let dir = test_dir("unfinished");
        let path = dir.join("out.bam");
        {
            let mut writer = OutputWriter::bam_from_path(&path, 2).unwrap();
            writer.write_header("@HD\tVN:1.6\n").unwrap();
            // Dropped without finish, as on an error.
        }
        assert!(!path.exists());
        assert!(!PendingFile::tmp_path_for(&path).exists());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_unfinished_output_keeps_existing_file() {
        let dir = test_dir("existing");
        let path = dir.join("out.sam");
        std::fs::write(&path, b"old\n").unwrap();
        {
            let mut writer =
                OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
            writer.write_bytes(b"new\n").unwrap();
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"old\n");

        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"new\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_symlinked_output_writes_through_link() {
        let dir = test_dir("symlink");
        let target = dir.join("target.sam");
        let link = dir.join("link.sam");
        std::fs::write(&target, b"old\n").unwrap();
        std::os::unix::fs::symlink(&target, &link).unwrap();

        let mut writer =
            OutputWriter::from_path_with_compression(&link, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();

        assert!(std::fs::symlink_metadata(&link).unwrap().is_symlink(), "link must survive");
        assert_eq!(std::fs::read(&target).unwrap(), b"new\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_existing_temp_path_is_not_followed_or_reused() {
        let dir = test_dir("tmp-exists");
        let path = dir.join("out.sam");
        let victim = dir.join("victim.sam");
        std::fs::write(&victim, b"victim\n").unwrap();
        let tmp = PendingFile::tmp_path_for(&path);
        std::os::unix::fs::symlink(&victim, &tmp).unwrap();

        let result = OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2);
        assert!(result.is_err(), "a file already at the temporary path is refused");

        assert_eq!(std::fs::read(&victim).unwrap(), b"victim\n");
        assert!(!path.exists());
        assert!(std::fs::symlink_metadata(&tmp).unwrap().is_symlink(), "not ours to remove");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_finished_bam_ends_with_bgzf_eof() {
        // The 28-byte empty BGZF block that marks a complete BAM.
        const BGZF_EOF: [u8; 28] = [
            0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
            0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];
        let dir = test_dir("bam");
        let path = dir.join("out.bam");
        let mut writer = OutputWriter::bam_from_path(&path, 2).unwrap();
        writer.write_header("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n").unwrap();
        writer.finish().unwrap();
        let bytes = std::fs::read(&path).unwrap();
        assert!(bytes.ends_with(&BGZF_EOF));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_finished_gzip_output_is_complete() {
        use std::io::Read;
        let dir = test_dir("gzip");
        let path = dir.join("out.sam.gz");
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::Gzip, 2).unwrap();
        // Several BGZF blocks, so a lost or reordered later block would show.
        let mut records = String::new();
        for i in 0..30_000 {
            std::fmt::Write::write_fmt(&mut records, format_args!("record {i}\n")).unwrap();
        }
        writer.write_bytes(records.as_bytes()).unwrap();
        writer.finish().unwrap();
        let mut decoded = String::new();
        flate2::read::MultiGzDecoder::new(File::open(&path).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert!(decoded == records, "decoded output is not what was written");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_finished_bzip2_output_is_complete() {
        use std::io::Read;
        let dir = test_dir("bzip2");
        let path = dir.join("out.sam.bz2");
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::Bzip2, 2).unwrap();
        writer.write_bytes(b"record\n").unwrap();
        writer.finish().unwrap();
        let mut decoded = String::new();
        bzip2::read::BzDecoder::new(File::open(&path).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(decoded, "record\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A writer that accepts `left` bytes and then fails, as on a full disk.
    struct FailAfter {
        left: usize,
    }

    impl Write for FailAfter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.left == 0 {
                return Err(io::Error::other("disk full"));
            }
            let n = buf.len().min(self.left);
            self.left -= n;
            Ok(n)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn test_text_sink_finish_reports_trailer_write_error() {
        // The disk fills 10 bytes in; the compressed data and trailer are written
        // by finish(), whose error must be reported, not dropped as an encoder's
        // Drop would.
        for mode in [CompressionMode::Gzip, CompressionMode::Bzip2] {
            let mut sink = TextSink::new(Box::new(FailAfter { left: 10 }), mode, 2);
            sink.write_all(b"record\n").unwrap();
            let err = sink.finish().unwrap_err();
            assert_eq!(err.to_string(), "disk full");
        }
        let mut plain = TextSink::new(Box::new(FailAfter { left: 10 }), CompressionMode::None, 2);
        plain.write_all(b"record\n").unwrap();
        plain.finish().unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn test_hard_linked_output_keeps_links() {
        let dir = test_dir("hardlink");
        let path = dir.join("out.sam");
        let other = dir.join("other.sam");
        std::fs::write(&path, b"old\n").unwrap();
        std::fs::hard_link(&path, &other).unwrap();
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&other).unwrap(), b"new\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_writing_to_dev_null_is_direct() {
        let mut writer = OutputWriter::from_path_with_compression(
            Path::new("/dev/null"),
            CompressionMode::None,
            2,
        )
        .unwrap();
        assert!(writer.pending.is_none());
        writer.write_bytes(b"record\n").unwrap();
        writer.finish().unwrap();
    }

    #[test]
    fn test_directory_output_path_is_rejected_up_front() {
        let dir = test_dir("dir-path");
        let err = OutputWriter::bam_from_path(&dir, 2).err().expect("directory must be rejected");
        assert!(err.to_string().contains("output path is a directory"), "{err}");
        let trailing = PathBuf::from(format!("{}/out.bam/", dir.display()));
        let err = OutputWriter::bam_from_path(&trailing, 2).err().expect("trailing slash rejected");
        assert!(err.to_string().contains("output path is a directory"), "{err}");
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0, "nothing may be created");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_replacing_output_keeps_its_permissions() {
        use std::os::unix::fs::PermissionsExt;
        let dir = test_dir("perms");
        let path = dir.join("out.sam");
        std::fs::write(&path, b"old\n").unwrap();
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o640)).unwrap();
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();
        let mode = std::fs::metadata(&path).unwrap().permissions().mode() & 0o777;
        assert_eq!(mode, 0o640);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_replacing_output_keeps_its_group() {
        use std::os::unix::fs::MetadataExt;
        let dir = test_dir("group");
        let path = dir.join("out.sam");
        std::fs::write(&path, b"old\n").unwrap();
        // Give the output another of this user's groups than a new file in `dir` gets.
        let new_file_gid = std::fs::metadata(&path).unwrap().gid();
        let Some(other_gid) = user_groups().into_iter().find(|&gid| gid != new_file_gid) else {
            eprintln!("skipping: this user belongs to only one group");
            return;
        };
        std::os::unix::fs::chown(&path, None, Some(other_gid)).unwrap();
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"new\n");
        assert_eq!(std::fs::metadata(&path).unwrap().gid(), other_gid);
        std::fs::remove_dir_all(&dir).ok();
    }

    /// The groups this process's user belongs to.
    #[cfg(unix)]
    fn user_groups() -> Vec<u32> {
        let mut groups = vec![0 as libc::gid_t; 256];
        // SAFETY: `groups` has room for the count passed.
        let n = unsafe { libc::getgroups(groups.len() as libc::c_int, groups.as_mut_ptr()) };
        groups.truncate(usize::try_from(n).unwrap_or(0));
        groups
    }

    #[cfg(unix)]
    #[test]
    fn test_read_only_output_file_is_rejected_up_front() {
        use std::os::unix::fs::PermissionsExt;
        let dir = test_dir("ro-file");
        let path = dir.join("out.sam");
        std::fs::write(&path, b"old\n").unwrap();
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o444)).unwrap();
        let result = OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2);
        // Root bypasses file permissions, so only check when they apply.
        if std::fs::OpenOptions::new().write(true).open(&path).is_err() {
            let err = result.err().expect("read-only output must be rejected");
            assert!(err.to_string().contains("cannot write"), "{err}");
            assert_eq!(std::fs::read(&path).unwrap(), b"old\n");
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_writable_file_in_read_only_directory_is_written() {
        use std::os::unix::fs::PermissionsExt;
        let dir = test_dir("ro-dir");
        let path = dir.join("out.sam");
        std::fs::write(&path, b"old\n").unwrap();
        std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o555)).unwrap();
        let result = (|| -> Result<()> {
            let mut writer =
                OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2)?;
            writer.write_bytes(b"new\n")?;
            writer.finish()
        })();
        std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o755)).unwrap();
        result.unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"new\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_failed_rename_keeps_completed_output() {
        let dir = test_dir("rename-fails");
        let path = dir.join("out.sam");
        let tmp = PendingFile::tmp_path_for(&path);
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        writer.write_bytes(b"record\n").unwrap();
        // A non-empty directory now occupies the final path, so the rename fails.
        std::fs::create_dir(&path).unwrap();
        std::fs::write(path.join("blocker"), b"").unwrap();
        let err = writer.finish().unwrap_err();
        assert!(format!("{err:#}").contains(&format!("left at {}", tmp.display())), "{err:#}");
        assert_eq!(std::fs::read(&tmp).unwrap(), b"record\n");
        assert!(!signal_cleanup::is_tracked(&tmp));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_name_too_long_for_temp_suffix_is_written_in_place() {
        let dir = test_dir("long-name");
        // Fits the usual 255-byte file name limit, but not with the temp suffix.
        let path = dir.join(format!("{}.sam", "x".repeat(246)));
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        assert!(writer.pending.is_none());
        writer.write_bytes(b"record\n").unwrap();
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"record\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_temp_file_is_tracked_for_signal_cleanup_until_finished() {
        let dir = test_dir("signal");
        let path = dir.join("out.sam");
        let tmp = PendingFile::tmp_path_for(&path);
        let writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None, 2).unwrap();
        assert!(signal_cleanup::is_tracked(&tmp));
        writer.finish().unwrap();
        assert!(!signal_cleanup::is_tracked(&tmp));

        let writer = OutputWriter::bam_from_path(&path, 2).unwrap();
        assert!(signal_cleanup::is_tracked(&tmp));
        drop(writer);
        assert!(!signal_cleanup::is_tracked(&tmp));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_parse_sq_lines() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n";
        let refs = parse_sq_lines(header);
        assert_eq!(refs, vec![("chr1", 248_956_422), ("chr2", 242_193_529)]);
    }

    #[test]
    fn test_parse_sq_lines_empty() {
        let header = "@HD\tVN:1.6\n@CO\tsome comment\n";
        let refs = parse_sq_lines(header);
        assert!(refs.is_empty());
    }

    #[test]
    fn test_build_ref_name_to_id() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:200\n";
        let map = build_ref_name_to_id(header);
        assert_eq!(map.get("chr1"), Some(&0));
        assert_eq!(map.get("chr2"), Some(&1));
        assert_eq!(map.get("chr3"), None);
    }

    #[test]
    fn test_write_bam_header() {
        let header = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n";
        let mut buf = Vec::new();
        write_bam_header(&mut buf, header).unwrap();

        // Check magic.
        assert_eq!(&buf[0..4], b"BAM\x01");

        // Check header text length.
        let text_len = i32::from_le_bytes(buf[4..8].try_into().unwrap());
        assert_eq!(text_len, header.len() as i32);

        // Check header text.
        let text_end = 8 + text_len as usize;
        assert_eq!(&buf[8..text_end], header.as_bytes());

        // Check number of references.
        let n_ref = i32::from_le_bytes(buf[text_end..text_end + 4].try_into().unwrap());
        assert_eq!(n_ref, 1);

        // Check first reference name.
        let name_len = i32::from_le_bytes(buf[text_end + 4..text_end + 8].try_into().unwrap());
        assert_eq!(name_len, 5); // "chr1\0"
        let name_end = text_end + 8 + name_len as usize;
        assert_eq!(&buf[text_end + 8..name_end - 1], b"chr1");
        assert_eq!(buf[name_end - 1], 0); // null terminator

        // Check reference length.
        let ref_len = i32::from_le_bytes(buf[name_end..name_end + 4].try_into().unwrap());
        assert_eq!(ref_len, 100);
    }
}
