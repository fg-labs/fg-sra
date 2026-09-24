//! Output dispatch for SAM, BAM, FASTA, and FASTQ formats.
//!
//! Manages the output pipeline including optional gzip/bzip2 compression
//! and BGZF for BAM output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bzip2::write::BzEncoder;
use flate2::write::GzEncoder;

/// Default buffer size for the text output writer (256 KB).
const OUTPUT_BUF_SIZE: usize = 256 * 1024;

/// Compression mode for text output.
#[derive(Clone, Copy)]
pub enum CompressionMode {
    None,
    Gzip,
    Bzip2,
}

/// Abstraction over output destinations (stdout or file, with optional compression).
///
/// For SAM/FASTA/FASTQ text output, uses a `BufWriter` with optional gzip/bzip2
/// compression. For BAM output, uses a BGZF writer that handles block compression.
pub struct OutputWriter {
    inner: WriterInner,
    /// For file output, the temporary file being written and its final path.
    pending: Option<PendingFile>,
}

/// An output file written under a temporary name in the destination directory
/// and renamed into place by [`OutputWriter::finish`], so an interrupted or
/// failed run never leaves a partial file at the requested path.
///
/// The temporary file is removed if the writer is dropped without finishing
/// (an error) or the process receives SIGINT, SIGTERM or SIGHUP (unless the
/// signal is ignored); a SIGKILL leaves it behind. If the final rename fails,
/// the completed temporary file is kept and named in the error.
struct PendingFile {
    tmp_path: PathBuf,
    final_path: PathBuf,
    /// Keep the temporary file when dropped (the rename failed after the output
    /// was complete).
    keep: bool,
}

/// How [`PendingFile::create`] writes a requested output path.
#[derive(Debug, PartialEq, Eq)]
enum WriteMode {
    /// Write a temporary sibling and rename it into place on success.
    TempThenRename,
    /// Open and write the path itself.
    Direct,
}

impl PendingFile {
    /// Open the file to write for `final_path`, returning the pending temporary
    /// file when the output will be renamed into place (see [`write_mode`]).
    ///
    /// Fails before any output is produced if `final_path` names a directory or
    /// is an existing file that is not writable. A temporary file copies the
    /// permissions of an existing output file. `final_path` is written directly
    /// instead when the temporary file cannot be created (e.g. a writable file
    /// in a read-only directory, a name too long for the suffix, or something
    /// already at the temporary path) or would have a different owner or group
    /// than the existing file (the rename would change them).
    fn create(final_path: &Path) -> Result<(File, Option<Self>)> {
        let is_dir = final_path.as_os_str().as_encoded_bytes().ends_with(b"/")
            || std::fs::metadata(final_path).is_ok_and(|m| m.is_dir());
        anyhow::ensure!(!is_dir, "output path is a directory: {}", final_path.display());

        let create_direct = || {
            File::create(final_path)
                .with_context(|| format!("failed to create {}", final_path.display()))
        };
        if write_mode(final_path) == WriteMode::Direct {
            return Ok((create_direct()?, None));
        }
        // A rename would replace an existing file even without write permission
        // on it; require it, as writing in place would, before doing any work.
        let existing = std::fs::metadata(final_path).ok();
        if existing.is_some() {
            std::fs::OpenOptions::new()
                .write(true)
                .open(final_path)
                .with_context(|| format!("cannot write {}", final_path.display()))?;
        }

        let tmp_path = Self::tmp_path_for(final_path);
        // `create_new` so a file or symlink already at the (predictable) temporary
        // path is never followed, reused, or later renamed into place.
        let mut options = std::fs::OpenOptions::new();
        options.write(true).create_new(true);
        // Start no more permissive than an existing output, so the temporary file
        // is never briefly readable by users who cannot read that output.
        #[cfg(unix)]
        if let Some(meta) = &existing {
            use std::os::unix::fs::{OpenOptionsExt, PermissionsExt};
            options.mode(meta.permissions().mode() & 0o777);
        }
        let file = match options.open(&tmp_path) {
            Ok(file) => file,
            Err(e) => {
                eprintln!(
                    "[output] cannot create temporary file {} ({e}); writing {} in place",
                    tmp_path.display(),
                    final_path.display()
                );
                return Ok((create_direct()?, None));
            }
        };
        let pending = Self { tmp_path, final_path: final_path.to_path_buf(), keep: false };
        signal_cleanup::register(&pending.tmp_path);
        if let Some(meta) = existing {
            if !same_owner(&meta, &file.metadata()?) {
                drop(pending); // Removes the temporary file.
                return Ok((create_direct()?, None));
            }
            file.set_permissions(meta.permissions()).with_context(|| {
                format!("failed to set permissions on {}", pending.tmp_path.display())
            })?;
        }
        Ok((file, Some(pending)))
    }

    /// The temporary path for `final_path`: a sibling (so the rename stays on one
    /// filesystem) named `<name>.<pid>.fg-sra-tmp`, so concurrent runs writing
    /// the same path do not share a temporary file.
    fn tmp_path_for(final_path: &Path) -> PathBuf {
        let mut name =
            final_path.file_name().map(std::ffi::OsStr::to_os_string).unwrap_or_default();
        name.push(format!(".{}.fg-sra-tmp", std::process::id()));
        final_path.with_file_name(name)
    }

    /// Flush the complete temporary file to disk and rename it to the final path.
    ///
    /// If the rename fails, the temporary file (the complete output) is kept and
    /// the error names it.
    fn commit(mut self) -> Result<()> {
        // A write handle: on Windows `sync_all` fails on a read-only one.
        std::fs::OpenOptions::new()
            .write(true)
            .open(&self.tmp_path)
            .and_then(|f| f.sync_all())
            .with_context(|| format!("failed to sync {}", self.tmp_path.display()))?;
        if let Err(e) = std::fs::rename(&self.tmp_path, &self.final_path) {
            self.keep = true;
            return Err(e).with_context(|| {
                format!(
                    "failed to rename {} to {}; the completed output was left at {}",
                    self.tmp_path.display(),
                    self.final_path.display(),
                    self.tmp_path.display()
                )
            });
        }
        // Persist the rename itself (best-effort: not every platform can sync a
        // directory).
        // A bare file name's parent is empty: the current directory.
        let dir = match self.final_path.parent() {
            Some(dir) if !dir.as_os_str().is_empty() => dir,
            _ => Path::new("."),
        };
        let _ = File::open(dir).and_then(|d| d.sync_all());
        Ok(())
        // `self` drops here; its temporary path no longer exists after a
        // successful rename, so the cleanup in `drop` is a no-op.
    }
}

impl Drop for PendingFile {
    fn drop(&mut self) {
        if !self.keep {
            let _ = std::fs::remove_file(&self.tmp_path);
        }
        signal_cleanup::unregister(&self.tmp_path);
    }
}

/// Whether two files have the same owner and group (always true off unix).
fn same_owner(a: &std::fs::Metadata, b: &std::fs::Metadata) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        a.uid() == b.uid() && a.gid() == b.gid()
    }
    #[cfg(not(unix))]
    {
        let _ = (a, b);
        true
    }
}

/// How to write `path`: write a temporary file and rename it into place, unless
/// a rename would change what the path refers to or cannot replace it.
///
/// Written directly: a symlink (the rename would replace the link, not its
/// target; this includes `/dev/stdout`), an existing non-regular file such as a
/// named pipe or device (a rename cannot replace it with a regular file), and a
/// regular file with several hard links (a rename would detach this path from
/// the others).
fn write_mode(path: &Path) -> WriteMode {
    let Ok(meta) = std::fs::symlink_metadata(path) else {
        return WriteMode::TempThenRename; // Does not exist yet.
    };
    let multiply_linked = {
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            meta.nlink() > 1
        }
        #[cfg(not(unix))]
        {
            false
        }
    };
    if meta.is_symlink() || !meta.is_file() || multiply_linked {
        WriteMode::Direct
    } else {
        WriteMode::TempThenRename
    }
}

/// Removal of pending temporary files when the process is terminated by
/// SIGINT, SIGTERM or SIGHUP, which would otherwise skip [`PendingFile`]'s
/// `Drop` and leave a multi-gigabyte temporary file behind.
mod signal_cleanup {
    use std::path::{Path, PathBuf};
    use std::sync::{Mutex, MutexGuard, Once};

    static PENDING: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());
    static INSTALL: Once = Once::new();

    fn pending() -> MutexGuard<'static, Vec<PathBuf>> {
        PENDING.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Track `path` for removal on termination, installing the handler once.
    pub(super) fn register(path: &Path) {
        INSTALL.call_once(install);
        pending().push(path.to_path_buf());
    }

    /// Stop tracking `path` (after it was renamed into place or removed).
    pub(super) fn unregister(path: &Path) {
        pending().retain(|p| p != path);
    }

    /// Whether `path` is tracked for removal.
    #[cfg(test)]
    pub(super) fn is_tracked(path: &Path) -> bool {
        pending().iter().any(|p| p == path)
    }

    /// Remove every tracked file (best-effort) and stop tracking it.
    pub(super) fn remove_all() {
        for path in pending().drain(..) {
            let _ = std::fs::remove_file(path);
        }
    }

    /// Spawn a thread that, on SIGINT/SIGTERM/SIGHUP, removes the tracked files
    /// and then terminates the process as the signal's default action would.
    ///
    /// Signals the process inherited as ignored (e.g. SIGHUP under `nohup`, or
    /// SIGINT for a background job of a non-interactive shell) stay ignored, so
    /// the handler never makes them fatal.
    #[cfg(unix)]
    fn install() {
        use signal_hook::consts::{SIGHUP, SIGINT, SIGTERM};
        let handled: Vec<i32> =
            [SIGINT, SIGTERM, SIGHUP].into_iter().filter(|&sig| !is_ignored(sig)).collect();
        if handled.is_empty() {
            return;
        }
        let Ok(mut signals) = signal_hook::iterator::Signals::new(handled) else {
            return; // Without the handler, termination just leaves the file behind.
        };
        std::thread::spawn(move || {
            if let Some(signal) = signals.forever().next() {
                remove_all();
                let _ = signal_hook::low_level::emulate_default_handler(signal);
                std::process::exit(128 + signal);
            }
        });
    }

    /// Whether `signal`'s current disposition is to ignore it.
    #[cfg(unix)]
    pub(super) fn is_ignored(signal: i32) -> bool {
        let mut current: libc::sigaction = unsafe { std::mem::zeroed() };
        // SAFETY: with a null new action, `sigaction` only reads the current
        // disposition of `signal` into `current`, a valid, writable, zeroed
        // `sigaction`; it changes no process state.
        let rc = unsafe { libc::sigaction(signal, std::ptr::null(), &raw mut current) };
        rc == 0 && current.sa_sigaction == libc::SIG_IGN
    }

    #[cfg(not(unix))]
    fn install() {}
}

/// A text output stream, optionally gzip/bzip2-compressed.
///
/// Kept as concrete encoder types (not a `Box<dyn Write>`) so [`Self::finish`]
/// can end a compressed stream and report errors writing its final block and
/// trailer, which the encoders' `Drop` would silently discard.
enum TextSink {
    Plain(Box<dyn Write>),
    Gzip(GzEncoder<Box<dyn Write>>),
    Bzip2(BzEncoder<Box<dyn Write>>),
}

impl TextSink {
    fn new(raw: Box<dyn Write>, mode: CompressionMode) -> Self {
        match mode {
            CompressionMode::None => Self::Plain(raw),
            CompressionMode::Gzip => {
                Self::Gzip(GzEncoder::new(raw, flate2::Compression::default()))
            }
            CompressionMode::Bzip2 => {
                Self::Bzip2(BzEncoder::new(raw, bzip2::Compression::default()))
            }
        }
    }

    /// End the stream (writing any compressed trailer) and flush it.
    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(w) => w.finish()?.flush(),
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
    Bgzf(noodles_bgzf::io::Writer<Box<dyn Write>>),
}

impl OutputWriter {
    /// Create a text writer to stdout with the given compression.
    pub fn stdout_with_compression(mode: CompressionMode) -> Self {
        let sink = TextSink::new(Box::new(io::stdout().lock()), mode);
        Self {
            inner: WriterInner::Text(BufWriter::with_capacity(OUTPUT_BUF_SIZE, sink)),
            pending: None,
        }
    }

    /// Create a text writer to a file with the given compression.
    pub fn from_path_with_compression(path: &Path, mode: CompressionMode) -> Result<Self> {
        let (file, pending) = PendingFile::create(path)?;
        let sink = TextSink::new(Box::new(file), mode);
        Ok(Self {
            inner: WriterInner::Text(BufWriter::with_capacity(OUTPUT_BUF_SIZE, sink)),
            pending,
        })
    }

    /// Create a BAM writer to stdout (BGZF-compressed).
    pub fn bam_stdout() -> Self {
        let raw: Box<dyn Write> = Box::new(io::stdout().lock());
        Self { inner: WriterInner::Bgzf(noodles_bgzf::io::Writer::new(raw)), pending: None }
    }

    /// Create a BAM writer to a file (BGZF-compressed).
    pub fn bam_from_path(path: &Path) -> Result<Self> {
        let (file, pending) = PendingFile::create(path)?;
        let raw: Box<dyn Write> = Box::new(file);
        Ok(Self { inner: WriterInner::Bgzf(noodles_bgzf::io::Writer::new(raw)), pending })
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
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
            let mut writer = OutputWriter::bam_from_path(&path).unwrap();
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
                OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
            writer.write_bytes(b"new\n").unwrap();
        }
        assert_eq!(std::fs::read(&path).unwrap(), b"old\n");

        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
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
            OutputWriter::from_path_with_compression(&link, CompressionMode::None).unwrap();
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

        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();

        assert_eq!(std::fs::read(&victim).unwrap(), b"victim\n");
        assert!(std::fs::symlink_metadata(&path).unwrap().is_file(), "output must not be a link");
        assert_eq!(std::fs::read(&path).unwrap(), b"new\n");
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
        let mut writer = OutputWriter::bam_from_path(&path).unwrap();
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::Gzip).unwrap();
        writer.write_bytes(b"record\n").unwrap();
        writer.finish().unwrap();
        let mut decoded = String::new();
        flate2::read::GzDecoder::new(File::open(&path).unwrap())
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(decoded, "record\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn test_finished_bzip2_output_is_complete() {
        use std::io::Read;
        let dir = test_dir("bzip2");
        let path = dir.join("out.sam.bz2");
        let mut writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::Bzip2).unwrap();
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
        // The disk fills after the 10-byte gzip header; the compressed data and
        // trailer are written in finish(), whose error must be reported, not
        // dropped as the encoders' Drop would.
        for mode in [CompressionMode::Gzip, CompressionMode::Bzip2] {
            let mut sink = TextSink::new(Box::new(FailAfter { left: 10 }), mode);
            sink.write_all(b"record\n").unwrap();
            let err = sink.finish().unwrap_err();
            assert_eq!(err.to_string(), "disk full");
        }
        let mut plain = TextSink::new(Box::new(FailAfter { left: 10 }), CompressionMode::None);
        plain.write_all(b"record\n").unwrap();
        plain.finish().unwrap();
    }

    #[test]
    fn test_write_mode() {
        let dir = test_dir("write-mode");
        let regular = dir.join("regular.sam");
        std::fs::write(&regular, b"x").unwrap();
        assert_eq!(write_mode(&dir.join("missing.sam")), WriteMode::TempThenRename);
        assert_eq!(write_mode(&regular), WriteMode::TempThenRename);
        #[cfg(unix)]
        {
            let link = dir.join("link.sam");
            std::os::unix::fs::symlink(&regular, &link).unwrap();
            assert_eq!(write_mode(&link), WriteMode::Direct, "symlink");
            assert_eq!(write_mode(Path::new("/dev/null")), WriteMode::Direct, "device");
            let hard = dir.join("hard.sam");
            std::fs::hard_link(&regular, &hard).unwrap();
            assert_eq!(write_mode(&regular), WriteMode::Direct, "hard-linked");
        }
        std::fs::remove_dir_all(&dir).ok();
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
        writer.write_bytes(b"new\n").unwrap();
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&other).unwrap(), b"new\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_writing_to_dev_null_is_direct() {
        let mut writer =
            OutputWriter::from_path_with_compression(Path::new("/dev/null"), CompressionMode::None)
                .unwrap();
        assert!(writer.pending.is_none());
        writer.write_bytes(b"record\n").unwrap();
        writer.finish().unwrap();
    }

    #[test]
    fn test_directory_output_path_is_rejected_up_front() {
        let dir = test_dir("dir-path");
        let err = OutputWriter::bam_from_path(&dir).err().expect("directory must be rejected");
        assert!(err.to_string().contains("output path is a directory"), "{err}");
        let trailing = PathBuf::from(format!("{}/out.bam/", dir.display()));
        let err = OutputWriter::bam_from_path(&trailing).err().expect("trailing slash rejected");
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
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
        let result = OutputWriter::from_path_with_compression(&path, CompressionMode::None);
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
                OutputWriter::from_path_with_compression(&path, CompressionMode::None)?;
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
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
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
        assert!(writer.pending.is_none());
        writer.write_bytes(b"record\n").unwrap();
        writer.finish().unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"record\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[cfg(unix)]
    #[test]
    fn test_is_ignored_reports_ignored_signals() {
        // SIGUSR2 is not used by the tests or the handler; ignore it, check, and
        // restore its default disposition.
        let signal = libc::SIGUSR2;
        assert!(!signal_cleanup::is_ignored(signal));
        // SAFETY: setting SIGUSR2's disposition to ignore/default is sound; nothing
        // in this process relies on it.
        unsafe { libc::signal(signal, libc::SIG_IGN) };
        assert!(signal_cleanup::is_ignored(signal));
        unsafe { libc::signal(signal, libc::SIG_DFL) };
        assert!(!signal_cleanup::is_ignored(signal));
    }

    #[test]
    fn test_temp_file_is_tracked_for_signal_cleanup_until_finished() {
        let dir = test_dir("signal");
        let path = dir.join("out.sam");
        let tmp = PendingFile::tmp_path_for(&path);
        let writer =
            OutputWriter::from_path_with_compression(&path, CompressionMode::None).unwrap();
        assert!(signal_cleanup::is_tracked(&tmp));
        writer.finish().unwrap();
        assert!(!signal_cleanup::is_tracked(&tmp));

        let writer = OutputWriter::bam_from_path(&path).unwrap();
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
