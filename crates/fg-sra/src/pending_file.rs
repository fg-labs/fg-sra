//! Output files written under a temporary name and renamed into place on
//! success, so an interrupted or failed run never leaves a partial file at the
//! requested path, nor destroys a file that was already there.

use std::fs::File;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

/// An output file written under a temporary name in the destination directory
/// and renamed into place by [`PendingFile::commit`], so an interrupted or
/// failed run never leaves a partial file at the requested path.
///
/// The temporary file is removed if the writer is dropped without finishing
/// (an error) or the process receives SIGINT, SIGTERM or SIGHUP (unless the
/// signal is ignored); a SIGKILL leaves it behind. If the final rename fails,
/// the completed temporary file is kept and named in the error.
pub(crate) struct PendingFile {
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
    pub(crate) fn create(final_path: &Path) -> Result<(File, Option<Self>)> {
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
    pub(crate) fn tmp_path_for(final_path: &Path) -> PathBuf {
        let mut name =
            final_path.file_name().map(std::ffi::OsStr::to_os_string).unwrap_or_default();
        name.push(format!(".{}.fg-sra-tmp", std::process::id()));
        final_path.with_file_name(name)
    }

    /// Flush the complete temporary file to disk and rename it to the final path.
    ///
    /// If the rename fails, the temporary file (the complete output) is kept and
    /// the error names it.
    pub(crate) fn commit(mut self) -> Result<()> {
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
pub(crate) mod signal_cleanup {
    use std::path::{Path, PathBuf};
    use std::sync::{Mutex, MutexGuard, Once};

    static PENDING: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());
    static INSTALL: Once = Once::new();

    fn pending() -> MutexGuard<'static, Vec<PathBuf>> {
        PENDING.lock().unwrap_or_else(std::sync::PoisonError::into_inner)
    }

    /// Track `path` for removal on termination, installing the handler once.
    pub(crate) fn register(path: &Path) {
        INSTALL.call_once(install);
        pending().push(path.to_path_buf());
    }

    /// Stop tracking `path` (after it was renamed into place or removed).
    pub(crate) fn unregister(path: &Path) {
        pending().retain(|p| p != path);
    }

    /// Whether `path` is tracked for removal.
    #[cfg(test)]
    pub(crate) fn is_tracked(path: &Path) -> bool {
        pending().iter().any(|p| p == path)
    }

    /// Remove every tracked file (best-effort) and stop tracking it.
    pub(crate) fn remove_all() {
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
    pub(crate) fn is_ignored(signal: i32) -> bool {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// A fresh, empty directory under the system temp dir for one test.
    fn test_dir(name: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("fg-sra-pending-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
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
}
