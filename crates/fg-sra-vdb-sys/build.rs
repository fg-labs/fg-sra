//! Build script for fg-sra-vdb-sys.
//!
//! Generates Rust FFI bindings to ncbi-vdb via bindgen.
//!
//! Library resolution (checked in order):
//! 1. `VDB_INCDIR` / `VDB_LIBDIR` env vars — use a pre-built ncbi-vdb
//! 2. `vendored` cargo feature — build ncbi-vdb from `vendor/ncbi-vdb` (a git submodule in the
//!    repository, and part of the published package) via cmake
//! 3. Neither — fail with a helpful error message
//!
//! The vendored source is copied to `OUT_DIR` and built there, never in place: the copy is
//! patched with `vendor/patches/ncbi-vdb/*.patch` (fixes not yet in an ncbi-vdb release), and
//! cmake's configure step writes into the source tree, which `cargo publish` rejects. With the
//! `zlib-ng` feature, the vendored library's bundled zlib is removed after the build so that
//! zlib-ng, from `libz-sys`, provides zlib instead.

use std::env;
use std::path::{Path, PathBuf};

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Determine include and library paths.
    let (inc_dir, lib_dir) = match (env::var("VDB_INCDIR"), env::var("VDB_LIBDIR")) {
        (Ok(inc), Ok(lib)) => {
            println!("cargo:warning=Using pre-built VDB: inc={inc}, lib={lib}");
            (PathBuf::from(inc), PathBuf::from(lib))
        }
        _ => build_vendored(&out_dir),
    };

    // Tell cargo where to find the libraries.
    println!("cargo:rustc-link-search=native={}", lib_dir.display());

    // Link against the ncbi-vdb static library.
    // ncbi-vdb is an uber-library that bundles all needed sub-libraries.
    println!("cargo:rustc-link-lib=static=ncbi-vdb");

    // System libraries required by ncbi-vdb.
    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=framework=Security");
        println!("cargo:rustc-link-lib=dylib=c++");
    } else {
        println!("cargo:rustc-link-lib=dylib=stdc++");
        println!("cargo:rustc-link-lib=dylib=dl");
        println!("cargo:rustc-link-lib=dylib=pthread");
    }

    // zlib: the system's, unless zlib-ng (from `libz-sys`) replaces the bundled one.
    if env::var_os("CARGO_FEATURE_ZLIB_NG").is_none() {
        println!("cargo:rustc-link-lib=dylib=z");
    }

    // Generate FFI bindings via bindgen.
    generate_bindings(&inc_dir, &out_dir);

    // Rerun if the wrapper header changes.
    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-env-changed=VDB_INCDIR");
    println!("cargo:rerun-if-env-changed=VDB_LIBDIR");
}

/// Build ncbi-vdb from the vendored submodule, or fail if vendored feature is disabled.
fn build_vendored(_out_dir: &Path) -> (PathBuf, PathBuf) {
    #[cfg(feature = "vendored")]
    {
        build_ncbi_vdb(_out_dir)
    }

    #[cfg(not(feature = "vendored"))]
    {
        panic!(
            "\n\
            ncbi-vdb not found. Either:\n\
            \n\
            1. Set VDB_INCDIR and VDB_LIBDIR to point at a pre-built ncbi-vdb, or\n\
            2. Enable the `vendored` feature (on by default) to build the vendored ncbi-vdb:\n\
            \n\
            cargo build --features fg-sra-vdb-sys/vendored\n\
            "
        );
    }
}

/// Build ncbi-vdb from the vendored source using cmake, in a patched copy under `out_dir`.
#[cfg(feature = "vendored")]
fn build_ncbi_vdb(out_dir: &Path) -> (PathBuf, PathBuf) {
    let vendor = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()).join("vendor");
    let vendored_src = vendor.join("ncbi-vdb");
    assert!(
        vendored_src.join("CMakeLists.txt").exists(),
        "{} not found; did you initialize the git submodule?\n\
         Run: git submodule update --init --recursive",
        vendored_src.display()
    );
    // Rebuild if the vendored source changes (the copy is made again).
    println!("cargo:rerun-if-changed={}", vendored_src.display());

    let patch_dir = vendor.join("patches/ncbi-vdb");
    // Rerun when a patch is added or removed, not only when one changes.
    println!("cargo:rerun-if-changed={}", patch_dir.display());

    let vdb_src = out_dir.join("ncbi-vdb");
    copy_source(&vendored_src, &vdb_src);
    apply_patches(&vdb_src, &patch_dir);
    // The copy keeps the vendored files' mtimes, so cmake rebuilds only what changed. A
    // patch that is dropped or narrowed, though, leaves files older than the objects built
    // from their patched versions: build from scratch whenever the set of patches changes.
    let stamp = out_dir.join("ncbi-vdb-patches.stamp");
    let patches = patch_set_stamp(&patch_dir);
    if std::fs::read_to_string(&stamp).ok().as_deref() != Some(patches.as_str()) {
        let build_dir = out_dir.join("build");
        if build_dir.exists() {
            std::fs::remove_dir_all(&build_dir)
                .unwrap_or_else(|e| panic!("{}: {e}", build_dir.display()));
        }
    }

    let dst =
        cmake::Config::new(&vdb_src).define("LIBS_ONLY", "ON").build_target("ncbi-vdb").build();

    let inc_dir = vdb_src.join("interfaces");
    let build_dir = dst.join("build");

    // cmake puts the static uber-library in lib/, and helper libs in ilib/.
    let lib_dir = build_dir.join("lib");
    let helper_lib_dir = build_dir.join("ilib");

    // Determine which directory contains the uber-library.
    let final_lib_dir = if lib_dir.join("libncbi-vdb.a").exists() {
        lib_dir
    } else if helper_lib_dir.join("libncbi-vdb.a").exists() {
        helper_lib_dir.clone()
    } else {
        panic!(
            "libncbi-vdb.a not found after cmake build. Searched:\n  {}\n  {}",
            lib_dir.join("libncbi-vdb.a").display(),
            helper_lib_dir.join("libncbi-vdb.a").display()
        );
    };

    #[cfg(feature = "zlib-ng")]
    remove_bundled_zlib(&vdb_src, &final_lib_dir.join("libncbi-vdb.a"));

    std::fs::write(&stamp, patches).unwrap_or_else(|e| panic!("{}: {e}", stamp.display()));

    // mbedcrypto is built as a separate static lib in ilib/.
    if helper_lib_dir.join("libmbedcrypto.a").exists() {
        println!("cargo:rustc-link-search=native={}", helper_lib_dir.display());
        println!("cargo:rustc-link-lib=static=mbedcrypto");
    }

    (inc_dir, final_lib_dir)
}

/// The `*.patch` files in `patch_dir`, in name order.
#[cfg(feature = "vendored")]
fn patch_files(patch_dir: &Path) -> Vec<PathBuf> {
    let mut patches: Vec<PathBuf> = std::fs::read_dir(patch_dir)
        .unwrap_or_else(|e| panic!("{}: {e}", patch_dir.display()))
        .filter_map(|entry| Some(entry.ok()?.path()))
        .filter(|path| path.extension().is_some_and(|ext| ext == "patch"))
        .collect();
    patches.sort();
    patches
}

/// A stamp of the patches in `patch_dir`: a hash of their names and contents.
#[cfg(feature = "vendored")]
fn patch_set_stamp(patch_dir: &Path) -> String {
    use std::hash::{DefaultHasher, Hash, Hasher};

    let mut hasher = DefaultHasher::new();
    for patch in patch_files(patch_dir) {
        patch.file_name().hash(&mut hasher);
        std::fs::read(&patch)
            .unwrap_or_else(|e| panic!("{}: {e}", patch.display()))
            .hash(&mut hasher);
    }
    format!("{:016x}", hasher.finish())
}

/// Copy the ncbi-vdb source at `src` to `dst`, replacing any earlier copy, with each file's
/// modification time. Its tests, Python bindings and git metadata are left out: the library
/// build needs none of them, and the published package leaves the same parts out (see
/// `exclude` in Cargo.toml, which must match).
#[cfg(feature = "vendored")]
fn copy_source(src: &Path, dst: &Path) {
    /// Top-level entries left out, besides anything named `.git*`.
    const SKIPPED: [&str; 2] = ["test", "py_vdb"];

    fn copy_dir(src: &Path, dst: &Path, skipped: &[&str]) {
        std::fs::create_dir_all(dst).unwrap_or_else(|e| panic!("{}: {e}", dst.display()));
        let entries = std::fs::read_dir(src).unwrap_or_else(|e| panic!("{}: {e}", src.display()));
        for entry in entries {
            let entry = entry.unwrap_or_else(|e| panic!("{}: {e}", src.display()));
            let file_name = entry.file_name();
            let top_level_skip = !skipped.is_empty()
                && (skipped.iter().any(|name| file_name == *name)
                    || file_name.to_string_lossy().starts_with(".git"));
            if top_level_skip {
                continue;
            }
            let (from, to) = (entry.path(), dst.join(entry.file_name()));
            let file_type = entry.file_type().unwrap_or_else(|e| panic!("{}: {e}", from.display()));
            if file_type.is_dir() {
                copy_dir(&from, &to, &[]);
            } else {
                std::fs::copy(&from, &to)
                    .unwrap_or_else(|e| panic!("copying {}: {e}", from.display()));
                // Keep the modification time (as `fs::copy` does only on some platforms), so
                // cmake rebuilds only the files that changed since the last copy.
                let modified = entry
                    .metadata()
                    .and_then(|metadata| metadata.modified())
                    .unwrap_or_else(|e| panic!("{}: {e}", from.display()));
                std::fs::File::open(&to)
                    .and_then(|file| file.set_modified(modified))
                    .unwrap_or_else(|e| panic!("{}: {e}", to.display()));
            }
        }
    }

    if dst.exists() {
        std::fs::remove_dir_all(dst).unwrap_or_else(|e| panic!("{}: {e}", dst.display()));
    }
    copy_dir(src, dst, &SKIPPED);
}

/// Apply each `*.patch` in `patch_dir` to the ncbi-vdb source at `vdb_src`, in name order,
/// skipping any already applied. The options used mean the same to GNU patch and to macOS's
/// BSD patch.
#[cfg(feature = "vendored")]
fn apply_patches(vdb_src: &Path, patch_dir: &Path) {
    use std::process::Command;

    for patch in &patch_files(patch_dir) {
        println!("cargo:rerun-if-changed={}", patch.display());
        let run = |extra: &[&str]| {
            Command::new("patch")
                .args(["-p1", "-f", "-s", "-F0", "-d"])
                .arg(vdb_src)
                .args(extra)
                .arg("-i")
                .arg(patch)
                .output()
                .expect("failed to run patch")
        };
        // A patch that reverses cleanly is already applied.
        if run(&["-R", "--dry-run"]).status.success() {
            continue;
        }
        let applied = run(&[]);
        assert!(
            applied.status.success(),
            "failed to apply {} to {}:\n{}",
            patch.display(),
            vdb_src.display(),
            String::from_utf8_lossy(&applied.stdout)
        );
    }
}

/// Delete the bundled zlib's objects from the uber-library, so that zlib-ng provides zlib.
///
/// Members are named after zlib's sources (`inflate.c.o`, …). `crc32.c.o` and `compress.c.o`
/// are left: other bundled libraries have members of the same names, and those two zlib objects
/// are referenced only by the zlib objects removed here, so they are never linked.
#[cfg(feature = "zlib-ng")]
fn remove_bundled_zlib(vdb_src: &Path, library: &Path) {
    use std::process::Command;

    const SHARED_NAMES: [&str; 2] = ["crc32.c.o", "compress.c.o"];
    // The uber-library is a symlink to a versioned file; `ar` must rewrite the file itself.
    let library = library.canonicalize().expect("libncbi-vdb.a");
    let zlib_objects: Vec<String> = std::fs::read_dir(vdb_src.join("libs/ext/zlib"))
        .expect("libs/ext/zlib")
        .filter_map(|entry| {
            let name = entry.ok()?.file_name().into_string().ok()?;
            name.ends_with(".c").then(|| format!("{name}.o"))
        })
        .filter(|member| !SHARED_NAMES.contains(&member.as_str()))
        .collect();
    let listing = Command::new("ar").arg("t").arg(&library).output().expect("ar t");
    let members = String::from_utf8_lossy(&listing.stdout).into_owned();
    let present: Vec<&String> =
        zlib_objects.iter().filter(|o| members.lines().any(|m| m == o.as_str())).collect();
    if present.is_empty() {
        return;
    }
    let status = Command::new("ar").arg("d").arg(&library).args(&present).status().expect("ar d");
    assert!(status.success(), "failed to remove zlib objects from {}", library.display());
    let status = Command::new("ranlib").arg(&library).status().expect("ranlib");
    assert!(status.success(), "ranlib failed on {}", library.display());
}

/// Returns the OS-specific include directory for ncbi-vdb headers.
fn os_include_dir(inc_dir: &Path) -> PathBuf {
    if cfg!(target_os = "macos") {
        inc_dir.join("os/mac")
    } else if cfg!(target_os = "linux") {
        inc_dir.join("os/linux")
    } else if cfg!(target_os = "windows") {
        inc_dir.join("os/win")
    } else {
        inc_dir.join("os/linux") // fallback
    }
}

/// Generate Rust FFI bindings from the VDB C headers.
fn generate_bindings(inc_dir: &Path, out_dir: &Path) {
    let wrapper_path = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()).join("wrapper.h");

    let bindings = bindgen::Builder::default()
        .header(wrapper_path.to_str().unwrap())
        // Include path for VDB headers.
        .clang_arg(format!("-I{}", inc_dir.display()))
        // OS-specific include path.
        .clang_arg(format!("-I{}", os_include_dir(inc_dir).display()))
        // Allowlist only the functions and types we need.
        // VDB Manager
        .allowlist_function("VDBManagerMakeRead")
        .allowlist_function("VDBManagerRelease")
        .allowlist_function("VDBManagerOpenDBRead")
        .allowlist_function("VDBManagerOpenTableRead")
        .allowlist_function("VDBManagerPathType")
        .allowlist_function("VDBManagerDisablePagemapThread")
        // VFS manager / resolver: process-wide remote-access control
        .allowlist_function("VFSManagerMake")
        .allowlist_function("VFSManagerRelease")
        .allowlist_function("VFSManagerGetResolver")
        .allowlist_function("VResolverRelease")
        .allowlist_function("VResolverRemoteEnable")
        .allowlist_item("vrAlwaysDisable")
        // VDatabase
        .allowlist_function("VDatabaseRelease")
        .allowlist_function("VDatabaseOpenTableRead")
        .allowlist_function("VDatabaseListTbl")
        .allowlist_function("VDatabaseOpenMetadataRead")
        // VDBDependencies
        .allowlist_function("VDatabaseListDependencies")
        .allowlist_function("VDBDependenciesRelease")
        .allowlist_function("VDBDependenciesCount")
        .allowlist_function("VDBDependenciesSeqId")
        .allowlist_function("VDBDependenciesLocal")
        // VTable
        .allowlist_function("VTableRelease")
        .allowlist_function("VTableCreateCursorRead")
        .allowlist_function("VTableCreateCachedCursorRead")
        .allowlist_function("VTableListReadableColumns")
        .allowlist_function("VTableListPhysColumns")
        .allowlist_function("VTableOpenMetadataRead")
        // VCursor
        .allowlist_function("VCursorAddColumn")
        .allowlist_function("VCursorOpen")
        .allowlist_function("VCursorCellDataDirect")
        .allowlist_function("VCursorIdRange")
        .allowlist_function("VCursorGetBlobDirect")
        // VBlob
        .allowlist_function("VBlobIdRange")
        .allowlist_function("VBlobCellData")
        .allowlist_function("VBlobRelease")
        .allowlist_function("VCursorRelease")
        // KMetadata / KMDataNode
        .allowlist_function("KMetadataRelease")
        .allowlist_function("KMetadataOpenNodeRead")
        .allowlist_function("KMDataNodeRelease")
        .allowlist_function("KMDataNodeRead")
        .allowlist_function("KMDataNodeReadAsU64")
        .allowlist_function("KMDataNodeReadAttr")
        .allowlist_function("KMDataNodeListChildren")
        // KNamelist
        .allowlist_function("KNamelistRelease")
        .allowlist_function("KNamelistCount")
        .allowlist_function("KNamelistGet")
        // ReferenceList / ReferenceObj
        .allowlist_function("ReferenceList_MakeDatabase")
        .allowlist_function("ReferenceList_Release")
        .allowlist_function("ReferenceList_Count")
        .allowlist_function("ReferenceList_Get")
        .allowlist_function("ReferenceList_Find")
        .allowlist_function("ReferenceObj_Name")
        .allowlist_function("ReferenceObj_SeqId")
        .allowlist_function("ReferenceObj_SeqLength")
        .allowlist_function("ReferenceObj_Idx")
        .allowlist_function("ReferenceObj_IdRange")
        .allowlist_function("ReferenceObj_Read")
        .allowlist_function("ReferenceObj_Circular")
        .allowlist_function("ReferenceObj_External")
        .allowlist_function("ReferenceObj_MakePlacementIterator")
        .allowlist_function("ReferenceObj_Release")
        // AlignMgr / PlacementSetIterator
        .allowlist_function("AlignMgrMakeRead")
        .allowlist_function("AlignMgrRelease")
        .allowlist_function("AlignMgrMakePlacementSetIterator")
        .allowlist_function("PlacementSetIteratorAddPlacementIterator")
        .allowlist_function("PlacementSetIteratorNextReference")
        .allowlist_function("PlacementSetIteratorNextWindow")
        .allowlist_function("PlacementSetIteratorNextAvailPos")
        .allowlist_function("PlacementSetIteratorNextRecordAt")
        .allowlist_function("PlacementSetIteratorRelease")
        .allowlist_function("PlacementIteratorRelease")
        // Types we need.
        .allowlist_type("rc_t")
        .allowlist_type("INSDC_coord_zero")
        .allowlist_type("INSDC_coord_one")
        .allowlist_type("INSDC_coord_len")
        .allowlist_type("INSDC_coord_val")
        .allowlist_type("PlacementRecord")
        .allowlist_type("PlacementRecordExtendFuncs")
        .allowlist_type("align_id_src")
        .allowlist_type("VDBDependencies")
        // KPathType / KDBPathType constants, for interpreting VDBManagerPathType.
        .allowlist_item("kpt.*")
        // rc.h constants for error decoding.
        .allowlist_var("rcDone")
        // Derive traits.
        .derive_debug(true)
        .derive_default(true)
        // Layout tests can be noisy; disable if needed.
        .layout_tests(false)
        // Generate bindings.
        .generate()
        .expect("failed to generate FFI bindings");

    bindings.write_to_file(out_dir.join("bindings.rs")).expect("failed to write bindings");
}
