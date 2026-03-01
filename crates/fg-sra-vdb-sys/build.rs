//! Build script for fg-sra-vdb-sys.
//!
//! Builds ncbi-vdb from the git submodule via cmake, then generates
//! Rust FFI bindings via bindgen.
//!
//! Supports env var override for pre-built VDB:
//!   VDB_INCDIR - path to ncbi-vdb interfaces/ directory
//!   VDB_LIBDIR - path to directory containing libncbi-vdb.a

use std::env;
use std::path::{Path, PathBuf};

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    // Determine include and library paths, either from env vars or by building.
    let (inc_dir, lib_dir) = match (env::var("VDB_INCDIR"), env::var("VDB_LIBDIR")) {
        (Ok(inc), Ok(lib)) => {
            println!("cargo:warning=Using pre-built VDB: inc={inc}, lib={lib}");
            (PathBuf::from(inc), PathBuf::from(lib))
        }
        _ => build_ncbi_vdb(&out_dir),
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

    // Additional system libraries.
    println!("cargo:rustc-link-lib=dylib=z");

    // Generate FFI bindings via bindgen.
    generate_bindings(&inc_dir, &out_dir);

    // Rerun if the wrapper header changes.
    println!("cargo:rerun-if-changed=wrapper.h");
    println!("cargo:rerun-if-env-changed=VDB_INCDIR");
    println!("cargo:rerun-if-env-changed=VDB_LIBDIR");
}

/// Build ncbi-vdb from the vendored submodule using cmake.
fn build_ncbi_vdb(_out_dir: &Path) -> (PathBuf, PathBuf) {
    let vdb_src = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap())
        .join("../../vendor/ncbi-vdb")
        .canonicalize()
        .expect("vendor/ncbi-vdb not found; did you initialize the git submodule?");

    let dst = cmake::Config::new(&vdb_src)
        .define("LIBS_ONLY", "ON")
        .build_target("ncbi-vdb")
        .build();

    let inc_dir = vdb_src.join("interfaces");
    let build_dir = dst.join("build");

    // cmake puts the static uber-library in lib/, and helper libs in ilib/.
    let lib_dir = build_dir.join("lib");
    let ilib_dir = build_dir.join("ilib");

    // Determine which directory contains the uber-library.
    let final_lib_dir = if lib_dir.join("libncbi-vdb.a").exists() {
        lib_dir
    } else if ilib_dir.join("libncbi-vdb.a").exists() {
        ilib_dir.clone()
    } else {
        panic!(
            "libncbi-vdb.a not found after cmake build. Searched:\n  {}\n  {}",
            lib_dir.join("libncbi-vdb.a").display(),
            ilib_dir.join("libncbi-vdb.a").display()
        );
    };

    // mbedcrypto is built as a separate static lib in ilib/.
    if ilib_dir.join("libmbedcrypto.a").exists() {
        println!("cargo:rustc-link-search=native={}", ilib_dir.display());
        println!("cargo:rustc-link-lib=static=mbedcrypto");
    }

    // Tell cargo to rebuild if the ncbi-vdb source changes.
    println!("cargo:rerun-if-changed={}", vdb_src.display());

    (inc_dir, final_lib_dir)
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
        .allowlist_function("VDBManagerDisablePagemapThread")
        // VDatabase
        .allowlist_function("VDatabaseRelease")
        .allowlist_function("VDatabaseOpenTableRead")
        .allowlist_function("VDatabaseListTbl")
        .allowlist_function("VDatabaseOpenMetadataRead")
        // VTable
        .allowlist_function("VTableRelease")
        .allowlist_function("VTableCreateCursorRead")
        .allowlist_function("VTableCreateCachedCursorRead")
        .allowlist_function("VTableListReadableColumns")
        // VCursor
        .allowlist_function("VCursorAddColumn")
        .allowlist_function("VCursorOpen")
        .allowlist_function("VCursorCellDataDirect")
        .allowlist_function("VCursorIdRange")
        .allowlist_function("VCursorRelease")
        // KMetadata / KMDataNode
        .allowlist_function("KMetadataRelease")
        .allowlist_function("KMetadataOpenNodeRead")
        .allowlist_function("KMDataNodeRelease")
        .allowlist_function("KMDataNodeRead")
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

    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("failed to write bindings");
}
