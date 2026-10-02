#!/usr/bin/env bash
#
# Build fg-sra-buildable-src-<version>.tar.gz from a checkout: every file git tracks, the
# ncbi-vdb submodule's included, under fg-sra-<version>/. GitHub's own tag archives leave
# submodules out and so don't build; this one does, and bioconda's recipe builds from it.
#
# Usage: scripts/buildable-src-tarball.sh <checkout> <output-dir>
#
# Prints the tarball's path. The version is the workspace's, from <checkout>/Cargo.toml.
# The tarball is reproducible: entries are owned by 0:0 and all stamped with the checkout's
# commit time. That single mtime also keeps ncbi-vdb's generated parsers no older than their
# grammars, so its build doesn't regenerate them with whatever bison is installed. Needs GNU
# tar (as `gtar` or `tar`) and a checkout with its submodules.

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <checkout> <output-dir>" >&2
    exit 2
fi
checkout=$1
output_dir=$2

tar=""
for candidate in gtar tar; do
    if command -v "$candidate" >/dev/null && "$candidate" --version 2>/dev/null | grep -q 'GNU tar'; then
        tar=$candidate
        break
    fi
done
if [ -z "$tar" ]; then
    echo "GNU tar is required (as gtar or tar)" >&2
    exit 1
fi

if git -C "$checkout" submodule status --recursive | grep -q '^-'; then
    echo "$checkout has uninitialized submodules: run git submodule update --init --recursive" >&2
    exit 1
fi

version=$(sed -n '/^\[workspace\.package\]/,/^\[/s/^version = "\(.*\)"$/\1/p' "$checkout/Cargo.toml")
if [ -z "$version" ]; then
    echo "no [workspace.package] version in $checkout/Cargo.toml" >&2
    exit 1
fi
prefix="fg-sra-$version"
mkdir -p "$output_dir"
tarball="$(cd "$output_dir" && pwd)/fg-sra-buildable-src-$version.tar.gz"
commit_time=$(git -C "$checkout" log -1 --format=%ct)

# `S` leaves symbolic link targets alone; only entry names get the prefix.
git -C "$checkout" ls-files --recurse-submodules -z |
    "$tar" -C "$checkout" --null --files-from=- --format=gnu \
        --transform="s,^,$prefix/,S" --mtime="@$commit_time" \
        --owner=0 --group=0 --numeric-owner --mode=a+rX,u+w,go-w -cf - |
    gzip -n > "$tarball"

tracked=$(git -C "$checkout" ls-files --recurse-submodules -z | tr -cd '\0' | wc -c)
archived=$("$tar" -tzf "$tarball" | wc -l)
if [ "$tracked" -ne "$archived" ]; then
    echo "$tarball holds $archived files but $checkout tracks $tracked" >&2
    exit 1
fi
if ! "$tar" -tzf "$tarball" "$prefix/crates/fg-sra-vdb-sys/vendor/ncbi-vdb/CMakeLists.txt" >/dev/null 2>&1; then
    echo "$tarball lacks the vendored ncbi-vdb" >&2
    exit 1
fi

echo "$tarball"
