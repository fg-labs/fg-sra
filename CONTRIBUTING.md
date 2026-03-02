# Contributing to fg-sra

## Development Setup

### Prerequisites

- Rust (stable toolchain, see `rust-toolchain.toml`)
- CMake (for building the vendored ncbi-vdb C library)

### Building

```bash
git clone --recurse-submodules https://github.com/fulcrumgenomics/fg-sra
cd fg-sra
cargo build --release
```

### Install Git Hooks

We use pre-commit hooks to ensure code quality. Install them after cloning:

```bash
./scripts/install-hooks
```

This installs hooks that run before each commit:
- `cargo ci-fmt` - Check code formatting
- `cargo ci-lint` - Run clippy lints

### Running Checks Manually

```bash
# Format check (fails if formatting differs)
cargo ci-fmt

# Lint check (fails on any warnings)
cargo ci-lint

# Run all tests
cargo ci-test
```

## Code Style

- Run `cargo fmt` before committing
- Fix all clippy warnings
- Add backticks around identifiers in doc comments (e.g., `` `read_name` ``)

## Testing

All new features should include tests. Run the full test suite with:

```bash
cargo ci-test
```

Integration tests that require network access to resolve SRA accessions are
marked `#[ignore]` and can be run with:

```bash
cargo test -- --ignored
```

## Pull Requests

1. Ensure all CI checks pass (`cargo ci-fmt`, `cargo ci-lint`, `cargo ci-test`)
2. Keep PRs focused and reasonably sized
3. Include tests for new functionality
