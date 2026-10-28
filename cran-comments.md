## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Rust

* The package compiles 'Rust' code through 'extendr'. All crate dependencies
  are vendored in `src/rust/vendor.tar.xz`, and CRAN builds run offline with
  `cargo build -j 2 --offline`.
* The cargo and rustc versions are printed during installation by
  `tools/msrv.R`; the minimum supported Rust version is 1.65.
* Authors and licenses of the vendored crates are listed in `inst/AUTHORS`.

## Test environments

<!-- TODO: fill in after running the checks in CRAN-checklist.md -->
* local Windows 11, R 4.6.1
* win-builder (devel)
* macOS builder (release)
* R-hub (linux)
