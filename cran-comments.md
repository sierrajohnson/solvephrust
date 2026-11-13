## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* "enthalpies" in DESCRIPTION is flagged as possibly misspelled; it is the
  correct plural of "enthalpy" (reaction enthalpies used for temperature
  correction of equilibrium constants).

## Rust

* The package compiles 'Rust' code through 'extendr'. All crate dependencies
  are vendored in `src/rust/vendor.tar.xz`, and CRAN builds run offline with
  `cargo build -j 2 --offline`.
* The cargo and rustc versions are printed during installation by
  `tools/msrv.R`; the minimum supported Rust version is 1.65.
* Authors and licenses of the vendored crates are listed in `inst/AUTHORS`.

## Test environments

* local Windows 11 x64, R 4.6.2: 0 errors | 0 warnings | 0 notes
* win-builder, Windows Server 2022 x64, R-devel (2026-11-12 r90881 ucrt):
  0 errors | 0 warnings | 1 note (new submission; "enthalpies" spelling)
* macOS builder, macOS 26.6 arm64 (Apple M1), R 4.6.1: 0 errors | 0 warnings |
  0 notes
* GitHub Actions, Ubuntu 26.04 LTS, R 4.6.2 (release) and R 4.5.3 (oldrel-1):
  0 errors | 0 warnings | 0 notes
