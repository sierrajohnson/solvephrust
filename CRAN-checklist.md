# CRAN submission checklist

Remaining steps before submitting solvephrust 0.1.0 to CRAN. Based on the
"Releasing to CRAN" chapter of [R Packages (2e)](https://r-pkgs.org/release.html)
and CRAN's [Using Rust in CRAN packages](https://cran.r-project.org/web/packages/using_rust.html)
policy. This file is excluded from the package build.

## Done (cran-prep branch)

- [x] Real maintainer and contributors in `Authors@R`; version `0.1.0`
- [x] `URL` and `BugReports` fields
- [x] Software names quoted in `Description` (`'Rust'`) and a reference for
      the default constants
- [x] Vendored crate authors and licenses in `inst/AUTHORS`, credited as `cph`
      in `Authors@R`
- [x] `@examples` for every exported function
- [x] `.Rbuildignore` entries for `data-raw/`, `.posit/`, `tools/benchmark.R`,
      and release files
- [x] `NEWS.md` and `cran-comments.md`

## To confirm

- [ ] Carter's family name (and optionally email/ORCID) in `Authors@R`;
      currently `person("Carter", role = "ctb")`
- [ ] Maintainer email: CRAN emails it for the life of the package, so it must
      be an address you will keep monitoring
- [ ] ISBN for Benjamin (2015) in `Description` (`9781478623083`, *Water
      Chemistry*, 2nd ed.)
- [ ] `LICENSE` year and copyright holder (`2025`, `solvephrust authors`)

## Local checks

- [ ] Install LaTeX so the PDF manual check runs, e.g. `tinytex::install_tinytex()`
- [ ] Install pandoc (or run from RStudio/Positron) so `README.md` and
      `NEWS.md` are checked
- [ ] `devtools::check(remote = TRUE, manual = TRUE)` with 0 errors, 0 warnings,
      and only the "New submission" note
- [ ] `urlchecker::url_check()`
- [ ] `spelling::spell_check_package()` (install `spelling` first); add a
      `inst/WORDLIST` for chemistry terms
- [ ] Build from a clean checkout so stale `src/*.o` and `src/*.dll` files are
      not picked up

## Other platforms

- [ ] `devtools::check_win_devel()`
- [ ] `devtools::check_win_release()`
- [ ] `devtools::check_mac_release()`
- [ ] R-hub Linux checks (`rhub::rhub_check()`), including a Fedora or Debian
      platform with the distribution's `rustc`
- [ ] Record results in `cran-comments.md` under "Test environments"

## Submit

- [ ] Optionally run `usethis::use_release_issue()` for the full checklist
- [ ] `devtools::submit_cran()` and confirm the email CRAN sends
- [ ] After acceptance: tag the release (`usethis::use_github_release()`) and
      bump to the development version (`usethis::use_dev_version()`)
