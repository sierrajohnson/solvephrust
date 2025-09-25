# Solve Ph Rust

<!-- badges: start -->
[![R-CMD-check](https://github.com/sierrajohnson/solvephrust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sierrajohnson/solvephrust/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This repository implements a solver to determine water Ph equilibrium.

The actual solver is implemented in Rust as is ~50x faster than the equivalent R implementation using uniroot.

This package is designed to work with [tidywater](https://github.com/BrownAndCaldwell-Public/tidywater), but can be used independently.

This is an isolated package to prevent developers on tidywater from requiring a local Rust installation. 
