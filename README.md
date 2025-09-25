# Solve Ph Rust

This repository implements a solver to determine water Ph equilibrium.

The actual solver is implemented in Rust as is ~50x faster than the equivalent R implementation using uniroot.

This package is designed to work with [tidywater](https://github.com/BrownAndCaldwell-Public/tidywater), but can be used independently.

This is an isolated package to prevent developers on tidywater from requiring a local Rust installation. 
