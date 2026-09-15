# solvephrust

<!-- badges: start -->
[![R-CMD-check](https://github.com/sierrajohnson/solvephrust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sierrajohnson/solvephrust/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This repository implements a configurable solver for water pH equilibrium.

The actual solver is implemented in Rust as is ~50x faster than the equivalent R implementation using uniroot.

This package is designed to work with [tidywater](https://github.com/BrownAndCaldwell-Public/tidywater), but can be used independently.

This is an isolated package to prevent developers on tidywater from requiring a local Rust installation.

```r
solver <- new_default_solver()

solve(
  solver,
  temp = 25,
  kw = 1e-14,
  ionic_strength = 0.01,
  ph_dependent = list(
    co3 = list(total = 0.0025, initial = list(0.0024, 0))
  )
)
```

Use `new_solver()` to configure any set of pH-dependent equilibrium compounds
and fixed-charge compounds. The package defaults are generated as internal data
from `data-raw/compile_data.R`; after editing that source file, regenerate them
from the package root with:

```sh
Rscript data-raw/compile_data.R
```

Solver objects are owning external pointers to parsed Rust configurations. They
avoid reparsing the compound catalog on each solve, but cannot be serialized and
restored across R sessions.

The default solver uses lowercase chemical-formula names: `co3`, `po4`, `ocl`,
`nh3`, `ch3coo`, `bo3`, and `sio4` for pH-dependent compounds, and `so4`, `na`,
`ca`, `mg`, `cl`, `mno4`, and `no3` for fixed-charge compounds.

Dissociation constants are ordered by distance from neutral charge: `K1`
describes the transition adjacent to neutral, followed by `K2` and `K3` at
increasing absolute charge for both negative and positive compounds.

## Benchmark

The dependency-free benchmark exercises the public R API with a cached default
solver and every semi-real-world case in
`tests/testthat/solve_ph_inputs_outputs.csv`. After installing the package, run
it from the package root:

```sh
Rscript tools/benchmark.R
```

Pass an optional number of complete dataset passes for shorter or longer runs,
for example `Rscript tools/benchmark.R 100`. The output reports case and solve
counts, elapsed time, solves per second, and average time per solve. A second
argument of `backend` bypasses S3 dispatch and the R method checks while using
the same Rust FFI, validation, and numerical core:

```sh
Rscript tools/benchmark.R 1000 backend
```

Comparing `public` (the default) with `backend` separates the small R dispatch
cost from the FFI-plus-solver cost without changing the workload.
