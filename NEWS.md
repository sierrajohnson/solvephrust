# solvephrust 0.1.0

* Initial CRAN submission.
* Configurable pH solver: `new_solver()` builds a solver from any set of
  pH-dependent and fixed-charge compounds, and `new_default_solver()` provides
  the built-in water chemistry catalog.
* `solve()` method for `solvephrust_solver` objects computes equilibrium pH
  with a safeguarded Newton solver implemented in 'Rust'.
