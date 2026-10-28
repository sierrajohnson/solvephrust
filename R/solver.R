#' Create a configurable pH solver
#'
#' @param ph_dependent A named list of equilibrium compound definitions. Each
#'   definition contains `constants`, an ordered list of one to
#'   three `k` and `delta_h` pairs, and `charge`, either `-1` or `1`. Constants
#'   use Benjamin's distance-from-neutral ordering: `K1` is adjacent to neutral,
#'   followed by `K2` and `K3` at increasing absolute charge.
#' @param ph_independent_charges A named list mapping fixed-charge compound names
#'   directly to nonzero integer charges.
#' @return A `solvephrust_solver` external pointer for use with [solve()]. The
#'   object owns its parsed Rust configuration and is intentionally not
#'   serializable across R sessions.
#' @examples
#' # A monoprotic weak acid (pKa 7) balanced by a fixed-charge cation
#' solver <- new_solver(
#'   ph_dependent = list(
#'     weak_acid = list(
#'       constants = list(list(k = 1e-7, delta_h = 0)),
#'       charge = -1
#'     )
#'   ),
#'   ph_independent_charges = list(cation = 1)
#' )
#'
#' solve(
#'   solver,
#'   temp = 25,
#'   kw = 1e-14,
#'   ph_dependent = list(weak_acid = list(total = 1e-3, initial = list(5e-4))),
#'   ph_independent = list(cation = 5e-4),
#'   h_i = 1e-7,
#'   oh_i = 1e-7
#' )
#' @export
new_solver <- function(ph_dependent = list(), ph_independent_charges = list()) {
  pointer <- create_solver(ph_dependent, ph_independent_charges)
  structure(pointer, class = "solvephrust_solver")
}

#' Create a solver with solvephrust's default compounds
#'
#' The catalog is generated from `data-raw/compile_data.R` and contains the
#' pH-dependent compounds `co3`, `po4`, `ocl`, `nh3`, `ch3coo`, `bo3`, and
#' `sio4`, and the fixed-charge compounds `so4`, `na`, `ca`, `mg`, `cl`, `mno4`,
#' and `no3`.
#' @return A new `solvephrust_solver` external pointer.
#' @examples
#' solver <- new_default_solver()
#'
#' # Water with 2.5 mM total carbonate, nearly all as bicarbonate
#' solve(
#'   solver,
#'   temp = 25,
#'   kw = 1e-14,
#'   ionic_strength = 0.01,
#'   ph_dependent = list(
#'     co3 = list(total = 0.0025, initial = list(0.0024, 0))
#'   )
#' )
#' @export
new_default_solver <- function() {
  new_solver(
    .solvephrust_default_compounds$ph_dependent,
    .solvephrust_default_compounds$ph_independent_charges
  )
}

#' Solve for equilibrium pH
#'
#' `temp` must be supplied by name. The `b` and `...` formals exist only because
#' this is a method for base R's [solve()] generic; positional or unused values
#' supplied through them are rejected.
#'
#' @param a A solver created by [new_solver()] or [new_default_solver()].
#' @param b Must be missing. Positional temperature arguments are not supported.
#' @param ... Must be empty; required by the base `solve()` generic.
#' @param temp Temperature in degrees Celsius, supplied by name.
#' @param kw Water dissociation constant.
#' @param ionic_strength Ionic strength in molar units, or `NULL`.
#' @param ph_dependent Named runtime values for configured equilibrium compounds.
#'   Each entry contains `total` and an `initial` list ordered by absolute charge.
#'   At least one compound is required.
#' @param ph_independent Named runtime values for configured fixed-charge compounds.
#'   Each value is that compound's dose.
#' @param h_i Initial hydrogen-ion concentration.
#' @param oh_i Initial hydroxide-ion concentration.
#' @return The equilibrium pH as a full-precision numeric scalar.
#' @examples
#' solver <- new_default_solver()
#'
#' # Carbonate and phosphate buffered water dosed with sodium and chloride
#' solve(
#'   solver,
#'   temp = 20,
#'   kw = 1e-14,
#'   ionic_strength = 0.005,
#'   ph_dependent = list(
#'     co3 = list(total = 2e-3, initial = list(1.9e-3, 0)),
#'     po4 = list(total = 1e-4, initial = list(5e-5, 5e-5, 0))
#'   ),
#'   ph_independent = list(na = 1e-3, cl = 1e-3),
#'   h_i = 1e-7,
#'   oh_i = 1e-7
#' )
#'
#' # A solver can be reused across many solves without reparsing its catalog
#' vapply(c(5, 15, 25), function(temp) {
#'   solve(
#'     solver,
#'     temp = temp,
#'     kw = 1e-14,
#'     ph_dependent = list(co3 = list(total = 2e-3, initial = list(1.9e-3, 0)))
#'   )
#' }, numeric(1))
#' @export
solve.solvephrust_solver <- function(a, b, ..., temp, kw,
                                     ionic_strength = NULL,
                                     ph_dependent = list(),
                                     ph_independent = list(),
                                     h_i = 0, oh_i = 0) {
  if (!missing(b)) {
    stop("`temp` must be supplied by name.", call. = FALSE)
  }
  if (...length()) {
    stop("Unused positional or named arguments are not supported.", call. = FALSE)
  }
  .Call(
    wrap__solve_generic,
    a, temp, ionic_strength, kw, ph_dependent, ph_independent, h_i, oh_i
  )
}
