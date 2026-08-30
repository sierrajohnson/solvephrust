.rust_value <- function(result) {
  if (!is.null(result$error)) {
    stop(result$error, call. = FALSE)
  }
  result$value
}

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
#' @export
new_solver <- function(ph_dependent = list(), ph_independent_charges = list()) {
  pointer <- .rust_value(create_solver(ph_dependent, ph_independent_charges))
  structure(pointer, class = "solvephrust_solver")
}

#' Create a solver with solvephrust's default compounds
#'
#' The catalog is generated from `data-raw/compile_data.R` and contains the
#' pH-dependent compounds `co3`, `po4`, `ocl`, `nh3`, `ch3coo`, `bo3`, and
#' `sio4`, and the fixed-charge compounds `so4`, `na`, `ca`, `mg`, `cl`, `mno4`,
#' and `no3`.
#' @return A new `solvephrust_solver` external pointer.
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
#' @export
solve.solvephrust_solver <- function(a, b, ..., temp, kw,
                                     ionic_strength = NULL,
                                     ph_dependent = list(),
                                     ph_independent = list(),
                                     h_i = 0, oh_i = 0) {
  if (!missing(b)) {
    stop("`temp` must be supplied by name.", call. = FALSE)
  }
  if (length(list(...))) {
    stop("Unused positional or named arguments are not supported.", call. = FALSE)
  }
  .rust_value(solve_generic(
    solver = a,
    temp = temp,
    ionic_strength = ionic_strength,
    kw = kw,
    dependent_compounds = ph_dependent,
    independent_compounds = ph_independent,
    h_i = h_i,
    oh_i = oh_i
  ))
}
