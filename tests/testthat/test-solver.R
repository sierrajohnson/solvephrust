row_to_inputs <- function(row) {
  list(
    ph_dependent = list(
      po4 = list(total = row$tot_po4,
                 initial = list(row$h2po4_i, row$hpo4_i, row$po4_i)),
      co3 = list(total = row$tot_co3,
                 initial = list(row$carbonate_alk_eq, 0)),
      ocl = list(total = row$tot_ocl, initial = list(row$ocl_i)),
      nh3 = list(total = row$tot_nh3, initial = list(row$nh4_i)),
      ch3coo = list(total = row$tot_ch3coo, initial = list(row$ch3coo_i))
    ),
    ph_independent = list(
      so4 = row$so4_dose,
      na = row$na_dose,
      ca = row$ca_dose,
      mg = row$mg_dose,
      cl = row$cl_dose,
      mno4 = row$mno4_dose,
      no3 = row$no3_dose
    )
  )
}

test_that("default solver matches legacy CSV results at legacy precision", {
  csv_path <- file.path("solve_ph_inputs_outputs.csv")
  if (!file.exists(csv_path)) csv_path <- file.path("tests", "testthat", csv_path)
  test_data <- read.csv(csv_path, stringsAsFactors = FALSE)
  solver <- new_default_solver()

  actual <- vapply(seq_len(nrow(test_data)), function(index) {
    row <- test_data[index, ]
    values <- row_to_inputs(row)
    solve(solver, temp = row$temp, kw = row$kw,
          ionic_strength = if (is.na(row$ionic_strength)) NULL else row$ionic_strength,
          ph_dependent = values$ph_dependent,
          ph_independent = values$ph_independent,
          h_i = row$h_i, oh_i = row$oh_i)
  }, numeric(1))

  legacy_full_precision <- c(
    7.4999999999999, 5.58148635197709, 10.4142492489473, 11.4955114381775,
    10.4142492489473, 5.58148635197709, 6.77544140531117, 6.39542366432055,
    11.1545001133204, 9.00000052785054, 6.37314914848028, 10.3121469819271,
    11.5663027225355, 10.3121469819271, 6.37314914848028, 7.39529830969527,
    6.90291149701913, 11.0961527849585, 6, 3.02629190060529,
    8.24525138898043, 10.9401586608698, 8.24525138898043, 3.02629190060529,
    5.60371860078235, 4.42362598330093, 10.3013033359448, 6.9999999999947,
    4.47660899143716, 10.3326078973603, 11.4841044400374, 10.3326078973603,
    4.47660899143716, 6.53045546692531, 6.16450357384107, 11.1251785827266,
    7.20000000000149, 4.49465527119014, 10.4544421253929, 11.5188809461145,
    10.4544421253929, 4.49465527119014, 6.61358563096248, 6.21706520063784,
    11.1935988837206, 7.80242351414112, 4.12935115098263, 10.6250055362978,
    11.56405129558, 10.6250055362978, 4.12935115098263, 6.72749963407485,
    6.26191603104397, 11.2767502347221
  )
  expect_equal(actual, legacy_full_precision, tolerance = 1e-10)
  expect_equal(round(actual, 2), test_data$ph_r_backend, tolerance = 0.01)
  expect_true(any(actual != round(actual, 2)))
})

test_that("constructors return fresh owning Rust pointers", {
  first <- new_default_solver()
  second <- new_default_solver()
  expect_s3_class(first, "solvephrust_solver")
  expect_type(first, "externalptr")
  expect_type(second, "externalptr")
  expect_false(identical(first, second))
})

test_that("an aliased pointer keeps the Rust configuration alive", {
  solver <- new_default_solver()
  alias <- solver
  rm(solver)
  gc()
  result <- solve(
    alias, temp = 25, kw = 1e-14,
    ph_dependent = list(co3 = list(total = 0, initial = list(0, 0))),
    h_i = 1e-7, oh_i = 1e-7
  )
  expect_equal(result, 7, tolerance = 1e-6)
})

test_that("solver pointers are intentionally not serializable", {
  restored <- unserialize(serialize(new_default_solver(), NULL))
  expect_error(
    solve(
      restored, temp = 25, kw = 1e-14,
      ph_dependent = list(co3 = list(total = 0, initial = list(0, 0)))
    ),
    "non-null pointer"
  )
})

test_that("custom compounds require no backend changes", {
  solver <- new_solver(
    ph_dependent = list(custom_acid = list(
      constants = list(list(k = 1e-7, delta_h = 0)),
      charge = -1
    )),
    ph_independent_charges = list(custom_cation = 1)
  )
  result <- solve(
    solver, temp = 25, kw = 1e-14,
    ph_dependent = list(custom_acid = list(total = 1e-3, initial = list(5e-4))),
    ph_independent = list(custom_cation = 5e-4),
    h_i = 1e-7, oh_i = 1e-7
  )
  expect_type(result, "double")
  expect_true(is.finite(result))
})

test_that("multivalent positive constants are ordered by distance from neutral", {
  solver <- new_solver(ph_dependent = list(dipositive = list(
    constants = list(
      list(k = 1e-5, delta_h = 0),  # +1 -> neutral
      list(k = 1e-8, delta_h = 0)   # +2 -> +1
    ),
    charge = 1
  )))
  weights <- c(1, 1e-6 / 1e-5, (1e-6 / 1e-5) * (1e-6 / 1e-8))
  alphas <- weights / sum(weights)
  result <- solve(
    solver, temp = 25, kw = 1e-14,
    ph_dependent = list(dipositive = list(
      total = 1,
      initial = list(alphas[[2]], alphas[[3]])
    )),
    h_i = 1e-6, oh_i = 1e-8
  )
  expect_equal(result, 6, tolerance = 1e-6)
})

test_that("a solve requires at least one dependent compound", {
  solver <- new_default_solver()
  expect_error(solve(solver, temp = 25, kw = 1e-14), "at least one compound")
  expect_error(solve(solver, 25, kw = 1e-14), "supplied by name")
})

test_that("borate and silicate defaults participate in equilibrium", {
  solver <- new_default_solver()
  result <- solve(
    solver, temp = 25, kw = 1e-14,
    ph_dependent = list(
      bo3 = list(total = 1e-4, initial = list(0)),
      sio4 = list(total = 1e-4, initial = list(0, 0))
    ),
    h_i = 1e-7, oh_i = 1e-7
  )
  expect_true(is.finite(result))
})

test_that("Rust validates and accepts R scalar representations", {
  solver <- new_default_solver()
  inputs <- list(co3 = list(total = 0, initial = list(0, 0)))
  expect_true(is.finite(solve(
    solver, temp = 25L, kw = 1e-14, ionic_strength = 0L,
    ph_dependent = inputs, h_i = 1e-7, oh_i = 1e-7
  )))
  expect_error(solve(
    solver, temp = NA_real_, kw = 1e-14, ph_dependent = inputs
  ), "temp")
  expect_error(solve(
    solver, temp = 25, kw = 1e-14, ionic_strength = -0.1,
    ph_dependent = inputs
  ), "ionic_strength")
})

test_that("constructor reports invalid schemas", {
  step <- list(k = 1e-7, delta_h = 0)
  expect_error(new_solver(ph_dependent = list(list(constants = list(step), charge = -1))),
               "non-empty name")
  expect_error(new_solver(ph_dependent = list(x = list(constants = list(), charge = -1))),
               "non-empty list")
  expect_error(new_solver(ph_dependent = list(x = list(constants = list(step), charge = 2))),
               "either -1 or 1")
  expect_error(new_solver(ph_independent_charges = list(x = 0)), "nonzero integer")
  expect_error(new_solver(ph_dependent = list(x = list(
    constants = rep(list(step), 4), charge = -1
  ))), "one to three")
  expect_error(new_solver(ph_dependent = list(x = list(constants = list(step), charge = -1)),
                          ph_independent_charges = list(x = -1)), "cannot be both")
})

test_that("solve reports invalid runtime values", {
  solver <- new_default_solver()
  expect_error(solve(solver, temp = 25, kw = 1e-14,
                     ph_dependent = list(unknown = list(total = 0, initial = list(0)))),
               "Unknown compound")
  expect_error(solve(solver, temp = 25, kw = 1e-14,
                     ph_dependent = list(co3 = list(total = 0, initial = list(0)))),
               "length 2")
  expect_error(solve(solver, temp = 25, kw = 1e-14,
                     ph_dependent = list(co3 = list(total = 0, initial = list(0, 0))),
                     ph_independent = list(na = -1)), ">= 0")
  expect_error(solve(solver, temp = -273.15, kw = 1e-14,
                     ph_dependent = list(co3 = list(total = 0, initial = list(0, 0)))), "temp")
  expect_error(solve(solver, temp = 25, kw = 0,
                     ph_dependent = list(co3 = list(total = 0, initial = list(0, 0)))), "kw")
  expect_error(solve(solver, temp = 25, kw = 1e-14,
                     ph_dependent = list(co3 = list(total = 0, initial = list(0, 0))),
                     ph_independent = list(na = 1000)),
               "failed to converge")
})
