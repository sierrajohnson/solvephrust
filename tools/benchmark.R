# Run from the package root after installing the package:
#   R CMD INSTALL .
#   Rscript tools/benchmark.R [passes]

suppressPackageStartupMessages(library(solvephrust))

args <- commandArgs(trailingOnly = TRUE)
passes <- if (length(args)) as.integer(args[[1]]) else 1000L
if (length(passes) != 1L || is.na(passes) || passes < 1L) {
  stop("`passes` must be one positive integer.", call. = FALSE)
}

csv_path <- file.path("tests", "testthat", "solve_ph_inputs_outputs.csv")
if (!file.exists(csv_path)) {
  stop("Run this benchmark from the package root; could not find `", csv_path, "`.",
       call. = FALSE)
}

test_data <- read.csv(csv_path, stringsAsFactors = FALSE)
solver <- new_default_solver()

# Convert data frames and compound lists outside the timed section so the
# benchmark measures the public solver API rather than CSV handling.
cases <- lapply(seq_len(nrow(test_data)), function(index) {
  row <- test_data[index, ]
  list(
    temp = row$temp,
    kw = row$kw,
    ionic_strength = if (is.na(row$ionic_strength)) NULL else row$ionic_strength,
    ph_dependent = list(
      po4 = list(
        total = row$tot_po4,
        initial = list(row$h2po4_i, row$hpo4_i, row$po4_i)
      ),
      co3 = list(
        total = row$tot_co3,
        initial = list(row$carbonate_alk_eq, 0)
      ),
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
    ),
    h_i = row$h_i,
    oh_i = row$oh_i
  )
})

run_case <- function(case) {
  solve(
    solver,
    temp = case$temp,
    kw = case$kw,
    ionic_strength = case$ionic_strength,
    ph_dependent = case$ph_dependent,
    ph_independent = case$ph_independent,
    h_i = case$h_i,
    oh_i = case$oh_i
  )
}

run_pass <- function() {
  vapply(cases, run_case, numeric(1))
}

# Warm the package and verify the cases before starting the clock.
results <- run_pass()
if (!isTRUE(all.equal(round(results, 2), test_data$ph_r_backend, tolerance = 0.01))) {
  stop("Benchmark results do not match the CSV's expected pH values.", call. = FALSE)
}

timing <- system.time({
  for (pass in seq_len(passes)) results <- run_pass()
})
elapsed <- unname(timing[["elapsed"]])
total_solves <- passes * length(cases)

cat(sprintf("cases per pass: %d\n", length(cases)))
cat(sprintf("passes:         %d\n", passes))
cat(sprintf("total solves:   %d\n", total_solves))
cat(sprintf("elapsed:        %.3f s\n", elapsed))
cat(sprintf("throughput:     %.0f solves/s\n", total_solves / elapsed))
cat(sprintf("time per solve: %.2f us\n", elapsed * 1e6 / total_solves))
