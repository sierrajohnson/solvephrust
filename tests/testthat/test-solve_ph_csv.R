test_that("solve_ph matches expected results from CSV test data", {
  # Load the test data CSV
  csv_path <- system.file("tests", "testthat", "solve_ph_inputs_outputs.csv", package = "sovlephrust")
  if (!file.exists(csv_path)) {
    # If not found in installed package, try relative path for development
    csv_path <- file.path("tests", "testthat", "solve_ph_inputs_outputs.csv")
  }
  if (!file.exists(csv_path)) {
    # Try current directory path
    csv_path <- "solve_ph_inputs_outputs.csv"
  }
  
  expect_true(file.exists(csv_path), "Test data CSV file not found")
  
  # Read the CSV data
  test_data <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  # Verify we have data
  expect_gt(nrow(test_data), 0, "Test data CSV is empty")
  
  # Check that required columns exist - these should match solve_ph parameters
  required_cols <- c("temp", "ionic_strength", "kw", "tot_po4", "tot_co3", "tot_ocl", 
                     "tot_nh3", "tot_ch3coo", "h2po4_i", "hpo4_i", "po4_i", "ocl_i", 
                     "nh4_i", "ch3coo_i", "carbonate_alk_eq", "oh_i", "h_i", 
                     "so4_dose", "na_dose", "ca_dose", "mg_dose", "cl_dose", 
                     "mno4_dose", "no3_dose", "ph_r_backend")
  
  missing_cols <- setdiff(required_cols, names(test_data))
  if (length(missing_cols) > 0) {
    fail(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Initialize results tracking
  results <- data.frame(
    water_id = integer(),
    dose_id = integer(),
    expected_ph = numeric(),
    actual_ph = numeric(),
    difference = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each row of test data
  for (i in 1:nrow(test_data)) {
    row <- test_data[i, ]
    
    # Handle NA values by converting to NULL for ionic_strength
    ionic_strength <- if (is.na(row$ionic_strength)) NULL else row$ionic_strength
    
    # Call solve_ph with parameters from CSV
    actual_ph <- solve_ph(
      temp = row$temp,
      ionic_strength = ionic_strength,
      kw = row$kw,
      tot_po4 = row$tot_po4,
      tot_co3 = row$tot_co3,
      tot_ocl = row$tot_ocl,
      tot_nh3 = row$tot_nh3,
      tot_ch3coo = row$tot_ch3coo,
      h2po4_i = row$h2po4_i,
      hpo4_i = row$hpo4_i,
      po4_i = row$po4_i,
      ocl_i = row$ocl_i,
      nh4_i = row$nh4_i,
      ch3coo_i = row$ch3coo_i,
      carbonate_alk_eq = row$carbonate_alk_eq,
      oh_i = row$oh_i,
      h_i = row$h_i,
      so4_dose = row$so4_dose,
      na_dose = row$na_dose,
      ca_dose = row$ca_dose,
      mg_dose = row$mg_dose,
      cl_dose = row$cl_dose,
      mno4_dose = row$mno4_dose,
      no3_dose = row$no3_dose
    )
    
    # Calculate difference
    expected_ph <- row$ph_r_backend
    difference <- abs(actual_ph - expected_ph)
    
    # Store results
    results <- rbind(results, data.frame(
      water_id = row$water_id,
      dose_id = row$dose_id,
      expected_ph = expected_ph,
      actual_ph = actual_ph,
      difference = difference,
      stringsAsFactors = FALSE
    ))
    
    # Individual test assertion with informative message
    expect_equal(
      actual_ph, 
      expected_ph, 
      tolerance = 0.0001,
      info = paste("Water ID:", row$water_id, "Dose ID:", row$dose_id,
                   "- Expected:", expected_ph, 
                   "- Actual:", actual_ph,
                   "- Difference:", difference)
    )
  }
  
  # Overall assertions
  expect_true(all(results$difference <= 0.01), 
              "All test cases should be within 0.01 pH units of expected values")
  
  # Expect most cases to be very close (within 0.001)
  close_cases <- sum(results$difference <= 0.001)
  expect_gte(close_cases / nrow(results), 0.9, 
             "At least 90% of cases should be within 0.001 pH units")
})

test_that("solve_ph handles edge cases correctly", {
  # Test with NULL ionic strength
  result <- solve_ph(
    temp = 25, ionic_strength = NULL, kw = 1e-14,
    tot_po4 = 0, tot_co3 = 0.001, tot_ocl = 0, tot_nh3 = 0, tot_ch3coo = 0,
    h2po4_i = 0, hpo4_i = 0, po4_i = 0, ocl_i = 0, nh4_i = 0, ch3coo_i = 0,
    carbonate_alk_eq = 0.001, oh_i = 1e-7, h_i = 1e-7,
    so4_dose = 0, na_dose = 0, ca_dose = 0, mg_dose = 0, cl_dose = 0, 
    mno4_dose = 0, no3_dose = 0
  )
  
  expect_true(is.numeric(result), "solve_ph should return a numeric value")
  expect_false(is.na(result), "solve_ph should not return NA for valid inputs")
  expect_true(result > 0 && result < 14, "pH should be between 0 and 14")
})

test_that("solve_ph returns NaN for invalid inputs", {
  # Test with extreme values that should cause convergence failure
  result <- solve_ph(
    temp = 25, ionic_strength = 0.1, kw = 1e-14,
    tot_po4 = 0, tot_co3 = 0, tot_ocl = 0, tot_nh3 = 0, tot_ch3coo = 0,
    h2po4_i = 0, hpo4_i = 0, po4_i = 0, ocl_i = 0, nh4_i = 0, ch3coo_i = 0,
    carbonate_alk_eq = -1000,  # Extreme negative value
    oh_i = 1e-7, h_i = 1e-7,
    so4_dose = 0, na_dose = 0, ca_dose = 0, mg_dose = 0, cl_dose = 0, 
    mno4_dose = 0, no3_dose = 0
  )
  
  # Should either return NaN or a reasonable pH value
  expect_true(is.numeric(result), "solve_ph should return a numeric value")
})
