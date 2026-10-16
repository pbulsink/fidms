test_that("process_sample runs end to end on the sample exports", {
  info <- sample_info(
    "Diesel-1",
    "Test sample",
    "MD",
    "pbulsink",
    density = 0.8519,
    methanol_wt_pct = 0,
    water_wt_pct = 0
  )
  # The sample excerpt contains no calibration targets, which is exactly the
  # macro's "no targets identified" branch, so it is quantified internally.
  res <- process_sample(
    test_path("fixtures", "FIDpeaks.csv"),
    test_path("fixtures", "MSDpeaks.csv"),
    test_path("fixtures", "MSDsearch.csv"),
    info = info,
    style = "internal",
    reference_peak = 5
  )
  expect_named(res$families, c("compound_type", "wt_pct"))
  expect_equal(nrow(res$families), 8)
  expect_equal(nrow(res$elements), 5)
  expect_equal(res$totals$element_sum, 100, tolerance = 1e-6)
  expect_gt(res$qc$total_peaks, 0)
  expect_equal(res$qc$aligned_peaks + res$qc$unique_peaks, res$qc$total_peaks)
  expect_false(any(is.na(res$peaks$family)))
})

test_that("external calibration errors when the sample has no targets", {
  info <- sample_info("Diesel-1", density = 0.85)
  expect_error(
    process_sample(
      test_path("fixtures", "FIDpeaks.csv"),
      test_path("fixtures", "MSDpeaks.csv"),
      test_path("fixtures", "MSDsearch.csv"),
      info = info,
      style = "slope"
    ),
    "No calibration targets"
  )
})

test_that("total_compound_table uses the macro's sentinel strings", {
  peaks <- data.frame(
    peak_label = c("1", "774 - 778"),
    name = c("Butane", "Unknown, Extrapolated Combined FID Peaks"),
    combined_peak = c(FALSE, TRUE),
    rt = c(1.2, 20.1),
    formula = c("C4H10", "C20H42"),
    mw = c(58.1, 282.5),
    area = c(100, 200),
    amount_ug_ml = c(1, 2),
    quant_basis = c("Target", "Relative"),
    score = c(99.2, NA),
    stringsAsFactors = FALSE
  )
  tbl <- total_compound_table(list(peaks = peaks))
  expect_equal(tbl[["Retention Time"]][2], "Combined FID Peaks")
  expect_type(tbl[["Peak No."]], "character")
  expect_equal(ncol(tbl), 9)
})

test_that("total_compound_table omits class-ion columns by default", {
  peaks <- data.frame(
    peak_label = c("1", "774 - 778"),
    name = c("Butane", "Unknown, Extrapolated Combined FID Peaks"),
    combined_peak = c(FALSE, TRUE),
    rt = c(1.2, 20.1),
    formula = c("C4H10", "C20H42"),
    mw = c(58.1, 282.5),
    area = c(100, 200),
    amount_ug_ml = c(1, 2),
    quant_basis = c("Target", "Relative"),
    score = c(99.2, NA),
    class_primary = c("Paraffin", NA),
    class_secondary = c(NA, NA),
    estimated_carbon_number = c(4.1, NA),
    assignment_tier = factor(c("C", NA), levels = c("A", "B", "C", "D", "M")),
    mixed_peak_flag = c(FALSE, NA),
    stringsAsFactors = FALSE
  )
  tbl <- total_compound_table(list(peaks = peaks))
  expect_equal(ncol(tbl), 9)
  expect_false("Hydrocarbon Class" %in% names(tbl))
})

test_that("total_compound_table appends class-ion columns when requested", {
  peaks <- data.frame(
    peak_label = c("1", "2"),
    name = c(NA, "Xylene"),
    combined_peak = c(FALSE, FALSE),
    rt = c(1.2, 5.4),
    formula = c(NA, "C8H10"),
    mw = c(NA, 106.2),
    area = c(100, 200),
    amount_ug_ml = c(1, 2),
    quant_basis = c("Relative", "Target"),
    score = c(NA, 90),
    class_primary = c("Paraffin", NA),
    class_secondary = c(NA, NA),
    estimated_carbon_number = c(4.14, NA),
    assignment_tier = factor(c("C", "A"), levels = c("A", "B", "C", "D", "M")),
    mixed_peak_flag = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  tbl <- total_compound_table(list(peaks = peaks), include_class_ion = TRUE)
  expect_equal(ncol(tbl), 14)
  expect_equal(tbl$`Hydrocarbon Class`, c("Paraffin", NA))
  expect_equal(tbl$`Est. Carbon No.`, c(4.1, NA))
  expect_equal(tbl$`Assignment Tier`, c("C", "A"))
})

test_that("total_compound_table degrades gracefully when class-ion columns are absent", {
  peaks <- data.frame(
    peak_label = "1",
    name = "Butane",
    combined_peak = FALSE,
    rt = 1.2,
    formula = "C4H10",
    mw = 58.1,
    area = 100,
    amount_ug_ml = 1,
    quant_basis = "Target",
    score = 99.2,
    stringsAsFactors = FALSE
  )
  tbl <- total_compound_table(list(peaks = peaks), include_class_ion = TRUE)
  expect_equal(ncol(tbl), 14)
  expect_true(is.na(tbl$`Hydrocarbon Class`))
})

test_that("excel report writes both sheets with the expected names", {
  skip_if_not_installed("openxlsx")
  info <- sample_info("Diesel-1", "Test", "MD", "pbulsink", density = 0.85)
  res <- process_sample(
    test_path("fixtures", "FIDpeaks.csv"),
    test_path("fixtures", "MSDpeaks.csv"),
    test_path("fixtures", "MSDsearch.csv"),
    info = info,
    style = "internal",
    reference_peak = 5
  )
  path <- withr::local_tempfile(fileext = ".xlsx")
  write_excel_report(res, path)
  expect_true(file.exists(path))
  expect_equal(
    openxlsx::getSheetNames(path),
    c("Diesel-1 Total Compound Report", "Diesel-1 Summary Report")
  )
})

test_that("excel report adds an MS Class Features sheet when the result has feature data", {
  skip_if_not_installed("openxlsx")
  info <- sample_info("Diesel-1", "Test", "MD", "pbulsink", density = 0.85)
  res <- process_sample(
    test_path("fixtures", "FIDpeaks.csv"),
    test_path("fixtures", "MSDpeaks.csv"),
    test_path("fixtures", "MSDsearch.csv"),
    info = info,
    style = "internal",
    reference_peak = 5
  )
  res$ms_class_features <- tibble::tibble(
    fid_peak_id = 1L,
    feature_id = 1L,
    rt_apex = 10,
    rt_start = 9.9,
    rt_end = 10.1,
    class = "Paraffin",
    class_fraction = 0.8,
    carbon_number = 8,
    residual = 0.01
  )
  path <- withr::local_tempfile(fileext = ".xlsx")
  write_excel_report(res, path)
  expect_equal(
    openxlsx::getSheetNames(path),
    c(
      "Diesel-1 Total Compound Report",
      "Diesel-1 Summary Report",
      "MS Class Features"
    )
  )
})
