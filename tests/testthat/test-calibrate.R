test_that(".fit_calibration recovers a known linear relationship", {
  pts <- data.frame(
    target = "fake", level = 1:5,
    resp = c(100, 200, 300, 400, 500),
    conc = c(1, 2, 3, 4, 5) * 2
  )
  f <- suppressWarnings(.fit_calibration(pts))
  expect_equal(f$slope, 0.02, tolerance = 1e-8)
  expect_equal(f$slope_b, 0.02, tolerance = 1e-8)
  expect_equal(f$intercept, 0, tolerance = 1e-8)
  expect_equal(f$r2_slope_b, 1, tolerance = 1e-8)
})

test_that("packaged calibration fits are sane", {
  f <- .fit_calibration()
  expect_equal(nrow(f), length(unique(calibration_points$target)))
  expect_true(all(f$slope > 0))
  expect_true(all(f$r2_slope_b > 0.9, na.rm = TRUE))
})

test_that(".match_targets matches names and synonyms", {
  peaks <- data.frame(name = c("Cyclohexane", "cyclohexane, methyl-", "Nonexistent"))
  m <- .match_targets(peaks)
  expect_equal(m[1], "cyclohexane")
  expect_equal(m[2], "methylcyclohexane")
  expect_true(is.na(m[3]))
})

make_peaks <- function() {
  p <- data.frame(
    name = c("cyclohexane", "unknown A", "phenol", "unknown B"),
    formula = c("C6H12", "C7H16", "C6H6O", "C3H8O"),
    area = c(1e6, 5e5, 8e5, 2e5),
    stringsAsFactors = FALSE
  )
  p <- cbind(p, .parse_formula(p$formula))
  p$mw <- .molecular_weight(p$formula)
  p
}

test_that(".assign_target_rows routes peaks to the nearest preceding target of type", {
  p <- make_peaks()
  m <- .match_targets(p)
  gov <- .assign_target_rows(p, m)
  # row 2 is a hydrocarbon, governed by cyclohexane at row 1
  expect_equal(gov[2], 1L)
  # row 4 is an oxygenate, governed by phenol at row 3
  expect_equal(gov[4], 3L)
})

test_that("oxygenates fall back to hydrocarbon targets when none are calibrated", {
  p <- make_peaks()[c(1, 4), ]
  m <- .match_targets(p)
  gov <- .assign_target_rows(p, m)
  expect_equal(gov[2], 1L)
})

test_that(".polyarc_ratio scales by area, carbon number and molar mass", {
  p <- make_peaks()
  r <- fidms:::.polyarc_ratio(p, 2L, 1L)
  expected <- (p$area[2] / p$area[1]) * (p$C[1] / p$C[2]) * (p$mw[2] / p$mw[1])
  expect_equal(r, expected)
  # a peak identical to its reference gets a ratio of exactly 1
  expect_equal(fidms:::.polyarc_ratio(p, 1L, 1L), 1)
})

test_that("relative peaks are scaled off their governing target", {
  p <- make_peaks()
  out <- calibrate(p, style = "slope")
  expect_equal(out$quant_basis, c("Target", "Relative", "Target", "Relative"))
  expect_equal(
    out$amount_ug_ml[2],
    out$amount_ug_ml[1] * fidms:::.polyarc_ratio(p, 2L, 1L)
  )
})

test_that("calibrated target amounts are never negative", {
  p <- make_peaks()
  p$area <- c(1, 1, 1, 1)
  out <- calibrate(p, style = "slope_intercept")
  expect_true(all(out$amount_ug_ml >= 0, na.rm = TRUE))
})

test_that("internal standardisation normalises to 100 percent", {
  p <- make_peaks()
  out <- calibrate(p, style = "internal", reference_peak = 1, density = 0.85)
  # amounts sum to 1e6 ug per mL-equivalent scaled by density, i.e. 100 wt%
  expect_equal(sum(out$amount_ug_ml) / 10000 / 0.85, 100)
  expect_equal(out$quant_basis, c("Target", "Relative", "Relative", "Relative"))
})

test_that("the internal reference choice does not change relative proportions", {
  p <- make_peaks()
  a <- calibrate(p, style = "internal", reference_peak = 1)$amount_ug_ml
  b <- calibrate(p, style = "internal", reference_peak = 3)$amount_ug_ml
  expect_equal(a, b)
})

test_that("calibrate errors without a reference peak or targets", {
  p <- make_peaks()
  expect_error(calibrate(p, style = "internal"), "reference_peak")
  no_targets <- p[c(2, 4), ]
  expect_error(calibrate(no_targets, style = "slope"), "No calibration targets")
})

test_that(".recovery_excluded drops methanol only when reported separately", {
  expect_false(any(grepl("Methanol", .recovery_excluded(FALSE))))
  expect_true(any(grepl("Methanol", .recovery_excluded(TRUE))))
  expect_true("Acetonitrile" %in% .recovery_excluded(FALSE))
})

test_that(".calibration_recovery reports every style", {
  p <- make_peaks()
  info <- .sample_info("X", density = 0.85, methanol_wt_pct = 2, water_wt_pct = 1)
  r <- .calibration_recovery(p, info, reference_peak = 1)
  expect_equal(r$style, c("slope", "slope_intercept", "quadratic", "internal"))
  expect_equal(r$total_recovery, r$gc_recovery + 3)
  expect_equal(r$gc_recovery[r$style == "internal"], 100)
})
