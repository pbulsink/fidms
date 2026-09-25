test_that("interpolation fills every unidentified peak", {
  peaks <- data.frame(
    peak = 1:6,
    peak_label = as.character(1:6),
    combined_peak = FALSE,
    rt = c(1, 2, 3, 4, 5, 6),
    area = 100,
    ri = c(600, 650, 700, 750, 800, 850),
    name = c("Hexane", NA, NA, "Octane", NA, NA),
    formula = c("C6H14", NA, NA, "C8H18", NA, NA),
    identified = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
    mw = NA_real_,
    stringsAsFactors = FALSE
  )
  peaks <- cbind(peaks, .parse_formula(peaks$formula))
  out <- .interpolate_unknowns(peaks)
  expect_false(any(is.na(out$formula)))
  expect_equal(
    out$origin,
    c(
      "msd",
      "interpolated",
      "interpolated",
      "msd",
      "extrapolated",
      "extrapolated"
    )
  )
  # interior carbon counts sit between the bracketing identifications
  expect_true(all(out$C[2:3] >= 6 & out$C[2:3] <= 8))
  expect_equal(out$name[2], "Unknown, Interpolated")
  expect_equal(out$name[5], "Unknown, Extrapolated")
  # mw of the identified rows comes from .join_msd_identity; here only the
  # newly-filled rows should gain a mass
  expect_false(any(is.na(out$mw[!peaks$identified])))
})

test_that(".build_formula omits zero counts and the implicit 1", {
  expect_equal(fidms:::.build_formula(1, 4, 0), "CH4")
  expect_equal(fidms:::.build_formula(2, 6, 1), "C2H6O")
})
