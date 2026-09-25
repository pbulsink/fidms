test_that(".parse_formula counts elements", {
  out <- .parse_formula(c("C6H12", "C2H6O", "CH3Cl", "C6H5NO2"))
  expect_equal(out$C, c(6L, 2L, 1L, 6L))
  expect_equal(out$H, c(12L, 6L, 3L, 5L))
  expect_equal(out$O, c(0L, 1L, 0L, 2L))
  expect_equal(out$Cl, c(0L, 0L, 1L, 0L))
  expect_equal(out$N, c(0L, 0L, 0L, 1L))
})

test_that(".parse_formula expands parentheses and handles NA", {
  expect_equal(.parse_formula("C(CH3)3")$H, 9L)
  expect_true(is.na(.parse_formula(NA_character_)$C))
})

test_that(".molecular_weight matches hand calculations", {
  expect_equal(.molecular_weight("C6H12"), 84.162, tolerance = 1e-3)
  expect_equal(.molecular_weight("H2O"), 18.015, tolerance = 1e-3)
  expect_equal(.molecular_weight("C2H6O"), 46.069, tolerance = 1e-3)
})

test_that("retention index model reproduces its anchors", {
  m <- .fit_retention_index()
  ref <- rical_reference_points
  pred <- .compute_retention_index(ref$time, m)
  expect_equal(cor(pred, ref$ri), 1, tolerance = 0.01)
  # C6 alkane anchor should land near 600
  expect_equal(.compute_retention_index(2.209, m), 600, tolerance = 40)
})
