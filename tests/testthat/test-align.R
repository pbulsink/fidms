test_that(".greedy_match pairs nearest within tolerance, one-to-one", {
  m <- fidms:::.greedy_match(c(1, 2, 3), c(1.01, 2.02, 9), tol = 0.05)
  expect_equal(m, c(1L, 2L, NA_integer_))
})

test_that(".find_rt_shift recovers a known synthetic shift", {
  set.seed(4711)
  fid_rt <- sort(runif(60, 1, 30))
  true_shift <- -0.02
  msd_rt <- fid_rt[1:50] + true_shift
  s <- .find_rt_shift(
    fid_rt,
    msd_rt,
    fidms_params(expected_shift = 0, shift_window = 0.1, shift_step = 0.001)
  )
  expect_equal(s$shift, true_shift, tolerance = 0.003)
})

test_that(".align_peaks keeps one row per FID peak and drops MSD-only peaks", {
  fid <- .read_fid_peaks(test_path("fixtures", "FIDpeaks.csv"))
  msd <- .read_msd_peaks(test_path("fixtures", "MSDpeaks.csv"))
  a <- .align_peaks(fid, msd)
  expect_equal(nrow(a$combined), nrow(fid))
  expect_lte(sum(!is.na(a$combined$msd_peak)), nrow(msd))
  expect_true(abs(a$shift) < 0.3)
})

test_that(".merge_close_peaks labels merged rows with a range", {
  fid <- .read_fid_peaks(test_path("fixtures", "FIDpeaks.csv"))
  msd <- .read_msd_peaks(test_path("fixtures", "MSDpeaks.csv"))
  search <- .read_msd_search(test_path("fixtures", "MSDsearch.csv"))
  p <- .merge_close_peaks(.join_msd_identity(
    .align_peaks(fid, msd)$combined,
    search
  ))
  expect_lte(nrow(p), nrow(fid))
  expect_true(all(grepl(" - ", p$peak_label[p$combined_peak])))
  expect_equal(sum(p$area), sum(fid$area), tolerance = 1e-6)
})
