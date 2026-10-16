test_that("the GCalignR engine returns the same structure as the grid engine", {
  skip_if_not_installed("GCalignR")
  fid <- .read_fid_peaks(test_path("fixtures", "FIDpeaks.csv"))
  msd <- .read_msd_peaks(test_path("fixtures", "MSDpeaks.csv"))

  g <- .align_peaks(fid, msd, engine = "gcalignr")
  expect_equal(nrow(g$combined), nrow(fid))
  expect_true(all(c("msd_peak", "msd_rt", "msd_area") %in% names(g$combined)))
  expect_lte(sum(!is.na(g$combined$msd_peak)), nrow(msd))
  expect_false(is.null(g$gcalignr))
})

test_that("both engines match a comparable number of peaks", {
  skip_if_not_installed("GCalignR")
  fid <- .read_fid_peaks(test_path("fixtures", "FIDpeaks.csv"))
  msd <- .read_msd_peaks(test_path("fixtures", "MSDpeaks.csv"))

  n_grid <- sum(
    !is.na(.align_peaks(fid, msd, engine = "grid")$combined$msd_peak)
  )
  n_gca <- sum(
    !is.na(.align_peaks(fid, msd, engine = "gcalignr")$combined$msd_peak)
  )
  # The two algorithms are not expected to agree exactly; this is a sanity
  # bound, not a parity check.
  expect_lte(abs(n_grid - n_gca), max(2, 0.3 * n_grid))
})

test_that("process_sample runs with the GCalignR engine", {
  skip_if_not_installed("GCalignR")
  info <- sample_info("Diesel-1", "Test", "MD", "pbulsink", density = 0.85)
  res <- process_sample(
    test_path("fixtures", "FIDpeaks.csv"),
    test_path("fixtures", "MSDpeaks.csv"),
    test_path("fixtures", "MSDsearch.csv"),
    info = info,
    engine = "gcalignr",
    style = "internal",
    reference_peak = 5
  )
  expect_gt(res$qc$total_peaks, 0)
  expect_equal(nrow(res$families), 8)
})
