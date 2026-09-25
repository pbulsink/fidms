fixture <- function(name) testthat::test_path("fixtures", name)

test_that("FID peaks parse, drop the banner and sort by RT", {
  fid <- .read_fid_peaks(fixture("FIDpeaks.csv"))
  expect_named(fid, c("peak", "rt_start", "rt", "rt_end", "area", "height"))
  expect_gt(nrow(fid), 5)
  expect_false(any(is.na(fid$rt)))
  expect_true(!is.unsorted(fid$rt))
  expect_type(fid$peak, "integer")
})

test_that("MSD peaks parse with the blank Name column", {
  msd <- .read_msd_peaks(fixture("MSDpeaks.csv"))
  expect_true(all(c("peak", "name", "rt", "area") %in% names(msd)))
  expect_true(all(is.na(msd$name)))
  expect_gt(nrow(msd), 5)
})

test_that("MSD search parses and drops padded blank rows", {
  s <- .read_msd_search(fixture("MSDsearch.csv"))
  expect_true(all(nzchar(s$formula)))
  expect_false(any(is.na(s$rt)))
  expect_true(all(s$score > 0, na.rm = TRUE))
})

test_that("missing columns error clearly", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c("#banner,,", "#,,", "Peak,Area,Height", "1,10,5"), tmp)
  expect_error(.read_fid_peaks(tmp), "missing required column")
})
