make_peaks <- function() {
  p <- data.frame(
    peak = 1:6,
    peak_label = c("1", "2", "3", "4", "5", "6 - 8"),
    name = c(
      "Hexane", "Toluene", "Unknown, Interpolated", "Octane",
      "Methane, dichloro-", "Unknown, Extrapolated Combined FID Peaks"
    ),
    formula = c("C6H14", "C7H8", "C3H8", "C8H18", "CH2Cl2", "C20H42"),
    area = c(1000, 10000, 50000, 4000, 30000, 2000),
    score = c(95, 99, NA, 80, 97, NA),
    identified = c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE),
    combined_peak = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  cbind(p, .parse_formula(p$formula))
}

test_that(".reference_peak_ok rejects interpolated and merged peaks", {
  p <- make_peaks()
  expect_equal(
    .reference_peak_ok(p),
    c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)
  )
})

test_that("halogen, B, P and Si disqualify a peak but O, N and S do not", {
  p <- make_peaks()[1, ]
  allowed <- c("C6H14O", "C6H15N", "C6H14S")
  for (f in allowed) {
    q <- cbind(p[, c("name", "area", "score", "combined_peak", "peak_label")],
               .parse_formula(f))
    expect_true(.reference_peak_ok(q), info = f)
  }
  blocked <- c("C6H13Cl", "C6H13F", "C6H13Br", "C6H13I", "C6H15B", "C6H15P", "C6H15Si")
  for (f in blocked) {
    q <- cbind(p[, c("name", "area", "score", "combined_peak", "peak_label")],
               .parse_formula(f))
    expect_false(.reference_peak_ok(q), info = f)
  }
})

test_that("a match score below 30 is rejected", {
  p <- make_peaks()
  p$score[1] <- 29.9
  expect_false(.reference_peak_ok(p)[1])
  p$score[1] <- 30
  expect_true(.reference_peak_ok(p)[1])
})

test_that("peaks below 0.001 percent of total area are rejected", {
  p <- make_peaks()
  # total area is 97000; 0.001% of that is 0.97
  p$area[1] <- 0.5
  expect_false(.reference_peak_ok(p)[1])
  p$area[1] <- 500
  expect_true(.reference_peak_ok(p)[1])
})

test_that("an unnamed peak is rejected even with a score", {
  p <- make_peaks()
  p$name[1] <- ""
  expect_false(.reference_peak_ok(p)[1])
  p$name[1] <- NA_character_
  expect_false(.reference_peak_ok(p)[1])
})

test_that("combined peaks are detected from the label when no flag is present", {
  p <- make_peaks()
  p$name[6] <- "Hexane"
  p$score[6] <- 95
  p$combined_peak <- NULL
  expect_false(.reference_peak_ok(p)[6])
})

test_that("suggestion falls back when the largest peak fails the check", {
  p <- make_peaks()
  s <- .suggest_reference_peak(p)
  # peak 3 is the largest but is interpolated, and peak 5 is the next largest
  # but is chlorinated, so toluene is proposed instead
  expect_equal(s$peak_index, 2L)
  expect_false(s$approved)
  expect_equal(s$name, "Toluene")
  expect_equal(s$n_peaks, 6L)
})

test_that("area percent is relative to total FID area", {
  p <- make_peaks()
  p$name[3] <- "Propane"
  p$score[3] <- 90
  p$identified[3] <- TRUE
  s <- .suggest_reference_peak(p)
  expect_true(s$approved)
  expect_equal(s$peak_index, 3L)
  expect_equal(s$area_percent, 50000 / 97000 * 100)
  expect_equal(s$total_area, 97000)
})

test_that("nothing is proposed when no peak qualifies", {
  p <- make_peaks()
  p$score <- NA_real_
  s <- .suggest_reference_peak(p)
  expect_true(is.na(s$peak_index))
  expect_false(s$approved)
})
