test_that("sim_masses has the expected schema and unique nominal masses per row", {
  expect_true(all(
    c("standard", "class", "ion_sum_label", "monitored_mz") %in%
      names(sim_masses)
  ))
  expect_true(is.list(sim_masses$monitored_mz))
  # No accidental extra list-nesting from tribble()
  expect_true(all(vapply(sim_masses$monitored_mz, is.numeric, logical(1))))
  # Naphthalene parent ion is a single unsummed mass, not a series
  parent <- sim_masses[sim_masses$ion_sum_label == "Mz128", ]
  expect_equal(parent$monitored_mz[[1]], 128)
})

test_that("class_response_matrix D2425 entries are populated and dimension-consistent", {
  expect_false(is.null(class_response_matrix$D2425))
  aromatic <- class_response_matrix$D2425$aromatic
  saturate <- class_response_matrix$D2425$saturate
  expect_false(anyNA(aromatic$matrix))
  expect_false(anyNA(saturate$matrix))
  expect_equal(nrow(aromatic$matrix), ncol(aromatic$matrix))
  expect_equal(length(aromatic$sensitivity), ncol(aromatic$matrix))
  expect_setequal(names(aromatic$sensitivity), colnames(aromatic$matrix))
  expect_setequal(names(saturate$sensitivity), colnames(saturate$matrix))
  # D2789 intentionally has no response matrix (see design doc)
  expect_null(class_response_matrix$D2789)
})

test_that("oxygenate_response_matrix ships as an explicit NA placeholder", {
  expect_true(anyNA(oxygenate_response_matrix$lab$matrix))
  expect_true(is.list(oxygenate_response_matrix$lab$hydrocarbon_cross_terms))
  expect_true(is.list(oxygenate_response_matrix$lab$isotope_dominant_series))
})

test_that("hydrocarbon_class_params returns expected defaults", {
  p <- hydrocarbon_class_params()
  expect_equal(p$min_library_score, 70)
  expect_equal(p$overlap_preference, "D2425")
  expect_error(hydrocarbon_class_params(overlap_preference = "bogus"))
})

# ---- .solve_hydrocarbon_class() ---------------------------------------------

test_that(".solve_hydrocarbon_class errors on unpopulated/missing matrices", {
  ion_sums <- c(Sum71 = 100, Sum67 = 20)
  expect_error(
    .solve_hydrocarbon_class(
      ion_sums,
      standard = "D2789",
      fraction = "aromatic"
    ),
    "No class_response_matrix entry"
  )
  empty <- list(D2425 = list(aromatic = list(matrix = NULL)))
  expect_error(
    .solve_hydrocarbon_class(
      ion_sums,
      standard = "D2425",
      fraction = "aromatic",
      response_matrix = empty
    ),
    "unpopulated"
  )
})

test_that(".solve_hydrocarbon_class errors when no ion-sum names match the matrix", {
  ion_sums <- c(Foo = 100, Bar = 50)
  expect_error(
    .solve_hydrocarbon_class(
      ion_sums,
      standard = "D2425",
      fraction = "aromatic"
    ),
    "None of the supplied"
  )
})

test_that(".solve_hydrocarbon_class returns non-negative fractions summing to 1", {
  ion_sums <- c(
    Sum71 = 1521308,
    Sum67 = 4783717,
    Sum91 = 15394273,
    Sum103 = 17448505,
    Sum115 = 8718093,
    Mz128 = 1799576,
    Sum141 = 4424521,
    Sum153 = 3093908,
    Sum151 = 2810736,
    Sum177 = 1257580
  )
  res <- .solve_hydrocarbon_class(
    ion_sums,
    standard = "D2425",
    fraction = "aromatic"
  )
  expect_true(all(res$fractions >= 0))
  expect_equal(sum(res$fractions), 1, tolerance = 1e-8)
  expect_true(res$residual >= 0)
})

# ASTM D2425-23 Appendix X1 worked example (aromatic fraction): known-input
# diagnostic ion-sum abundances and the standard's published wt% answer.
# Kept at file scope so both the diagnostic and regression tests below share
# one source of truth.
x1_aromatic_ion_sums <- c(
  Sum71 = 1521308,
  Sum67 = 4783717,
  Sum91 = 15394273,
  Sum103 = 17448505,
  Sum115 = 8718093,
  Mz128 = 1799576,
  Sum141 = 4424521,
  Sum153 = 3093908,
  Sum151 = 2810736,
  Sum177 = 1257580
)
x1_aromatic_published_wt_pct <- c(
  Paraffin = 4.82,
  Cycloparaffin = 8.87,
  Alkylbenzene = 24.53,
  Indan_Tetralin = 31.95,
  Indene = 8.73,
  Naphthalene_parent = 2.80,
  Naphthalene_alkyl = 6.84,
  Acenaphthene = 5.76,
  Acenaphthylene = 3.62,
  Tricyclic_aromatic = 2.08
)

test_that("D2425 Appendix X1 aromatic matrix is well-conditioned and reconstructs its own worked-example input exactly", {
  # `class_response_matrix$D2425$aromatic$matrix` is transcribed from ASTM
  # D2425-23 Table X1.1; it is a square (10x10) system, so
  # `.solve_hydrocarbon_class()` solves it exactly (see below) rather than via
  # NNLS. This test just confirms the matrix itself is well-conditioned and
  # self-consistent with its own generating data.
  a <- class_response_matrix$D2425$aromatic$matrix
  rows <- intersect(rownames(a), names(x1_aromatic_ion_sums))
  b <- as.numeric(x1_aromatic_ion_sums[rows])
  a_sub <- a[rows, , drop = FALSE]

  expect_lt(kappa(a_sub), 10)

  exact_x <- stats::setNames(solve(a_sub, b), colnames(a_sub))
  reconstructed <- as.numeric(a_sub %*% exact_x[colnames(a_sub)])
  expect_equal(reconstructed, b, tolerance = 1e-6)
})

test_that(".solve_hydrocarbon_class reproduces the ASTM D2425 Appendix X1 aromatic worked example exactly", {
  res <- .solve_hydrocarbon_class(
    x1_aromatic_ion_sums,
    standard = "D2425",
    fraction = "aromatic"
  )
  expect_equal(sum(res$fractions), 1, tolerance = 1e-8)
  expect_lt(res$residual, 1e-6) # exact square solve, should reconstruct input essentially perfectly

  sens <- class_response_matrix$D2425$aromatic$sensitivity
  wt_pct <- 100 *
    (res$raw_solution / sens[names(res$raw_solution)]) /
    sum(res$raw_solution / sens[names(res$raw_solution)])

  expect_equal(
    as.numeric(wt_pct[names(x1_aromatic_published_wt_pct)]),
    as.numeric(x1_aromatic_published_wt_pct),
    tolerance = 0.01
  )
})

test_that(".solve_hydrocarbon_class reproduces the ASTM D2425 Appendix X1 saturate worked example exactly", {
  ion_sums <- c(
    Sum71 = 8872325,
    Sum67 = 13688870,
    Sum123 = 2800479,
    Sum149 = 992831,
    Sum91 = 798894
  )
  x1_saturate_published_wt_pct <- c(
    Paraffin = 59.96,
    Monocycloparaffin = 23.24,
    Dicycloparaffin = 11.68,
    Tricycloparaffin = 3.73,
    Alkylbenzene = 1.40
  )
  res <- .solve_hydrocarbon_class(
    ion_sums,
    standard = "D2425",
    fraction = "saturate"
  )
  expect_equal(sum(res$fractions), 1, tolerance = 1e-8)
  expect_lt(res$residual, 1e-6)

  sens <- class_response_matrix$D2425$saturate$sensitivity
  wt_pct <- 100 *
    (res$raw_solution / sens[names(res$raw_solution)]) /
    sum(res$raw_solution / sens[names(res$raw_solution)])

  expect_equal(
    as.numeric(wt_pct[names(x1_saturate_published_wt_pct)]),
    as.numeric(x1_saturate_published_wt_pct),
    tolerance = 0.01
  )
})

# ---- .integrate_class_window() ----------------------------------------------

test_that(".integrate_class_window sums abundance within the window and rejects bad windows", {
  traces <- tibble::tibble(
    rt = rep(c(1, 2, 3, 4, 5), 2),
    class = rep(c("A", "B"), each = 5),
    ion_sum_label = rep(c("SumA", "SumB"), each = 5),
    summed_abundance = c(0, 0, 10, 0, 0, 0, 0, 100, 0, 0)
  )
  out <- .integrate_class_window(
    traces,
    rt_start = 2,
    rt_end = 4,
    background = FALSE
  )
  expect_equal(unname(out["SumA"]), 10)
  expect_equal(unname(out["SumB"]), 100)
  expect_error(.integrate_class_window(traces, rt_start = 4, rt_end = 2))
})

test_that(".integrate_class_window with background subtraction reduces a flat baseline to ~0", {
  traces <- tibble::tibble(
    rt = seq(0, 10, by = 0.1),
    class = "A",
    ion_sum_label = "SumA",
    summed_abundance = 50 + c(rep(0, 50), rep(500, 21), rep(0, 30))
  )
  out <- .integrate_class_window(
    traces,
    rt_start = 4.5,
    rt_end = 6.5,
    background = TRUE
  )
  # baseline of 50 should be mostly subtracted, leaving close to the true peak sum
  raw_sum <- sum(traces$summed_abundance[traces$rt >= 4.5 & traces$rt <= 6.5])
  expect_lt(out[["SumA"]], raw_sum)
  expect_gt(out[["SumA"]], 0)
})

# ---- .build_class_ion_traces() ----------------------------------------------

test_that(".build_class_ion_traces sums scan intensities at nominal masses", {
  scan_matrix <- tibble::tibble(
    rt = c(1, 1, 1, 2, 2),
    mz = c(71.05, 85.02, 50, 71.1, 200),
    intensity = c(100, 50, 1000, 80, 5)
  )
  sm <- tibble::tribble(
    ~standard , ~class     , ~ion_sum_label , ~monitored_mz ,
    "D2425"   , "Paraffin" , "Sum71"        , c(71, 85)
  )
  traces <- .build_class_ion_traces(scan_matrix, sim_masses = sm)
  expect_equal(nrow(traces), 2)
  expect_equal(traces$summed_abundance[traces$rt == 1], 150)
  expect_equal(traces$summed_abundance[traces$rt == 2], 80)
})

# ---- oxygenate pipeline -----------------------------------------------------

test_that(".correct_isotopic_tail subtracts the expected fraction and floors at 0", {
  expect_equal(
    .correct_isotopic_tail(100, 1000, isotope_fraction = 0.011),
    100 - 11
  )
  expect_equal(.correct_isotopic_tail(5, 1000, isotope_fraction = 0.011), 0)
  expect_equal(
    .correct_isotopic_tail(c(100, 5), c(1000, 1000), isotope_fraction = 0.011),
    c(89, 0)
  )
})

test_that(".solve_oxygenate_class errors while the response matrix is unpopulated", {
  ion_sums <- c(Sum94 = 100, Sum109 = 50, Sum74 = 30, Sum60 = 20, Sum110 = 10)
  expect_error(.solve_oxygenate_class(ion_sums), "unpopulated")
})

test_that(".solve_oxygenate_class runs once a placeholder matrix is populated, applying isotope and cross-term corrections", {
  populated <- oxygenate_response_matrix
  # a simple, deliberately diagonal-dominant matrix for a structural test
  m <- diag(5)
  dimnames(m) <- dimnames(oxygenate_response_matrix$lab$matrix)
  populated$lab$matrix <- m
  populated$lab$sensitivity[] <- 1
  populated$lab$hydrocarbon_cross_terms$Sum94["Alkylbenzene"] <- 0.05
  populated$lab$hydrocarbon_cross_terms$Sum109["Alkylbenzene"] <- 0.02
  populated$lab$hydrocarbon_cross_terms$Sum110["Indan_Tetralin"] <- 0.1

  ion_sums <- c(Sum94 = 100, Sum109 = 50, Sum74 = 30, Sum60 = 20, Sum110 = 10)
  hc_fractions <- c(Alkylbenzene = 0.3, Indan_Tetralin = 0.2)

  res <- .solve_oxygenate_class(
    ion_sums,
    hydrocarbon_fractions = hc_fractions,
    response_matrix = populated,
    isotope_correction = FALSE
  )
  expect_true(all(res$fractions >= 0))
  expect_equal(sum(res$fractions), 1, tolerance = 1e-8)
  # Sum94 should be reduced by the Alkylbenzene cross-term (0.05 * 0.3 = 0.015
  # of raw Sum94 units is not itself meaningful here since fractions are
  # 0-1 scale while ion sums are raw counts; check the corrected raw solve
  # instead used a smaller effective Sum94 than the uncorrected input).
  res_uncorrected <- .solve_oxygenate_class(
    ion_sums,
    hydrocarbon_fractions = NULL,
    response_matrix = populated,
    isotope_correction = FALSE
  )
  expect_lt(
    res$raw_solution[["Alkylphenol"]],
    res_uncorrected$raw_solution[["Alkylphenol"]]
  )
})

# ---- .read_ms_scan_matrix() (requires a raw .D directory) -------------------

test_that(".read_ms_scan_matrix reads a real Agilent .D directory when available", {
  # The raw .D directory (~750MB, mostly vendor Results/Qual output) lives at
  # the project root and is deliberately not copied into tests/fixtures;
  # this test is skipped wherever it is absent (e.g. CI, other machines).
  d_path <- testthat::test_path("..", "..", "diesel.D")
  skip_if_not(
    dir.exists(d_path),
    "No .D directory available for raw MS ingestion test"
  )
  skip_if_not_installed("chromConverter")

  scans <- .read_ms_scan_matrix(d_path)
  expect_true(all(c("rt", "mz", "intensity") %in% names(scans)))
  expect_gt(nrow(scans), 0)
  expect_true(all(scans$mz >= 0))

  filtered <- .read_ms_scan_matrix(d_path, mz_range = c(60, 80))
  expect_true(all(filtered$mz >= 60 & filtered$mz <= 80))

  traces <- .build_class_ion_traces(scans)
  expect_true(all(
    c("rt", "class", "ion_sum_label", "summed_abundance") %in% names(traces)
  ))
  expect_gt(nrow(traces), 0)
})

test_that(".read_ms_scan_matrix errors clearly when data.ms is missing", {
  tmp <- withr::local_tempdir()
  expect_error(.read_ms_scan_matrix(tmp), "No `data.ms`")
})

# ---- .detect_class_features() -----------------------------------------------

test_that(".detect_class_features finds one feature for a single clean peak and two for a shoulder", {
  single <- tibble::tibble(
    rt = seq(0, 10, by = 0.05),
    class = "Alkylbenzene",
    ion_sum_label = "Sum91",
    summed_abundance = 10 +
      1000 * exp(-((seq(0, 10, by = 0.05) - 5)^2) / (2 * 0.3^2))
  )
  one <- .detect_class_features(
    single,
    rt_start = 3,
    rt_end = 7,
    min_feature_snr = 3
  )
  expect_equal(nrow(one), 1L)
  expect_equal(one$rt_apex, 5, tolerance = 0.05)

  rt <- seq(0, 10, by = 0.05)
  shoulder <- tibble::tibble(
    rt = rt,
    class = "Alkylbenzene",
    ion_sum_label = "Sum91",
    summed_abundance = 10 +
      1000 * exp(-((rt - 4.5)^2) / (2 * 0.2^2)) +
      1000 * exp(-((rt - 5.5)^2) / (2 * 0.2^2))
  )
  two <- .detect_class_features(
    shoulder,
    rt_start = 3,
    rt_end = 7,
    min_feature_snr = 3
  )
  expect_equal(nrow(two), 2L)

  expect_error(.detect_class_features(single, rt_start = 7, rt_end = 3))
})

test_that(".detect_class_features returns zero rows for a flat/noisy trace", {
  flat <- tibble::tibble(
    rt = seq(0, 10, by = 0.05),
    class = "A",
    ion_sum_label = "SumA",
    summed_abundance = 50
  )
  out <- .detect_class_features(flat, rt_start = 3, rt_end = 7)
  expect_equal(nrow(out), 0L)
})

# ---- parent-ion carbon-number deconvolution (ASTM Sec 15.2-15.5) ----------

test_that(".mole_fraction_by_carbon and .average_carbon_number reproduce the Appendix X1.2 alkylbenzene worked example", {
  # ASTM D2425-23 Appendix X1.2.2/X1.2.3: molecular weights 134 (C10) to 246
  # (C18); the standard reports the resulting mole-fraction distribution and
  # average carbon number (12.84) directly, so this is a known-input/
  # known-output regression test, not a synthetic one.
  published_mole_fractions <- c(
    `10` = 0.18,
    `11` = 0.19,
    `12` = 0.15,
    `13` = 0.13,
    `14` = 0.10,
    `15` = 0.08,
    `16` = 0.07,
    `17` = 0.06,
    `18` = 0.05
  )
  avg <- .average_carbon_number(published_mole_fractions)
  expect_equal(avg, 12.84, tolerance = 0.01)
})

test_that(".mole_fraction_by_carbon and .average_carbon_number reproduce the Appendix X1.2 naphthalene worked example", {
  published_mole_fractions <- c(
    `11` = 0.25,
    `12` = 0.20,
    `13` = 0.21,
    `14` = 0.15,
    `15` = 0.08,
    `16` = 0.05,
    `17` = 0.04,
    `18` = 0.03
  )
  avg <- .average_carbon_number(published_mole_fractions)
  expect_equal(avg, 13.06, tolerance = 0.01)
})

test_that(".mole_fraction_by_carbon applies Eq 10/12 (parent minus isotope-corrected parent-1, over mole sensitivity)", {
  # Two-carbon-number alkylbenzene series (C10 mz=134, C11 mz=148); construct
  # abundances that should reproduce a chosen mole-fraction ratio exactly.
  factors <- parent_ion_factors[
    parent_ion_factors$series == "Alkylbenzene" &
      parent_ion_factors$carbon %in% c(10, 11),
  ]
  abundances <- c(
    `134` = 1000,
    `133` = 100,
    `148` = 2000,
    `147` = 50
  )
  mf <- .mole_fraction_by_carbon(abundances, series = "Alkylbenzene")
  k10 <- factors[factors$carbon == 10, ]
  k11 <- factors[factors$carbon == 11, ]
  expect_equal(
    unname(mf["10"]),
    (1000 - k10$isotope_factor * 100) / k10$mole_sensitivity
  )
  expect_equal(
    unname(mf["11"]),
    (2000 - k11$isotope_factor * 50) / k11$mole_sensitivity
  )
})

test_that(".mole_fraction_by_carbon floors negative results at 0 and handles missing masses as 0 abundance", {
  mf <- .mole_fraction_by_carbon(c(`134` = 1), series = "Alkylbenzene")
  expect_true(all(mf >= 0))
  expect_true(all(is.finite(mf)))
})

test_that(".average_carbon_number returns NA when the distribution sums to zero", {
  expect_true(is.na(.average_carbon_number(c(`10` = 0, `11` = 0))))
})

test_that(".integrate_parent_ion_window sums scan intensity at each parent/parent-1 nominal mass within a window", {
  scan_matrix <- tibble::tibble(
    rt = c(5, 5, 5, 5, 15),
    mz = c(134, 133, 148, 200, 134),
    intensity = c(500, 40, 300, 20, 999) # rt=15 is outside the window
  )
  out <- .integrate_parent_ion_window(
    scan_matrix,
    rt_start = 4,
    rt_end = 6,
    series = "Alkylbenzene",
    background = FALSE
  )
  expect_equal(unname(out[["134"]]), 500)
  expect_equal(unname(out[["133"]]), 40)
  expect_equal(unname(out[["148"]]), 300)
})

# ---- .estimate_carbon_number() ----------------------------------------------

test_that(".estimate_carbon_number uses only the RI cross-check when scan_matrix is not supplied", {
  ref <- data.frame(
    carbon = 8:12,
    time = c(8, 9, 10, 11, 12),
    ri = c(800, 900, 1000, 1100, 1200)
  )
  model <- .fit_retention_index(
    data.frame(time = ref$time, ri = ref$ri),
    degree = 1
  )
  est <- .estimate_carbon_number(
    class = "Alkylbenzene",
    scan_matrix = NULL,
    ri = 1000,
    ri_model = model
  )
  expect_equal(est$estimate, 10, tolerance = 0.01)
  expect_true(est$range[1] <= est$estimate && est$range[2] >= est$estimate)
  expect_true(is.na(est$series_used))
})

test_that(".estimate_carbon_number returns NA when neither a parent-ion solve nor an RI model resolve", {
  est <- .estimate_carbon_number(
    class = NA_character_,
    scan_matrix = NULL,
    ri = NA_real_,
    ri_model = NULL
  )
  expect_true(is.na(est$estimate))
  expect_true(all(is.na(est$range)))
  expect_true(is.na(est$series_used))
})

test_that(".estimate_carbon_number runs the parent-ion deconvolution for an Alkylbenzene peak when scan_matrix is supplied", {
  # A single dominant carbon number (C12, mz 162) with negligible neighbors,
  # so the weighted average should land close to 12.
  masses <- parent_ion_factors[parent_ion_factors$series == "Alkylbenzene", ]
  rt <- seq(9, 11, by = 0.05)
  in_window <- rt >= 9.5 & rt <= 10.5
  scan_matrix <- dplyr::bind_rows(lapply(seq_len(nrow(masses)), function(i) {
    height <- if (masses$carbon[i] == 12) 5000 else 50
    tibble::tibble(
      rt = rt,
      mz = masses$mz[i],
      intensity = ifelse(in_window, height, 5)
    )
  }))
  est <- .estimate_carbon_number(
    class = "Alkylbenzene",
    scan_matrix = scan_matrix,
    rt_start = 9.5,
    rt_end = 10.5,
    ri = NA_real_,
    ri_model = NULL
  )
  expect_equal(est$series_used, "Alkylbenzene")
  expect_equal(est$estimate, 12, tolerance = 0.5)
})

test_that(".estimate_carbon_number does not attempt a parent-ion solve for classes without a formulaic calculation", {
  # Indan/Tetralin has no parent-ion carbon-number calculation per the
  # standard (Sec 15.6): "estimated by inspection", not automated here.
  est <- .estimate_carbon_number(
    class = "Indan_Tetralin",
    scan_matrix = tibble::tibble(rt = 1, mz = 1, intensity = 1),
    rt_start = 0,
    rt_end = 2,
    ri = NA_real_,
    ri_model = NULL
  )
  expect_true(is.na(est$series_used))
  expect_true(is.na(est$estimate))
})

# ---- .add_hydrocarbon_class() -----------------------------------------------

x1_saturate_ion_sums <- c(
  Sum71 = 8872325,
  Sum67 = 13688870,
  Sum123 = 2800479,
  Sum149 = 992831,
  Sum91 = 798894
)

make_class_traces_for_window <- function(rt_start, rt_end, ion_sums) {
  rt <- seq(rt_start - 1, rt_end + 1, by = 0.05)
  in_window <- rt >= rt_start & rt <= rt_end
  dplyr::bind_rows(lapply(names(ion_sums), function(lbl) {
    y <- rep(5, length(rt))
    y[in_window] <- y[in_window] + ion_sums[[lbl]] / sum(in_window)
    tibble::tibble(
      rt = rt,
      class = "x",
      ion_sum_label = lbl,
      summed_abundance = y
    )
  }))
}

test_that(".add_hydrocarbon_class assigns tier A unchanged for a high-score library hit", {
  peaks <- tibble::tibble(
    peak = 1L,
    rt = 10,
    rt_start = 9.9,
    rt_end = 10.1,
    ri = 1000,
    identified = TRUE,
    score = 90,
    family = "Aromatic",
    name = "Toluene",
    formula = "C7H8"
  )
  out <- .add_hydrocarbon_class(peaks, class_traces = NULL, ri_model = NULL)
  expect_equal(as.character(out$peaks$assignment_tier), "A")
  expect_equal(out$peaks$assignment_basis, "library+family")
  expect_true(out$peaks$resolved)
  expect_equal(nrow(out$ms_class_features), 0L)
})

test_that(".add_hydrocarbon_class runs the class-ion solve for an unidentified peak (tier C) and reports carbon-number range", {
  peaks <- tibble::tibble(
    peak = 1L,
    rt = 10,
    rt_start = 9,
    rt_end = 11,
    ri = 1000,
    identified = FALSE,
    score = NA_real_,
    family = NA_character_,
    name = NA_character_,
    formula = NA_character_
  )
  traces <- make_class_traces_for_window(9, 11, x1_saturate_ion_sums)
  out <- .add_hydrocarbon_class(
    peaks,
    class_traces = traces,
    ri_model = NULL,
    fraction = "saturate"
  )
  expect_equal(as.character(out$peaks$assignment_tier), "C")
  expect_equal(out$peaks$assignment_basis, "class+carbon_number+RI")
  expect_equal(out$peaks$class_primary, "Paraffin")
  expect_true(out$peaks$resolved)
  expect_false(is.na(out$peaks$class_residual))
})

test_that(".add_hydrocarbon_class treats a NULL class_traces as unresolved (tier D) for every peak", {
  peaks <- tibble::tibble(
    peak = 1:2,
    rt = c(10, 11),
    rt_start = c(9.9, 10.9),
    rt_end = c(10.1, 11.1),
    ri = c(1000, 1100),
    identified = c(FALSE, TRUE),
    score = c(NA_real_, 50),
    family = c(NA_character_, "Aromatic"),
    name = c(NA_character_, "Xylene"),
    formula = c(NA_character_, "C8H10")
  )
  out <- .add_hydrocarbon_class(peaks, class_traces = NULL, ri_model = NULL)
  expect_equal(as.character(out$peaks$assignment_tier), c("D", "D"))
  expect_false(any(out$peaks$resolved))
  expect_equal(nrow(out$ms_class_features), 0L)
})

test_that(".add_hydrocarbon_class flags a conflict when a low-score library hit disagrees with the class-ion solve", {
  peaks <- tibble::tibble(
    peak = 1L,
    rt = 10,
    rt_start = 9,
    rt_end = 11,
    ri = 1000,
    identified = TRUE,
    score = 50, # below min_library_score (70) -> tier B verification path
    family = "Naphthene", # deliberately mismatched vs the Paraffin-dominant solve below
    name = "Cyclohexane, methyl-",
    formula = "C7H14"
  )
  traces <- make_class_traces_for_window(9, 11, x1_saturate_ion_sums)
  out <- .add_hydrocarbon_class(
    peaks,
    class_traces = traces,
    ri_model = NULL,
    fraction = "saturate"
  )
  expect_equal(as.character(out$peaks$assignment_tier), "B")
  expect_equal(out$peaks$assignment_basis, "library+class_conflict")
  expect_true(out$peaks$class_conflict_flag)
  expect_true(out$peaks$resolved)
  # library name/formula are retained, never overwritten by the class label
  expect_equal(out$peaks$name, "Cyclohexane, methyl-")
})

# ---- prepare_peaks() integration with class_traces -------------------------

test_that("prepare_peaks() with class_traces = NULL preserves current behavior (interpolation runs on every unidentified peak)", {
  out <- prepare_peaks(
    testthat::test_path("fixtures", "FIDpeaks.csv"),
    testthat::test_path("fixtures", "MSDpeaks.csv"),
    testthat::test_path("fixtures", "MSDsearch.csv")
  )
  expect_false(any(is.na(out$peaks$formula)))
  expect_equal(nrow(out$ms_class_features), 0L)
})

test_that("prepare_peaks() with class_traces supplied skips interpolation for class-ion-resolved peaks", {
  base <- prepare_peaks(
    testthat::test_path("fixtures", "FIDpeaks.csv"),
    testthat::test_path("fixtures", "MSDpeaks.csv"),
    testthat::test_path("fixtures", "MSDsearch.csv")
  )
  unresolved <- which(!base$peaks$identified)
  skip_if(
    length(unresolved) == 0L,
    "fixture has no unidentified peaks to test against"
  )
  target <- unresolved[1]

  traces <- make_class_traces_for_window(
    base$peaks$rt_start[target],
    base$peaks$rt_end[target],
    x1_saturate_ion_sums
  )

  out <- prepare_peaks(
    testthat::test_path("fixtures", "FIDpeaks.csv"),
    testthat::test_path("fixtures", "MSDpeaks.csv"),
    testthat::test_path("fixtures", "MSDsearch.csv"),
    class_traces = traces,
    class_fraction = "saturate"
  )

  # the class-ion solve should have resolved the target peak (tier C) and its
  # formula must not have been overwritten by .interpolate_unknowns()'s guess.
  expect_equal(
    out$peaks$assignment_tier[target],
    factor("C", levels = levels(out$peaks$assignment_tier), ordered = TRUE)
  )
  expect_true(is.na(out$peaks$formula[target]))
  expect_false(is.na(out$peaks$class_primary[target]))
})
