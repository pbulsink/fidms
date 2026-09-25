#' Extract scan-level ion abundance matrix from a raw MS file
#'
#' Reads the full-scan MS data from an Agilent `.D` directory. Raw data is
#' full-scan, as acquired; diagnostic class-ion sums are derived post-hoc from
#' this matrix so library search / non-target capability on the same raw file
#' is preserved.
#'
#' @param path Path to an Agilent `.D` directory (or its `data.ms` file
#'   directly).
#' @param mz_range Optional length-2 numeric vector; when supplied, scans are
#'   filtered to `mz >= mz_range[1] & mz <= mz_range[2]` before returning.
#' @return A tibble with columns `rt`, `mz`, `intensity` (one row per observed
#'   ion in each scan; long format).
#' @keywords internal
.read_ms_scan_matrix <- function(path, mz_range = NULL) {
  if (!requireNamespace("chromConverter", quietly = TRUE)) {
    stop(
      "Package 'chromConverter' is required to read raw MS data.",
      call. = FALSE
    )
  }
  data_ms <- if (grepl("\\.ms$", path, ignore.case = TRUE)) {
    path
  } else {
    file.path(path, "data.ms")
  }
  if (!file.exists(data_ms)) {
    stop("No `data.ms` file found at '", path, "'.", call. = FALSE)
  }
  scans <- chromConverter::read_chemstation_ms(
    data_ms,
    what = "MS1",
    format_out = "data.frame",
    data_format = "long"
  )
  scans <- tibble::as_tibble(scans[, c("rt", "mz", "intensity")])
  if (!is.null(mz_range)) {
    scans <- scans[
      scans$mz >= mz_range[1] & scans$mz <= mz_range[2],
      ,
      drop = FALSE
    ]
  }
  scans
}

#' Build diagnostic class-ion sum traces from a scan matrix
#'
#' Sums scan-level ion abundances into the diagnostic ion-sum series defined
#' by [sim_masses] (or [sim_masses_oxygenate]), one summed trace per scan
#' per class/family. Matching to `monitored_mz` is done at the nearest
#' integer nominal mass, since low-resolution quad MS data does not resolve
#' fractional masses.
#'
#' @param scan_matrix Long tibble from [.read_ms_scan_matrix()] (`rt`, `mz`,
#'   `intensity`).
#' @param sim_masses Ion-sum definition table, default [sim_masses]. Must
#'   have columns `class` (or `family`), `ion_sum_label`, `monitored_mz`
#'   (list-column of numeric vectors).
#' @return A tibble with columns `rt`, `class`, `ion_sum_label`,
#'   `summed_abundance` (one row per scan x family).
#' @keywords internal
.build_class_ion_traces <- function(
  scan_matrix,
  sim_masses = fidms::sim_masses
) {
  stopifnot(all(c("rt", "mz", "intensity") %in% names(scan_matrix)))
  class_col <- if ("class" %in% names(sim_masses)) "class" else "family"

  nominal_mz <- round(scan_matrix$mz)
  traces <- vector("list", nrow(sim_masses))
  for (i in seq_len(nrow(sim_masses))) {
    targets <- round(sim_masses$monitored_mz[[i]])
    keep <- nominal_mz %in% targets
    if (!any(keep)) {
      traces[[i]] <- tibble::tibble(
        rt = numeric(0),
        class = character(0),
        ion_sum_label = character(0),
        summed_abundance = numeric(0)
      )
      next
    }
    summed <- stats::aggregate(
      intensity ~ rt,
      data = scan_matrix[keep, , drop = FALSE],
      FUN = sum
    )
    traces[[i]] <- tibble::tibble(
      rt = summed$rt,
      class = sim_masses[[class_col]][i],
      ion_sum_label = sim_masses$ion_sum_label[i],
      summed_abundance = summed$intensity
    )
  }
  dplyr::bind_rows(traces)
}

#' Integrate class-ion traces over an FID peak's retention window
#'
#' Sums each diagnostic ion-family trace within `[rt_start, rt_end]`,
#' optionally subtracting a local linear background estimated from the
#' window edges.
#'
#' @param class_traces Long tibble from [.build_class_ion_traces()].
#' @param rt_start,rt_end Retention-time window bounds.
#' @param background Subtract a local linear baseline before integrating.
#' @return A named numeric vector, one entry per `ion_sum_label`, of
#'   background-subtracted ion-family sums for this window.
#' @keywords internal
.integrate_class_window <- function(
  class_traces,
  rt_start,
  rt_end,
  background = TRUE
) {
  if (rt_end <= rt_start) {
    stop("`rt_end` must be greater than `rt_start`.", call. = FALSE)
  }
  labels <- unique(class_traces$ion_sum_label)
  out <- stats::setNames(numeric(length(labels)), labels)
  for (lbl in labels) {
    tr <- class_traces[class_traces$ion_sum_label == lbl, , drop = FALSE]
    tr <- tr[order(tr$rt), , drop = FALSE]
    win <- tr[tr$rt >= rt_start & tr$rt <= rt_end, , drop = FALSE]
    if (!nrow(win)) {
      out[lbl] <- 0
      next
    }
    y <- win$summed_abundance
    if (background) {
      base <- .baseline_estimate(tr, rt_start, rt_end)
      y <- pmax(y - base, 0)
    }
    out[lbl] <- sum(y)
  }
  out
}

# Local linear baseline evaluated at each in-window rt, anchored on the
# median abundance just outside the window edges (falls back to 0 when no
# edge data is available, e.g. the window touches the trace boundary).
.baseline_estimate <- function(trace, rt_start, rt_end, edge_width = NULL) {
  if (is.null(edge_width)) {
    span <- diff(range(trace$rt))
    edge_width <- if (is.finite(span) && span > 0) span * 0.02 else 0
  }
  left <- trace[
    trace$rt >= rt_start - edge_width & trace$rt < rt_start,
    ,
    drop = FALSE
  ]
  right <- trace[
    trace$rt > rt_end & trace$rt <= rt_end + edge_width,
    ,
    drop = FALSE
  ]
  left_val <- if (nrow(left)) stats::median(left$summed_abundance) else NA_real_
  right_val <- if (nrow(right)) {
    stats::median(right$summed_abundance)
  } else {
    NA_real_
  }
  if (is.na(left_val) && is.na(right_val)) {
    return(0)
  }
  if (is.na(left_val)) {
    left_val <- right_val
  }
  if (is.na(right_val)) {
    right_val <- left_val
  }
  win <- trace[trace$rt >= rt_start & trace$rt <= rt_end, , drop = FALSE]
  frac <- if (rt_end > rt_start) {
    (win$rt - rt_start) / (rt_end - rt_start)
  } else {
    0
  }
  left_val + frac * (right_val - left_val)
}

#' Solve non-negative class composition for one ion-sum vector
#'
#' When the supplied ion sums exactly match the response matrix's rows
#' (a square system, as in ASTM D2425 Appendix X1's worked examples), solves
#' `ion_sums = A %*% fractions` exactly via [solve()] -- this matches the
#' standard's own worked method and is required to reproduce its published
#' results exactly (a non-negativity-constrained solve, e.g. NNLS, generally
#' does not: it can zero out small legitimately-solvable classes and bias the
#' remaining fractions). When more ion sums are supplied than classes (an
#' over-determined system, e.g. combining aromatic and saturate ion sums into
#' one solve, or future non-worked-example data), falls back to
#' [nnls::nnls()] to get a non-negative least-squares fit.
#'
#' @param ion_sums Named numeric vector of background-subtracted ion-sum
#'   abundances, as returned by [.integrate_class_window()]. Names must match
#'   (a subset of) the row names of the chosen response matrix.
#' @param standard Which standard's response matrix to use; currently only
#'   `"D2425"` has populated matrices.
#' @param fraction Which sub-matrix within the standard to solve against:
#'   `"aromatic"` or `"saturate"`.
#' @param response_matrix Nested list of response matrices, default
#'   [class_response_matrix]. Each `response_matrix[[standard]][[fraction]]`
#'   entry must be a list with a numeric `matrix` (rows = ion sums, columns =
#'   classes) and a named `sensitivity` vector (one entry per class).
#' @return A list with `fractions` (named numeric vector, one entry per
#'   class, normalized to sum to 1 -- non-negative when solved via NNLS, but
#'   an exact square solve can legitimately return small negative entries,
#'   matching the standard's own worked method), `residual` (normalized
#'   residual of the fit) and `raw_solution` (unnormalized solution).
#' @keywords internal
.solve_hydrocarbon_class <- function(
  ion_sums,
  standard = "D2425",
  fraction = c("aromatic", "saturate"),
  response_matrix = fidms::class_response_matrix
) {
  fraction <- match.arg(fraction)
  if (is.null(response_matrix[[standard]])) {
    stop(
      "No class_response_matrix entry for standard '",
      standard,
      "'.",
      call. = FALSE
    )
  }
  entry <- response_matrix[[standard]][[fraction]]
  if (is.null(entry) || is.null(entry$matrix) || anyNA(entry$matrix)) {
    stop(
      "class_response_matrix[['",
      standard,
      "']][['",
      fraction,
      "']] is ",
      "NULL or unpopulated; cannot solve without transcribed coefficients.",
      call. = FALSE
    )
  }
  a <- entry$matrix
  rows <- intersect(rownames(a), names(ion_sums))
  if (!length(rows)) {
    stop(
      "None of the supplied `ion_sums` names match rows of the '",
      standard,
      "'/'",
      fraction,
      "' response matrix.",
      call. = FALSE
    )
  }
  b <- as.numeric(ion_sums[rows])
  a_sub <- a[rows, , drop = FALSE]

  if (nrow(a_sub) == ncol(a_sub)) {
    raw <- stats::setNames(as.numeric(solve(a_sub, b)), colnames(a_sub))
    residuals <- b - as.numeric(a_sub %*% raw)
  } else {
    if (!requireNamespace("nnls", quietly = TRUE)) {
      stop(
        "Package 'nnls' is required to solve an over-determined ",
        "class-composition model.",
        call. = FALSE
      )
    }
    fit <- nnls::nnls(a_sub, b)
    raw <- stats::setNames(fit$x, colnames(a_sub))
    residuals <- fit$residuals
  }
  total <- sum(raw)
  fractions <- if (total > 0) raw / total else raw
  resid_norm <- if (sum(b^2) > 0) {
    sqrt(sum(residuals^2) / sum(b^2))
  } else {
    NA_real_
  }

  list(fractions = fractions, residual = resid_norm, raw_solution = raw)
}

#' Solve non-negative oxygenate-family composition
#'
#' Corrects for isotopic-tail and hydrocarbon-class cross-term interference
#' before solving the residual oxygenate-only ion-sum vector against the
#' oxygenate response matrix. This is a second, optional solve stage after
#' [.solve_hydrocarbon_class()], not merged into the hydrocarbon matrix, and is
#' an unverified proposed extension pending lab calibration data (see
#' `design/hydrocarbon-type-design.md`).
#'
#' @param ion_sums Named numeric vector of oxygenate diagnostic ion-sum
#'   abundances (see [sim_masses_oxygenate]).
#' @param hydrocarbon_fractions Named numeric vector, the `fractions` output
#'   of [.solve_hydrocarbon_class()], used to subtract hydrocarbon cross-term
#'   contributions before solving.
#' @param response_matrix Nested list, default [oxygenate_response_matrix].
#' @param isotope_correction Whether to apply [.correct_isotopic_tail()] first.
#' @param key Which entry of `response_matrix` to use (default `"lab"`, the
#'   only entry currently defined).
#' @return A list with `fractions`, `residual` and `raw_solution`, matching
#'   the shape of [.solve_hydrocarbon_class()]'s return value.
#' @keywords internal
.solve_oxygenate_class <- function(
  ion_sums,
  hydrocarbon_fractions = NULL,
  response_matrix = fidms::oxygenate_response_matrix,
  isotope_correction = TRUE,
  key = "lab"
) {
  if (!requireNamespace("nnls", quietly = TRUE)) {
    stop(
      "Package 'nnls' is required to solve the class-composition model.",
      call. = FALSE
    )
  }
  entry <- response_matrix[[key]]
  if (is.null(entry) || is.null(entry$matrix) || anyNA(entry$matrix)) {
    stop(
      "oxygenate_response_matrix[['",
      key,
      "']] is NULL or unpopulated; ",
      "requires lab calibration coefficients before it can be solved. ",
      "See `design/hydrocarbon-type-design.md`.",
      call. = FALSE
    )
  }
  a <- entry$matrix
  rows <- intersect(rownames(a), names(ion_sums))
  if (!length(rows)) {
    stop(
      "None of the supplied `ion_sums` names match rows of the oxygenate ",
      "response matrix.",
      call. = FALSE
    )
  }
  b <- stats::setNames(as.numeric(ion_sums[rows]), rows)

  if (isotope_correction && !is.null(entry$isotope_dominant_series)) {
    for (lbl in names(entry$isotope_dominant_series)) {
      if (
        lbl %in%
          names(b) &&
          entry$isotope_dominant_series[[lbl]] %in% names(ion_sums)
      ) {
        b[lbl] <- .correct_isotopic_tail(
          b[lbl],
          ion_sums[[entry$isotope_dominant_series[[lbl]]]]
        )
      }
    }
  }

  if (
    !is.null(hydrocarbon_fractions) && !is.null(entry$hydrocarbon_cross_terms)
  ) {
    for (lbl in names(entry$hydrocarbon_cross_terms)) {
      if (!lbl %in% names(b)) {
        next
      }
      cross <- entry$hydrocarbon_cross_terms[[lbl]]
      for (hc_class in names(cross)) {
        coef <- cross[[hc_class]]
        if (is.na(coef) || !hc_class %in% names(hydrocarbon_fractions)) {
          next
        }
        b[lbl] <- max(0, b[lbl] - coef * hydrocarbon_fractions[[hc_class]])
      }
    }
  }

  a_sub <- a[rows, , drop = FALSE]
  fit <- nnls::nnls(a_sub, as.numeric(b))
  raw <- stats::setNames(fit$x, colnames(a_sub))
  total <- sum(raw)
  fractions <- if (total > 0) raw / total else raw
  resid_norm <- if (sum(b^2) > 0) {
    sqrt(sum(fit$residuals^2) / sum(b^2))
  } else {
    NA_real_
  }

  list(fractions = fractions, residual = resid_norm, raw_solution = raw)
}

#' Subtract isotopic-tail contribution of a dominant series from an
#' oxygenate-family channel before NNLS solving
#'
#' Natural-abundance \eqn{^{13}}C isotope correction (not a co-fragmentation
#' interference): a fixed fraction of a strongly-responding series' summed
#' abundance is expected to appear one or two mass units higher, contaminating
#' an oxygenate channel at the same nominal mass. This is applied *before*
#' [.solve_oxygenate_class()]'s NNLS step, not as a matrix coefficient.
#'
#' @param oxygenate_ion_sum Single numeric value (or vector), the raw
#'   oxygenate-channel abundance.
#' @param dominant_series_sum Single numeric value (or vector, same length),
#'   the co-eluting dominant series' summed abundance (e.g. `Sum91`
#'   alkylbenzene for the `Sum94` alkylphenol channel).
#' @param isotope_fraction Fraction of `dominant_series_sum` expected to leak
#'   into the oxygenate channel. Defaults to the natural-abundance
#'   \eqn{^{13}}C estimate for a single carbon shift (~1.1%); pass an
#'   empirically-fit value once validated against oxygenate-free samples.
#' @return `oxygenate_ion_sum`, background-subtracted (floored at 0).
#' @keywords internal
.correct_isotopic_tail <- function(
  oxygenate_ion_sum,
  dominant_series_sum,
  isotope_fraction = 0.011
) {
  pmax(oxygenate_ion_sum - isotope_fraction * dominant_series_sum, 0)
}

#' Integrate individual parent/parent-1 ion masses over a retention window
#'
#' Sums raw scan-level abundance at each exact nominal mass in
#' `parent_ion_factors$mz` (and its polyisotopic `mz - 1` companion) within
#' `[rt_start, rt_end]`, optionally background-subtracted -- this is the
#' per-mass equivalent of [.integrate_class_window()], which instead sums
#' whole ion-sum *families*. Needed because the parent-ion carbon-number
#' calculation (ASTM D2425-23 Sec 15.2-15.5 / Appendix X1.2.2-X1.2.5,
#' [.average_carbon_number()]) requires the individual `P(m)` and `P(m-1)`
#' peak heights, not a pre-summed series.
#'
#' @param scan_matrix Long tibble from [.read_ms_scan_matrix()] (`rt`, `mz`,
#'   `intensity`).
#' @param rt_start,rt_end Retention-time window bounds.
#' @param series Which `parent_ion_factors$series` to integrate masses for
#'   (`"Alkylbenzene"` or `"Naphthalene"`).
#' @param parent_ion_factors Table of per-carbon-number parent masses,
#'   default [parent_ion_factors].
#' @param background Subtract a local linear baseline before integrating.
#' @return A named numeric vector keyed by nominal mass (as character, e.g.
#'   `"148"`), one entry per parent mass and its `mz - 1` companion.
#' @keywords internal
.integrate_parent_ion_window <- function(
  scan_matrix,
  rt_start,
  rt_end,
  series = c("Alkylbenzene", "Naphthalene"),
  parent_ion_factors = fidms::parent_ion_factors,
  background = TRUE
) {
  series <- match.arg(series)
  rows <- parent_ion_factors[
    parent_ion_factors$series == series,
    ,
    drop = FALSE
  ]
  masses <- sort(unique(c(rows$mz, rows$mz - 1)))
  # Reuse the tested .build_class_ion_traces()/.integrate_class_window()
  # machinery by treating each individual mass as its own one-mass "family".
  one_mass_sim <- tibble::tibble(
    class = paste0("mz", masses),
    ion_sum_label = paste0("mz", masses),
    monitored_mz = lapply(masses, identity)
  )
  traces <- .build_class_ion_traces(scan_matrix, sim_masses = one_mass_sim)
  out <- .integrate_class_window(
    traces,
    rt_start,
    rt_end,
    background = background
  )
  stats::setNames(out, sub("^mz", "", names(out)))
}

#' Mole fraction at each carbon number from parent-ion abundances
#'
#' Implements ASTM D2425-23 Eq 10/12 (Sec 15.2/15.4, Appendix X1.2.2/X1.2.4):
#' `x_n = [P(m) - K1 * P(m-1)] / K2`, where `P(m)` is the parent-ion
#' abundance at carbon number `n`, `P(m-1)` its polyisotopic companion one
#' mass unit lower, `K1` the isotope correction factor and `K2` the mole
#' sensitivity (both from [parent_ion_factors]/ASTM Table 1).
#'
#' @param parent_abundances Named numeric vector keyed by nominal mass (as
#'   from [.integrate_parent_ion_window()]).
#' @param series Which series to compute (`"Alkylbenzene"` or
#'   `"Naphthalene"`).
#' @param parent_ion_factors Table of per-carbon-number factors, default
#'   [parent_ion_factors].
#' @return A named numeric vector of mole fractions, one per carbon number,
#'   names are carbon numbers as character. Floored at 0 (a mole fraction
#'   cannot be negative; small negative values arise from measurement noise
#'   when a carbon number is essentially absent).
#' @keywords internal
.mole_fraction_by_carbon <- function(
  parent_abundances,
  series = c("Alkylbenzene", "Naphthalene"),
  parent_ion_factors = fidms::parent_ion_factors
) {
  series <- match.arg(series)
  rows <- parent_ion_factors[
    parent_ion_factors$series == series,
    ,
    drop = FALSE
  ]
  rows <- rows[order(rows$carbon), , drop = FALSE]
  p_m <- parent_abundances[as.character(rows$mz)]
  p_m1 <- parent_abundances[as.character(rows$mz - 1)]
  p_m[is.na(p_m)] <- 0
  p_m1[is.na(p_m1)] <- 0
  xn <- (p_m - rows$isotope_factor * p_m1) / rows$mole_sensitivity
  xn <- pmax(xn, 0)
  stats::setNames(xn, as.character(rows$carbon))
}

#' Average carbon number from a mole-fraction-by-carbon-number distribution
#'
#' Implements ASTM D2425-23 Eq 11/13 (Sec 15.3/15.5, Appendix
#' X1.2.3/X1.2.5): the abundance-weighted mean carbon number,
#' `A = sum(n * x_n) / sum(x_n)`.
#'
#' @param mole_fractions Named numeric vector from [.mole_fraction_by_carbon()]
#'   (names are carbon numbers, coercible to numeric).
#' @return A single numeric average carbon number, or `NA` when
#'   `mole_fractions` sums to zero (no usable signal).
#' @keywords internal
.average_carbon_number <- function(mole_fractions) {
  total <- sum(mole_fractions)
  if (!is.finite(total) || total <= 0) {
    return(NA_real_)
  }
  carbons <- as.numeric(names(mole_fractions))
  sum(carbons * mole_fractions) / total
}

#' Detect candidate MS class features on the class-ion traces
#'
#' A lightweight local-maxima peak picker run independently on each ion-sum
#' family's trace, restricted to one retention window (typically an FID
#' peak's `[rt_start, rt_end]`). This is intentionally simple (not a new
#' dependency) -- its job is only to decide whether a window hosts one or
#' several distinguishable MS class sub-features (e.g. a shoulder), not to do
#' general-purpose chromatographic peak detection.
#'
#' @param class_traces Long tibble from [.build_class_ion_traces()].
#' @param rt_start,rt_end Retention-time window to search within.
#' @param min_feature_snr Minimum ratio of a local maximum to the window's
#'   local baseline (see [.integrate_class_window()]) before it counts as a
#'   distinct feature.
#' @return A tibble with one row per detected feature: `rt_apex`, `rt_start`,
#'   `rt_end` (the half-window around the apex, clipped to the input window),
#'   and `snr`. Zero rows when nothing clears `min_feature_snr`.
#' @keywords internal
.detect_class_features <- function(
  class_traces,
  rt_start,
  rt_end,
  min_feature_snr = 3
) {
  if (rt_end <= rt_start) {
    stop("`rt_end` must be greater than `rt_start`.", call. = FALSE)
  }
  # Sum across ion-sum families at each scan time so a feature only needs to
  # be prominent in *some* family to be detected; per-family composition is
  # resolved later by .solve_hydrocarbon_class(), not here.
  totals <- stats::aggregate(
    summed_abundance ~ rt,
    data = class_traces[
      class_traces$rt >= rt_start & class_traces$rt <= rt_end,
      ,
      drop = FALSE
    ],
    FUN = sum
  )
  totals <- totals[order(totals$rt), , drop = FALSE]
  if (nrow(totals) < 3L) {
    return(.empty_class_features())
  }

  base <- .baseline_estimate(
    data.frame(rt = totals$rt, summed_abundance = totals$summed_abundance),
    rt_start,
    rt_end
  )
  y <- totals$summed_abundance
  noise <- stats::mad(y - base, na.rm = TRUE)
  if (!is.finite(noise) || noise <= 0) {
    noise <- max(1, stats::sd(y - base, na.rm = TRUE), na.rm = TRUE)
  }

  n <- length(y)
  is_max <- rep(FALSE, n)
  if (n >= 3L) {
    is_max[2:(n - 1L)] <- y[2:(n - 1L)] > y[1:(n - 2L)] &
      y[2:(n - 1L)] >= y[3:n]
  }
  snr <- (y - base) / noise
  is_max <- is_max & snr >= min_feature_snr
  if (!any(is_max)) {
    return(.empty_class_features())
  }

  apex_idx <- which(is_max)
  half_width <- max(diff(range(totals$rt)) / (2 * length(apex_idx)), 0)
  tibble::tibble(
    rt_apex = totals$rt[apex_idx],
    rt_start = pmax(rt_start, totals$rt[apex_idx] - half_width),
    rt_end = pmin(rt_end, totals$rt[apex_idx] + half_width),
    snr = snr[apex_idx]
  )
}

.empty_class_features <- function() {
  tibble::tibble(
    rt_apex = numeric(0),
    rt_start = numeric(0),
    rt_end = numeric(0),
    snr = numeric(0)
  )
}

#' Estimate carbon number from parent-ion peak deconvolution + retention index
#'
#' For the two series ASTM D2425-23 actually provides a mass-spectral
#' carbon-number calculation for (Sec 15.2-15.5, Appendix X1.2.2-X1.2.5):
#' alkylbenzenes and naphthalenes. When `class` maps to one of these and
#' `scan_matrix` is supplied, this integrates the parent/parent-1 peak
#' heights over the window ([.integrate_parent_ion_window()]), converts them
#' to a mole-fraction-by-carbon-number distribution
#' ([.mole_fraction_by_carbon()]), and takes the abundance-weighted mean
#' ([.average_carbon_number()]) -- i.e. it reproduces the standard's own
#' worked calculation, not a proxy. Per the standard's own note (Sec 15.6),
#' the other, lower-concentration classes (`'s 103, 115, 153, 151`, i.e.
#' Indan/Tetralin, Indene, Acenaphthene(-ylene), tricyclic aromatics, and all
#' saturate-fraction classes) have no such formulaic parent-ion calculation
#' -- their average carbon number is "estimated by inspection" in the
#' standard, which this function does not attempt to automate; for those
#' classes (or when `scan_matrix` is not supplied), the estimate falls back
#' to the retention-index cross-check alone.
#'
#' @param class Dominant class name (e.g. `"Alkylbenzene"`,
#'   `"Naphthalene_alkyl"`), used to select which parent-ion series (if any)
#'   applies.
#' @param scan_matrix Optional long tibble from [.read_ms_scan_matrix()] (`rt`,
#'   `mz`, `intensity`). Required for the parent-ion calculation; when `NULL`,
#'   only the RI cross-check is used.
#' @param rt_start,rt_end Retention-time window to integrate parent-ion
#'   abundances over (e.g. an FID peak's window). Required alongside
#'   `scan_matrix`.
#' @param ri Retention index for the same window (e.g. from
#'   [.compute_retention_index()]).
#' @param ri_model A retention-index model from [.fit_retention_index()]. Used
#'   only to confirm an RI value is meaningful (n-paraffin RI is
#'   ~100 * carbon number by the Kovats convention; used directly rather than
#'   inverting the fitted `time`/`ri` polynomial, which only maps
#'   `time -> ri`).
#' @param parent_ion_factors Table of per-carbon-number factors, default
#'   [parent_ion_factors].
#' @return A list with `estimate` (double), `range` (`c(low, high)`, spanning
#'   both available point estimates plus a +/-1 pad rather than reporting
#'   false precision), and `series_used` (`"Alkylbenzene"`, `"Naphthalene"`
#'   or `NA` when only the RI cross-check contributed).
#' @keywords internal
.estimate_carbon_number <- function(
  class = NA_character_,
  scan_matrix = NULL,
  rt_start = NA_real_,
  rt_end = NA_real_,
  ri = NA_real_,
  ri_model = NULL,
  parent_ion_factors = fidms::parent_ion_factors
) {
  series <- .parent_ion_series_for_class(class)
  series_est <- NA_real_
  if (
    !is.na(series) &&
      !is.null(scan_matrix) &&
      !is.na(rt_start) &&
      !is.na(rt_end)
  ) {
    parent_abundances <- .integrate_parent_ion_window(
      scan_matrix,
      rt_start,
      rt_end,
      series = series,
      parent_ion_factors = parent_ion_factors
    )
    mole_fractions <- .mole_fraction_by_carbon(
      parent_abundances,
      series = series,
      parent_ion_factors = parent_ion_factors
    )
    series_est <- .average_carbon_number(mole_fractions)
  }

  ri_est <- NA_real_
  if (
    !is.null(ri_model) && !is.na(ri) && inherits(ri_model, "fidms_ri_model")
  ) {
    # RI for n-paraffins is ~100 * carbon number by construction (Kovats
    # convention); use that directly rather than inverting the polynomial,
    # which is only fit for time -> RI, not RI -> carbon number.
    ri_est <- ri / 100
  }

  candidates <- c(series_est, ri_est)
  candidates <- candidates[!is.na(candidates)]
  if (!length(candidates)) {
    return(
      list(
        estimate = NA_real_,
        range = c(NA_real_, NA_real_),
        series_used = NA_character_
      )
    )
  }
  estimate <- mean(candidates)
  low <- min(candidates) - 1
  high <- max(candidates) + 1
  list(
    estimate = estimate,
    range = c(low, high),
    series_used = if (!is.na(series_est)) series else NA_character_
  )
}

# Maps a `class_primary` label onto the parent-ion series it corresponds to,
# per ASTM D2425-23 Sec 15.2-15.5 (only alkylbenzenes and naphthalenes have a
# formulaic parent-ion carbon-number calculation; see .estimate_carbon_number()
# docs for the rest).
.parent_ion_series_for_class <- function(class) {
  if (is.na(class)) {
    return(NA_character_)
  }
  if (identical(class, "Alkylbenzene")) {
    return("Alkylbenzene")
  }
  if (grepl("^Naphthalene", class)) {
    return("Naphthalene")
  }
  NA_character_
}

#' Attach hydrocarbon-class evidence to an aligned/identified peak table
#'
#' Implements the evidence-tier precedence from
#' `design/hydrocarbon-type-design.md`: a library hit with match `score`
#' at/above `params$min_library_score` is authoritative and passed through
#' unchanged (tier A). A lower-scoring hit triggers the class-ion solve as a
#' verification check (tier B): agreement raises confidence, disagreement is
#' flagged via `class_conflict_flag` without discarding the library name.
#' Peaks with no library hit at all attempt the class-ion solve as their
#' primary path (tier C); only peaks where that solve fails thresholds are
#' left unresolved (tier D candidates) for [.interpolate_unknowns()].
#'
#' @param peaks Peak table with `rt`, `rt_start`, `rt_end`, `identified`,
#'   `score`, `family` (from [.classify_family()]) and `ri` columns.
#' @param class_traces Long tibble from [.build_class_ion_traces()]. May be
#'   `NULL`, in which case every peak is returned unresolved (tier D
#'   candidates) and `ms_class_features` is empty -- this preserves the
#'   current CSV-only behavior where [.interpolate_unknowns()] runs on every
#'   unidentified peak.
#' @param scan_matrix Optional long tibble from [.read_ms_scan_matrix()] (`rt`,
#'   `mz`, `intensity`), the raw per-mass scan data underlying
#'   `class_traces`. When supplied, [.estimate_carbon_number()] can run the
#'   real ASTM parent-ion deconvolution for alkylbenzene/naphthalene classes
#'   instead of relying on the retention-index cross-check alone.
#' @param ri_model Retention-index model, see [.fit_retention_index()].
#' @param params Parameters from [hydrocarbon_class_params()].
#' @param standard,fraction Passed to [.solve_hydrocarbon_class()].
#' @return A list with `peaks` (augmented with `class_primary`,
#'   `class_secondary`, `class_fractions`, `class_residual`,
#'   `estimated_carbon_number`, `carbon_number_series`, `carbon_number_range`,
#'   `ms_feature_count`, `mixed_peak_flag`, `assignment_tier`,
#'   `assignment_basis`, `class_conflict_flag` and `resolved`) and
#'   `ms_class_features` (long detail table, tier B/C peaks only).
#' @keywords internal
.add_hydrocarbon_class <- function(
  peaks,
  class_traces = NULL,
  scan_matrix = NULL,
  ri_model = NULL,
  params = hydrocarbon_class_params(),
  standard = "D2425",
  fraction = c("aromatic", "saturate")
) {
  fraction <- match.arg(fraction)
  n <- nrow(peaks)
  peaks$class_primary <- NA_character_
  peaks$class_secondary <- NA_character_
  peaks$class_fractions <- vector("list", n)
  peaks$class_residual <- NA_real_
  peaks$estimated_carbon_number <- NA_real_
  peaks$carbon_number_series <- NA_character_
  peaks$carbon_number_range <- vector("list", n)
  peaks$ms_feature_count <- 0L
  peaks$mixed_peak_flag <- FALSE
  peaks$assignment_tier <- factor(
    NA_character_,
    levels = c("A", "B", "C", "D", "M"),
    ordered = TRUE
  )
  peaks$assignment_basis <- NA_character_
  peaks$class_conflict_flag <- FALSE
  # `resolved`: TRUE once a peak has authoritative (tier A) or class-ion
  # (tier B/C) evidence; only FALSE rows should proceed to
  # .interpolate_unknowns() as a tier D fallback.
  peaks$resolved <- FALSE

  ms_class_features <- vector("list", n)

  min_score <- params$min_library_score
  identified <- isTRUE(peaks$identified) |
    (!is.na(peaks$identified) & peaks$identified)
  identified[is.na(identified)] <- FALSE
  score <- peaks$score
  score[is.na(score)] <- -Inf # unscored library hits are treated as low-confidence, not skipped

  tier_a <- identified & score >= min_score
  peaks$assignment_tier[tier_a] <- "A"
  peaks$assignment_basis[tier_a] <- "library+family"
  peaks$resolved[tier_a] <- TRUE

  needs_class_ion <- which(!tier_a)

  for (i in needs_class_ion) {
    is_tier_b <- identified[i] # score < min_score here, since tier_a already excluded
    if (
      is.null(class_traces) ||
        is.na(peaks$rt_start[i]) ||
        is.na(peaks$rt_end[i])
    ) {
      peaks$assignment_tier[i] <- "D"
      peaks$assignment_basis[i] <- "interpolated_formula+family"
      next
    }

    solved <- tryCatch(
      {
        ion_sums <- .integrate_class_window(
          class_traces,
          peaks$rt_start[i],
          peaks$rt_end[i],
          background = TRUE
        )
        .solve_hydrocarbon_class(
          ion_sums,
          standard = standard,
          fraction = fraction
        )
      },
      error = function(e) NULL
    )

    features <- tryCatch(
      .detect_class_features(
        class_traces,
        peaks$rt_start[i],
        peaks$rt_end[i],
        min_feature_snr = params$min_feature_snr
      ),
      error = function(e) .empty_class_features()
    )
    peaks$ms_feature_count[i] <- nrow(features)
    peaks$mixed_peak_flag[i] <- nrow(features) > 1L

    if (
      is.null(solved) ||
        is.na(solved$residual) ||
        solved$residual > params$max_residual
    ) {
      peaks$assignment_tier[i] <- if (is_tier_b) "B" else "D"
      if (is_tier_b) {
        peaks$assignment_basis[i] <- "library+class_conflict"
        peaks$class_conflict_flag[i] <- TRUE
      } else {
        peaks$assignment_basis[i] <- "interpolated_formula+family"
      }
      next
    }

    ordered_fractions <- sort(solved$fractions, decreasing = TRUE)
    primary <- names(ordered_fractions)[1]
    secondary <- if (length(ordered_fractions) > 1L) {
      names(ordered_fractions)[2]
    } else {
      NA_character_
    }

    carbon <- .estimate_carbon_number(
      class = primary,
      scan_matrix = scan_matrix,
      rt_start = peaks$rt_start[i],
      rt_end = peaks$rt_end[i],
      ri = peaks$ri[i],
      ri_model = ri_model
    )

    peaks$class_primary[i] <- primary
    peaks$class_secondary[i] <- secondary
    peaks$class_fractions[[i]] <- solved$fractions
    peaks$class_residual[i] <- solved$residual
    peaks$estimated_carbon_number[i] <- carbon$estimate
    peaks$carbon_number_series[i] <- carbon$series_used
    peaks$carbon_number_range[[i]] <- carbon$range

    if (is_tier_b) {
      agrees <- !is.na(peaks$family[i]) &&
        identical(tolower(peaks$family[i]), tolower(gsub("_", "", primary)))
      # family names use spaces/no separators (e.g. "Naphthene"), class
      # names use underscores for multi-word classes (e.g.
      # "Naphthalene_alkyl") -- compare loosely rather than requiring exact
      # string equality between the two vocabularies.
      agrees <- agrees ||
        (!is.na(peaks$family[i]) &&
          grepl(
            gsub("_", ".*", tolower(primary)),
            tolower(peaks$family[i])
          ))
      if (agrees) {
        peaks$assignment_tier[i] <- "B"
        peaks$assignment_basis[i] <- "library+class_verified"
      } else {
        peaks$assignment_tier[i] <- "B"
        peaks$assignment_basis[i] <- "library+class_conflict"
        peaks$class_conflict_flag[i] <- TRUE
      }
      # tier B keeps the library name/formula as the reported value; class
      # evidence is corroboration only, not a replacement.
      peaks$resolved[i] <- TRUE
    } else {
      peaks$assignment_tier[i] <- "C"
      peaks$assignment_basis[i] <- "class+carbon_number+RI"
      peaks$resolved[i] <- TRUE
    }

    if (nrow(features)) {
      ms_class_features[[i]] <- tibble::tibble(
        fid_peak_id = peaks$peak[i],
        feature_id = seq_len(nrow(features)),
        rt_apex = features$rt_apex,
        rt_start = features$rt_start,
        rt_end = features$rt_end,
        class = primary,
        class_fraction = ordered_fractions[[1]],
        carbon_number = carbon$estimate,
        residual = solved$residual
      )
    }
  }

  peaks$assignment_tier[peaks$mixed_peak_flag] <- "M"

  list(
    peaks = peaks,
    ms_class_features = dplyr::bind_rows(ms_class_features)
  )
}

#' Default parameters for the class-ion pipeline
#'
#' @param background_window Minutes of edge data used for local baseline
#'   subtraction in [.integrate_class_window()].
#' @param min_feature_snr Minimum signal-to-noise (peak sum vs. local
#'   baseline) required before a class-ion result is treated as resolvable.
#' @param min_library_score Match-score threshold (see `peaks$score` from
#'   `.join_msd_identity()`) below which the class-ion solve is run as a
#'   verification check rather than skipped (tier A/B split).
#' @param max_residual Maximum acceptable normalized NNLS residual before a
#'   class-ion solve is treated as having failed thresholds (tier C -> D).
#' @param overlap_preference Which standard to prefer at the D2425/D2789
#'   carbon-number overlap (~C10); currently only `"D2425"` has populated
#'   response matrices, so this is a placeholder for a future gasoline-range
#'   project.
#' @return A named list of parameters.
#' @export
hydrocarbon_class_params <- function(
  background_window = 0.05,
  min_feature_snr = 3,
  min_library_score = 70,
  max_residual = 0.25,
  overlap_preference = c("D2425", "D2789")
) {
  overlap_preference <- match.arg(overlap_preference)
  list(
    background_window = background_window,
    min_feature_snr = min_feature_snr,
    min_library_score = min_library_score,
    max_residual = max_residual,
    overlap_preference = overlap_preference
  )
}
