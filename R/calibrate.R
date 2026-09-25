#' Fit calibration curves from raw standard-curve points
#'
#' Recomputes the three curve forms (slope only, slope+intercept, quadratic)
#' per target from the raw concentration/response points rather than using the
#' coefficients stored in the old workbook.
#'
#' @param points Standard-curve points with columns `target`, `conc`, `resp`.
#'   Defaults to the packaged [calibration_points].
#' @return A tibble with one row per target and columns `slope`, `slope_b`,
#'   `intercept`, `a2`, `b1`, `c0`, `r2_slope`, `r2_slope_b`, `r2_quad`.
#' @keywords internal
.fit_calibration <- function(points = calibration_points) {
  split(points, points$target) |>
    lapply(.fit_calibration_one) |>
    dplyr::bind_rows()
}

.fit_calibration_one <- function(d) {
  d <- d[!is.na(d$conc) & !is.na(d$resp), , drop = FALSE]
  r2 <- function(fit) if (is.null(fit)) NA_real_ else summary(fit)$r.squared
  safe_lm <- function(f) tryCatch(stats::lm(f, data = d), error = function(e) NULL)
  coefn <- function(fit, i) if (is.null(fit)) NA_real_ else unname(stats::coef(fit)[i])

  m1 <- safe_lm(conc ~ resp + 0)
  m2 <- safe_lm(conc ~ resp)
  m3 <- if (nrow(d) >= 3L) safe_lm(conc ~ resp + I(resp^2)) else NULL

  tibble::tibble(
    target = d$target[1],
    n_points = nrow(d),
    slope = coefn(m1, 1),
    intercept = coefn(m2, 1),
    slope_b = coefn(m2, 2),
    c0 = coefn(m3, 1),
    b1 = coefn(m3, 2),
    a2 = coefn(m3, 3),
    r2_slope = r2(m1),
    r2_slope_b = r2(m2),
    r2_quad = r2(m3)
  )
}

#' Match peaks to calibration targets
#'
#' Targets are matched on compound name or any of their synonyms
#' (case-insensitive, punctuation-insensitive).
#'
#' @param peaks Peak table with a `name` column.
#' @param targets Calibration target table, default [calibration_targets].
#' @return Character vector of matched target names (`NA` where no match).
#' @keywords internal
.match_targets <- function(peaks, targets = calibration_targets) {
  key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  lookup <- stats::setNames(
    rep(targets$target, 2L),
    c(key(targets$target), key(targets$synonyms))
  )
  lookup <- lookup[nzchar(names(lookup)) & names(lookup) != "na"]
  unname(lookup[key(peaks$name)])
}

#' Assign each peak to its governing calibration target
#'
#' Reproduces the workbook's `Closest Target Row` column, which is an
#' approximate-match `VLOOKUP` down a type-specific list of target row numbers:
#' each peak is governed by the nearest *preceding* matched target of its own
#' type (oxygenate when the formula contains oxygen, hydrocarbon otherwise),
#' falling back to the first target when a peak elutes before any of them.
#'
#' @param peaks Peak table with `O` and a `target` column of matched target
#'   names (as produced by [.match_targets()]).
#' @param matched Character vector of matched target names, one per peak.
#' @param targets Calibration target table, default [calibration_targets].
#' @return Integer vector of row indices into `peaks` giving the governing
#'   target peak for each row (`NA` when no target of any type was found).
#' @keywords internal
.assign_target_rows <- function(peaks, matched, targets = calibration_targets) {
  type <- targets$type[match(matched, targets$target)]
  is_oxy <- !is.na(peaks$O) & peaks$O > 0

  hc_rows <- which(!is.na(matched) & type == "hydrocarbon")
  ox_rows <- which(!is.na(matched) & type == "oxygenate")

  # An oxygenate with no calibrated oxygenates present falls back to the
  # hydrocarbon list, mirroring the workbook's nested IFNA.
  pick <- function(i) {
    pool <- if (is_oxy[i] && length(ox_rows)) ox_rows else hc_rows
    if (!length(pool)) pool <- ox_rows
    if (!length(pool)) return(NA_integer_)
    preceding <- pool[pool <= i]
    if (length(preceding)) max(preceding) else min(pool)
  }
  vapply(seq_len(nrow(peaks)), pick, integer(1))
}

#' Quantify peaks against a calibration
#'
#' Port of the VBA `Calibrate`. Peaks matching a calibration target are
#' quantified directly from that target's curve (negatives clamped to zero, as
#' the macro did). Every other peak is scaled off its governing target using
#' the Polyarc relative response
#' \deqn{A_i = A_t \frac{F_i}{F_t} \frac{C_t}{C_i} \frac{M_i}{M_t}}
#' where \eqn{F} is FID area, \eqn{C} carbon count and \eqn{M} molar mass.
#' This follows from the Polyarc reactor's uniform carbon response: all carbon
#' is converted to methane before the flame, so peak area is proportional to
#' moles of carbon, and mass is therefore proportional to area times molar mass
#' divided by carbon number.
#'
#' `style = "internal"` reproduces the macro's no-external-calibration path: a
#' relative molar response factor is taken against a reference peak, and the
#' sample is normalised so the peaks sum to 100% of the sample mass. It
#' therefore assumes complete recovery and needs no reference concentration.
#'
#' @param peaks Peak table from [.interpolate_unknowns()].
#' @param style One of `"slope"`, `"slope_intercept"`, `"quadratic"` or
#'   `"internal"`.
#' @param targets Calibration target table, default [calibration_targets].
#' @param points Calibration standard points, default [calibration_points].
#' @param reference_peak Row index of the internal-standard reference peak
#'   (required when `style = "internal"`, see [.suggest_reference_peak()]).
#' @param dilution Dilution factor applied to calibrated target amounts.
#' @param density Sample density in g/mL, used by the internal normalisation.
#' @return `peaks` with `amount_ug_ml` and `quant_basis` (`"Target"` or
#'   `"Relative"`) columns added.
#' @export
calibrate <- function(peaks,
                      style = c("slope", "slope_intercept", "quadratic", "internal"),
                      targets = calibration_targets,
                      points = calibration_points,
                      reference_peak = NULL,
                      dilution = 1,
                      density = 1) {
  style <- match.arg(style)

  if (style == "internal") {
    if (is.null(reference_peak) || is.na(reference_peak)) {
      stop("`reference_peak` is required for internal standardisation.", call. = FALSE)
    }
    peaks$amount_ug_ml <- .internal_amounts(peaks, reference_peak, density)
    peaks$quant_basis <- "Relative"
    peaks$quant_basis[reference_peak] <- "Target"
    peaks$calibration_target <- NA_character_
    return(peaks)
  }

  curves <- dplyr::left_join(.fit_calibration(points), targets, by = "target")
  matched <- .match_targets(peaks, targets)

  if (all(is.na(matched))) {
    stop(
      "No calibration targets were identified in this sample; ",
      "use style = \"internal\" to quantify against a reference peak.",
      call. = FALSE
    )
  }

  # Calibrated amounts for the target peaks themselves.
  amount <- rep(NA_real_, nrow(peaks))
  is_target <- !is.na(matched)
  row <- match(matched[is_target], curves$target)
  amount[is_target] <- pmax(
    .apply_curve(peaks$area[is_target], curves[row, , drop = FALSE], style), 0
  ) * dilution

  # Everything else is scaled off its governing target.
  gov <- .assign_target_rows(peaks, matched, targets)
  rel <- !is_target & !is.na(gov)
  amount[rel] <- amount[gov[rel]] * .polyarc_ratio(peaks, which(rel), gov[rel])

  peaks$amount_ug_ml <- amount
  peaks$quant_basis <- ifelse(is_target, "Target", "Relative")
  peaks$calibration_target <- ifelse(is_target, matched, matched[gov])
  peaks
}

#' Polyarc relative response ratio between peaks and their reference
#'
#' @param peaks Peak table.
#' @param i Integer row indices of the peaks being quantified.
#' @param ref Integer row indices of the corresponding reference peaks.
#' @return Numeric vector of response ratios.
#' @keywords internal
.polyarc_ratio <- function(peaks, i, ref) {
  (peaks$area[i] / peaks$area[ref]) *
    (peaks$C[ref] / peaks$C[i]) *
    (peaks$mw[i] / peaks$mw[ref])
}

.internal_amounts <- function(peaks, reference_peak, density) {
  # RMRF_A * Area, the workbook's column Z.
  rmrf_area <- (peaks$mw / peaks$mw[reference_peak]) *
    (peaks$C[reference_peak] / peaks$C) *
    peaks$area
  total <- sum(rmrf_area, na.rm = TRUE)
  # Mass fraction -> ppm by mass -> ug/mL via density.
  rmrf_area / total * 1e6 * density
}

.apply_curve <- function(area, curve, style) {
  switch(
    style,
    slope = area * curve$slope,
    slope_intercept = curve$intercept + area * curve$slope_b,
    quadratic = curve$c0 + area * curve$b1 + area^2 * curve$a2
  )
}

#' Compounds excluded from the GC recovery total
#'
#' The macro's recovery `SUMIFS` excluded the extraction solvents, and also
#' excluded methanol whenever a methanol content was entered separately (so it
#' is not counted twice).
#'
#' @param exclude_methanol Whether to also exclude methanol.
#' @return A character vector of compound names.
#' @keywords internal
.recovery_excluded <- function(exclude_methanol = FALSE) {
  base <- c("Methylene chloride", "Methane, dichloro-", "Acetonitrile", "Solvent")
  if (exclude_methanol) c(base, "Methyl Alcohol", "Methanol") else base
}

#' GC recovery for each calibration style
#'
#' Reproduces the recovery percentages `frmCalType` showed so the analyst can
#' compare curve forms before committing to one. The internal style is
#' normalised to 100% by construction.
#'
#' @inheritParams calibrate
#' @param info Sample metadata from [.sample_info()].
#' @return A tibble with one row per style and columns `style`, `gc_recovery`
#'   and `total_recovery` (GC plus the non-GC methanol and water contributions).
#' @keywords internal
.calibration_recovery <- function(peaks, info,
                                 targets = calibration_targets,
                                 points = calibration_points,
                                 reference_peak = NULL) {
  non_gc <- info$methanol_wt_pct + info$water_wt_pct
  drop <- .recovery_excluded(info$methanol_wt_pct > 0)
  key <- function(x) gsub("[^a-z0-9]", "", tolower(ifelse(is.na(x), "", x)))

  one <- function(style) {
    q <- tryCatch(
      calibrate(peaks, style = style, targets = targets, points = points,
                reference_peak = reference_peak, dilution = info$dilution,
                density = info$density),
      error = function(e) NULL
    )
    if (is.null(q)) return(NA_real_)
    keep <- !key(q$name) %in% key(drop)
    sum(q$amount_ug_ml[keep], na.rm = TRUE) / 10000 / info$density
  }

  styles <- c("slope", "slope_intercept", "quadratic", "internal")
  gc <- vapply(styles, one, numeric(1))
  tibble::tibble(
    style = styles,
    gc_recovery = unname(gc),
    total_recovery = unname(gc) + non_gc
  )
}
