#' Interpolate formulae for unidentified peaks
#'
#' Port of the VBA `Interpolate`. Interior gaps between identified peaks get a
#' linear interpolation of carbon/hydrogen/oxygen counts by position; peaks
#' eluting after the last identified peak are extrapolated from a linear
#' regression of carbon number on retention index.
#'
#' @param peaks Peak table from [.merge_close_peaks()] or [.join_msd_identity()].
#' @param skip Optional logical vector, same length as `nrow(peaks)`. `TRUE`
#'   marks peaks that already have evidence from
#'   [.add_hydrocarbon_class()] (tier A/B/C) and must be left untouched --
#'   neither used as an interpolation anchor (their formula/element columns
#'   may be `NA` when the evidence is class-only, not a library formula) nor
#'   overwritten with a formula guess. Defaults to `!peaks$identified` being
#'   the only fill criterion, i.e. today's unconditional behavior when
#'   `class_traces` is not supplied to [prepare_peaks()].
#' @return The table with `formula`, `name`, `mw` and element counts filled in
#'   for previously unidentified, non-`skip`ped peaks, plus an `origin` column
#'   taking values `"msd"`, `"interpolated"` or `"extrapolated"` (`NA` for
#'   `skip`ped rows).
#' @keywords internal
.interpolate_unknowns <- function(peaks, skip = NULL) {
  if (is.null(skip)) {
    skip <- rep(FALSE, nrow(peaks))
  }
  peaks$origin <- ifelse(peaks$identified, "msd", NA_character_) # ry: ignore[RY106]
  # Anchors for interpolation/extrapolation are genuine library
  # identifications only; `skip`ped (class-ion-resolved but formula-less)
  # peaks must not serve as a bracketing reference point.
  known <- which(peaks$identified & !skip)
  if (length(known) == 0L) {
    warning("No identified peaks; nothing to interpolate from.", call. = FALSE)
    return(peaks)
  }

  # interior gaps -----------------------------------------------------------
  for (k in seq_len(length(known) - 1L)) {
    lo <- known[k]
    hi <- known[k + 1L]
    gap <- seq_len(hi - lo - 1L) + lo
    gap <- gap[!skip[gap]]
    if (!length(gap)) {
      next
    }
    frac <- (gap - lo) / (hi - lo)
    for (el in c("C", "H", "O")) {
      peaks[[el]][gap] <- round(
        peaks[[el]][lo] + frac * (peaks[[el]][hi] - peaks[[el]][lo])
      )
    }
    peaks$origin[gap] <- "interpolated"
  }

  # tail extrapolation ------------------------------------------------------
  last <- max(known)
  tail_idx <- seq_len(nrow(peaks))[-seq_len(last)]
  tail_idx <- tail_idx[!skip[tail_idx]]
  if (length(tail_idx)) {
    ref <- peaks[known, , drop = FALSE]
    fit_c <- stats::lm(C ~ ri, data = ref)
    c_hat <- pmax(1, round(stats::predict(fit_c, newdata = peaks[tail_idx, ])))
    peaks$C[tail_idx] <- c_hat
    peaks$H[tail_idx] <- 2L * c_hat + 2L
    peaks$O[tail_idx] <- 0L
    peaks$origin[tail_idx] <- "extrapolated"
  }

  # leading peaks before the first identification ---------------------------
  first <- min(known)
  if (first > 1L) {
    head_idx <- seq_len(first - 1L)
    head_idx <- head_idx[!skip[head_idx]]
    for (el in c("C", "H", "O")) {
      peaks[[el]][head_idx] <- peaks[[el]][first]
    }
    peaks$origin[head_idx] <- "extrapolated"
  }

  filled <- !peaks$identified & !skip
  # interpolated formulae are C/H/O only; the remaining tracked elements are
  # known-absent rather than unknown, so zero them before computing mass.
  for (el in setdiff(.tracked_elements(), c("C", "H", "O"))) {
    peaks[[el]][filled] <- 0L
  }
  peaks$formula[filled] <- .build_formula(
    peaks$C[filled],
    peaks$H[filled],
    peaks$O[filled]
  )
  peaks$name[filled] <- ifelse(
    peaks$combined_peak[filled],
    ifelse(
      peaks$origin[filled] == "interpolated",
      "Unknown, Interpolated Combined FID Peaks",
      "Unknown, Extrapolated Combined FID Peaks"
    ),
    ifelse(
      peaks$origin[filled] == "interpolated",
      "Unknown, Interpolated",
      "Unknown, Extrapolated"
    )
  )
  peaks$mw[filled] <- .molecular_weight(peaks$formula[filled])
  peaks
}

.build_formula <- function(c_n, h_n, o_n) {
  part <- function(sym, n) ifelse(n > 0, paste0(sym, ifelse(n == 1, "", n)), "")
  trimws(paste(part("C", c_n), part("H", h_n), part("O", o_n)))
}
