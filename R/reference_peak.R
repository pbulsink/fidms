#' Suggest an internal-standard reference peak
#'
#' Port of `frmSuggested`'s initialisation: the largest-area peak is proposed,
#' reported alongside its area percent of the total FID area, its MS match
#' score and its compound name. If the proposal fails the quality check the
#' analyst browses the full peak list instead (`frmBrowse`), which the Shiny
#' peak-browser table replaces.
#'
#' @param peaks Peak table from [.interpolate_unknowns()] or
#'   [.join_msd_identity()].
#' @param min_score Minimum MS match score, see [.reference_peak_ok()].
#' @param min_area_percent Minimum share of total FID area, see
#'   [.reference_peak_ok()].
#' @return A list with the proposed `peak_index` (row of `peaks`, `NA` when
#'   nothing qualifies), `n_peaks`, `total_area`, `area_percent`, `score` and
#'   `name`, plus `approved` indicating whether the proposal passed
#'   [.reference_peak_ok()].
#' @keywords internal
.suggest_reference_peak <- function(
  peaks,
  min_score = 30,
  min_area_percent = 0.001
) {
  total_area <- sum(peaks$area, na.rm = TRUE)
  ok <- .reference_peak_ok(
    peaks,
    min_score = min_score,
    min_area_percent = min_area_percent
  )

  # Propose the largest peak overall, as the macro did; fall back to the
  # largest *qualifying* peak so the app has something to show when the
  # biggest peak fails the check.
  idx <- which.max(peaks$area)
  approved <- length(idx) == 1L && isTRUE(ok[idx])
  if (!approved) {
    qualifying <- which(ok)
    idx <- if (length(qualifying)) {
      qualifying[which.max(peaks$area[qualifying])]
    } else {
      NA_integer_
    }
  }

  list(
    peak_index = if (length(idx)) idx else NA_integer_,
    n_peaks = nrow(peaks),
    total_area = total_area,
    area_percent = if (is.na(idx)) {
      NA_real_
    } else {
      peaks$area[idx] / total_area * 100
    },
    score = if (is.na(idx)) NA_real_ else peaks$score[idx],
    name = if (is.na(idx)) NA_character_ else peaks$name[idx],
    approved = approved
  )
}

#' Quality check for internal-standard reference peaks
#'
#' Port of the VBA `PeakApprove()`. A peak is usable as an internal standard
#' only if it:
#'
#' * carries a compound name from the MS library search;
#' * has an MS match score (so it was not interpolated or extrapolated) of at
#'   least `min_score`;
#' * is not a combined/collapsed FID peak;
#' * contains no halogen, boron, phosphorus or silicon — oxygen, nitrogen and
#'   sulfur are allowed; and
#' * accounts for at least `min_area_percent` of the total FID area.
#'
#' @param peaks Peak table.
#' @param min_score Minimum MS match score.
#' @param min_area_percent Minimum share of total FID area, in percent.
#' @return A logical vector, one element per peak.
#' @keywords internal
.reference_peak_ok <- function(
  peaks,
  min_score = 30,
  min_area_percent = 0.001
) {
  named <- !is.na(peaks$name) & nzchar(peaks$name)
  scored <- !is.na(peaks$score) & peaks$score >= min_score

  combined <- if (!is.null(peaks$combined_peak)) {
    peaks$combined_peak %in% TRUE
  } else if (!is.null(peaks$peak_label)) {
    grepl(" - ", peaks$peak_label, fixed = TRUE)
  } else {
    FALSE
  }

  # Halogens, B, P and Si disqualify a peak; O, N and S are permitted.
  disqualifying <- c("F", "Cl", "Br", "I", "B", "P", "Si")
  hetero <- rowSums(peaks[, disqualifying, drop = FALSE], na.rm = TRUE) > 0

  total_area <- sum(peaks$area, na.rm = TRUE)
  big_enough <- !is.na(peaks$area) &
    total_area > 0 &
    (peaks$area / total_area) * 100 >= min_area_percent

  named & scored & !combined & !hetero & big_enough
}
