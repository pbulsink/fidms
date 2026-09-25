#' Default processing parameters
#'
#' Replaces the workbook's very-hidden `PARAM` sheet.
#'
#' @param rt_offset Fixed retention-time offset applied to FID times (`PARAM!B1`).
#' @param expected_shift Starting guess for the MS/FID shift (`PARAM!T4`).
#' @param shift_window Half-width of the shift search window (`PARAM!B4`).
#' @param shift_step Step size of the shift search (`PARAM!D4`).
#' @param partial_threshold Maximum |dRT| accepted as a match (`PARAM!J4`).
#' @param merge_max_diff Maximum |dRT| for merging split FID peaks (`PARAM!L4`).
#' @return A named list of parameters.
#' @export
fidms_params <- function(
  rt_offset = 0,
  expected_shift = -0.017,
  shift_window = 0.2,
  shift_step = 0.001,
  partial_threshold = 0.03,
  merge_max_diff = 0.05
) {
  list(
    rt_offset = rt_offset,
    expected_shift = expected_shift,
    shift_window = shift_window,
    shift_step = shift_step,
    partial_threshold = partial_threshold,
    merge_max_diff = merge_max_diff
  )
}

#' Find the retention-time shift between FID and MSD peak lists
#'
#' Grid search over candidate linear shifts, scoring each by the total
#' nearest-neighbour distance from every MSD peak to the closest shifted FID
#' peak. This is the same objective as `GCalignR`'s linear-shift step, kept
#' in-package so the two-chromatogram case needs no extra dependency.
#'
#' @param fid_rt,msd_rt Numeric retention-time vectors.
#' @param params Parameter list from [fidms_params()].
#' @return A list with the best `shift` and its `score`.
#' @keywords internal
.find_rt_shift <- function(fid_rt, msd_rt, params = fidms_params()) {
  candidates <- seq(
    params$expected_shift - params$shift_window,
    params$expected_shift + params$shift_window,
    by = params$shift_step
  )
  fid_sorted <- sort(fid_rt)
  scores <- vapply(
    candidates,
    function(d) {
      nn <- .nearest_distance(msd_rt, fid_sorted + d)
      sum(pmin(nn, params$partial_threshold * 3))
    },
    numeric(1)
  )
  best <- which.min(scores)
  list(
    shift = candidates[best],
    score = scores[best],
    candidates = candidates,
    scores = scores
  )
}

.nearest_distance <- function(x, sorted_ref) {
  idx <- findInterval(x, sorted_ref, all.inside = TRUE)
  lo <- sorted_ref[idx]
  hi <- sorted_ref[pmin(idx + 1L, length(sorted_ref))]
  pmin(abs(x - lo), abs(x - hi))
}

#' Align FID and MSD peak lists
#'
#' Two engines are available. `"grid"` (the default) applies the best linear
#' shift from [.find_rt_shift()] and then greedily matches each MSD peak to the
#' nearest unused FID peak within `partial_threshold`; it is a direct port of
#' the VBA `Lookup`/`Bump_Signal_Down` pair and is fast for two peak lists.
#' `"gcalignr"` delegates to `GCalignR::align_chromatograms()`, which performs
#' the same linear-shift step but groups peaks against a running mean; prefer it
#' as peak lists grow or when more than two detectors/samples must be aligned
#' together.
#'
#' @param fid FID peaks from [.read_fid_peaks()].
#' @param msd MSD peaks from [.read_msd_peaks()].
#' @param params Parameter list from [fidms_params()].
#' @param engine Alignment engine, `"grid"` or `"gcalignr"`.
#' @return A list with the `shift` used and a `combined` tibble: one row per FID
#'   peak, with matched MSD columns (`msd_peak`, `msd_rt`, `msd_area`) or `NA`.
#'   MSD-only peaks are dropped, matching the macro. The `"gcalignr"` engine
#'   additionally returns the raw `gcalignr` object.
#' @keywords internal
.align_peaks <- function(
  fid,
  msd,
  params = fidms_params(),
  engine = c("grid", "gcalignr")
) {
  engine <- match.arg(engine)
  if (engine == "gcalignr") {
    return(.align_peaks_gcalignr(fid, msd, params))
  }

  fid$rt <- fid$rt + params$rt_offset
  shift <- .find_rt_shift(fid$rt, msd$rt, params)$shift
  matches <- .greedy_match(msd$rt, fid$rt + shift, params$partial_threshold)

  combined <- fid
  combined$msd_peak <- NA_integer_
  combined$msd_rt <- NA_real_
  combined$msd_area <- NA_real_
  matched <- !is.na(matches)
  combined$msd_peak[matches[matched]] <- msd$peak[matched]
  combined$msd_rt[matches[matched]] <- msd$rt[matched]
  combined$msd_area[matches[matched]] <- msd$area[matched]

  list(
    shift = shift,
    n_msd_unmatched = sum(!matched),
    combined = combined
  )
}

#' Greedy nearest-neighbour matching
#'
#' @param x Values to match (MSD retention times).
#' @param y Candidate targets (shifted FID retention times).
#' @param tol Maximum absolute difference accepted.
#' @return Integer vector, same length as `x`, of indices into `y` (`NA` when
#'   no candidate is within `tol`).
#' @keywords internal
.greedy_match <- function(x, y, tol) {
  out <- rep(NA_integer_, length(x))
  ok_x <- which(!is.na(x))
  ok_y <- which(!is.na(y))
  if (!length(ok_x) || !length(ok_y)) {
    return(out)
  }
  pairs <- expand.grid(i = ok_x, j = ok_y)
  pairs$d <- abs(x[pairs$i] - y[pairs$j])
  pairs <- pairs[!is.na(pairs$d) & pairs$d <= tol, , drop = FALSE]
  pairs <- pairs[order(pairs$d), , drop = FALSE]
  used_y <- logical(length(y))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs$i[k]
    j <- pairs$j[k]
    if (is.na(out[i]) && !used_y[j]) {
      out[i] <- j
      used_y[j] <- TRUE
    }
  }
  out
}

#' Attach MSD library identifications to aligned peaks
#'
#' @param combined The `combined` tibble from [.align_peaks()].
#' @param search MSD search results from [.read_msd_search()].
#' @param tol Retention-time tolerance for the join.
#' @param ri_model Retention-index model, see [.fit_retention_index()].
#' @return `combined` with `name`, `formula`, `score`, element counts, `mw` and
#'   `ri` columns added.
#' @keywords internal
.join_msd_identity <- function(
  combined,
  search,
  tol = 0.02,
  ri_model = .fit_retention_index()
) {
  idx <- .greedy_match(combined$msd_rt, search$rt, tol)
  combined$name <- NA_character_
  combined$formula <- NA_character_
  combined$score <- NA_real_
  ok <- !is.na(idx)
  combined$name[ok] <- search$name[idx[ok]]
  combined$formula[ok] <- search$formula[idx[ok]]
  combined$score[ok] <- search$score[idx[ok]]

  counts <- .parse_formula(combined$formula)
  combined <- dplyr::bind_cols(combined, counts)
  combined$mw <- .molecular_weight(combined$formula)
  combined$ri <- .compute_retention_index(combined$rt, ri_model)
  combined$identified <- !is.na(combined$formula)
  combined
}

#' Merge split FID peaks that share a single MSD identification
#'
#' Port of the VBA `Collapse_Rows`: consecutive unidentified FID peaks closer
#' together than `merge_max_diff` are summed into one reported row, which is
#' flagged as a combined peak.
#'
#' @param combined Identified peak table.
#' @param params Parameter list from [fidms_params()].
#' @return A tibble with combined rows merged and a `combined_peak` flag plus a
#'   `peak_label` character column (`"774 - 778"` style for merged rows).
#' @keywords internal
.merge_close_peaks <- function(combined, params = fidms_params()) {
  n <- nrow(combined)
  group <- integer(n)
  g <- 0L
  for (i in seq_len(n)) {
    if (i == 1L) {
      g <- 1L
    } else {
      mergeable <- !combined$identified[i] &&
        !combined$identified[i - 1L] &&
        (combined$rt[i] - combined$rt[i - 1L]) <= params$merge_max_diff
      if (!mergeable) g <- g + 1L
    }
    group[i] <- g
  }
  combined$.group <- group

  out <- dplyr::summarise(
    dplyr::group_by(combined, .data$.group),
    peak = dplyr::first(.data$peak),
    peak_label = if (dplyr::n() > 1L) {
      paste(dplyr::first(.data$peak), "-", dplyr::last(.data$peak))
    } else {
      as.character(dplyr::first(.data$peak))
    },
    combined_peak = dplyr::n() > 1L,
    rt = stats::weighted.mean(.data$rt, .data$area),
    rt_start = min(.data$rt_start, na.rm = TRUE),
    rt_end = max(.data$rt_end, na.rm = TRUE),
    area = sum(.data$area, na.rm = TRUE),
    height = max(.data$height, na.rm = TRUE),
    msd_rt = dplyr::first(.data$msd_rt),
    name = dplyr::first(.data$name),
    formula = dplyr::first(.data$formula),
    score = dplyr::first(.data$score),
    mw = dplyr::first(.data$mw),
    ri = dplyr::first(.data$ri),
    identified = any(.data$identified),
    .groups = "drop"
  )
  out$.group <- NULL
  counts <- .parse_formula(out$formula)
  dplyr::bind_cols(out, counts)
}
