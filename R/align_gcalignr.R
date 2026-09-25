#' Align FID and MSD peak lists with GCalignR
#'
#' Alternative engine for [.align_peaks()]. Passes the two peak lists to
#' `GCalignR::align_chromatograms()` as two "samples", which performs the same
#' full-chromatogram linear shift followed by row-wise peak grouping against a
#' running mean. Prefer this engine when peak lists grow or when additional
#' detectors/samples need aligning together, since GCalignR's grouping handles
#' split and merged peaks more gracefully than a single nearest-neighbour pass.
#'
#' @inheritParams .align_peaks
#' @return The same structure as [.align_peaks()].
#' @keywords internal
.align_peaks_gcalignr <- function(fid, msd, params = fidms_params()) {
  rlang::check_installed("GCalignR")

  input <- list(
    FID = as.data.frame(fid[, c("rt", "area", "height")]),
    MSD = as.data.frame(msd[, c("rt", "area", "height")])
  )

  aligned <- suppressMessages(utils::capture.output(
    res <- GCalignR::align_chromatograms(
      data = input,
      rt_col_name = "rt",
      reference = "FID",
      max_linear_shift = params$shift_window,
      max_diff_peak2mean = params$partial_threshold,
      min_diff_peak2peak = params$merge_max_diff,
      permute = FALSE,
      write_output = NULL
    )
  ))
  aligned <- res

  rt <- aligned$aligned$rt
  rt <- rt[, c("FID", "MSD"), drop = FALSE]

  # GCalignR reports the *original* retention times in the aligned table, with
  # the linear shift recorded separately, so the peak lists are joined on the
  # unshifted values and the shift is only reported as the MS/FID offset time.
  shift <- .gcalignr_shift(aligned)

  # A zero entry means "no peak from this sample in this row".
  matched_rows <- rt$FID > 0 & rt$MSD > 0
  fid_idx <- match(round(rt$FID[matched_rows], 6), round(fid$rt, 6))
  msd_idx <- match(round(rt$MSD[matched_rows], 6), round(msd$rt, 6))
  keep <- !is.na(fid_idx) & !is.na(msd_idx)
  fid_idx <- fid_idx[keep]
  msd_idx <- msd_idx[keep]

  combined <- fid
  combined$rt <- combined$rt + params$rt_offset
  combined$msd_peak <- NA_integer_
  combined$msd_rt <- NA_real_
  combined$msd_area <- NA_real_
  combined$msd_peak[fid_idx] <- msd$peak[msd_idx]
  combined$msd_rt[fid_idx] <- msd$rt[msd_idx]
  combined$msd_area[fid_idx] <- msd$area[msd_idx]

  list(
    shift = shift,
    n_msd_unmatched = nrow(msd) - length(msd_idx),
    combined = combined,
    gcalignr = aligned
  )
}

.gcalignr_shift <- function(aligned) {
  ls <- tryCatch(aligned$Logfile$LinearShift, error = function(e) NULL)
  if (is.null(ls) || !"MSD" %in% ls$sample) {
    return(NA_real_)
  }
  # GCalignR records the shift added to MSD to bring it onto FID; the macro and
  # the grid engine report the offset in the opposite direction (msd ~ fid + d).
  -as.numeric(ls$shift[ls$sample == "MSD"])
}
