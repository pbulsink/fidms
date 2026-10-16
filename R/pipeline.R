#' Run the full FID/MSD processing pipeline
#'
#' Convenience wrapper mirroring the macro's `Process()` entry point: parse,
#' align, identify, merge, interpolate, quantify and summarise.
#'
#' @param fid_path,msd_path,search_path Paths to the three CSV exports.
#' @param info Sample metadata from [sample_info()].
#' @param params Processing parameters from [fidms_params()].
#' @param style Calibration style, see [calibrate()].
#' @param include_solvents Whether solvents count toward the summary.
#' @param reference_peak Internal-standard reference peak passed to
#'   [calibrate()], see [.suggest_reference_peak()].
#' @param targets,points Calibration tables passed to [calibrate()].
#' @param ri_reference RI reference points, see [.fit_retention_index()].
#' @param engine Alignment engine passed to [.align_peaks()].
#' @return The list returned by [summarize_report()], with the alignment
#'   `shift` attached.
#' @export
process_sample <- function(
  fid_path,
  msd_path,
  search_path,
  info,
  params = fidms_params(),
  style = "slope",
  include_solvents = FALSE,
  reference_peak = NULL,
  targets = calibration_targets,
  points = calibration_points,
  ri_reference = rical_reference_points,
  engine = c("grid", "gcalignr")
) {
  prepared <- prepare_peaks(
    fid_path,
    msd_path,
    search_path,
    params = params,
    ri_reference = ri_reference,
    engine = engine
  )

  peaks <- calibrate(
    prepared$peaks,
    style = style,
    targets = targets,
    points = points,
    reference_peak = reference_peak,
    dilution = info$dilution,
    density = info$density
  )

  out <- summarize_report(peaks, info, include_solvents, shift = prepared$shift)
  out$shift <- prepared$shift
  out$ms_class_features <- prepared$ms_class_features
  out
}

#' Parse, align, identify and interpolate a sample
#'
#' The first half of [process_sample()], stopping before quantitation. Split
#' out so an internal-standard reference peak can be suggested and reviewed
#' (see [.suggest_reference_peak()]) before a calibration style is applied.
#'
#' @inheritParams process_sample
#' @param class_traces Optional long tibble from [.build_class_ion_traces()].
#'   When supplied, [.add_hydrocarbon_class()] runs between identification and
#'   interpolation (see `design/hydrocarbon-type-design.md`), and
#'   [.interpolate_unknowns()] only fills peaks left unresolved by both the
#'   library-hit check and the class-ion attempt. When `NULL` (the default,
#'   e.g. CSV-only samples with no raw `.D` scan data), behavior is unchanged:
#'   [.interpolate_unknowns()] runs on every unidentified peak.
#' @param scan_matrix Optional long tibble from [.read_ms_scan_matrix()], the
#'   raw per-mass scan data underlying `class_traces`. When supplied
#'   alongside `class_traces`, enables the ASTM parent-ion carbon-number
#'   deconvolution in [.estimate_carbon_number()] for alkylbenzene/naphthalene
#'   classes; without it, the carbon-number estimate falls back to the
#'   retention-index cross-check alone.
#' @param class_params Parameters from [hydrocarbon_class_params()], used
#'   only when `class_traces` is supplied.
#' @param class_standard,class_fraction Passed to [.add_hydrocarbon_class()]
#'   (which forwards them to [.solve_hydrocarbon_class()]); used only when
#'   `class_traces` is supplied.
#' @return A list with the interpolated `peaks` table, the alignment `shift`,
#'   and `ms_class_features` (empty tibble when `class_traces` is `NULL`).
#' @export
prepare_peaks <- function(
  fid_path,
  msd_path,
  search_path,
  params = fidms_params(),
  ri_reference = rical_reference_points,
  engine = c("grid", "gcalignr"),
  class_traces = NULL,
  scan_matrix = NULL,
  class_params = hydrocarbon_class_params(),
  class_standard = "D2425",
  class_fraction = c("aromatic", "saturate")
) {
  engine <- match.arg(engine)
  class_fraction <- match.arg(class_fraction)
  fid <- .read_fid_peaks(fid_path)
  msd <- .read_msd_peaks(msd_path)
  search <- .read_msd_search(search_path)

  aligned <- .align_peaks(fid, msd, params, engine = engine)
  ri_model <- .fit_retention_index(ri_reference)
  peaks <- .join_msd_identity(aligned$combined, search, ri_model = ri_model)
  peaks <- .merge_close_peaks(peaks, params)

  ms_class_features <- dplyr::bind_rows(list())
  skip <- NULL
  if (!is.null(class_traces)) {
    peaks$family <- .classify_family(peaks)
    hc <- .add_hydrocarbon_class(
      peaks,
      class_traces = class_traces,
      scan_matrix = scan_matrix,
      ri_model = ri_model,
      params = class_params,
      standard = class_standard,
      fraction = class_fraction
    )
    peaks <- hc$peaks
    ms_class_features <- hc$ms_class_features
    # .interpolate_unknowns() must only fill peaks left unresolved by both the
    # library-hit check and the class-ion attempt (tier D candidates); tier
    # B/C peaks are `resolved` but not `identified` (no library formula), so
    # they must be `skip`ped rather than treated as interpolation anchors.
    skip <- peaks$resolved & !peaks$identified
  }

  peaks <- .interpolate_unknowns(peaks, skip = skip)

  list(
    peaks = peaks,
    shift = aligned$shift,
    ms_class_features = ms_class_features
  )
}

#' Build the Total Compound Report table
#'
#' @param result Output of [summarize_report()] or [process_sample()].
#' @param include_class_ion Whether to append the hydrocarbon-class/class-ion
#'   evidence columns (`Hydrocarbon Class`, `Secondary Class`, `Est. Carbon
#'   No.`, `Assignment Tier`, `Mixed Peak`) described in
#'   `design/hydrocarbon-type-design.md`. Defaults to `FALSE` so existing
#'   golden-file/snapshot tests and report layouts are unaffected; set to
#'   `TRUE` when `result$peaks` was produced with `class_traces` supplied to
#'   [prepare_peaks()]/[process_sample()]. Columns are filled with `NA` for
#'   peaks that never went through the class-ion solve (e.g. tier A/D), and
#'   the new columns are appended after the original nine so existing column
#'   positions/indices are unaffected.
#' @return A tibble with the nine client-facing columns (plus five more when
#'   `include_class_ion = TRUE`). Collapsed rows use the macro's sentinel
#'   strings, so `Peak No.` and `Retention Time` are character.
#' @export
total_compound_table <- function(result, include_class_ion = FALSE) {
  p <- result$peaks
  tbl <- tibble::tibble(
    `Peak No.` = p$peak_label,
    `Compound Name` = p$name,
    `Retention Time` = ifelse(
      p$combined_peak,
      "Combined FID Peaks",
      format(round(p$rt, 3), trim = TRUE)
    ),
    `Compound Formula` = p$formula,
    `Compound Molar Mass` = round(p$mw, 3),
    `Area (FID)` = p$area,
    "Amount (ug/mL)" = round(p$amount_ug_ml, 3),
    `Target or Relative` = p$quant_basis,
    `MS Match Score` = p$score
  )
  if (!include_class_ion) {
    return(tbl)
  }
  # Class-ion columns are only present in `peaks` when `prepare_peaks()`/
  # `process_sample()` was run with `class_traces` supplied; fill with NA
  # otherwise so `include_class_ion = TRUE` degrades gracefully rather than
  # erroring on CSV-only samples.
  na_col <- function(x, default) if (is.null(x)) rep(default, nrow(p)) else x
  tbl$`Hydrocarbon Class` <- na_col(p$class_primary, NA_character_)
  tbl$`Secondary Class` <- na_col(p$class_secondary, NA_character_)
  tbl$`Est. Carbon No.` <- round(na_col(p$estimated_carbon_number, NA_real_), 1)
  tbl$`Assignment Tier` <- as.character(na_col(
    p$assignment_tier,
    NA_character_
  ))
  tbl$`Mixed Peak` <- na_col(p$mixed_peak_flag, NA)
  tbl
}

#' Static NOTES boilerplate printed on the Summary Report
#' @return Character vector of note lines.
#' @keywords internal
.summary_notes <- function() {
  c(
    "Compound types are programmatically determined and may not be entirely inclusive of compound type.",
    "Additional unmatched compounds are summed in unknown heteroatomic or unknown hydrocarbon type.",
    "Oxygenate includes all compounds with an oxygen atom, including methanol",
    "Solvents are not included in summary.",
    "Compound types reported here is not a replacement for standard PIONA analysis"
  )
}

#' Note shown beside the elemental table
#' @return A single string.
#' @keywords internal
.element_note <- function() {
  "Elemental values are wt% sum of all GC determined compounds and does not include H2O"
}
