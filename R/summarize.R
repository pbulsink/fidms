#' Sample metadata
#'
#' Replaces the VBA input boxes and `frmWater` user form.
#'
#' @param sample_id Client sample number (e.g. `"D319"`).
#' @param sample_name Sample description.
#' @param project_code Project code used in the output file name.
#' @param analyst Analyst initials; multiple analysts are slash-joined.
#' @param density Sample density in g/mL.
#' @param methanol_wt_pct Methanol content, weight percent.
#' @param water_wt_pct Water content, weight percent.
#' @param dilution Dilution factor.
#' @param report_date Report date.
#' @return A list of class `fidms_sample`.
#' @keywords internal
.sample_info <- function(
  sample_id,
  sample_name = "",
  project_code = "",
  analyst = "",
  density = 1,
  methanol_wt_pct = 0,
  water_wt_pct = 0,
  dilution = 1,
  report_date = Sys.Date()
) {
  structure(
    list(
      sample_id = sample_id,
      sample_name = sample_name,
      project_code = project_code,
      analyst = paste(analyst, collapse = "/"),
      density = density,
      methanol_wt_pct = methanol_wt_pct,
      water_wt_pct = water_wt_pct,
      dilution = dilution,
      report_date = as.Date(report_date)
    ),
    class = "fidms_sample"
  )
}

#' Summarise a quantified peak table
#'
#' Computes the chemical-family weight-percent breakdown, the CHNSO elemental
#' weight percents, GC recovery and the QC peak counts shown on both reports.
#'
#' @param peaks Quantified peak table from [calibrate()].
#' @param info Sample metadata from [.sample_info()].
#' @param include_solvents Whether solvent peaks contribute to the summary.
#' @param shift MS/FID offset time (minutes) reported as QC.
#' @return A list with `peaks` (annotated), `families`, `elements`, `totals`
#'   and `qc` components.
#' @export
summarize_report <- function(
  peaks,
  info,
  include_solvents = FALSE,
  shift = NA_real_
) {
  stopifnot(inherits(info, "fidms_sample"))
  peaks$family <- .classify_family(peaks)
  peaks$solvent <- .is_solvent(peaks)

  # wt% = (ug/mL) / (g/mL) = ug/g, i.e. ppm by mass; /10000 converts to
  # percent. This matches the macro's `SUMIFS(...)/10000/density`.
  peaks$wt_pct <- peaks$amount_ug_ml / 10000 / info$density

  used <- if (include_solvents) rep(TRUE, nrow(peaks)) else !peaks$solvent

  # The macro's recovery SUMIFS excluded the extraction solvents always, and
  # methanol whenever it was reported separately, so it is not double-counted.
  key <- function(x) gsub("[^a-z0-9]", "", tolower(ifelse(is.na(x), "", x)))
  in_recovery <- used &
    !key(peaks$name) %in% key(.recovery_excluded(info$methanol_wt_pct > 0))

  family_levels <- c(
    "Oxygenate",
    "Paraffin",
    "Isoparaffin",
    "Naphthene",
    "Aromatic",
    "Unsaturate",
    "Unknown Heteroatomic",
    "Unknown Hydrocarbon"
  )
  families <- tibble::tibble(
    compound_type = family_levels,
    wt_pct = vapply(
      family_levels,
      function(f) {
        sum(peaks$wt_pct[used & peaks$family == f], na.rm = TRUE)
      },
      numeric(1)
    )
  )

  gc_recovery <- sum(peaks$wt_pct[in_recovery], na.rm = TRUE)

  elements <- tibble::tibble(
    element = c("Carbon", "Hydrogen", "Nitrogen", "Oxygen", "Sulfur"),
    symbol = c("C", "H", "N", "O", "S")
  )
  elements$wt_pct <- vapply(
    elements$symbol,
    function(sym) {
      mass_frac <- peaks[[sym]] *
        atomic_weights$mass[atomic_weights$symbol == sym] /
        peaks$mw
      sum(peaks$wt_pct[used] * mass_frac[used], na.rm = TRUE)
    },
    numeric(1)
  )
  # elemental table is normalised to the GC-determined material only
  if (sum(elements$wt_pct) > 0) {
    elements$wt_pct <- elements$wt_pct / sum(elements$wt_pct) * 100
  }
  elements$symbol <- NULL

  totals <- list(
    gc_recovery = gc_recovery,
    methanol = info$methanol_wt_pct,
    water = info$water_wt_pct,
    total = gc_recovery + info$methanol_wt_pct + info$water_wt_pct,
    element_sum = sum(elements$wt_pct)
  )

  qc <- list(
    aligned_peaks = sum(peaks$origin == "msd", na.rm = TRUE),
    unique_peaks = sum(peaks$origin != "msd", na.rm = TRUE),
    total_peaks = nrow(peaks),
    interpolated_peaks = sum(
      peaks$origin %in% c("interpolated", "extrapolated")
    ),
    ms_fid_offset = shift
  )

  list(
    peaks = peaks,
    families = families,
    elements = elements,
    totals = totals,
    qc = qc,
    info = info,
    include_solvents = include_solvents
  )
}
