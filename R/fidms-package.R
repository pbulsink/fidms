#' @keywords internal
"_PACKAGE"

#' @importFrom rlang .data
NULL

#' Atomic weights
#'
#' Conventional IUPAC atomic weights for the elements tracked by
#' [.parse_formula()]. Replaces the workbook's `PeriodicTable` sheet.
#'
#' @format A tibble with columns `symbol`, `element`, `mass` (g/mol).
"atomic_weights"

#' Calibration targets
#'
#' Target compounds, synonyms, type and reference RT/RI exported from the
#' workbook's `Calibration` sheet.
#'
#' @format A tibble with columns `target`, `synonyms`, `type`, `rt`, `ri`,
#'   `carbon`, `hydrogen`, `oxygen`.
"calibration_targets"

#' Calibration standard-curve points
#'
#' Raw concentration/response pairs per target. Levels with a zero response are
#' below the detection limit and carry no calibration-curve entry, so they are
#' excluded rather than fit as an observed zero.
#'
#' @format A tibble with columns `target`, `level`, `conc`, `resp`.
"calibration_points"

#' Retention-index reference points
#'
#' n-Alkane and target-compound RT/RI anchors from the workbook's `RICal`
#' sheet. The stored degree-6 polynomial is deliberately not carried over;
#' [.fit_retention_index()] refits from these points.
#'
#' @format A tibble with columns `carbon`, `hydrogen`, `oxygen`, `time`, `ri`,
#'   `type`.
"rical_reference_points"

#' Chemical family name fragments
#'
#' Name-fragment lookup replacing the workbook's `ChemicalFamily` sheet.
#'
#' @format A tibble with columns `family`, `term`.
"chemical_family_terms"

#' Known solvent names
#'
#' Compounds excluded from the summary report unless the analyst opts in.
#'
#' @format A character vector.
"solvent_names"

#' Diagnostic hydrocarbon-class ion-sum definitions
#'
#' Characteristic-ion sums per ASTM D2425-23 (middle distillates, C10-C18)
#' and D2789-95 (gasoline range, C6-C10, ion-sum definitions only; no
#' response matrix is shipped for this standard, see [class_response_matrix]).
#' Used by [.build_class_ion_traces()] to sum scan-level ion abundances into
#' class-diagnostic traces.
#'
#' @format A tibble with columns `standard`, `class`, `ion_sum_label`,
#'   `monitored_mz` (list-column of numeric mass vectors).
"sim_masses"

#' Diagnostic oxygenate-family ion-sum definitions (proposed extension)
#'
#' A user-proposed extension to the hydrocarbon-class deconvolution
#' framework, not transcribed from ASTM D2425/D2789 (neither standard's
#' applicability extends to oxygenate-containing samples). Tracked
#' separately from [sim_masses] for that reason. Unsaturated FAMEs are
#' molecular ions of named compounds rather than a fragment-sum series, so
#' those rows use `compound`/`mz` instead of `ion_sum_label`/`monitored_mz`.
#'
#' @format A tibble with columns `family`, `ion_sum_label`, `monitored_mz`
#'   (list-column), `hydrocarbon_interference`, `compound`, `mz`.
"sim_masses_oxygenate"

#' Hydrocarbon class-composition response matrices
#'
#' Carbon-number-region-specific sensitivity/pattern matrices `A(C)` used by
#' [.solve_hydrocarbon_class()], transcribed from ASTM D2425-23 Appendix X1
#' (the worked example: Table X1.1 for the aromatic fraction, Table 5 for the
#' saturate fraction). Coefficients are starting values requiring local
#' validation, not directly-applicable constants. No D2789 (gasoline-range)
#' entries are populated; the lab's samples are middle-distillate-oriented.
#'
#' @format A nested list keyed by standard (currently only `"D2425"`), then
#'   by fraction (`"aromatic"`, `"saturate"`), each a list with a numeric
#'   `matrix` (rows = ion sums, columns = classes), a named `sensitivity`
#'   vector and a `carbon_number_range`.
"class_response_matrix"

#' Oxygenate-family response matrix (schema only, unpopulated)
#'
#' No ASTM standard provides oxygenate-family deconvolution coefficients;
#' this ships as an explicit `NA`-filled placeholder so [.solve_oxygenate_class()]
#' and [.correct_isotopic_tail()] are structurally testable ahead of real lab
#' calibration data. [.solve_oxygenate_class()] errors explicitly rather than
#' silently proceeding when the required matrix entries are `NA`.
#'
#' @format A list keyed `"lab"`, with a numeric `matrix` (rows = oxygenate ion
#'   sums, columns = oxygenate families), a named `sensitivity` vector, a
#'   `hydrocarbon_cross_terms` list (per-channel named coefficient vectors)
#'   and an `isotope_dominant_series` list (oxygenate channel -> dominant
#'   hydrocarbon ion-sum label).
"oxygenate_response_matrix"

#' Parent-ion isotope factors and mole sensitivities
#'
#' Per-carbon-number parent-ion (`mz`) isotope correction factor and mole
#' sensitivity for the alkylbenzene and naphthalene series, transcribed from
#' ASTM D2425-23 Table 1. Supports the standard's own mass-spectral average
#' carbon-number calculation (Sec 15.2-15.5, Appendix X1.2.2-X1.2.5), used by
#' [.mole_fraction_by_carbon()], [.average_carbon_number()] and
#' [.estimate_carbon_number()]. These are the only two series the standard
#' provides a formulaic carbon-number calculation for; other classes'
#' average carbon numbers are "estimated by inspection" per the standard and
#' are not automated here.
#'
#' @format A tibble with columns `series` (`"Alkylbenzene"` or
#'   `"Naphthalene"`), `carbon`, `mz` (parent-ion nominal mass),
#'   `isotope_factor`, `mole_sensitivity`.
"parent_ion_factors"
