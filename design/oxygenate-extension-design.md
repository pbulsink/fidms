# Design: Oxygenate-family deconvolution (extension beyond D2425/D2789)

Status: hydrocarbon-class core pipeline (Phases 1-4) is complete and shipped — see git history /
`R/hydrocarbon_type.R` for `read_ms_scan_matrix()`, `build_class_ion_traces()`,
`solve_hydrocarbon_class()`, `add_hydrocarbon_class()`, `estimate_carbon_number()`, and their
`prepare_peaks()`/report/app integration. This document covers only the **not-yet-validated
oxygenate extension** that plugs into that pipeline as a second, optional solve stage. It
supersedes and replaces the old `design/hydrocarbon-type-design.md`, keeping only the
still-relevant oxygenate material from it.

This is a joint coding + wet-lab document: several items below are blocked on lab calibration
work, not code, and are called out as such throughout.

## Provenance caveat

Neither D2425 nor D2789 covers oxygenates — D2425 explicitly limits scope to hydrocarbon types
(§1.2 excludes established applicability to oxygenate-containing samples), and D2789 states
applicability to oxygenate-containing gasolines "has not been established" (§1.2). Everything in
this document is a **user-proposed extension** to the existing matrix-deconvolution framework, not
transcribed from either standard. It must be tracked and validated as a distinct add-on module —
same architecture (summed diagnostic ions, NNLS-style matrix solve) but its own calibration and
evidence tier, never presented with the same authority as the D2425/D2789-derived hydrocarbon
matrices.

## Architectural placement

The oxygenate layer plugs in as a **second, optional solve stage** after the hydrocarbon
class-ion solve (`solve_hydrocarbon_class()`), consuming its output (particularly the resolved
`Alkylbenzene` and `Indan/Tetralin` fractions) as inputs to its own cross-term-corrected NNLS
solve, rather than being merged into one combined matrix.

```r
#' Solve non-negative oxygenate-family composition, correcting for hydrocarbon
#' and isotopic interference
#' @export
solve_oxygenate_class(
  ion_sums,
  hydrocarbon_fractions,   # named vector, output of solve_hydrocarbon_class()
  carbon_number,
  response_matrix = fidms::oxygenate_response_matrix,
  isotope_correction = TRUE
)
# 1. correct_isotopic_tail() on channels with a known dominant-series overlap
# 2. subtract hydrocarbon_fractions-scaled cross-terms per family (e.g. Alkylphenol
#    minus Alkylbenzene contribution) using response_matrix's cross-term columns
# 3. nnls::nnls() on the residual oxygenate-only ion-sum vector against the
#    oxygenate-family sub-matrix (including the acid/FAME 2x2+ sub-block)
```

Both `solve_oxygenate_class()` and `correct_isotopic_tail()` are implemented in
`R/hydrocarbon_type.R` and structurally tested in `tests/testthat/test-hydrocarbon-type.R` against
a synthetic, fully-populated placeholder matrix. What remains is: (a) populating
`oxygenate_response_matrix` with real coefficients (lab work), and (b) wiring the gating decision
into `prepare_peaks()` (code work) — neither is done yet.

Oxygenate assignment should be treated the same way as hydrocarbon-class assignment: a
chromatographically localized, evidence-tiered, matrix-deconvolved estimate, never presented with
more certainty than a library ID, and gated behind the same "only for peaks `classify_family()`
left as `Unknown Heteroatomic`/oxygen-flagged-but-unidentified" trigger condition as the
hydrocarbon path. Peaks already resolved as `"Oxygenate"` by `classify_family()` via a real library
hit or interpolated-formula oxygen count remain tier A/D per the existing precedence rules; this
new solver only applies where no name/formula evidence exists at all.

## `sim_masses_oxygenate` (implemented)

One row per family/ion-sum, plus one row per named compound for the Unsaturated FAME family
(confirmed molecular ions, not a fragment series):

| family | ion_sum_label | monitored_mz | hydrocarbon_interference |
|---|---|---|---|
| Alkylphenol | Sum94 | 94, 108, 122, 136, 150, 164 | Alkylbenzene |
| Alkoxyphenol (Guaiacol) | Sum109 | 109, 123, 137, 151, 165 | Alkylbenzene (C3-C5 benzylic region) |
| Saturated FAME | Sum74 | 74, 87, 143 | Branched alkane |
| Carboxylic acid | Sum60 | 60, 74, 88 | Alkane/alkene |
| Polyol (aromatic diol) | Sum110 | 110, 124, 138, 152 | Indan/Tetralin |
| Unsaturated FAME | (per-compound `mz`) | Methyl Palmitoleate 268, Methyl Oleate 296, Methyl Linoleate 294, Methyl Linolenate 292 | Diaromatic/naphthenoaromatic (molecular-ion-region) |

Notes on cross-terms baked into this table:

- `Alkylphenol` and `Alkoxyphenol` both overlap alkylbenzene fragment intensity (already summed as
  `Sum91` in `sim_masses`) at the same nominal masses — each needs its own cross-term coefficient,
  not a shared one, since the alkoxyphenol overlap is the narrower C3-C5 benzylic region rather
  than the general alkylbenzene series.
- `Saturated FAME` (Sum74) and `Carboxylic acid` (Sum60) interfere with **each other**, not just
  with a hydrocarbon class: free acids contribute strongly to m/z 60 and weakly to m/z 74; FAMEs
  contribute heavily to m/z 74 and negligibly to m/z 60. This needs a small 2x2 (or larger, if
  ethyl/other esters are added) sub-block within the oxygenate response matrix, independent of the
  hydrocarbon cross-terms; Sum74's minor branched-alkane interference is a third term.
- `Polyol (aromatic diol)` overlaps indane/tetralin intensity (`Sum103` in `sim_masses`) — same
  cross-term pattern as the alkylphenol/alkylbenzene case.
- `Unsaturated FAME` masses are confirmed `[M+•]` molecular ions of named compounds, not a
  homologous fragment series, so they are monitored individually (background-subtracted locally)
  rather than folded into the generic oxygenate class-fraction/NNLS output — closer to
  `calibration_targets`-style named-target detection than to the summed-series deconvolution used
  for the other five families.

## Isotopic-tail correction ("Benzylic Shift Calibration")

Alkylbenzene's dominant `Sum91`/105/119 series contributes an M+1/M+2 isotopic tail into the
alkylphenol `Sum94`/108/122 channels (and similarly into `Sum109`). This is a **different
correction mechanism** than the class-vs-class cross-term above — a natural-abundance ¹³C isotope
correction, not a co-fragmentation interference — modeled as a fixed (or carbon-number-scaled)
fractional subtraction from the alkylbenzene channel's own summed intensity, computed *before* the
NNLS solve:

```r
#' Subtract isotopic-tail contribution of a dominant series from an
#' oxygenate-family channel before NNLS solving
#' @export
correct_isotopic_tail(oxygenate_ion_sum, dominant_series_sum, isotope_fraction = 0.011)
# isotope_fraction: theoretical/empirical M+1 (~1.1% x #C per ion) and M+2 fraction,
# scaled by dominant_series_sum's magnitude; subtracted from oxygenate_ion_sum
# prior to solve_oxygenate_class().
```

Implemented in `R/hydrocarbon_type.R` with the ~1.1%-per-carbon default. **Not yet validated**:
the predicted tail needs to be checked against known hydrocarbon-only samples (no oxygenate
present) to see whether it matches the observed `Sum94`/`Sum109`-region background — a testable,
falsifiable lab step, required before trusting any `Alkylphenol`/`Alkoxyphenol` result derived from
the default fraction.

## `oxygenate_response_matrix`: schema shipped, coefficients unpopulated

```r
oxygenate_response_matrix <- list(
  lab = list(
    matrix = matrix(NA_real_, nrow = 5, ncol = 5,
      dimnames = list(c("Sum94", "Sum109", "Sum74", "Sum60", "Sum110"),
                       c("Alkylphenol", "Alkoxyphenol_guaiacol", "Saturated_FAME",
                         "Carboxylic_acid", "Polyol_aromatic_diol"))),
    sensitivity = c(...),  # NA placeholder, one entry per class
    hydrocarbon_cross_terms = list(
      Sum94  = c(Alkylbenzene = NA_real_),
      Sum109 = c(Alkylbenzene = NA_real_),
      Sum110 = c(Indan_Tetralin = NA_real_)
    ),
    isotope_dominant_series = list(Sum94 = "Sum91", Sum109 = "Sum91")
  )
)
```

This exists today in `data-raw/hydrocarbon_type.R` / `data/oxygenate_response_matrix.rda` as an
explicit `NA`-filled placeholder (confirmed by
`test_that("oxygenate_response_matrix ships as an explicit NA placeholder", ...)`), so
`solve_oxygenate_class()` raises a clear error rather than silently proceeding when called with the
shipped default. No standard provides these coefficients (unlike the hydrocarbon matrices, which
at least have D2425/D2789 as an authoritative source) — **this entire matrix must be built from the
lab's own calibration standards** (spiked known phenol/guaiacol/FAME/acid/polyol mixtures in a
representative hydrocarbon matrix), not transcribed from any external document.

## Remaining work

### Coding

1. **Gating decision (open).** Decide precisely which `classify_family()` outcomes route a peak
   into `solve_oxygenate_class()` — proposed: `"Unknown Heteroatomic"`, or `"Oxygenate"` only when
   derived from a low-confidence interpolated formula rather than a library hit — and wire it into
   `prepare_peaks()`/`add_hydrocarbon_class()`. Currently `solve_oxygenate_class()` is callable but
   nothing in the pipeline calls it.
2. **Sensitivity analysis harness (no lab data required).** Before committing lab time to
   calibrating every coefficient, build a small synthetic-perturbation script (`data-raw/` or a
   `dev/` scratch file, not exported package code) that:
   - perturbs each matrix entry, each hydrocarbon cross-term, and the isotope fraction across a
     plausible range against a fixed synthetic `ion_sums` vector,
   - reports how much resolved `fractions`/`residual` move per perturbed parameter,
   - is used to prioritize which coefficients are worth the lab effort to calibrate precisely
     versus which can start from a rough default.
3. Once coefficients exist (see lab work below), populate `oxygenate_response_matrix$lab` in
   `data-raw/hydrocarbon_type.R` and re-run the existing structural tests plus new
   accuracy/regression tests analogous to the Appendix X1 worked-example tests used for the
   hydrocarbon matrices.
4. Extend report/app surfacing (`total_compound_table()`, Excel writer, Shiny app) to display
   oxygenate-family results once (1)-(3) are done — mirroring the `include_class_ion` pattern
   already used for hydrocarbon-class columns, likely as an `include_oxygenate` parameter. Not
   started; low priority until the solve itself is validated.

### Lab

1. **Calibration standards.** Spiked known phenol/guaiacol/FAME/acid/polyol mixtures in a
   representative hydrocarbon matrix, run individually (family alone, no hydrocarbon background) to
   get diagonal response coefficients, then spiked into a representative hydrocarbon matrix to get
   the cross-terms empirically rather than assumed.
2. **Isotope-tail validation.** Run hydrocarbon-only (oxygenate-free) samples and compare observed
   `Sum94`/`Sum109`-region background against the ~1.1%-per-carbon default prediction from
   `correct_isotopic_tail()`; fit and supply an empirical `isotope_fraction` if the default is off.
3. Sequence lab work using the coding-side sensitivity analysis (item 2 above) to prioritize which
   coefficients most affect the resolved family split, rather than calibrating all ~15+ entries
   (5x5 matrix + 3 cross-terms + isotope fraction) with equal priority.

## Explicit non-goals / cautions

- The oxygenate-family deconvolution is an **unverified proposed extension**, not standard-derived
  like the hydrocarbon matrices — its response-matrix coefficients require original lab calibration
  work and should not be treated as ready-to-implement in the same sense as the D2425/D2789
  hydrocarbon tables.
- `assignment_tier`/`assignment_basis` must never be collapsed into the same visual weight as a
  library-confirmed `name` in reports — keep as distinct columns, not a replacement. This applies
  equally to any future oxygenate-family tier/basis columns.
- No FID-area partitioning by oxygenate-family fraction in this phase.
