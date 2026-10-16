# fidms - Agent Guidelines

R replacement for the "Automated Polyarc" Excel/VBA workflow that parses Agilent GC-FID and GC-MSD peak tables, quantifies compounds against a calibration table, classifies them into chemical families, and writes the client-facing Total Compound Report and Summary Report as Excel or PDF.

- **Package:** `fidms` v`0.0.0.9000`
- **Authors:** Philip Bulsink (`pbulsink@example.com`, cre/aut)
- **License:** MIT + file LICENSE
- **Testthat edition:** 3

## Package Overview

Replaces the "Automated Polyarc" Excel/VBA workflow. Parses Agilent FID and MSD peak tables plus MSD library-search exports, aligns the FID and MSD retention-time axes, identifies and interpolates compounds, quantifies them against a calibration table, classifies them into chemical families, and writes the client-facing Total Compound and Summary reports as Excel or PDF. Includes a local Shiny front end and a real-sample walkthrough vignette.

## Project Structure

### R/ (Source Code)
- `align.R` - FID/MSD retention-time alignment: shift search, greedy matching, MSD identity join, split-peak merging; holds processing params
- `align_gcalignr.R` - alternative GCalignR-based alignment engine for growing peak lists
- `app.R` - Shiny front-end launcher (local only by design)
- `calibrate.R` - calibration-curve fitting and Polyarc relative-response quantitation, recovery comparison
- `chemistry.R` - molecular formula parsing, molecular weight, retention-index model fit/predict
- `classify.R` - chemical-family classification and solvent flagging
- `fidms-package.R` - package-level docs: roxygen dataset documentation for all 11 packaged data objects; `.data` import from rlang
- `hydrocarbon_type.R` - ASTM D2425/D2789 class-ion deconvolution from raw scan matrices (largest module)
- `interpolate.R` - formula interpolation/extrapolation for unidentified peaks
- `parse.R` - Agilent CSV readers (FID/MSD peak tables, library-search exports)
- `pipeline.R` - end-to-end orchestration and Total Compound table construction
- `reference_peak.R` - internal-standard reference peak suggestion and quality check
- `report.R` - Excel and Quarto-PDF client report writers
- `summarize.R` - sample metadata object; family/element/recovery summary

### data/ (Packaged Data)
11 `.rda` datasets, all documented in `R/fidms-package.R`:
- `atomic_weights`, `calibration_targets`, `calibration_points`, `chemical_family_terms`, `solvent_names` - lookups rebuilt from the old workbook by `data-raw/lookup_tables.R` / `data-raw/atomic_weights.R`
- `rical_reference_points` - RT/RI anchors; the RI polynomial is refit at runtime, not shipped
- `sim_masses`, `class_response_matrix`, `parent_ion_factors` - ASTM D2425-23 transcriptions (Appendix X1, Tables 1/5) for class-ion deconvolution
- `sim_masses_oxygenate`, `oxygenate_response_matrix` - user-proposed oxygenate extension; the response matrix is an explicit NA-filled placeholder pending lab calibration data

### data-raw/ (Data Rebuild Scripts)
- `lookup_tables.R` - rebuilds packaged lookup tables from `caltable.csv` and `rical.csv` (both now live in `data-raw/`, not `reference/`)
- `atomic_weights.R`, `hydrocarbon_type.R` - source the atomic-weight table and the ASTM-derived class-ion matrices
- `capture_app_screenshots.R` - dev-only shinytest2 script that regenerates the app screenshots in `vignettes/images/`; run from the package root, needs local `diesel.D`
- `render_pdf_page.R` - dev-only; renders the example diesel PDF report and saves page 1 to `vignettes/images/07_pdf_report_page1.png` via `pdftools::pdf_convert()` (needs Quarto + LaTeX)

### inst/ (App and Report Template)
- `app/app.R` - local Shiny app (serves via `run_app()`); client data never leaves the machine, do not add deployment targets
- `report/report.qmd` - Quarto template rendered by `write_pdf_report()`

### vignettes/ (Vignette)
- `real_sample_workflow.Rmd` - "Processing a Real GC-MS/FID Sample": end-to-end walkthrough of the example diesel sample in code and through the Shiny app. Built with knitr (`VignetteBuilder: knitr`, html_vignette output); pre-built outputs are committed alongside it (`real_sample_workflow.html`, `Processing a Real GC-MS_FID Sample.pdf`)
- Bundles the **full** example diesel peak-table CSVs (`FIDpeaks.csv`, `MSDpeaks.csv`, `MSDsearch.csv`) — complete sample runs, not the 10-peak test fixtures in `tests/testthat/fixtures/`
- `images/` - app screenshots (01–06) and page 1 of the generated PDF report (07)

Building the vignette:
- Rebuild with `devtools::build_vignettes()` (or any R CMD build; `knitr` + `rmarkdown` are in Suggests). The class-ion section only evaluates when a example diesel `.D` directory is found (`FIDMS_diesel` env var, or `../diesel.D`, or `diesel.D`)
- App screenshots: run `data-raw/capture_app_screenshots.R` from the package root — it drives a throwaway launcher (working-tree app via `pkgload::load_all()`) with shinytest2's `AppDriver` against local `diesel.D`, and rewrites all of `vignettes/images/01–06`. Rerun after any UI change
- PDF first page: run `data-raw/render_pdf_page.R` from the package root after `devtools::load_all()` — needs Quarto + LaTeX; uses dev-only (undeclared) dependency `pdftools` to rasterize report page 1 at dpi 110

Note: roxygen comments in `R/hydrocarbon_type.R` and `R/pipeline.R` cite `design/hydrocarbon-type-design.md`, which is not currently present in the repo; treat those citations as historical.

### tests/testthat/ (Unit Tests)
- `test-align.R`, `test-align-gcalignr.R` - alignment engines against fixture CSVs
- `test-calibrate.R` - curve fitting, target assignment, quantitation styles, recovery
- `test-chemistry.R` - formula parsing, molecular weight, retention-index model
- `test-classify.R` - chemical-family classification
- `test-hydrocarbon-type.R` - class-ion traces, solves, carbon-number estimation (largest test file)
- `test-interpolate.R`, `test-parse.R`, `test-reference-peak.R` - per-module coverage
- `test-pipeline.R` - end-to-end runs on the fixtures plus report-table layout assertions
- `fixtures/` - 10-peak sample excerpts: `FIDpeaks.csv`, `MSDpeaks.csv`, `MSDsearch.csv` (copied from `reference/`; these are excerpts, not full sample runs)

### Local-only files (gitignored — present on this machine, absent in fresh clones or CI)
- `reference/` - ASTM D2425-23/D2789-95/D8276 PDFs with `.png` screenshots of every table, the original `Automated Polyarc.xlsm`, 10-peak sample CSV excerpts (`FIDpeaks.csv`, etc.), `MD_diesel_GC-MS-FID.xlsx` (the macro's real output), and other reference material
- `diesel.D/` - a local Agilent `.D` raw-scan directory (diesel fuel sample) usable to exercise the class-ion pipeline via `.read_ms_scan_matrix()`. The original `diesel.D/` is no longer on this machine but is required by the vignette's class-ion section and by `data-raw/capture_app_screenshots.R`; restore it (or set `FIDMS_diesel`) to regenerate those
- `dev/` - local-only development area: `oxygenate-extension-design.md` + `oxygenate-extension-record.md` (moved out of the former tracked `design/` dir), and `oxygenate_insilico_overlap.R`, a scratch analysis over `data-raw/combined_library.msp` (also gitignored)

## Architecture Summary

The package is a linear pipeline mirroring the VBA macro: parse CSVs, align FID/MSD retention times, attach MSD library identities, merge split peaks, interpolate unknown formulae, quantify against calibration (or an internal reference peak), summarise into family/element/recovery totals, and render Excel or PDF reports. The hydrocarbon-class module is a side path that runs between identification and interpolation when raw scan data (`class_traces`) is supplied. Most functions are internal (dot-prefixed, `@keywords internal`); the public API is deliberately small.

| File | Role | Key Functions |
|------|------|---------------|
| `align.R` | FID/MSD alignment: shift search, greedy matching, identity join, split-peak merging | `fidms_params()`, `.align_peaks()`, `.merge_close_peaks()` |
| `align_gcalignr.R` | GCalignR-based alternative alignment engine (guarded; used only via `engine = "gcalignr"`) | `.align_peaks_gcalignr()` |
| `app.R` | Shiny front-end launcher, local-only by design | `run_app()` |
| `calibrate.R` | Calibration-curve fitting and Polyarc relative-response quantitation with recovery comparison | `calibrate()`, `.fit_calibration()`, `.assign_target_rows()` |
| `chemistry.R` | Formula parsing, molecular weight, retention-index model (refit per run) | `.parse_formula()`, `.molecular_weight()`, `.fit_retention_index()` |
| `classify.R` | Chemical-family classification and solvent flagging | `.classify_family()`, `.is_solvent()` |
| `fidms-package.R` | Dataset documentation for the 11 packaged data objects | (data only) |
| `hydrocarbon_type.R` | ASTM class-ion deconvolution from raw scan matrices; parent-ion carbon-number estimation | `hydrocarbon_class_params()`, `.add_hydrocarbon_class()`, `.solve_hydrocarbon_class()` |
| `interpolate.R` | Formula interpolation/extrapolation for unresolved peaks, with class-evidence skip logic | `.interpolate_unknowns()` |
| `parse.R` | Agilent CSV readers for FID/MSD peak tables and library-search exports | `.read_fid_peaks()`, `.read_msd_peaks()`, `.read_msd_search()` |
| `pipeline.R` | End-to-end orchestration; Total Compound table with optional class-ion columns | `process_sample()`, `prepare_peaks()`, `total_compound_table()` |
| `reference_peak.R` | Internal-standard reference peak suggestion and quality gate | `.suggest_reference_peak()`, `.reference_peak_ok()` |
| `report.R` | Two-sheet Excel writer; Quarto-rendered PDF writer | `write_excel_report()`, `write_pdf_report()` |
| `summarize.R` | Sample metadata object; family/element/recovery summary and QC counts | `summarize_report()`, `sample_info()` |

**Dependencies:**
- **Imports:** chromConverter, dplyr, GCalignR, nnls, openxlsx, rlang, stats, stringr, tibble, tidyr, utils
- **Suggested (optional features):** bslib, DT, knitr, quarto, rmarkdown, shiny, testthat (>= 3.0.0), withr

**Dependency guards in code:** `GCalignR` (`rlang::check_installed`, only when `engine = "gcalignr"`), `chromConverter` and `nnls` (`requireNamespace()`, class-ion path), `shiny`/`bslib`/`DT` (app) and `quarto` (PDF reports). Tests of the Excel writer use `skip_if_not_installed("openxlsx")`; `withr` is used for temp files.

## Standard Workflow

This section describes the typical development and testing workflow for this package. Agents should follow these steps when making changes:

### 1. Understand the Code Path
- Read roxygen documentation in `R/*.R`; internal functions are dot-prefixed with `@keywords internal`, exported ones are listed in NAMESPACE
- Trace data flow through the pipeline (parse → align → identify → merge → interpolate → calibrate → summarize → write) or through the class-ion side path (`dev/oxygenate-extension-design.md`, local-only, and the roxygen docs in `R/hydrocarbon_type.R` describe the evidence-tier A/B/C/D/M logic)
- Identify which files are affected by the proposed change

### 2. Make Changes
- Modify only the relevant file(s) in `R/`. Do not create new `.R` files unless adding a completely new capability that does not fit within existing modules
- Keep internal helpers dot-prefixed and marked `@keywords internal`; follow the existing style (tidyverse verbs with explicit namespace, base R otherwise)
- Use plain `stop("message", call. = FALSE)` for error messages, matching the codebase
- Preserve the roxygen cross-references (e.g. `[.align_peaks()]`) when editing docs

### 3. Add Tests
- Create or update tests in `tests/testthat/test-{module}.R`, matching the source file under test (`test-align-gcalignr.R` covers `R/align_gcalignr.R`)
- Place new tests adjacent to existing tests for the same function, or at the end of the file if no tests exist yet
- Reuse the 10-peak fixtures in `tests/testthat/fixtures/` via `test_path("fixtures", ...)`; build synthetic tibbles when a fixture cannot express the case (see `test-pipeline.R`)
- Run `devtools::test()` to verify all tests pass
- Ensure new exported functions have usage examples in roxygen `@examples` blocks

### 4. Update Documentation
- Regenerate man pages and NAMESPACE: `devtools::document()`
- If a packaged dataset changes, rebuild it with the matching `data-raw/` script rather than editing `.rda` files by hand, then re-run `devtools::document()` for its dataset docs in `R/fidms-package.R`
- Update the README pipeline diagram or VBA mapping table if user-facing behavior changes
- If app UI or report layout changes, regenerate `vignettes/images/` with `data-raw/capture_app_screenshots.R` and/or `data-raw/render_pdf_page.R`, then rebuild the vignette (`devtools::build_vignettes()`)

### 5. Verify Build
- Full package check: `devtools::check()` (or `R CMD check`)
- Review test coverage for edited files; aim for high coverage on new code (`cobertura.xml` at repo root is a prior covr export)
- Ensure no regressions in existing tests and that all file references remain valid

## Testing Philosophy

- Test outputs and behavior, not just function existence: expect exact values (sheet names, column counts, sentinel strings like `"Combined FID Peaks"`, recovery sums to 100 within tolerance) rather than mere non-error
- Fixtures are deliberately small 10-peak excerpts; end-to-end tests run the full pipeline on them with `style = "internal"` because the excerpt contains no calibration targets (exercising the macro's "no targets identified" branch)
- Failing tests are valuable — address failures rather than suppressing them

## Project Memory: PDF table extraction (ASTM D2425 / D2789 standards)

The ASTM standard PDFs in the repo reference folder (`D2425-23.pdf`, `D2789-95.pdf`, etc.)
contain complex tables (overlapping/merged cells, multi-row headers) that
`pdftotext -layout` extracts heavily garbled for — column/row alignment is not
reliable for these specific tables, even though it works fine for prose
sections. When a future task needs to re-verify or extend transcribed
table data from these PDFs (e.g. `class_response_matrix`, `sim_masses`),
don't rely solely on `pdftotext -layout`. Steps tried, and their outcomes,
when investigating this on this machine:

- `pdftotext -layout` on the table pages: garbled column alignment for wide
  multi-column tables (confirmed unreliable for Table 3/4/5 and Appendix X1
  matrices specifically).
- Checked for PDF-to-image tools to visually inspect the raw page layout
  instead of relying on text extraction:
  - `pdftoppm`, `pdftocairo`, `mutool`: **not installed** on this machine.
  - `python3 -c "import fitz"` (PyMuPDF): **not available**.
  - ImageMagick (`magick`) **is installed**, but has **no Ghostscript
    delegate**, so `magick -density 200 "file.pdf[N]" out.png` fails --
    rendering PDF pages to images is not currently possible on this machine
    without installing Ghostscript (or PyMuPDF) first.
- Fallback approach that worked reasonably well: cross-reference plain
  `pdftotext` (no `-layout`, natural reading order) against the `-layout`
  output for the same page range (`pdftotext -f N -l N file.pdf -`), using
  the reading-order stream to disambiguate which numbers belong to which
  column when `-layout`'s alignment is ambiguous. This is still
  best-effort, not a substitute for visual page inspection.
- There are `.png` files present in the references folder containing screenshots
  of all tables and table-like structures in the D2425 pdf document. These 
  can be reviewed if required.

**Recommendation for next time:** if a task requires high-confidence
transcription of one of these tables, ask the user whether Ghostscript (or
`poppler-utils`' `pdftoppm`/`pdftocairo`) can be installed, or have them
supply page images/screenshots directly, rather than spending effort on
further text-extraction workarounds.

## Deliberate Design Decisions and Open Items

From the README — do not "fix" these without confirmation:

- The stored degree-6 RI polynomial and stored calibration coefficients/R² from the workbook are **not** ported; both are refit at runtime from reference points
- Calibration levels with a zero response are below-detection points, excluded from fits rather than treated as observed zeros (confirmed with the analyst)
- `engine = "grid"` vs `"gcalignr"`: the two engines are not expected to agree exactly and do not on the 10-peak excerpt; the production default is still an open item pending comparison on a full real sample
- Golden-file regression against the macro's output for a full real sample is **not yet done** (the repo CSVs are excerpts, not the full run behind `MD_diesel_GC-MS-FID.xlsx`)
- Deferred: calibration-report import and local persistence
