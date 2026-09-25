#' Write the client-facing Excel report
#'
#' Reproduces the two-sheet workbook the macro saved: `{id} Total Compound
#' Report` and `{id} Summary Report`.
#'
#' @param result Output of [process_sample()] / [summarize_report()].
#' @param path Output file path. Defaults to
#'   `{ProjectCode}_{SampleID}_GC-MS-FID.xlsx` in the working directory.
#' @return The path written, invisibly.
#' @export
write_excel_report <- function(result, path = NULL) {
  info <- result$info
  if (is.null(path)) {
    path <- sprintf("%s_%s_GC-MS-FID.xlsx", info$project_code, info$sample_id)
  }
  wb <- openxlsx::createWorkbook()
  date_style <- openxlsx::createStyle(numFmt = "yyyy-mm-dd")
  bold <- openxlsx::createStyle(textDecoration = "bold")

  # ---- Total Compound Report ------------------------------------------------
  s1 <- paste(info$sample_id, "Total Compound Report")
  openxlsx::addWorksheet(wb, s1)
  t <- result$totals
  q <- result$qc
  meta <- data.frame(
    A = c("Sample ID", "Sample Name", "Project Code", "Report Date", "Analyst"),
    B = c(
      info$sample_id,
      info$sample_name,
      info$project_code,
      NA,
      info$analyst
    ),
    C = c("Density", "Methanol", "Water", "GC Recovery", "Total:"),
    D = c(
      info$density,
      info$methanol_wt_pct,
      info$water_wt_pct,
      round(t$gc_recovery, 1),
      round(t$total, 1)
    ),
    E = c("g/mL", "wt. %", "wt. %", "wt. %", "wt. %"),
    F = c(
      "Aligned Peaks",
      "Unique Peaks",
      "Total Peaks",
      "Interpolated Peaks",
      "MS/FID Offset Time"
    ),
    G = c(
      q$aligned_peaks,
      q$unique_peaks,
      q$total_peaks,
      q$interpolated_peaks,
      round(q$ms_fid_offset, 3)
    ),
    H = c(NA, NA, NA, NA, "min"),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, s1, meta, colNames = FALSE)
  openxlsx::writeData(wb, s1, info$report_date, startRow = 4, startCol = 2)
  openxlsx::addStyle(wb, s1, date_style, rows = 4, cols = 2)
  openxlsx::writeData(wb, s1, total_compound_table(result), startRow = 6)
  openxlsx::addStyle(wb, s1, bold, rows = 6, cols = 1:9, gridExpand = TRUE)
  openxlsx::setColWidths(wb, s1, cols = 1:9, widths = "auto")

  # ---- Summary Report -------------------------------------------------------
  s2 <- paste(info$sample_id, "Summary Report")
  openxlsx::addWorksheet(wb, s2)
  head_block <- data.frame(
    A = c(
      "Sample ID",
      "Sample Name",
      "Project Code",
      "Report Date",
      "Analyst",
      "Density (g/mL)"
    ),
    B = c(
      info$sample_id,
      info$sample_name,
      info$project_code,
      NA,
      info$analyst,
      info$density
    ),
    stringsAsFactors = FALSE
  )
  openxlsx::writeData(wb, s2, head_block, colNames = FALSE)
  openxlsx::writeData(wb, s2, info$report_date, startRow = 4, startCol = 2)
  openxlsx::addStyle(wb, s2, date_style, rows = 4, cols = 2)

  fam <- rbind(
    data.frame(
      a = result$families$compound_type,
      b = round(result$families$wt_pct, 1)
    ),
    data.frame(
      a = c("Methanol (external)", "Water", "TOTAL"),
      b = round(c(t$methanol, t$water, t$total), 1)
    )
  )
  names(fam) <- c("Compound Type", "Total Amount (wt. %)")
  openxlsx::writeData(wb, s2, fam, startRow = 7)

  el <- data.frame(
    Element = c(result$elements$element, "Sum"),
    `Calculated wt. %` = round(c(result$elements$wt_pct, t$element_sum), 2),
    check.names = FALSE
  )
  openxlsx::writeData(wb, s2, el, startRow = 7, startCol = 3)
  openxlsx::writeData(
    wb,
    s2,
    c("NOTE:", .element_note()),
    startRow = 7,
    startCol = 6,
    colNames = FALSE
  )
  openxlsx::addStyle(wb, s2, bold, rows = 7, cols = c(1:4), gridExpand = TRUE)

  notes_row <- 7 + nrow(fam) + 2
  openxlsx::writeData(wb, s2, "NOTES:", startRow = notes_row, startCol = 1)
  openxlsx::writeData(
    wb,
    s2,
    paste0(seq_along(.summary_notes()), ". ", .summary_notes()),
    startRow = notes_row + 1,
    startCol = 1,
    colNames = FALSE
  )
  openxlsx::setColWidths(wb, s2, cols = 1:6, widths = "auto")

  # ---- MS Class Features (optional) -----------------------------------------
  # Only present when `result` came from a run with `class_traces` supplied
  # (see `design/hydrocarbon-type-design.md`); an auditable detail table, kept
  # off the main Total Compound sheet.
  feats <- result$ms_class_features
  if (!is.null(feats) && nrow(feats) > 0) {
    s3 <- "MS Class Features"
    openxlsx::addWorksheet(wb, s3)
    openxlsx::writeData(wb, s3, feats)
    openxlsx::addStyle(
      wb,
      s3,
      bold,
      rows = 1,
      cols = seq_len(ncol(feats)),
      gridExpand = TRUE
    )
    openxlsx::setColWidths(wb, s3, cols = seq_len(ncol(feats)), widths = "auto")
  }

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(path)
}

#' Render the client-facing PDF report
#'
#' Renders the packaged Quarto template with the processed result.
#'
#' @param result Output of [process_sample()].
#' @param path Output PDF path.
#' @return The path written, invisibly.
#' @export
write_pdf_report <- function(result, path = NULL) {
  rlang::check_installed("quarto")
  info <- result$info
  if (is.null(path)) {
    path <- sprintf("%s_%s_GC-MS-FID.pdf", info$project_code, info$sample_id)
  }
  template <- system.file("report", "report.qmd", package = "fidms")
  if (!nzchar(template)) {
    template <- "inst/report/report.qmd"
  }
  rds <- tempfile(fileext = ".rds")
  saveRDS(result, rds)
  quarto::quarto_render(
    input = template,
    output_file = basename(path),
    execute_params = list(result_rds = rds)
  )
  rendered <- file.path(dirname(template), basename(path))
  if (rendered != path) {
    file.rename(rendered, path)
  }
  invisible(path)
}
