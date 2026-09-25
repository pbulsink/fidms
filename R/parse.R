#' Read an Agilent peak or search table CSV
#'
#' The exports carry a two-line banner (`#Peaks: ...` / `#`) before the real
#' header row and may be padded with hundreds of empty trailing rows.
#'
#' @param path Path to the CSV file.
#' @return A tibble with the banner and padding removed.
#' @keywords internal
.read_agilent_csv <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^\ufeff", "", lines)
  header <- which(!startsWith(trimws(lines), "#") & nzchar(trimws(gsub(",", "", lines))))[1]
  if (is.na(header)) {
    stop("No header row found in ", path, call. = FALSE)
  }
  out <- utils::read.csv(
    text = paste(lines[header:length(lines)], collapse = "\n"),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  out <- out[, !is.na(names(out)) & nzchar(names(out)), drop = FALSE]
  .drop_blank_rows(tibble::as_tibble(out))
}

.drop_blank_rows <- function(df) {
  if (nrow(df) == 0) return(df)
  blank <- apply(df, 1, function(r) all(is.na(r) | !nzchar(trimws(as.character(r)))))
  df[!blank, , drop = FALSE]
}

#' Read an FID peak table
#'
#' @param path Path to an `FIDpeaks.csv` export.
#' @return A tibble with columns `peak`, `rt_start`, `rt`, `rt_end`, `area`,
#'   `height`, sorted by retention time.
#' @keywords internal
.read_fid_peaks <- function(path) {
  raw <- .read_agilent_csv(path)
  required <- c("Peak", "Start X", "Center X", "End X", "Area", "Height")
  .check_columns(raw, required, "FID peak table")
  out <- tibble::tibble(
    peak = as.integer(raw[["Peak"]]),
    rt_start = as.numeric(raw[["Start X"]]),
    rt = as.numeric(raw[["Center X"]]),
    rt_end = as.numeric(raw[["End X"]]),
    area = as.numeric(raw[["Area"]]),
    height = as.numeric(raw[["Height"]])
  )
  out <- out[!is.na(out$rt), , drop = FALSE]
  out[order(out$rt), , drop = FALSE]
}

#' Read an MSD peak table
#'
#' @param path Path to an `MSDpeaks.csv` export.
#' @return A tibble with columns `peak`, `name`, `rt_start`, `rt`, `rt_end`,
#'   `area`, `height`, sorted by retention time.
#' @keywords internal
.read_msd_peaks <- function(path) {
  raw <- .read_agilent_csv(path)
  required <- c("Peak", "Start X", "Center X", "End X", "Area", "Height")
  .check_columns(raw, required, "MSD peak table")
  out <- tibble::tibble(
    peak = as.integer(raw[["Peak"]]),
    name = if ("Name" %in% names(raw)) as.character(raw[["Name"]]) else NA_character_,
    rt_start = as.numeric(raw[["Start X"]]),
    rt = as.numeric(raw[["Center X"]]),
    rt_end = as.numeric(raw[["End X"]]),
    area = as.numeric(raw[["Area"]]),
    height = as.numeric(raw[["Height"]])
  )
  out <- out[!is.na(out$rt), , drop = FALSE]
  out[order(out$rt), , drop = FALSE]
}

#' Read an MSD library search export
#'
#' @param path Path to an `MSDsearch.csv` export.
#' @return A tibble with columns `cpd`, `name`, `formula`, `rt`, `area`,
#'   `height`, `ri_lib`, `score`, `cas`.
#' @keywords internal
.read_msd_search <- function(path) {
  raw <- .read_agilent_csv(path)
  .check_columns(raw, c("Formula", "RT", "Score"), "MSD search table")
  pick <- function(nm) if (nm %in% names(raw)) raw[[nm]] else NA
  out <- tibble::tibble(
    cpd = suppressWarnings(as.integer(pick("Cpd"))),
    name = trimws(as.character(pick("Name"))),
    formula = trimws(as.character(raw[["Formula"]])),
    rt = as.numeric(raw[["RT"]]),
    area = suppressWarnings(as.numeric(pick("Area"))),
    height = suppressWarnings(as.numeric(pick("Height"))),
    ri_lib = suppressWarnings(as.numeric(pick("RI"))),
    score = suppressWarnings(as.numeric(raw[["Score"]])),
    cas = as.character(pick("CAS ID"))
  )
  out <- out[!is.na(out$rt) & !is.na(out$formula) & nzchar(out$formula), , drop = FALSE]
  out[order(out$rt), , drop = FALSE]
}

.check_columns <- function(df, required, what) {
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(
      what, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
