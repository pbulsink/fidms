#' Parse a molecular formula into element counts
#'
#' Replaces the VBA `ChemRegex` user-defined function. Handles Hill-notation
#' formulae such as `C6H12`, `C2H6O`, `CH3Cl`, and bracketed/hydrated forms are
#' flattened by simple expansion of a single level of parentheses.
#'
#' @param formula Character vector of molecular formulae.
#' @return A tibble with one row per input formula and one integer column per
#'   element tracked (`C`, `H`, `N`, `O`, `F`, `Cl`, `Br`, `I`, `P`, `S`, `B`,
#'   `Si`).
#' @keywords internal
.parse_formula <- function(formula) {
  elements <- .tracked_elements()
  out <- matrix(
    0L,
    nrow = length(formula), ncol = length(elements),
    dimnames = list(NULL, elements)
  )
  for (i in seq_along(formula)) {
    f <- formula[i]
    if (is.na(f) || !nzchar(trimws(f))) {
      out[i, ] <- NA_integer_
      next
    }
    counts <- .parse_formula_one(trimws(f))
    known <- intersect(names(counts), elements)
    out[i, known] <- as.integer(counts[known])
  }
  tibble::as_tibble(out)
}

.parse_formula_one <- function(f) {
  f <- .expand_groups(f)
  m <- stringr::str_match_all(f, "([A-Z][a-z]?)(\\d*)")[[1]]
  if (nrow(m) == 0) return(stats::setNames(integer(0), character(0)))
  el <- m[, 2]
  n <- ifelse(m[, 3] == "", 1L, suppressWarnings(as.integer(m[, 3])))
  tapply(n, el, sum)
}

.expand_groups <- function(f) {
  pattern <- "\\(([^()]*)\\)(\\d*)"
  while (stringr::str_detect(f, pattern)) {
    m <- stringr::str_match(f, pattern)
    mult <- if (m[3] == "") 1L else as.integer(m[3])
    inner <- stringr::str_match_all(m[2], "([A-Z][a-z]?)(\\d*)")[[1]]
    expanded <- paste0(
      inner[, 2],
      ifelse(inner[, 3] == "", 1L, as.integer(inner[, 3])) * mult,
      collapse = ""
    )
    f <- sub(pattern, expanded, f)
  }
  f
}

.tracked_elements <- function() {
  c("C", "H", "N", "O", "F", "Cl", "Br", "I", "P", "S", "B", "Si")
}

#' Molecular weight from a formula
#'
#' Replaces the VBA `molweight()` function and its `PeriodicTable` sheet.
#'
#' @param formula Character vector of molecular formulae.
#' @return Numeric vector of molar masses in g/mol; `NA` for unparseable input.
#' @keywords internal
.molecular_weight <- function(formula) {
  counts <- .parse_formula(formula)
  w <- atomic_weights[match(names(counts), atomic_weights$symbol), "mass", drop = TRUE]
  as.numeric(as.matrix(counts) %*% w)
}

#' Fit a retention-time to retention-index model
#'
#' The workbook stored a degree-6 polynomial; per the port design the fit is
#' recomputed from the reference points each run.
#'
#' @param reference Data frame of reference points with numeric columns `time`
#'   and `ri`. Defaults to the packaged [rical_reference_points].
#' @param degree Polynomial degree (default 6).
#' @return An object of class `fidms_ri_model`.
#' @keywords internal
.fit_retention_index <- function(reference = rical_reference_points, degree = 6L) {
  ref <- reference[!is.na(reference$time) & !is.na(reference$ri), , drop = FALSE]
  if (nrow(ref) < degree + 1L) {
    degree <- max(1L, nrow(ref) - 1L)
  }
  fit <- stats::lm(ri ~ stats::poly(time, degree, raw = TRUE), data = ref)
  structure(
    list(fit = fit, degree = degree, range = range(ref$time)),
    class = "fidms_ri_model"
  )
}

#' Compute retention indices for observed retention times
#'
#' @param rt Numeric vector of retention times (minutes).
#' @param model A model from [.fit_retention_index()].
#' @return Numeric vector of retention indices.
#' @keywords internal
.compute_retention_index <- function(rt, model = .fit_retention_index()) {
  stopifnot(inherits(model, "fidms_ri_model"))
  as.numeric(stats::predict(model$fit, newdata = data.frame(time = rt)))
}
