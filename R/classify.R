#' Classify peaks into chemical families
#'
#' Port of the VBA `ChemicalFamily` logic: name-fragment lookups first, then
#' formula-based degree-of-unsaturation heuristics, falling back to
#' `"Unknown Hydrocarbon"` / `"Unknown Heteroatomic"`.
#'
#' @param peaks Peak table with `name`, `formula` and element-count columns.
#' @param terms Name-fragment lookup table, default [chemical_family_terms].
#' @return Character vector of family labels, one per peak.
#' @keywords internal
.classify_family <- function(peaks, terms = chemical_family_terms) {
  nm <- tolower(ifelse(is.na(peaks$name), "", peaks$name))
  hetero <- rowSums(
    peaks[, c("N", "F", "Cl", "Br", "I", "P", "S", "B", "Si")],
    na.rm = TRUE
  ) >
    0
  oxy <- !is.na(peaks$O) & peaks$O > 0
  dbe <- peaks$C - peaks$H / 2 + 1 # rings + double bonds

  family <- rep(NA_character_, nrow(peaks))

  name_family <- function(f) {
    pat <- paste0(terms$term[terms$family == f], collapse = "|")
    grepl(pat, nm)
  }
  family[is.na(family) & name_family("Aromatic")] <- "Aromatic"
  family[is.na(family) & name_family("Naphthene")] <- "Naphthene"

  # formula heuristics for anything still unlabelled
  unk <- is.na(family)
  family[unk & oxy] <- "Oxygenate"
  unk <- is.na(family)
  family[unk & hetero] <- "Unknown Heteroatomic"
  unk <- is.na(family)
  family[unk & !is.na(dbe) & dbe >= 4] <- "Aromatic"
  unk <- is.na(family)
  family[unk & !is.na(dbe) & dbe >= 1] <- "Unsaturate"
  unk <- is.na(family)
  family[unk & !is.na(dbe) & dbe <= 0] <- ifelse(
    grepl("methyl|ethyl|propyl|iso|,", nm[unk & !is.na(dbe) & dbe <= 0]),
    "Isoparaffin",
    "Paraffin"
  )
  family[is.na(family)] <- "Unknown Hydrocarbon"

  # oxygenate always wins per the report note "includes all compounds with an
  # oxygen atom, including methanol"
  family[oxy] <- "Oxygenate"
  family
}

#' Flag known solvent peaks
#'
#' @param peaks Peak table with a `name` column.
#' @param solvents Character vector of solvent names, default [solvent_names].
#' @return Logical vector.
#' @keywords internal
.is_solvent <- function(peaks, solvents = solvent_names) {
  key <- function(x) gsub("[^a-z0-9]", "", tolower(ifelse(is.na(x), "", x)))
  key(peaks$name) %in% key(solvents)
}
