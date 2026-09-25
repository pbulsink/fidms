## Atomic weights replacing the workbook's hidden `PeriodicTable` sheet.
## IUPAC 2021 conventional atomic weights, g/mol.
library(tibble)

atomic_weights <- tibble::tribble(
  ~symbol, ~element,   ~mass,
  "H",     "Hydrogen", 1.008,
  "B",     "Boron",    10.81,
  "C",     "Carbon",   12.011,
  "N",     "Nitrogen", 14.007,
  "O",     "Oxygen",   15.999,
  "F",     "Fluorine", 18.998403,
  "Si",    "Silicon",  28.085,
  "P",     "Phosphorus", 30.973762,
  "S",     "Sulfur",   32.06,
  "Cl",    "Chlorine", 35.45,
  "Br",    "Bromine",  79.904,
  "I",     "Iodine",   126.90447
)

usethis::use_data(atomic_weights, overwrite = TRUE)
