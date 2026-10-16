## Builds packaged lookup tables from the staged CSV exports of the workbook's
## `Calibration`, `RICal` and `ChemicalFamily` sheets.
library(dplyr)
library(tibble)

source("data-raw/atomic_weights.R")

# ---- calibration targets ----------------------------------------------------
raw_cal <- readr::read_csv("caltable.csv", show_col_types = FALSE)
raw_cal <- raw_cal[!is.na(raw_cal$Target) & nzchar(raw_cal$Target), ]

levels_ <- c(
  "10000",
  "5000",
  "2500",
  "1000",
  "500",
  "250",
  "100",
  "50",
  "25",
  "10"
)

# Long table of standard-curve points. A zero response is a below-detection
# level with no calibration-curve entry (confirmed with the analyst), not an
# observed blank, so those points are dropped rather than fit as (x, 0).
points <- purrr::map_dfr(levels_, function(lv) {
  tibble::tibble(
    target = raw_cal$Target,
    level = as.numeric(lv),
    conc = suppressWarnings(as.numeric(raw_cal[[paste0(lv, "Conc")]])),
    resp = suppressWarnings(as.numeric(raw_cal[[paste0(lv, "Resp")]]))
  )
})
blank <- tibble::tibble(
  target = raw_cal$Target,
  level = 0,
  conc = suppressWarnings(as.numeric(raw_cal[["0 Point"]])),
  resp = suppressWarnings(as.numeric(raw_cal[["0R Point"]]))
)
calibration_points <- dplyr::bind_rows(points, blank) |>
  dplyr::filter(!is.na(conc), !is.na(resp)) |>
  dplyr::filter(resp > 0, conc > 0) |>
  dplyr::arrange(target, level)

calibration_targets <- tibble::tibble(
  target = raw_cal$Target,
  synonyms = raw_cal$Synonyms,
  type = tolower(trimws(raw_cal$Type)),
  rt = suppressWarnings(as.numeric(raw_cal$RT)),
  ri = suppressWarnings(as.numeric(raw_cal$RI)),
  carbon = suppressWarnings(as.integer(raw_cal$Carbon)),
  hydrogen = suppressWarnings(as.integer(raw_cal$Hydrogen)),
  oxygen = tidyr::replace_na(suppressWarnings(as.integer(raw_cal$Oxygen)), 0L)
)

# ---- RI reference points ----------------------------------------------------
raw_ri <- readr::read_csv("rical.csv", show_col_types = FALSE)
rical_reference_points <- tibble::tibble(
  carbon = suppressWarnings(as.integer(raw_ri$Carbon)),
  hydrogen = suppressWarnings(as.integer(raw_ri$Hydrogen)),
  oxygen = tidyr::replace_na(suppressWarnings(as.integer(raw_ri$Oxygen)), 0L),
  time = suppressWarnings(as.numeric(raw_ri$Time)),
  ri = suppressWarnings(as.numeric(raw_ri$RI)),
  type = tolower(trimws(raw_ri$Type))
) |>
  dplyr::filter(!is.na(time), !is.na(ri))

# ---- chemical family name fragments ----------------------------------------
# Saturates are picked up by CnH2n+2 so no need to list them here.
chemical_family_terms <- tibble::tribble(
  ~family      , ~term          ,
  "Aromatic"   , "benzene"      ,
  "Aromatic"   , "toluene"      ,
  "Aromatic"   , "xylene"       ,
  "Aromatic"   , "phenyl"       ,
  "Aromatic"   , "naphthalene"  ,
  "Aromatic"   , "indene"       ,
  "Aromatic"   , "styrene"      ,
  "Aromatic"   , "phenol"       ,
  "Aromatic"   , "cresol"       ,
  "Aromatic"   , "anisole"      ,
  "Aromatic"   , "thiophene"    ,
  "Aromatic"   , "pyridine"     ,
  "Naphthene"  , "cyclo"        ,
  "Naphthene"  , "cyclopropane" ,
  "Naphthene"  , "cyclobutane"  ,
  "Naphthene"  , "cyclopentane" ,
  "Naphthene"  , "cyclohexane"  ,
  "Naphthene"  , "cycloheptane" ,
  "Naphthene"  , "cyclooctane"  ,
  "Naphthene"  , "decalin"      ,
  "Naphthene"  , "bicyclo"      ,
  "Unsaturate" , "ene"          ,
  "Unsaturate" , "yne"          ,
  "Unsaturate" , "diene"        ,
  "Oxygenate"  , "ol"           ,
  "Oxygenate"  , "one"          ,
  "Oxygenate"  , "al"           ,
  "Oxygenate"  , "acid"         ,
  "Oxygenate"  , "ether"        ,
  "Oxygenate"  , "ester"        ,
  "Oxygenate"  , "furan"
)

solvent_names <- c(
  "methanol",
  "methyl alcohol",
  "dichloromethane",
  "methylene chloride",
  "acetone",
  "carbon disulfide",
  "water",
  "acetonitrile"
)

usethis::use_data(
  calibration_targets,
  calibration_points,
  rical_reference_points,
  chemical_family_terms,
  solvent_names,
  overwrite = TRUE
)
