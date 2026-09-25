## Builds packaged data for the hydrocarbon-type (SIM/class-ion) pipeline:
## `sim_masses`, `sim_masses_oxygenate`, `class_response_matrix`,
## `oxygenate_response_matrix` and `parent_ion_factors`. See
## design/hydrocarbon-type-design.md.
##
## Hydrocarbon-side numbers are transcribed from ASTM D2425-23 (middle
## distillates) and D2789-95 (gasoline range, ion-sum definitions only).
## `class_response_matrix$D2425` is built from the Appendix X1 worked example
## (the cleanest, self-validating source per the design doc): Table X1.1
## ("Aromatics Matrix") for the aromatic fraction and Table 5 ("Saturate
## Concentration Matrix") for the saturate fraction. Both are single
## representative-carbon-number matrices (not a full per-carbon-number
## interpolation grid); `.solve_hydrocarbon_class()` errors explicitly rather
## than extrapolating silently outside of what is transcribed here.
## `parent_ion_factors` (Table 1) supports the separate parent-ion average
## carbon-number calculation (Sec 15.2-15.5) used by
## `.average_carbon_number()`/`.estimate_carbon_number()`.
library(tibble)

# ---- sim_masses: diagnostic ion-sum definitions -----------------------------
# `monitored_mz` is a flat numeric vector of individual masses to sum; per
# design discussion, ordering/pairing within a class carries no meaning here
# (carbon-number attribution is a Phase 3 concern), so the D2425 "doublet"
# series are flattened rather than kept as carbon-number-keyed pairs.
sim_masses <- tibble::tribble(
  ~standard , ~class               , ~ion_sum_label , ~monitored_mz                                                                 ,
  # --- D2789-95 (gasoline range, C6-C10), Eq 1-6 ---
  "D2789"   , "Paraffin"           , "Sum43"        , list(c(43, 57, 71, 85, 99))                                                   ,
  "D2789"   , "Monocycloparaffin"  , "Sum41"        , list(c(41, 55, 69, 83, 97))                                                   ,
  "D2789"   , "Dicycloparaffin"    , "Sum67"        , list(c(67, 68, 81, 82, 95, 96))                                               ,
  "D2789"   , "Alkylbenzene"       , "Sum77"        , list(c(77, 78, 79, 91, 92, 105, 106, 119, 120, 133, 134, 147, 148, 161, 162)) ,
  "D2789"   , "Indan_Tetralin"     , "Sum103"       , list(c(103, 104, 117, 118, 131, 132, 145, 146, 159, 160))                     ,
  "D2789"   , "Naphthalene"        , "Sum128"       , list(c(128, 141, 142, 155, 156))                                              ,
  # --- D2425-23 (middle distillate, C10-C18), Sec 4.2 / Table 3 ---
  "D2425"   , "Paraffin"           , "Sum71"        , list(c(71, 85))                                                               ,
  "D2425"   , "Cycloparaffin"      , "Sum67"        , list(c(67, 68, 69, 81, 82, 83, 96, 97))                                       ,
  "D2425"   , "Dicycloparaffin"    , "Sum123"       , list(seq(123, 249, by = 14), seq(124, 250, by = 14))                          ,
  "D2425"   , "Tricycloparaffin"   , "Sum149"       , list(seq(149, 247, by = 14), seq(150, 248, by = 14))                          ,
  "D2425"   , "Alkylbenzene"       , "Sum91"        , list(seq(91, 175, by = 14), seq(92, 176, by = 14))                            ,
  "D2425"   , "Indan_Tetralin"     , "Sum103"       , list(seq(103, 187, by = 14), seq(104, 188, by = 14))                          ,
  "D2425"   , "Indene"             , "Sum115"       , list(seq(115, 185, by = 14), seq(116, 186, by = 14))                          ,
  "D2425"   , "Naphthalene_parent" , "Mz128"        , list(128)                                                                     ,
  "D2425"   , "Naphthalene_alkyl"  , "Sum141"       , list(seq(141, 239, by = 14), seq(142, 240, by = 14))                          ,
  "D2425"   , "Acenaphthene"       , "Sum153"       , list(seq(153, 251, by = 14), seq(154, 252, by = 14))                          ,
  "D2425"   , "Acenaphthylene"     , "Sum151"       , list(seq(151, 249, by = 14), seq(152, 250, by = 14))                          ,
  "D2425"   , "Tricyclic_aromatic" , "Sum177"       , list(seq(177, 247, by = 14), seq(178, 248, by = 14))
)
# tribble() double-wraps a `list(...)` cell value; unwrap the extra level so
# `monitored_mz` is a flat list-column of numeric vectors.
sim_masses$monitored_mz <- lapply(sim_masses$monitored_mz, unlist)

# ---- class_response_matrix: D2425 Appendix X1 worked example ---------------
# Rows = diagnostic ion sums (matching `sim_masses$ion_sum_label` for D2425),
# columns = hydrocarbon classes; cell = pattern-coefficient contribution of
# that class to that ion sum. Transcribed from ASTM D2425-23 Appendix X1
# (Table X1.1, aromatic fraction) and Table 5 (saturate fraction; the X1.5
# worked-example matrix itself was not cleanly OCR-legible, so the equivalent
# clean Table 5 matrix is used, per the extraction subagent's note).
ion_sums_aromatic <- c(
  "Sum71",
  "Sum67",
  "Sum91",
  "Sum103",
  "Sum115",
  "Mz128",
  "Sum141",
  "Sum153",
  "Sum151",
  "Sum177"
)
classes_aromatic <- c(
  "Paraffin",
  "Cycloparaffin",
  "Alkylbenzene",
  "Indan_Tetralin",
  "Indene",
  "Naphthalene_parent",
  "Naphthalene_alkyl",
  "Acenaphthene",
  "Acenaphthylene",
  "Tricyclic_aromatic"
)
# Transcribed directly from ASTM D2425-23 Appendix X1, Table X1.1 "Aromatics
# Matrix" (reference/ASTM D2425 table X1-1.png). Row order matches the table's
# printed row order (Sum71, Sum67, Sum91, Sum103, Sum115, Mz128, Sum141,
# Sum153, Sum151, Sum177); column order matches `classes_aromatic` (Paraffin,
# Cycloparaffin, Alkylbenzene, Indan_Tetralin, Indene, Naphthalene_parent,
# Naphthalene_alkyl, Acenaphthene, Acenaphthylene, Tricyclic_aromatic), which
# is the same left-to-right order as the table's ten "13/14.5" columns.
aromatic_matrix <- matrix(
  c(
    # Sum71
    100,
    6,
    0.4,
    1,
    1.7,
    0.5,
    2,
    1,
    1,
    0.6,
    # Sum67
    23,
    100,
    2,
    2,
    6,
    0.8,
    2,
    2,
    5,
    0.7,
    # Sum91
    0.4,
    0.2,
    100,
    15,
    6.2,
    0.1,
    1,
    5,
    3,
    18,
    # Sum103
    0,
    0,
    10,
    100,
    20.3,
    0.6,
    0.1,
    3,
    3,
    1.5,
    # Sum115
    0,
    1,
    5,
    25,
    100,
    11.4,
    18,
    0.8,
    2.7,
    1,
    # Mz128
    0,
    0,
    1,
    0,
    13,
    100,
    5.6,
    0.7,
    0.1,
    0.8,
    # Sum141
    10,
    2,
    0,
    0,
    28,
    0,
    100,
    10,
    0,
    0.3,
    # Sum153
    0,
    2,
    0,
    0,
    6.1,
    0,
    10,
    100,
    15,
    3.5,
    # Sum151
    0,
    7,
    0,
    0,
    4.5,
    0,
    7,
    20,
    100,
    30,
    # Sum177
    0,
    2,
    0,
    0,
    0.6,
    0,
    0,
    4,
    15,
    100
  ),
  nrow = 10,
  byrow = TRUE,
  dimnames = list(ion_sums_aromatic, classes_aromatic)
)
# Mass-sensitivity per class, same column order (ASTM D2425-23 Table 3/X1.3).
aromatic_sensitivity <- c(
  Paraffin = 97,
  Cycloparaffin = 204,
  Alkylbenzene = 256,
  Indan_Tetralin = 241,
  Indene = 200,
  Naphthalene_parent = 184,
  Naphthalene_alkyl = 224,
  Acenaphthene = 196,
  Acenaphthylene = 205,
  Tricyclic_aromatic = 205
)

ion_sums_saturate <- c("Sum71", "Sum67", "Sum123", "Sum149", "Sum91")
classes_saturate <- c(
  "Paraffin",
  "Monocycloparaffin",
  "Dicycloparaffin",
  "Tricycloparaffin",
  "Alkylbenzene"
)
# Transcribed directly from ASTM D2425-23 Appendix X1, Table X1.5 "Saturates
# Matrix" (reference/ASTM D2425 table X1-5.png) -- the worked-example matrix,
# not the general Table 5 "Saturate Concentration Matrix" (whose cell values
# differ slightly, e.g. Sum71 = 100,6,1.1,1,0.4 there vs 100,6,1.5,2,0.5 here).
# Using the Appendix X1 matrix keeps this self-validating against the
# published X1.6/X1.7/X1.8 solved results. Row order: Sum71, Sum67, Sum123,
# Sum149, Sum91; column order matches `classes_saturate` (Paraffin,
# Monocycloparaffin, Dicycloparaffin, Tricycloparaffin, Alkylbenzene).
saturate_matrix <- matrix(
  c(
    # Sum71
    100,
    6,
    1.1,
    1,
    0.4,
    # Sum67
    23,
    100,
    130,
    170,
    2,
    # Sum123
    0.1,
    1,
    100,
    10,
    0.2,
    # Sum149
    0,
    0,
    5,
    100,
    1.5,
    # Sum91
    0.4,
    0.2,
    4,
    15,
    100
  ),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(ion_sums_saturate, classes_saturate)
)
saturate_sensitivity <- c(
  Paraffin = 97,
  Monocycloparaffin = 204,
  Dicycloparaffin = 156,
  Tricycloparaffin = 158,
  Alkylbenzene = 256
)

class_response_matrix <- list(
  D2425 = list(
    aromatic = list(
      matrix = aromatic_matrix,
      sensitivity = aromatic_sensitivity,
      carbon_number_range = c(11, 14.5)
    ),
    saturate = list(
      matrix = saturate_matrix,
      sensitivity = saturate_sensitivity,
      carbon_number_range = c(13, 15.5)
    )
  )
  # D2789 intentionally excluded: lab prefers middle-distillate (D2425)
  # matrices at the carbon-number overlap; D2789's own gasoline range
  # (C6-C10) sensitivity/inverse matrices are out of scope for this project.
)

# ---- parent_ion_factors: ASTM D2425-23 Table 1 ------------------------------
# Per-carbon-number parent-ion isotope correction factor (K1/L1) and mole
# sensitivity (K2/L2) for the two series the standard actually deconvolutes
# to an *average carbon number* by mass-spectral calculation (Sec 15.2-15.5,
# Appendix X1.2.2-X1.2.5) -- alkylbenzenes (parent ion, m/e = 14n + 4) and
# naphthalenes (parent ion, m/e = 14n - 12). This is a genuinely different
# calculation from the class-composition NNLS/exact solve above: it walks
# each series' *parent* (M) and *parent-1* (M-1) polyisotopic peak heights,
# not the low-mass fragment sums in `sim_masses`, and produces one weighted
# average carbon number per fraction (not a full mole-fraction-by-class
# breakdown). Transcribed directly from ASTM D2425-23 Table 1
# (reference/ASTM D2425 table 1.png).
parent_ion_factors <- tibble::tribble(
  ~series        , ~carbon , ~mz , ~isotope_factor , ~mole_sensitivity ,
  "Alkylbenzene" ,      10 , 134 , 0.1101          ,                85 ,
  "Alkylbenzene" ,      11 , 148 , 0.1212          ,                63 ,
  "Alkylbenzene" ,      12 , 162 , 0.1323          ,                60 ,
  "Alkylbenzene" ,      13 , 176 , 0.1434          ,                57 ,
  "Alkylbenzene" ,      14 , 190 , 0.1545          ,                54 ,
  "Alkylbenzene" ,      15 , 204 , 0.1656          ,                51 ,
  "Alkylbenzene" ,      16 , 218 , 0.1767          ,                48 ,
  "Alkylbenzene" ,      17 , 232 , 0.1878          ,                45 ,
  "Alkylbenzene" ,      18 , 246 , 0.1989          ,                42 ,
  "Naphthalene"  ,      11 , 142 , 0.1201          ,               194 ,
  "Naphthalene"  ,      12 , 156 , 0.1314          ,               166 ,
  "Naphthalene"  ,      13 , 170 , 0.1425          ,               150 ,
  "Naphthalene"  ,      14 , 184 , 0.1536          ,               150 ,
  "Naphthalene"  ,      15 , 198 , 0.1647          ,               150 ,
  "Naphthalene"  ,      16 , 212 , 0.1758          ,               150 ,
  "Naphthalene"  ,      17 , 226 , 0.1871          ,               150 ,
  "Naphthalene"  ,      18 , 240 , 0.1982          ,               150
)

# ---- sim_masses_oxygenate: proposed extension, unverified ------------------
# Not transcribed from any ASTM standard (neither covers oxygenates); this is
# the user-proposed extension from the design doc, tracked separately so it
# is never confused with the D2425/D2789-derived hydrocarbon tables.
sim_masses_oxygenate <- tibble::tribble(
  ~family                 , ~ion_sum_label , ~monitored_mz                        , ~hydrocarbon_interference      , ~compound                     , ~mz      ,
  "Alkylphenol"           , "Sum94"        , list(c(94, 108, 122, 136, 150, 164)) , "Alkylbenzene"                 , NA_character_                 , NA_real_ ,
  "Alkoxyphenol_guaiacol" , "Sum109"       , list(c(109, 123, 137, 151, 165))     , "Alkylbenzene"                 , NA_character_                 , NA_real_ ,
  "Saturated_FAME"        , "Sum74"        , list(c(74, 87, 143))                 , "Branched alkane"              , NA_character_                 , NA_real_ ,
  "Carboxylic_acid"       , "Sum60"        , list(c(60, 74, 88))                  , "Alkane/alkene"                , NA_character_                 , NA_real_ ,
  "Polyol_aromatic_diol"  , "Sum110"       , list(c(110, 124, 138, 152))          , "Indan/Tetralin"               , NA_character_                 , NA_real_ ,
  # Unsaturated FAMEs are confirmed molecular ions ([M+.]) of named compounds,
  # not a homologous fragment series -- one row per compound, monitored
  # individually rather than summed (see design doc).
  "Unsaturated_FAME"      , NA_character_  , list(268)                            , "Diaromatic/naphthenoaromatic" , "Methyl Palmitoleate (C16:1)" ,      268 ,
  "Unsaturated_FAME"      , NA_character_  , list(296)                            , "Diaromatic/naphthenoaromatic" , "Methyl Oleate (C18:1)"       ,      296 ,
  "Unsaturated_FAME"      , NA_character_  , list(294)                            , "Diaromatic/naphthenoaromatic" , "Methyl Linoleate (C18:2)"    ,      294 ,
  "Unsaturated_FAME"      , NA_character_  , list(292)                            , "Diaromatic/naphthenoaromatic" , "Methyl Linolenate (C18:3)"   ,      292
)
sim_masses_oxygenate$monitored_mz <- lapply(
  sim_masses_oxygenate$monitored_mz,
  unlist
)

# ---- oxygenate_response_matrix: schema only, unpopulated --------------------
# No standard provides these coefficients; must come from the lab's own
# calibration standards (see design doc). Shipped as an explicit NA-filled
# placeholder so the solve/interference-correction code path exists and is
# structurally testable ahead of real calibration data.
oxygenate_classes <- c(
  "Alkylphenol",
  "Alkoxyphenol_guaiacol",
  "Saturated_FAME",
  "Carboxylic_acid",
  "Polyol_aromatic_diol"
)
oxygenate_ion_sums <- c("Sum94", "Sum109", "Sum74", "Sum60", "Sum110")
oxygenate_response_matrix <- list(
  lab = list(
    matrix = matrix(
      NA_real_,
      nrow = length(oxygenate_ion_sums),
      ncol = length(oxygenate_classes),
      dimnames = list(oxygenate_ion_sums, oxygenate_classes)
    ),
    sensitivity = stats::setNames(
      rep(NA_real_, length(oxygenate_classes)),
      oxygenate_classes
    ),
    # cross-terms: fraction of a hydrocarbon class's resolved abundance
    # expected to appear in a given oxygenate ion-sum channel, subtracted
    # before the oxygenate NNLS solve (see .solve_oxygenate_class()).
    # keyed by ion_sum_label (row name), matching `.solve_oxygenate_class()`'s
    # lookup against the ion-sum vector it is correcting.
    hydrocarbon_cross_terms = list(
      Sum94 = c(Alkylbenzene = NA_real_),
      Sum109 = c(Alkylbenzene = NA_real_),
      Sum110 = c(Indan_Tetralin = NA_real_)
    ),
    # isotopic-tail correction: which dominant ion-sum channel's natural
    # 13C tail contaminates each oxygenate channel at the same nominal mass
    # (see .correct_isotopic_tail()). Alkylbenzene's Sum91 series is the
    # dominant contributor into both Sum94 (Alkylphenol) and Sum109
    # (Alkoxyphenol/Guaiacol).
    isotope_dominant_series = list(
      Sum94 = "Sum91",
      Sum109 = "Sum91"
    )
  )
)

usethis::use_data(
  sim_masses,
  sim_masses_oxygenate,
  class_response_matrix,
  oxygenate_response_matrix,
  parent_ion_factors,
  overwrite = TRUE
)
