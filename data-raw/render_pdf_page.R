# Dev-only: renders the example diesel PDF report and saves page 1 for the vignette.
# Needs Quarto + LaTeX. Run from the package root after devtools::load_all().
v <- "vignettes"
info <- sample_info(
  "example diesel",
  "Demo fuel sample",
  "DEMO",
  "PB",
  density = 0.8519,
  methanol_wt_pct = 0,
  water_wt_pct = 53 / 1e4
)
prep <- prepare_peaks(
  file.path(v, "FIDpeaks.csv"),
  file.path(v, "MSDpeaks.csv"),
  file.path(v, "MSDsearch.csv")
)
res <- summarize_report(
  calibrate(prep$peaks, style = "slope", dilution = 1, density = info$density),
  info,
  shift = prep$shift
)
pdf <- write_pdf_report(res, file.path(tempdir(), "report_diesel.pdf"))
pdftools::pdf_convert(
  pdf,
  pages = 1,
  dpi = 110,
  verbose = FALSE,
  filename = file.path(v, "images", "07_pdf_report_page1.png")
)
