# Dev-only: captures the screenshots embedded in vignettes/real_sample_workflow.Rmd.
# Run from the package root; needs the local (gitignored) D319.D directory.
library(shinytest2)

pkg <- normalizePath(".", winslash = "/")
d_dir <- normalizePath("D319.D", winslash = "/")
out <- file.path(pkg, "vignettes", "images")
unlink(out, recursive = TRUE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# The app keeps no persistent state, so a throwaway launcher that loads the
# working-tree package is enough to keep the run isolated.
app_dir <- file.path(tempfile("fidms_app"))
dir.create(app_dir)
writeLines(
  c(
    sprintf("pkgload::load_all('%s', quiet = TRUE)", pkg),
    sprintf(
      "source('%s', local = TRUE)$value",
      file.path(pkg, "inst/app/app.R")
    )
  ),
  file.path(app_dir, "app.R")
)

app <- AppDriver$new(app_dir, width = 1400, height = 900, load_timeout = 60000)
on.exit(app$stop(), add = TRUE)
shot <- function(name) app$get_screenshot(file.path(out, name))

shot("01_app_start.png")

app$set_inputs(
  d_path = d_dir,
  sample_id = "D319",
  sample_name = "Demo fuel sample",
  project_code = "DEMO",
  analyst = "ab/cd",
  density = 0.85,
  water = 0.5,
  dilution = 1,
  wait_ = FALSE
)
app$run_js("document.querySelectorAll('.accordion-button')[0].click()")
Sys.sleep(1)
app$run_js("document.querySelectorAll('.accordion-button')[1].click()")
Sys.sleep(1)
shot("02_sample_inputs.png")

app$run_js("document.querySelectorAll('.accordion-button')[1].click()")
Sys.sleep(1)
app$run_js("document.querySelectorAll('.accordion-button')[2].click()")
Sys.sleep(1)
Sys.sleep(1)
shot("02b_processing_options.png")
app$click("run", timeout_ = 120000)
app$wait_for_idle(timeout = 120000)
Sys.sleep(3)
shot("03_peaks.png")

app$run_js("document.querySelector(\"a[data-value='Summary']\").click()")
app$wait_for_idle(timeout = 60000)
Sys.sleep(2)
shot("04_summary.png")

app$run_js(
  "document.querySelector(\"a[data-value='Retention index points']\").click()"
)
Sys.sleep(1)
shot("05_ri_points.png")

app$run_js(
  "document.querySelector(\"a[data-value='MS Class Features']\").click()"
)
Sys.sleep(1)
shot("06_class_features.png")
