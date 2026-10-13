#' Launch the local FID/MS processing app
#'
#' Client data never leaves the machine: the app is intended to be run locally
#' from Positron, not deployed to a server.
#'
#' @param ... Passed to [shiny::runApp()].
#' @return Invisibly, the result of [shiny::runApp()].
#' @export
run_app <- function(...) {
  rlang::check_installed(c("shiny", "bslib", "DT"))
  app_dir <- system.file("app", package = "fidms")
  if (!nzchar(app_dir)) {
    app_dir <- "inst/app"
  }
  shiny::runApp(app_dir, ...)
}
