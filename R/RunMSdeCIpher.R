#' Start MSdeCIpher
#'
#'Opens MSdeCIpher shiny user interface in the system's default browser.
#'@export
RunMSdeCIpher<- function() {
  appDir <- system.file("ShinyApp", "MSdeCIpher", package = "MSdeCIpher")
  if (appDir == "") {
    stop("Could not find package. Try re-installing.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
