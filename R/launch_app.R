#' Estimate sample size for phase I trial via simulation interactively
#'
#' This function launch a shiny-based web app that allows the user calculate sample size for multiple arm multiple endpoint for phase I trial
#' The users are allowed to customize the parameters for simulation
#'
#' @examples
#' \dontrun{
#'    launch_app()
#'
#' }
#'
#' @export

launch_app <- function(){
  appDir <- system.file("app", package = "simsamplesize")
  shiny::runApp(appDir, display.mode = "normal")
}
