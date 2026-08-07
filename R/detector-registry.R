#' Names of the available detectors
#'
#' @return A character vector of detector names.
#' @export
available_methods <- function() {
  c("rcma", "ottawa", "barrowman", "sufficiency", "simulation")
}
