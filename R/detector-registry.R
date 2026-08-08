#' Names of the available detectors
#'
#' @return A character vector of detector names.
#' @examples
#' available_methods()
#'
#' # The names accepted by check_currency() and backtest().
#' identical(available_methods(), c("rcma", "ottawa", "barrowman",
#'                                  "sufficiency", "simulation"))
#' @export
available_methods <- function() {
  c("rcma", "ottawa", "barrowman", "sufficiency", "simulation")
}
