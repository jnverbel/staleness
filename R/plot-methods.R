#' Plot detector calibration from a backtest
#'
#' Base graphics on purpose: a methods package meant to last should not depend
#' on a plotting stack it does not control.
#'
#' @param x A `staleness_backtest`.
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`.
#' @param ... Passed to [graphics::barplot()].
#' @return `x`, invisibly.
#' @examples
#' library(metafor)
#' bcg <- data.frame(
#'   yi   = c(-0.89, -1.59, -1.35, -1.44, -0.22, -0.79, -1.62,
#'             0.01, -0.47, -1.37, -0.34,  0.45, -0.02),
#'   vi   = c(0.326, 0.195, 0.415, 0.020, 0.051, 0.007, 0.223,
#'            0.004, 0.056, 0.073, 0.012, 0.533, 0.071),
#'   year = c(1948, 1949, 1960, 1977, 1973, 1953, 1973,
#'            1980, 1968, 1961, 1974, 1969, 1976)
#' )
#' bt <- backtest(evidence_stream(rma(yi, vi, data = bcg, measure = "RR"),
#'                                date = bcg$year))
#'
#' op <- graphics::par(no.readonly = TRUE)
#' plot(bt)
#' plot(bt, truth = "conclusion")
#' graphics::par(op)
#' @export
plot.staleness_backtest <- function(x, truth = "shift", ...) {
  cal <- calibration(x, truth = truth)
  # calibration() now returns a row for every requested method, including ones
  # with no eligible rows and therefore NA metrics (see ?calibration). A run
  # where *every* method is in that state has nothing to draw, and must say so
  # rather than hand an all-NA matrix to barplot() and produce an empty pair of
  # axes that looks like a result.
  if (is.null(cal) || nrow(cal) == 0 ||
      all(is.na(cal$sensitivity) & is.na(cal$specificity))) {
    stop("nothing to plot for truth = \"", truth, "\": every row was ",
         "not_applicable, censored, or had an undeterminable truth value",
         call. = FALSE)
  }
  vals <- rbind(sensitivity = cal$sensitivity, specificity = cal$specificity)
  colnames(vals) <- cal$method

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)

  # Defaults, but only where the caller has not spoken. Setting them here AND
  # forwarding `...` meant that supplying main, ylab or ylim -- the three most
  # obvious things to want to change -- failed on R's argument matching.
  dots <- list(...)
  args <- list(height = vals, beside = TRUE, legend.text = rownames(vals))
  if (is.null(dots$ylim)) args$ylim <- c(0, 1)
  if (is.null(dots$ylab)) args$ylab <- "rate"
  if (is.null(dots$main)) {
    args$main <- paste0("Detector calibration (truth: ", truth, ")")
  }
  do.call(graphics::barplot, c(args, dots))
  if (any(cal$contaminated)) {
    graphics::mtext(
      paste("contaminated (shares logic with this truth):",
            paste(cal$method[cal$contaminated], collapse = ", ")),
      side = 1, line = 3, cex = 0.8
    )
  }
  invisible(x)
}
