#' Plot detector calibration from a backtest
#'
#' Base graphics on purpose: a methods package meant to last should not depend
#' on a plotting stack it does not control.
#'
#' @param x A `staleness_backtest`.
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`.
#' @param ... Passed to [graphics::barplot()].
#' @return `x`, invisibly.
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

  graphics::barplot(vals, beside = TRUE, ylim = c(0, 1),
                    ylab = "rate", legend.text = rownames(vals),
                    main = paste0("Detector calibration (truth: ", truth, ")"),
                    ...)
  if (any(cal$contaminated)) {
    graphics::mtext(
      paste("contaminated (shares logic with this truth):",
            paste(cal$method[cal$contaminated], collapse = ", ")),
      side = 1, line = 3, cex = 0.8
    )
  }
  invisible(x)
}
