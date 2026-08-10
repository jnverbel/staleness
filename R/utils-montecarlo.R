# Monte Carlo error in the two detectors that simulate.
#
# simulation() estimates a power from B draws and sufficiency_changepoint()
# estimates a p-value from n_perm permutations. Both report a number to three
# decimals, and both would return a different number under a different seed.
# That variability is not a defect -- it is what a finite simulation costs --
# but it is invisible in the verdict, and it matters most exactly where the
# estimate sits beside the threshold that decides the verdict.
#
# Nothing here changes a verdict. The verdict stays the one the published rule
# gives on the point estimate; what is added is the answer to "would another
# seed have said something else?".

#' Wilson score interval for a simulated proportion
#'
#' @param successes,draws Counts from the simulation.
#' @param conf Coverage.
#' @return A length-two numeric vector, or two `NA`s when there is nothing to
#'   interval.
#' @keywords internal
#' @noRd
mc_interval <- function(successes, draws, conf = 0.95) {
  na <- c(NA_real_, NA_real_)
  if (!is.finite(successes) || !is.finite(draws) || draws < 1) return(na)
  # Wilson, not Wald. The Wald interval is `p +- z * sqrt(p(1-p)/n)`, whose
  # width is exactly zero at p = 0 and p = 1 -- so a power of 1.000 from 200
  # replicates would be reported as certain, which is the one place a reader
  # most needs to be told it is not. Wilson keeps a width there and stays
  # inside [0, 1] everywhere.
  z  <- stats::qnorm(1 - (1 - conf) / 2)
  ph <- successes / draws
  d  <- 1 + z^2 / draws
  centre <- ph + z^2 / (2 * draws)
  halfw  <- z * sqrt(ph * (1 - ph) / draws + z^2 / (4 * draws^2))
  c(max(0, (centre - halfw) / d), min(1, (centre + halfw) / d))
}

#' Monte Carlo standard error of a simulated proportion
#'
#' @inheritParams mc_interval
#' @return A single number, or `NA_real_`.
#' @keywords internal
#' @noRd
mc_se <- function(successes, draws) {
  if (!is.finite(successes) || !is.finite(draws) || draws < 1) return(NA_real_)
  ph <- successes / draws
  sqrt(ph * (1 - ph) / draws)
}

#' Whether the deciding threshold falls inside the Monte Carlo interval
#'
#' The practical reading of the interval: `TRUE` means a rerun with a different
#' seed could plausibly have returned the other verdict, and the estimate
#' should be treated as undecided rather than as the number it happens to be.
#' `FALSE` means the simulation is resolved on this question -- not that the
#' detector is right.
#'
#' @param interval As returned by [mc_interval()].
#' @param threshold The value the verdict compares against.
#' @return A single logical, `FALSE` when no interval was computed.
#' @keywords internal
#' @noRd
mc_near_threshold <- function(interval, threshold) {
  if (anyNA(interval) || !is.finite(threshold)) return(FALSE)
  threshold >= interval[1] && threshold <= interval[2]
}
