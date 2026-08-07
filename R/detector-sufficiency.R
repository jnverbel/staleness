#' Rosenthal's fail-safe N
#'
#' Included because the sufficiency and stability method specifies it. Note that
#' Rosenthal's fail-safe N has been discredited as a measure of publication bias
#' since Becker (2005). It is implemented faithfully so that the backtesting
#' engine can settle the question with data rather than opinion.
#'
#' @param yi,vi Effect sizes and their variances.
#' @param z_crit One-sided critical value, 1.645 for alpha = 0.05.
#' @return The fail-safe N.
#' @keywords internal
failsafe_n <- function(yi, vi, z_crit = 1.645) {
  z <- yi / sqrt(vi)
  (sum(z)^2) / (z_crit^2) - length(z)
}

#' Sufficiency and stability detector
#'
#' Sufficiency is the fail-safe N scaled by `5k + 10`; a review is sufficient
#' when this index exceeds 1, Rosenthal's own rule of thumb for a pooled effect
#' being robust to unpublished null studies. Stability is the slope of the
#' cumulative pooled effect regressed on accumulated information; a review is
#' stable when that slope does not differ significantly from zero.
#'
#' Per the primary source (Mullen, Muellerleile & Bryant 2001, as applied by
#' Pattanittum et al. 2012), an out-of-date review is one that is BOTH
#' sufficient and unstable: enough evidence has accumulated to be confident the
#' effect is real, but the pooled estimate is still drifting, so its magnitude
#' is not yet settled. Insufficient evidence alone is never grounds for
#' "out of date" — see the design doc (section 4.4) for how this was confirmed
#' against the primary source, which is the opposite combination from what a
#' first reading of the secondary source suggests.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param min_k Minimum number of studies. Below this the stability regression
#'   is meaningless.
#' @param alpha_slope Significance level for the stability slope test.
#' @return A `staleness_verdict`.
#' @export
sufficiency <- function(prev, new_ma, min_k = 5, alpha_slope = 0.05) {
  yi <- as.numeric(new_ma$yi)
  vi <- as.numeric(new_ma$vi)
  k  <- length(yi)
  if (k < min_k) {
    return(verdict_na("sufficiency",
      paste0("needs at least ", min_k, " studies; found ", k)))
  }

  index <- failsafe_n(yi, vi) / (5 * k + 10)
  sufficient <- index > 1

  # cumulative fixed-effect estimate after each study, in input order
  w      <- 1 / vi
  cum_wy <- cumsum(w * yi)
  cum_w  <- cumsum(w)
  cum_theta <- cum_wy / cum_w

  idx <- seq_len(k)
  fit <- stats::lm(cum_theta[-1] ~ idx[-1])
  # A near-constant or near-perfectly-linear cumulative effect (identical
  # studies, or a very clean drift) triggers summary.lm()'s "essentially
  # perfect fit" warning. That is a floating-point artifact of the fit, not a
  # sign of a broken model, so it is suppressed rather than left to leak into
  # every caller's console.
  p_slope <- suppressWarnings(summary(fit)$coefficients[2, 4])
  stable  <- p_slope >= alpha_slope

  # Out-of-date requires sufficiency AND instability together. Insufficient
  # evidence (index <= 1) never triggers "out_of_date" by itself, regardless
  # of stability.
  out <- if (sufficient && !stable) "out_of_date" else "current"
  new_verdict("sufficiency", out, signal = index,
              detail = list(index = index, sufficient = sufficient,
                            stable = stable, slope = unname(stats::coef(fit)[2]),
                            p_slope = p_slope, k = k))
}
