#' Prospective-power simulation detector
#'
#' Simulates new studies with the parameters of the recent literature, adds them
#' to the prior meta-analysis, and reports how often the pooled result would
#' reach significance. Signals when that power crosses the threshold.
#'
#' In the published comparison this method flagged none of 80 reviews. It is
#' either extremely conservative or misspecified; the backtesting engine exists
#' to find out which.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_evidence A list with `yi`, `vi` and `k`, as returned by
#'   [window_between()].
#' @param B Number of simulation replicates.
#' @param alpha Significance level for each simulated meta-analysis.
#' @param power_threshold Signal when simulated power reaches this value.
#' @param seed Optional integer seed. Results are not reproducible without it.
#' @return A `staleness_verdict`.
#' @export
simulation <- function(prev, new_evidence, B = 10000, alpha = 0.05,
                       power_threshold = 0.80, seed = NULL) {
  if (prev$pval < alpha) {
    return(verdict_na("simulation",
      "prior meta-analysis was already significant; the method does not apply"))
  }
  if (is.null(new_evidence$k) || new_evidence$k < 1) {
    return(verdict_na("simulation",
      "no recent studies to estimate simulation parameters from"))
  }
  if (!is.null(seed)) set.seed(seed)

  k_new  <- new_evidence$k
  vi_new <- mean(new_evidence$vi)
  theta  <- as.numeric(prev$beta)
  tau2   <- prev$tau2
  sd_new <- sqrt(vi_new + tau2)

  yi_prev <- as.numeric(prev$yi)
  vi_prev <- as.numeric(prev$vi)

  hits <- 0L
  for (b in seq_len(B)) {
    yi_sim <- stats::rnorm(k_new, mean = theta, sd = sd_new)
    yi_all <- c(yi_prev, yi_sim)
    vi_all <- c(vi_prev, rep(vi_new, k_new))
    # Fixed-effect pooling: fast and adequate for a significance count.
    w      <- 1 / vi_all
    est    <- sum(w * yi_all) / sum(w)
    se     <- sqrt(1 / sum(w))
    p      <- 2 * stats::pnorm(-abs(est / se))
    if (p < alpha) hits <- hits + 1L
  }

  power <- hits / B
  new_verdict("simulation",
              if (power >= power_threshold) "out_of_date" else "current",
              signal = power,
              detail = list(B = B, k_new = k_new, vi_new = vi_new,
                            power_threshold = power_threshold, seed = seed))
}
