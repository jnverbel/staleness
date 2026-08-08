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
#'   Either way the caller's own random stream is saved before the simulation
#'   and restored afterwards, so calling this detector — directly or through
#'   [check_currency()] or [backtest()] — never changes what a caller's
#'   subsequent `runif()`, `sample()` or `rnorm()` returns.
#' @return A `staleness_verdict`.
#' @examples
#' library(metafor)
#' # Like barrowman(), this one applies only to an inconclusive prior review.
#' prev <- rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
#'             vi = c(0.16, 0.20, 0.18, 0.15, 0.22), measure = "RR")
#' new <- list(yi = c(-0.45, -0.38, -0.52, -0.29, -0.41),
#'             vi = c(0.05, 0.04, 0.06, 0.05, 0.04), k = 5)
#'
#' # The signal is simulated power, so it needs a seed to be reproducible.
#' simulation(prev, new, B = 2000, seed = 1)
#'
#' # The caller's random stream is restored afterwards: running the detector
#' # does not move anyone else's simulation along.
#' set.seed(42); before <- runif(1)
#' set.seed(42); invisible(simulation(prev, new, B = 500, seed = 1))
#' identical(before, runif(1))
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

  k_new  <- new_evidence$k
  vi_new <- mean(new_evidence$vi)
  theta  <- as.numeric(prev$beta)
  tau2   <- prev$tau2
  sd_new <- sqrt(vi_new + tau2)

  yi_prev <- as.numeric(prev$yi)
  vi_prev <- as.numeric(prev$vi)

  # The whole simulation runs inside with_preserved_seed(): both the optional
  # set.seed() and the B * k_new draws it consumes are invisible to the
  # caller's stream, which is put back exactly as it was found.
  hits <- with_preserved_seed(seed = seed, {
    h <- 0L
    for (b in seq_len(B)) {
      yi_sim <- stats::rnorm(k_new, mean = theta, sd = sd_new)
      yi_all <- c(yi_prev, yi_sim)
      vi_all <- c(vi_prev, rep(vi_new, k_new))
      # Fixed-effect pooling: fast and adequate for a significance count.
      w      <- 1 / vi_all
      est    <- sum(w * yi_all) / sum(w)
      se     <- sqrt(1 / sum(w))
      p      <- 2 * stats::pnorm(-abs(est / se))
      if (p < alpha) h <- h + 1L
    }
    h
  })

  power <- hits / B
  new_verdict("simulation",
              if (power >= power_threshold) "out_of_date" else "current",
              signal = power,
              detail = list(B = B, k_new = k_new, vi_new = vi_new,
                            power_threshold = power_threshold, seed = seed))
}
