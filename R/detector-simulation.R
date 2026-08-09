#' Prospective-power simulation detector
#'
#' Simulates new studies with the parameters of the recent literature, adds them
#' to the prior meta-analysis, and reports how often the pooled result would
#' reach significance. Signals when that power crosses the threshold.
#'
#' In the published comparison this method flagged none of 80 reviews. It is
#' either extremely conservative or misspecified; the backtesting engine exists
#' to find out which. The highest simulated power in that cohort was 63.4%,
#' against a threshold of 80%.
#'
#' @section How faithfully this follows the published procedure:
#' The method is set out step by step in Pattanittum et al. (2012), Appendix
#' S1, and this implementation follows it: the new study's effect is drawn from
#' a **t** distribution with the prior meta-analysis's pooled effect and
#' variance; **one** study is simulated, carrying the combined precision of the
#' recent studies rather than one study per recent study; the updated pooled
#' effect is tested at 5%; the draw is repeated `B` times; and power **strictly
#' above** the threshold is the signal, as the source's "Power >80%" reads. The
#' degrees of freedom are the one parameter the source leaves implicit; they
#' are taken as the prior study count minus one.
#'
#' One deviation is unavoidable and is declared rather than hidden. The source
#' simulates at the level of **participants** — binomial draws for dichotomous
#' outcomes, normal draws for continuous ones, from event rates and means
#' recovered per arm — and then computes an effect from them. This package
#' never sees participant-level data. It works from `yi` and `vi`, so it
#' simulates the effect directly, with variance `vi_new + tau^2`. For a study
#' of any reasonable size the two coincide closely, since the effect estimate
#' is asymptotically normal about the same mean with the same variance; they
#' can differ for very small or very sparse studies, where the participant-level
#' distribution is skewed and the effect-level one is not. An implementation
#' that wanted the source's exact procedure would need the 2x2 tables, which is
#' a different package from this one.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_evidence A list with `yi`, `vi` and `k`, as returned by
#'   [window_between()]. The three must agree: `yi` and `vi` of the same
#'   length, and `k` equal to that length. `k` decides whether there is any
#'   new evidence while `yi` and `vi` are what gets used, so a mismatch is
#'   answered from one and reported from the other; it is refused rather than
#'   resolved. The values must also be usable: `yi` finite, `vi` finite and
#'   strictly positive.
#' @param B Number of simulation replicates.
#' @param alpha Significance level for each simulated meta-analysis.
#' @param power_threshold Signal when simulated power *exceeds* this value.
#'   The comparison is strict, following the source's "Power >80%", so a
#'   simulated power of exactly `power_threshold` reads as `"current"`.
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
  check_rma_uni(prev, "prev")
  check_count(B, "B")
  check_probability(alpha, "alpha")
  # Closed here and open for alpha: a threshold of 0 means every power fires
  # and 1 means none does, and both are coherent requests.
  check_probability(power_threshold, "power_threshold", closed = TRUE)
  # Checked here rather than only at with_preserved_seed(), which sits past
  # the early returns below: a not_applicable verdict must not swallow a
  # malformed seed.
  check_seed(seed)
  # Same contract as check_currency(): this function is exported too, and its
  # `new_evidence` reaches it directly as often as through that one.
  check_new_evidence(new_evidence, "new_evidence")
  if (prev$pval < alpha) {
    return(verdict_na("simulation",
      "prior meta-analysis was already significant; the method does not apply"))
  }
  if (is.null(new_evidence$k) || new_evidence$k < 1) {
    return(verdict_na("simulation",
      "no recent studies to estimate simulation parameters from"))
  }

  k_new  <- new_evidence$k
  yi_prev <- as.numeric(prev$yi)
  vi_prev <- as.numeric(prev$vi)

  # The published procedure simulates ONE new study carrying the combined size
  # of the recent ones, not one study per recent study, so its variance is the
  # summed precision rather than the mean variance.
  vi_new <- 1 / sum(1 / new_evidence$vi)
  theta  <- as.numeric(prev$beta)
  tau2   <- prev$tau2
  sd_new <- sqrt(vi_new + tau2)
  # Drawn from a t-distribution with the prior meta-analysis's parameters, as
  # the source specifies; degrees of freedom from the prior study count, which
  # the source leaves implicit.
  df <- max(1L, length(yi_prev) - 1L)

  # Each replicate must be tested under the model the caller actually fitted.
  # It used to be pooled fixed-effect with a normal p-value regardless, which
  # made simulation() the one detector that ignored what evidence_stream() and
  # check_currency() go to some trouble to carry: a prior fitted with
  # test = "knha" got the same power as one fitted with test = "z", though
  # their p-values differ by a factor of two.
  #
  # It was also incoherent with itself. The replicate is drawn with
  # sd_new = sqrt(vi_new + tau2) -- under heterogeneity -- and was then pooled
  # without it. Simulated in one world, tested in another.
  prev_test     <- if (is.null(prev$test)) "z" else prev$test
  prev_weighted <- if (is.null(prev$weighted)) TRUE else isTRUE(prev$weighted)
  # tau^2 is treated as KNOWN, at the value the prior meta-analysis estimated.
  # That is not a shortcut, it is the only reading coherent with the draw: the
  # replicate is generated with sd = sqrt(vi_new + tau2), so tau^2 is already
  # a parameter of the simulated world. Re-estimating it from each replicate
  # would test in a different world from the one drawn.
  #
  # With tau^2 held, inverse-variance pooling and a z test have a closed form
  # that IS rma()'s answer -- verified equal to twelve digits against
  # rma(tau2 = ), and against rma(method = "FE") in the tau^2 = 0 case. So the
  # fast path is exact wherever it is used, and it covers metafor's own
  # default (REML with test = "z"). Only a non-z test needs a refit.
  closed_form_exact <- identical(prev_test, "z") && prev_weighted

  out <- with_preserved_seed(seed = seed, {
    h <- 0L; nonconv <- 0L
    for (b in seq_len(B)) {
      yi_sim <- theta + stats::rt(1L, df = df) * sd_new
      yi_all <- c(yi_prev, yi_sim)
      vi_all <- c(vi_prev, vi_new)
      if (closed_form_exact) {
        # tau^2 enters the weights, which is what the fixed-effect pooling
        # used to leave out: drawn under heterogeneity, pooled without it.
        w   <- 1 / (vi_all + tau2)
        est <- sum(w * yi_all) / sum(w)
        se  <- sqrt(1 / sum(w))
        p   <- 2 * stats::pnorm(-abs(est / se))
      } else {
        fit <- tryCatch(suppressWarnings(metafor::rma(
          yi       = yi_all,
          vi       = vi_all,
          measure  = prev$measure,
          test     = prev_test,
          weighted = prev_weighted,
          tau2     = tau2
        )), error = function(e) NULL)
        # REML does not always converge on a simulated draw -- measured at
        # roughly one replicate in ten. Counted rather than swept up: dividing
        # by B would score every failure as "not significant" and bias power
        # downwards, so the denominator is the replicates that could be
        # evaluated, and the count travels in the verdict.
        if (is.null(fit) || !is.finite(fit$pval)) { nonconv <- nonconv + 1L; next }
        p <- fit$pval
      }
      if (p < alpha) h <- h + 1L
    }
    list(hits = h, nonconv = nonconv)
  })

  evaluable <- B - out$nonconv
  # A power estimated from a small remnant is not an estimate. One replicate in
  # five failing means the model cannot be fitted to this kind of draw, which
  # is a fact about the analysis, not a number to report.
  if (evaluable < 1L || out$nonconv > B / 5) {
    return(verdict_na("simulation", paste0(
      "the prior model could not be refitted on ", out$nonconv, " of ", B,
      " simulated replicates; the power estimate would rest on too few")))
  }
  power <- out$hits / evaluable
  new_verdict("simulation",
              # Strictly above: the source reads "Power >80%".
              if (power > power_threshold) "out_of_date" else "current",
              signal = power,
              detail = list(B = B, k_new = k_new, k_simulated = 1L,
                            vi_new = vi_new, df = df,
                            power_threshold = power_threshold, seed = seed,
                            # What each replicate was tested under, and how
                            # many could not be. A caller comparing two runs
                            # needs to see that these differed.
                            model = prev$method, test = prev_test,
                            closed_form = closed_form_exact,
                            n_evaluable = evaluable,
                            n_nonconverged = out$nonconv))
}
