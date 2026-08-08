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

#' Cumulative fixed-effect estimate after each study
#'
#' @param yi,vi Effect sizes and their variances, in the order studies are to
#'   be accumulated.
#' @return Numeric vector of length `length(yi)`.
#' @keywords internal
cumulative_effect <- function(yi, vi) {
  w <- 1 / vi
  cumsum(w * yi) / cumsum(w)
}

#' Slope of a cumulative series against study index
#'
#' The first element of a cumulative series is just the first study, carrying
#' no information about accumulation, so it is dropped before the slope is
#' taken. Computed in closed form rather than through [stats::lm()]: the same
#' number, but cheap enough to repeat a thousand times inside a permutation
#' loop.
#'
#' @param cum_theta Numeric vector, as returned by [cumulative_effect()].
#' @return The ordinary-least-squares slope, or `NA_real_` if the series is
#'   too short to have one.
#' @keywords internal
cum_drift_slope <- function(cum_theta) {
  y <- cum_theta[-1]
  if (length(y) < 2) return(NA_real_)
  x  <- seq_along(y)
  xc <- x - mean(x)
  sum(xc * (y - mean(y))) / sum(xc^2)
}

#' Sufficiency and stability detector
#'
#' Sufficiency is the fail-safe N scaled by `5k + 10`; a review is sufficient
#' when this index exceeds 1, Rosenthal's own rule of thumb for a pooled effect
#' being robust to unpublished null studies. Stability is the slope of the
#' cumulative pooled effect regressed on study index; a review is stable when
#' that slope is no larger than the same studies would produce in a random
#' order.
#'
#' The two indicators deliberately read from different bodies of evidence.
#' Per the primary source that operationalises this method as a two-snapshot
#' comparison (Pattanittum et al. 2012, Table 1), sufficiency — both the
#' fail-safe N and the `k` in `5k + 10` — is computed on the meta-analysis
#' **as previously published** (`prev`), while stability is computed on the
#' cumulative series of the **updated** meta-analysis (`new_ma`), which is
#' where new studies show up as drift or lack of it. Confirmed against a
#' second, independent secondary source; the original method paper (Mullen,
#' Muellerleile & Bryant 2001) was not reachable in full text, and does not
#' itself define a two-snapshot variant to compare against — see the design
#' doc (section 4.4) for the full trail of evidence.
#'
#' Per the same source, an out-of-date review is one that is BOTH sufficient
#' and unstable: enough evidence had already accumulated, as of the prior
#' review, to be confident the effect is real, but the pooled estimate
#' (including what came after) is still drifting, so its magnitude is not yet
#' settled. Insufficient evidence alone is never grounds for "out of date" —
#' the opposite combination from what a first reading of the secondary source
#' suggests, and also confirmed in the design doc.
#'
#' @section How stability is tested, and why not as published:
#' The source states the instability criterion as an *"absolute slope of the
#' linear regression >0"*. Taken literally that rule is degenerate: on
#' continuous data the slope of a cumulative series is never exactly zero, so
#' every review with `index > 1` would be flagged out of date and the
#' indicator would carry no information. Some significance rule has to stand
#' in for it.
#'
#' The obvious substitute — the t-test on the slope reported by
#' [stats::lm()] — is invalid here, and measurably so. A cumulative mean is
#' near-perfectly autocorrelated by construction and converges on the pooled
#' effect by the law of large numbers, so an OLS test on it detects
#' *convergence* and reports *instability*: over 300 samples of genuinely
#' unchanging evidence (20 prior plus 10 new studies from one distribution) it
#' returned `out_of_date` 209 times, where `rcma()` and `ottawa()` returned it
#' none.
#'
#' This package therefore substitutes a **permutation test over study order**.
#' The observed slope is compared against the slopes obtained by reshuffling
#' the same studies `n_perm` times and recomputing the cumulative series each
#' time; `p_slope` is the two-sided proportion of permuted slopes at least as
#' large in absolute value, using the `(1 + count) / (n_perm + 1)` estimator
#' so that the p-value is never zero. This assumes nothing about the
#' independence of a series that is autocorrelated by design, and asks the
#' scientifically meaningful question directly: is the drift in the cumulative
#' effect larger than what *these same studies* would produce in a random
#' order? Under evidence that is not changing, study order is exchangeable and
#' the test is exact, so the no-change invariant holds by construction rather
#' than by luck. The same 300-sample experiment run against this
#' implementation returns `out_of_date` 12 times — the nominal 5% rate, within
#' sampling error, instead of 70%. See `vignette("methods")` for the
#' declaration of this adaptation.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param min_k Minimum number of studies in `new_ma`. Below this the
#'   stability regression is meaningless.
#' @param alpha_slope Cutoff for the stability permutation p-value: the review
#'   is unstable when `p_slope < alpha_slope`.
#' @param n_perm Number of order permutations used to build the null
#'   distribution of the slope.
#' @param seed Integer seed for the permutation draw. Fixed by default so
#'   verdicts are reproducible; the caller's own random stream is saved and
#'   restored around the draw, so calling this detector never changes it.
#' @return A `staleness_verdict`.
#' @export
sufficiency <- function(prev, new_ma, min_k = 5, alpha_slope = 0.05,
                        n_perm = 999, seed = 20260807) {
  yi <- as.numeric(new_ma$yi)
  vi <- as.numeric(new_ma$vi)
  k_new <- length(yi)
  if (k_new < min_k) {
    return(verdict_na("sufficiency",
      paste0("needs at least ", min_k, " studies; found ", k_new)))
  }

  # Sufficiency is read off the meta-analysis AS PREVIOUSLY PUBLISHED, not the
  # updated one (see the roxygen note above and design doc section 4.4).
  yi_prev <- as.numeric(prev$yi)
  vi_prev <- as.numeric(prev$vi)
  k_prev  <- length(yi_prev)
  index <- failsafe_n(yi_prev, vi_prev) / (5 * k_prev + 10)
  sufficient <- index > 1

  # Stability, by contrast, is read off the updated evidence: the cumulative
  # fixed-effect estimate after each study, in input order.
  cum_theta <- cumulative_effect(yi, vi)

  if (!all(is.finite(cum_theta))) {
    return(verdict_na("sufficiency",
      "the cumulative effect is not finite (a study variance of zero?); stability cannot be assessed"))
  }

  # Degenerate short circuit, BEFORE any fitting. Byte-identical studies give a
  # cumulative series whose entire spread is rounding noise -- 12 identical
  # studies span 2.2e-16 -- and fitting a model to that returns a slope of
  # 1e-17 with a "significant" p-value. Such a series is not drifting; it is
  # constant, which is maximal stability. Never fit a model to rounding noise.
  fitted_series <- cum_theta[-1]
  scale <- max(abs(fitted_series))
  tol <- 4 * k_new * .Machine$double.eps * max(scale, .Machine$double.eps)
  if (diff(range(fitted_series)) <= tol) {
    slope   <- 0
    p_slope <- 1
    stable  <- TRUE
  } else {
    slope <- cum_drift_slope(cum_theta)
    # Permutation test over study order (see the roxygen section above). The
    # seed is applied through with_preserved_seed() so that the caller's own
    # random stream is left exactly as it was found.
    perm <- with_preserved_seed(seed = seed, {
      vapply(seq_len(n_perm), function(i) {
        o <- sample.int(k_new)
        cum_drift_slope(cumulative_effect(yi[o], vi[o]))
      }, numeric(1))
    })
    p_slope <- (1 + sum(abs(perm) >= abs(slope), na.rm = TRUE)) / (n_perm + 1)
    # Defensive: a p-value that is not finite must be resolved deliberately,
    # never left to turn `if (sufficient && !stable)` into `if (NA)`. An
    # undeterminable drift is not evidence of drift, so it reads as stable.
    stable <- if (is.finite(p_slope)) p_slope >= alpha_slope else TRUE
  }

  # Out-of-date requires sufficiency AND instability together. Insufficient
  # evidence (index <= 1) never triggers "out_of_date" by itself, regardless
  # of stability.
  out <- if (sufficient && !stable) "out_of_date" else "current"
  new_verdict("sufficiency", out, signal = index,
              detail = list(index = index, sufficient = sufficient,
                            stable = stable, slope = slope,
                            p_slope = p_slope, k = k_prev, k_new = k_new))
}
