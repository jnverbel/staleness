#' Barrowman participant-ratio detector
#'
#' Compares the participants contributed by new studies against the number that
#' would be needed for the pooled result to reach significance.
#'
#' Applies only when the prior meta-analysis was not statistically significant.
#' When it does not apply the verdict is `"not_applicable"`, never `"current"`:
#' treating inapplicability as evidence of currency would bias any calibration
#' study in this detector's favour.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param n_prev Numeric, total participants in the prior meta-analysis.
#' @param n_new Numeric, participants contributed by the new studies.
#' Omitting a sample size altogether is a malformed call and is an error.
#' Supplying one that is present but not a usable number — `NA`, `NaN`,
#' infinite — is a fact about the data, not about the call, so it yields
#' `"not_applicable"` with the reason attached, the same treatment every other
#' un-answerable case gets.
#'
#' @param alpha Significance level used to decide applicability.
#' @param z_crit Critical value, 1.96 in the published method.
#' @return A `staleness_verdict`.
#' @examples
#' library(metafor)
#' # This detector only speaks when the prior review was inconclusive, so the
#' # example needs a prior that genuinely failed to reach significance.
#' prev <- rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
#'             vi = c(0.16, 0.20, 0.18, 0.15, 0.22), measure = "RR")
#' prev$pval > 0.05
#'
#' # 2265 new participants against 555 before: enough to matter.
#' barrowman(prev, n_prev = 555, n_new = 2265)
#'
#' # On a prior that was already significant it declines to answer. That is
#' # deliberately not "current": scoring inapplicability as a correct call
#' # would flatter this detector in any calibration.
#' decided <- rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4), measure = "RR")
#' barrowman(decided, n_prev = 555, n_new = 2265)
#' @export
barrowman <- function(prev, n_prev, n_new, alpha = 0.05, z_crit = 1.96) {
  if (is.null(n_prev) || is.null(n_new)) {
    stop("`barrowman()` needs the sample size of both the prior meta-analysis ",
         "and the new studies", call. = FALSE)
  }
  if (!is.finite(n_prev) || !is.finite(n_new)) {
    return(verdict_na("barrowman",
      "sample size is missing or not a finite number; the participant ratio cannot be computed"))
  }
  if (!is.finite(prev$pval)) {
    return(verdict_na("barrowman",
      "the prior meta-analysis has no usable p-value; applicability cannot be decided"))
  }
  if (prev$pval < alpha) {
    return(verdict_na("barrowman",
      "prior meta-analysis was already significant; the method does not apply"))
  }
  z <- prev$zval
  if (!is.finite(z) || abs(z) < 1e-6) {
    return(verdict_na("barrowman",
      "prior test statistic is too close to zero; required sample size diverges"))
  }
  n_required <- (z_crit^2 * n_prev) / z^2
  q <- n_new / n_required
  new_verdict("barrowman",
              if (q > 1) "out_of_date" else "current",
              signal = q,
              detail = list(n_required = n_required, n_new = n_new, z_prev = z))
}
