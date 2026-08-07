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
#' @param alpha Significance level used to decide applicability.
#' @param z_crit Critical value, 1.96 in the published method.
#' @return A `staleness_verdict`.
#' @export
barrowman <- function(prev, n_prev, n_new, alpha = 0.05, z_crit = 1.96) {
  if (is.null(n_prev) || is.null(n_new) ||
      !is.finite(n_prev) || !is.finite(n_new)) {
    stop("`barrowman()` needs the sample size of both the prior meta-analysis ",
         "and the new studies", call. = FALSE)
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
