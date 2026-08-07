#' Recursive cumulative meta-analysis detector
#'
#' Signals when the pooled effect, on its natural scale, changes by at least
#' half relative to the prior estimate.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param lower,upper Signal thresholds on the ratio scale.
#' @return A `staleness_verdict`.
#' @export
rcma <- function(prev, new_ma, lower = 0.5, upper = 1.5) {
  r <- effect_ratio(
    theta_new  = as.numeric(new_ma$beta),
    theta_prev = as.numeric(prev$beta),
    measure    = prev$measure,
    se_prev    = prev$se
  )
  if (is.na(r$ratio)) return(verdict_na("rcma", r$reason))

  out <- if (r$ratio <= lower || r$ratio >= upper) "out_of_date" else "current"
  new_verdict("rcma", out, signal = r$ratio,
              detail = list(lower = lower, upper = upper))
}
