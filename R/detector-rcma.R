#' Recursive cumulative meta-analysis detector
#'
#' Signals when the pooled effect, on its natural scale, changes by at least
#' half relative to the prior estimate.
#'
#' @section Relationship to `ottawa()`:
#' These two rules look like the same rule and are not, which earlier versions
#' of this documentation asserted wrongly. Both compare an updated effect
#' against a prior one at thresholds of 0.5 and 1.5, but they compare
#' different quantities. `rcma()` takes the ratio of the effects. [ottawa()]
#' takes the ratio of the relative risk *reductions*,
#' `(1 - RR_new) / (1 - RR_prev)` -- the definition its published application
#' uses (Pattanittum et al. 2012, Table 1), verified in this package's test
#' suite against the ten reviews in that study's Appendix S3.
#'
#' On difference measures the two do coincide, because the source defers to
#' the rCMA rule there. On ratio measures they can disagree completely: of
#' those ten published reviews, every one fires on the RRR ratio and **not
#' one** fires on the ratio of the risk ratios. A review moving from RR 0.995
#' to RR 0.848 barely moves the risk ratio, while its risk reduction goes from
#' 0.5\% to 15\%.
#'
#' So there is no containment on ratio measures, and an `rcma` firing is not a
#' subset of the `ottawa` firings. This section stays because the relationship
#' is easy to assume and wrong in both directions: anyone recomputing the
#' inter-method agreement this package exists to re-ask should treat these as
#' two criteria, not as one criterion counted twice.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param lower,upper Signal thresholds on the ratio scale.
#' @return A `staleness_verdict`.
#' @examples
#' library(metafor)
#' # A published review of four trials, all agreeing on a halving of risk.
#' prev <- rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4), measure = "RR")
#'
#' # Four newer trials find no benefit, and the pooled effect moves 53%.
#' updated <- rma(yi = c(rep(log(0.50), 4), rep(log(1.10), 4)),
#'                vi = c(rep(0.05, 4), rep(0.02, 4)), measure = "RR")
#' rcma(prev, updated)
#'
#' # Containment in ottawa() is arithmetic, not an empirical agreement:
#' # whenever rcma fires, ottawa fires on the very same ratio.
#' c(rcma = rcma(prev, updated)$verdict, ottawa = ottawa(prev, updated)$verdict)
#'
#' # A smaller shift stays below the 1.5 threshold and reads as current.
#' mild <- rma(yi = c(rep(log(0.50), 4), rep(log(0.90), 3)),
#'             vi = c(rep(0.05, 4), rep(0.02, 3)), measure = "RR")
#' rcma(prev, mild)
#' @export
rcma <- function(prev, new_ma, lower = 0.5, upper = 1.5) {
  r <- effect_ratio(
    theta_new  = as.numeric(new_ma$beta),
    theta_prev = as.numeric(prev$beta),
    measure    = prev$measure,
    se_prev    = prev$se
  )
  if (is.na(r$ratio)) return(verdict_na("rcma", r$reason))
  # exp(theta_new - theta_prev) overflows once the prior effect underflows to
  # zero on the natural scale. Inf survives is.na(), and comparing it against
  # the thresholds yields out_of_date with a signal nobody can act on.
  if (!is.finite(r$ratio)) {
    return(verdict_na("rcma",
      "the effect ratio is not finite; the prior effect underflowed on its natural scale"))
  }

  out <- if (r$ratio <= lower || r$ratio >= upper) "out_of_date" else "current"
  new_verdict("rcma", out, signal = r$ratio,
              detail = list(lower = lower, upper = upper))
}
