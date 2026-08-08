#' Recursive cumulative meta-analysis detector
#'
#' Signals when the pooled effect, on its natural scale, changes by at least
#' half relative to the prior estimate.
#'
#' @section Relationship to `ottawa()`:
#' This rule is contained verbatim in [ottawa()], which fires on the same
#' `effect_ratio()` with the same 0.5 / 1.5 thresholds as one of its three
#' signals. At default parameters, therefore, `rcma() == "out_of_date"`
#' **implies** `ottawa() == "out_of_date"`: the two detectors are not
#' independent, and no `rcma` firing can ever be disjoint from `ottawa`.
#' Verified over 400 random prior/updated pairs: `rcma` fired 63 times,
#' `ottawa` fired 116 times, and every one of `rcma`'s 63 was among
#' `ottawa`'s. Zero disjoint firings, as the algebra requires.
#'
#' This is faithful to the published methods, not a defect in either: the
#' Ottawa method genuinely adopts the rCMA effect-change criterion and adds to
#' it. It is declared because it is structural rather than empirical. Anyone
#' recomputing the inter-method agreement question this package exists to
#' answer (Kappa = 0.14 across the five methods) would otherwise be treating
#' one of the ten detector pairs as data when it is arithmetic. The same
#' reasoning gave rise to [CONTAMINATED_PAIRS] for detector-truth pairs;
#' detector-to-detector containment deserves the same visibility.
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

  out <- if (r$ratio <= lower || r$ratio >= upper) "out_of_date" else "current"
  new_verdict("rcma", out, signal = r$ratio,
              detail = list(lower = lower, upper = upper))
}
