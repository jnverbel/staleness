#' Ottawa method detector
#'
#' Two quantitative signals: a change in statistical significance, and a change
#' of at least 50 percent in the pooled effect.
#'
#' The published method also lists four qualitative signals: opposing findings,
#' substantial harm, a superior new therapy, and a change in effectiveness. This
#' function does not infer them. They are supplied by the analyst through
#' `qualitative` and reported separately. An algorithm that claimed to judge
#' "substantial harm" on its own would be overselling what the data support.
#'
#' @section Relationship to `rcma()`:
#' The effect-change half of this method is [rcma()]'s entire rule: the same
#' `effect_ratio()`, the same 0.5 / 1.5 thresholds. At default parameters
#' `rcma() == "out_of_date"` therefore **implies** `ottawa() == "out_of_date"`
#' by construction, so the two detectors cannot disagree in that direction and
#' their agreement is not an empirical quantity. Faithful to the published
#' methods; declared so that nobody reads a structural identity as a finding.
#' See [rcma()] for the measured demonstration.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param alpha Significance threshold. The published method uses 0.04, not 0.05.
#' @param sig_change `"gain"` counts only non-significant becoming significant,
#'   as in the original. `"any"` also counts the loss of significance.
#' @param qualitative Character vector of qualitative signals declared by the
#'   analyst. Any non-empty entry triggers a signal.
#' @return A `staleness_verdict`.
#' @examples
#' library(metafor)
#' prev <- rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4), measure = "RR")
#' updated <- rma(yi = c(rep(log(0.50), 4), rep(log(1.10), 4)),
#'                vi = c(rep(0.05, 4), rep(0.02, 4)), measure = "RR")
#' ottawa(prev, updated)
#'
#' # The four qualitative signals are never inferred from the data. The
#' # analyst declares them, and any non-empty entry is enough on its own:
#' # here the effect has barely moved, yet the verdict changes.
#' mild <- rma(yi = c(rep(log(0.50), 4), rep(log(0.90), 3)),
#'             vi = c(rep(0.05, 4), rep(0.02, 3)), measure = "RR")
#' ottawa(prev, mild)
#' ottawa(prev, mild, qualitative = "a superior therapy has been licensed")
#' @export
ottawa <- function(prev, new_ma, alpha = 0.04,
                   sig_change = c("gain", "any"),
                   qualitative = character()) {
  sig_change <- match.arg(sig_change)

  was_sig <- prev$pval   < alpha
  is_sig  <- new_ma$pval < alpha
  sig_signal <- if (sig_change == "gain") (!was_sig && is_sig) else (was_sig != is_sig)

  # The Ottawa effect criterion is the ratio of relative risk REDUCTIONS, not
  # of the effects -- see rrr_ratio() and Pattanittum et al. (2012), Table 1.
  r <- rrr_ratio(
    theta_new  = as.numeric(new_ma$beta),
    theta_prev = as.numeric(prev$beta),
    measure    = prev$measure,
    se_prev    = prev$se
  )
  eff_signal <- !is.na(r$ratio) && (r$ratio <= 0.5 || r$ratio >= 1.5)
  qual_signal <- length(qualitative) > 0 && any(nzchar(qualitative))

  fired <- sig_signal || eff_signal || qual_signal
  new_verdict(
    "ottawa",
    if (fired) "out_of_date" else "current",
    signal = r$ratio,
    reason = r$reason,
    detail = list(
      signal_significance = sig_signal,
      signal_effect       = eff_signal,
      signal_qualitative  = qual_signal,
      qualitative         = qualitative,
      p_prev              = prev$pval,
      p_new               = new_ma$pval
    )
  )
}
