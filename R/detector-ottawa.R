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
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param alpha Significance threshold. The published method uses 0.04, not 0.05.
#' @param sig_change `"gain"` counts only non-significant becoming significant,
#'   as in the original. `"any"` also counts the loss of significance.
#' @param qualitative Character vector of qualitative signals declared by the
#'   analyst. Any non-empty entry triggers a signal.
#' @return A `staleness_verdict`.
#' @export
ottawa <- function(prev, new_ma, alpha = 0.04,
                   sig_change = c("gain", "any"),
                   qualitative = character()) {
  sig_change <- match.arg(sig_change)

  was_sig <- prev$pval   < alpha
  is_sig  <- new_ma$pval < alpha
  sig_signal <- if (sig_change == "gain") (!was_sig && is_sig) else (was_sig != is_sig)

  r <- effect_ratio(
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
