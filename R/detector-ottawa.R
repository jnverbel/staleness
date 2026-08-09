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
#' @section Which reading of "50% change" this implements, and why:
#' The original description (Shojania et al. 2007) states the quantitative
#' signals as "changes in statistical significance or relative changes in
#' effect magnitude of at least 50%". That phrasing does not say *of what*,
#' and the full text is behind a subscription, so it could not be consulted
#' here. Declared rather than glossed over.
#'
#' The only published application of the method resolves it explicitly, and on
#' data that can be checked: Pattanittum et al. (2012), Table 1, computes the
#' change on relative risk **reductions** for ratio measures, and defers to the
#' rCMA rule for mean differences. This package follows that reading, because
#' it is the one that reproduces the ten worked reviews in that study's
#' Appendix S3: all ten fire under it, none fires under the ratio of the
#' effects themselves. A reader who takes Shojania's sentence to mean the ratio
#' of effects would get a materially different detector, and should know the
#' choice was made deliberately and on evidence.
#'
#' A second, independent source settles it the same way. Kuhnisch et al.
#' (2013), applying the method and citing Shojania for it, state signal B2 as
#' "a change in relative effect size of at least 50%" and work two examples:
#' RR 2.10 to 1.51, and RR 2.61 to 1.66, both declared to meet the criterion.
#' Neither meets it on the ratio of the effects (0.719 and 0.636, inside the
#' band) nor as a percentage change in the estimate (28% and 36%, under the
#' 50% bar). Both meet it on the ratio of risk reductions (0.464 and 0.410,
#' i.e. changes of 54% and 59%). Two papers, four worked examples, one
#' reading.
#'
#' @section The effect criterion is on risk REDUCTIONS, and is unstable:
#' For ratio measures the comparison is
#' `(1 - RR_new) / (1 - RR_prev)`, not the ratio of the effects. For mean
#' differences it defers to [rcma()], so **there** the two coincide -- but on
#' ratio measures they do not, and an `rcma` firing is not a subset of the
#' `ottawa` firings. Earlier versions of this documentation declared such a
#' containment; it was an artefact of implementing both criteria with the same
#' internal function.
#'
#' The consequence is worth stating plainly, because it is unflattering and
#' measurable. `1 - RR_prev` goes to zero as the prior effect approaches no
#' effect, so the ratio explodes exactly where the Ottawa method is meant to be
#' used -- on meta-analyses whose result is not yet significant. Measured on
#' evidence containing **no change at all**: the effect signal fires on 64% of
#' samples under a null effect, and on 0% once the effect is real and precise.
#' That is not an implementation artefact; it follows from the criterion.
#'
#' It also explains the published comparison. Ottawa flagged 34 of 80 reviews
#' where recursive CMA and Barrowman each flagged 7 -- and all 80 were null
#' meta-analyses by inclusion criterion. Confirmed outside simulation on
#' `metadat::dat.laopaiboon2015`, a null review where this detector's
#' specificity falls to 0.14 while the other four hold at 1.00.
#'
#' **This is not corrected here, on purpose.** Correcting it would mean
#' implementing a different method, and the reason this package exists is to
#' find out how the published ones behave. What is corrected is the silence:
#' `detail$effect_unstable` marks the verdicts where the denominator is near
#' zero, and `detail$rrr_prev` reports the denominator itself, so nobody
#' receives a ratio of -19 without being told it came from dividing by
#' -0.005. Same principle as [CONTAMINATED_PAIRS]: flag the result, do not
#' quietly remove or repair it.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#'   It must be on the same `measure` as `prev`: comparing a risk ratio
#'   against a mean difference does not define a ratio, and used to return one
#'   anyway, with a verdict attached.
#' @param alpha Significance threshold. The published method uses 0.04, not 0.05.
#' @param sig_change `"gain"` counts only non-significant becoming significant,
#'   as in the original. `"any"` also counts the loss of significance.
#' @param qualitative Character vector of qualitative signals declared by the
#'   analyst. Any non-empty entry triggers a signal. `NA` is an error: an
#'   unknown qualitative signal is not a present one, and `nzchar(NA)` is
#'   `TRUE`, so accepting it returned `"out_of_date"` for "I don't know". Pass
#'   `""` for a signal that was checked and found absent.
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
  # An NA alpha made both significance comparisons NA, which the code below
  # read as a change: the detector returned out_of_date without a warning.
  check_rma_uni(prev, "prev")
  check_rma_uni(new_ma, "new_ma")
  check_probability(alpha, "alpha")
  check_same_measure(prev, new_ma)
  # nzchar(NA_character_) is TRUE, so an analyst recording "unknown" for a
  # qualitative signal used to get out_of_date. Unknown is not present, and
  # this is an argument rather than a datum, so it is refused.
  # nzchar() coerces, so any non-empty value of any type counted as a declared
  # signal: qualitative = 0 fired, and so did a list. This argument carries an
  # analyst's judgement in words.
  if (!is.character(qualitative)) {
    stop("`qualitative` must be a character vector of declared signals; got ",
         class(qualitative)[1], call. = FALSE)
  }
  if (anyNA(qualitative)) {
    stop("`qualitative` has missing values; an unknown qualitative signal is ",
         "not a signal. Leave it out, or pass \"\" for one that was checked ",
         "and found absent", call. = FALSE)
  }
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
      p_new               = new_ma$pval,
      # The effect criterion divides by 1 - RR_prev. When the prior effect is
      # indistinguishable from no effect that denominator sits near zero and
      # the ratio is arbitrarily large. The verdict is left exactly as the
      # method dictates -- correcting it would be implementing a different
      # method -- and the instability is reported alongside it instead.
      rrr_prev            = r$rrr_prev,
      effect_unstable     = isTRUE(r$unstable),
      # Whether this verdict reproduces the published procedure's arithmetic.
      #
      # Shojania et al. (2007), AHRQ Technical Review 16, section Methods:
      # "we performed the updated meta-analyses by combining the original
      # pooled result with the individual results of eligible new trials" --
      # the prior pooled estimate enters as ONE point, not as its constituent
      # studies. They did that for practical reasons (the original reviews
      # often did not report individual trials) and noted it is equivalent
      # under fixed effects.
      #
      # Verified here rather than taken on faith. With method = "FE" the two
      # give identical estimates and p-values to eight decimals. With REML
      # they do not: on one worked example, beta -0.2607 against -0.2167 and
      # p 0.00072 against 0.02824, a factor of 39.
      #
      # check_currency() refits over all studies, which is the better estimate
      # and the caller's own model. Under a non-fixed-effect fit that is no
      # longer the published arithmetic, so it is flagged rather than
      # silently substituted -- the same principle as effect_unstable.
      reproduces_published = identical(prev$method, "FE")
    )
  )
}
