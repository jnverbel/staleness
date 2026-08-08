# Measures that metafor stores on the log scale. Verified against metafor 4.x:
# rma() on measure = "RR" returns beta = log(RR).
RATIO_MEASURES <- c("RR", "OR", "PETO", "IRR", "HR", "ROM", "PLO", "IRLN", "MPRR", "MPOR")

#' Is a metafor effect measure stored on the log scale?
#'
#' @param measure Character, the `measure` field of an `rma` object.
#' @return `TRUE` if the measure is a ratio stored as a logarithm.
#' @keywords internal
is_ratio_measure <- function(measure) {
  if (is.null(measure) || !nzchar(measure)) return(FALSE)
  measure %in% RATIO_MEASURES
}

#' Convert a model coefficient to its natural scale
#'
#' `metafor` stores ratio measures as logarithms. Computing a ratio of effects
#' on the log scale silently produces wrong but plausible numbers, so every
#' ratio in this package goes through here first.
#'
#' @param theta Numeric, a coefficient from an `rma` object.
#' @param measure Character, the `measure` field of that object.
#' @return The effect on its natural scale.
#' @keywords internal
to_natural <- function(theta, measure) {
  if (is_ratio_measure(measure)) exp(theta) else theta
}

#' Ratio of an updated effect to a prior effect, on the natural scale
#'
#' For difference measures the ratio is undefined when the prior effect is
#' indistinguishable from zero: it explodes and produces spurious signals. The
#' original formulations of the rCMA and Ottawa methods do not address this.
#'
#' @param theta_new,theta_prev Numeric coefficients as stored by `metafor`.
#' @param measure Character, the `measure` field.
#' @param se_prev Numeric, standard error of the prior effect.
#' @param min_z Numeric, the prior effect must be at least this many standard
#'   errors from zero for a difference ratio to be computed.
#' @return List with `ratio` (possibly `NA_real_`) and `reason`.
#' @keywords internal
effect_ratio <- function(theta_new, theta_prev, measure, se_prev, min_z = 2) {
  if (is_ratio_measure(measure)) {
    # exp(a)/exp(b) == exp(a - b), numerically safer
    return(list(ratio = exp(theta_new - theta_prev), reason = ""))
  }
  if (!is.finite(se_prev) || se_prev <= 0) {
    return(list(ratio = NA_real_, reason = "prior standard error is not usable"))
  }
  if (abs(theta_prev) < min_z * se_prev) {
    return(list(
      ratio  = NA_real_,
      reason = "prior effect is indistinguishable from zero; the ratio is unstable"
    ))
  }
  list(ratio = theta_new / theta_prev, reason = "")
}

#' Ratio of relative risk reductions, the Ottawa effect criterion
#'
#' The Ottawa method's "change in effect size of at least 50%" is **not**
#' computed on the effect itself. Pattanittum et al. (2012), Table 1, states
#' it explicitly: for treatment effects measured as a relative ratio (RR, OR)
#' the comparison is between the *relative risk reduction* of the updated
#' meta-analysis and the RRR of the previous one, where `RRR = 1 - RR`. For
#' mean differences the relative change is computed as `rcma()` computes it.
#'
#' The distinction is not cosmetic. Across the ten reviews with the largest
#' indicator in that study's Appendix S3, all ten fire on the RRR ratio and
#' **none** of them fires on the ratio of the risk ratios: a review moving
#' from RR 0.995 to RR 0.848 barely moves the risk ratio (0.85) while its risk
#' reduction goes from 0.5% to 15%, a factor of 30.
#'
#' That sensitivity is also the criterion's weakness. `1 - RR` is tiny
#' whenever the prior effect is near no-effect, so the ratio explodes on
#' rounding alone, and the same guard the difference-measure branch uses
#' applies here: a prior effect indistinguishable from no effect yields
#' `NA_real_` rather than a large and meaningless number.
#'
#' @inheritParams effect_ratio
#' @return List with `ratio` (possibly `NA_real_`) and `reason`.
#' @keywords internal
rrr_ratio <- function(theta_new, theta_prev, measure, se_prev, min_z = 2) {
  # Mean differences: the source defers to the rCMA rule.
  if (!is_ratio_measure(measure)) {
    return(effect_ratio(theta_new, theta_prev, measure, se_prev, min_z))
  }
  # NOT guarded by min_z, deliberately. The difference-measure branch refuses
  # a prior effect within min_z standard errors of zero, and transplanting
  # that rule here would disable the detector on precisely the reviews it was
  # designed for: the Ottawa method targets meta-analyses whose result is not
  # significant, and in the published application every one of the eighty
  # reviews was a null meta-analysis by inclusion criterion. A guard that
  # rejects all of them answers a different question from the method's.
  #
  # Only the genuine singularity is refused: RR exactly 1 makes RRR exactly 0
  # and the ratio undefined.
  rrr_prev <- 1 - exp(theta_prev)
  if (!is.finite(rrr_prev) || rrr_prev == 0) {
    return(list(
      ratio  = NA_real_,
      reason = "prior effect is exactly no effect; the risk-reduction ratio is undefined"
    ))
  }
  list(ratio = (1 - exp(theta_new)) / rrr_prev, reason = "")
}
