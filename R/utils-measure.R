# Two questions that a single list used to answer, wrongly, as though they were
# one: HOW a measure is stored, and WHETHER it is a comparative effect at all.
#
# The old RATIO_MEASURES was documented as "measures that metafor stores on the
# log scale" and then used to decide which measures the detectors could act on.
# It contained PLO and IRLN, which are single-group summaries: PLO is a
# proportion on the logit scale, IRLN an incidence rate on the log scale.
# Neither compares a treated arm against a control, so `1 - exp(theta)` -- the
# Ottawa criterion's relative risk reduction -- is not a risk reduction of
# anything. On a pooled proportion of 0.24 it returned 0.68 and Ottawa issued a
# verdict on it.
#
# PLO was wrong on the storage question too: logit is inverted by plogis(), not
# exp(), so exp(beta) returned the ODDS (0.32) where the pooled proportion was
# 0.24.

# Stored as a natural logarithm; exp() recovers the natural scale.
LOG_SCALE_MEASURES <- c("RR", "OR", "PETO", "IRR", "HR", "ROM", "MPRR", "MPOR",
                        "IRLN", "MNLN", "CVLN", "SDLN", "PLN")

# Stored as a logit; plogis() recovers the natural scale.
LOGIT_SCALE_MEASURES <- c("PLO")

# Comparative effects of one arm against another, on a ratio scale. These, and
# only these, are what the rCMA and Ottawa effect criteria are defined on: both
# ask how much a TREATMENT EFFECT moved, and a single-group summary has no
# treatment effect to move. Every measure outside this list and outside the
# difference measures is refused with a reason rather than run through a
# formula that does not apply to it.
COMPARATIVE_RATIO_MEASURES <- c("RR", "OR", "PETO", "IRR", "HR", "ROM",
                                "MPRR", "MPOR")

# Narrower still, and for a different reason. The Ottawa effect criterion is
# `(1 - RR_new) / (1 - RR_prev)`: a ratio of relative risk REDUCTIONS, so
# `1 - exp(theta)` has to be a risk reduction. That holds for a risk ratio and,
# by the published application's own usage, for an odds ratio. It does not hold
# for the rest of the comparative measures:
#
#   ROM  is a ratio of MEANS on a continuous outcome. `1 - ROM` is a percentage
#        change in a mean, not a reduction in anybody's risk.
#   HR   `1 - HR` is a relative reduction in HAZARD, which is a rate over time,
#        not a risk.
#   IRR  same, for an incidence rate.
#
# The ten worked reviews in Pattanittum et al. (2012), Appendix S3, against
# which this criterion is verified, are all risk ratios. Extending the formula
# past that is arithmetic without an argument, so the effect signal declines on
# those measures. B1 (the significance change) and the qualitative signals do
# not depend on the RRR and keep working: ottawa() still answers, it just does
# not answer on effect.
OTTAWA_RRR_MEASURES <- c("RR", "OR", "PETO", "MPRR", "MPOR")

#' Is the Ottawa relative-risk-reduction criterion defined for this measure?
#'
#' @param measure Character, the `measure` field of an `rma` object.
#' @return `TRUE` when `1 - exp(theta)` is a risk reduction.
#' @keywords internal
is_rrr_measure <- function(measure) {
  if (is.null(measure) || !nzchar(measure)) return(FALSE)
  measure %in% OTTAWA_RRR_MEASURES
}

# Single-group summaries: a proportion, a rate, a mean, a correlation. metafor
# offers many; naming them explicitly is better than inferring, because a
# measure this package has never heard of should also be refused.
SINGLE_GROUP_MEASURES <- c("PLO", "PR", "PAS", "PFT", "PLN",
                           "IR", "IRLN", "IRS", "IRFT",
                           "MN", "MNLN", "CVLN", "SDLN",
                           "COR", "UCOR", "ZCOR")

#' Is a metafor effect measure stored on the log scale?
#'
#' @param measure Character, the `measure` field of an `rma` object.
#' @return `TRUE` if `exp()` recovers the natural scale.
#' @keywords internal
is_log_scale_measure <- function(measure) {
  if (is.null(measure) || !nzchar(measure)) return(FALSE)
  measure %in% LOG_SCALE_MEASURES
}

#' Is this a comparative effect on a ratio scale?
#'
#' The question the detectors actually need answered, and the one the old
#' `is_ratio_measure()` was silently standing in for.
#'
#' @param measure Character, the `measure` field of an `rma` object.
#' @return `TRUE` for a two-arm ratio measure.
#' @keywords internal
is_comparative_ratio <- function(measure) {
  if (is.null(measure) || !nzchar(measure)) return(FALSE)
  measure %in% COMPARATIVE_RATIO_MEASURES
}

#' Is this a summary of one group rather than a comparison?
#'
#' @param measure Character, the `measure` field of an `rma` object.
#' @return `TRUE` for a proportion, rate, mean or correlation.
#' @keywords internal
is_single_group <- function(measure) {
  if (is.null(measure) || !nzchar(measure)) return(FALSE)
  measure %in% SINGLE_GROUP_MEASURES
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
  if (is_log_scale_measure(measure)) return(exp(theta))
  # A logit is inverted by plogis(), not exp(). PLO used to take the exp()
  # branch and return the odds where the caller asked for the proportion.
  if (!is.null(measure) && nzchar(measure) && measure %in% LOGIT_SCALE_MEASURES) {
    return(stats::plogis(theta))
  }
  theta
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
  # Refused before either branch: a proportion, a rate, a mean or a correlation
  # has no treatment effect for these criteria to compare, and running one
  # through the difference branch would be as meaningless as the ratio branch.
  if (is_single_group(measure)) {
    return(list(
      ratio  = NA_real_,
      reason = paste0("`", measure, "` summarises one group rather than ",
                      "comparing two; this criterion is defined on a ",
                      "treatment effect")
    ))
  }
  if (is_comparative_ratio(measure)) {
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
#' rounding alone -- and near no-effect is precisely where the Ottawa method
#' is meant to be applied.
#'
#' The difference-measure branch refuses such a case outright. This one does
#' **not**, and the asymmetry is deliberate: every one of the eighty reviews
#' in the method's published application was a null meta-analysis by inclusion
#' criterion, so a guard that rejected them would answer a different question
#' from the method's. The ratio is returned as the method defines it, with
#' `rrr_prev` and `unstable` attached so the caller can see what it was
#' divided by. Only `RR == 1` exactly, where the ratio is undefined rather
#' than merely large, yields `NA_real_`.
#'
#' @inheritParams effect_ratio
#' @return List with `ratio` (possibly `NA_real_`), `reason`, `rrr_prev` (the
#'   denominator `1 - RR_prev`, `NA` on difference measures) and `unstable`,
#'   which is `TRUE` when the prior effect sits within `min_z` standard errors
#'   of no effect and the ratio is therefore arbitrarily large. `unstable`
#'   never changes the ratio: it reports it.
#' @keywords internal
rrr_ratio <- function(theta_new, theta_prev, measure, se_prev, min_z = 2) {
  # Mean differences: the source defers to the rCMA rule. No RRR exists on
  # that scale, so there is nothing for it to be unstable about -- but the
  # fields are still filled in, so no caller has to test for their absence.
  # Anything that is not a two-arm ratio measure -- a mean difference, or a
  # single-group summary -- goes to effect_ratio(), which defers to the rCMA
  # rule for the former and refuses the latter with a reason.
  if (!is_comparative_ratio(measure)) {
    r <- effect_ratio(theta_new, theta_prev, measure, se_prev, min_z)
    r$rrr_prev <- NA_real_
    r$unstable <- FALSE
    return(r)
  }
  # A comparative ratio, but not one whose complement is a risk reduction:
  # a hazard ratio, an incidence rate ratio, a ratio of means. The criterion
  # is not defined there, so it declines rather than computing a number the
  # formula cannot mean.
  if (!is_rrr_measure(measure)) {
    return(list(
      ratio    = NA_real_,
      reason   = paste0("the Ottawa effect criterion is a ratio of relative ",
                        "risk reductions, and `1 - exp(theta)` is not a risk ",
                        "reduction for `", measure, "`"),
      rrr_prev = NA_real_,
      unstable = FALSE
    ))
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
      ratio    = NA_real_,
      reason   = "prior effect is exactly no effect; the risk-reduction ratio is undefined",
      rrr_prev = rrr_prev,
      unstable = TRUE
    ))
  }
  # The ratio is returned whatever the denominator, because that is the
  # published method. But a caller handed a ratio of -19 deserves to know it
  # came from dividing by -0.005, so the denominator travels with it, together
  # with a flag using the same test the difference branch uses to refuse
  # outright: a prior effect within min_z standard errors of no effect.
  list(ratio    = (1 - exp(theta_new)) / rrr_prev,
       reason   = "",
       rrr_prev = rrr_prev,
       unstable = is.finite(se_prev) && se_prev > 0 &&
                  abs(theta_prev) < min_z * se_prev)
}
