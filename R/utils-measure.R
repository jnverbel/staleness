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
