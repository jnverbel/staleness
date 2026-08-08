#' Ground-truth definitions for backtesting
#'
#' Three independent answers to "was this meta-analysis really out of date at
#' time t?". Independence matters: defining truth with the same rule a detector
#' uses makes that detector correct by construction, which would turn any
#' calibration result into an artefact.
#'
#' @param theta_t Pooled effect at the cut point.
#' @param se_t Standard error at the cut point.
#' @param theta_final Pooled effect using all evidence.
#' @param se_final Standard error using all evidence.
#' @param p_t,p_final Two-sided p-values at the cut point and at the end.
#' @param threshold Number of standard errors for `truth_shift`.
#' @param alpha Significance level for `truth_conclusion`.
#' @name truth
NULL

#' @rdname truth
#' @export
truth_shift <- function(theta_t, theta_final, se_final, threshold = 1.96) {
  abs(theta_final - theta_t) / se_final > threshold
}

#' @rdname truth
#' @export
truth_surprise <- function(theta_t, se_t, theta_final, threshold = 1.96) {
  abs(theta_final - theta_t) / se_t > threshold
}

#' @rdname truth
#' @export
truth_conclusion <- function(theta_t, p_t, theta_final, p_final, alpha = 0.05) {
  sign_flip <- sign(theta_t) != sign(theta_final)
  sig_flip  <- (p_t < alpha) != (p_final < alpha)
  sign_flip || sig_flip
}

#' @rdname truth
#' @export
available_truths <- function() c("shift", "surprise", "conclusion")

#' Detector-truth pairs that share logic
#'
#' `ottawa` signals on a change of significance, which is exactly what
#' `truth_conclusion` measures. Any metric computed for that pair is circular.
#'
#' Such pairs are **flagged, not dropped**: [calibration()] returns a
#' `contaminated` column on every row, looked up from this table and sitting
#' next to the numbers it warns about. Dropping a row makes a reader notice an
#' absence; flagging it puts the warning in front of them. It also leaves the
#' circular comparison available to anyone who actually wants it — "does
#' `ottawa` reproduce its own rule correctly" is a legitimate question, just a
#' different one from "does `ottawa` predict the future" — and filtering on
#' `contaminated` takes one line for anyone who does not.
#'
#' @export
CONTAMINATED_PAIRS <- data.frame(
  method = "ottawa",
  truth  = "conclusion",
  stringsAsFactors = FALSE
)
