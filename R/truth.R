#' Ground-truth definitions for backtesting
#'
#' Three independent answers to "was this meta-analysis really out of date at
#' time t?". Independence matters: defining truth with the same rule a detector
#' uses makes that detector correct by construction, which would turn any
#' calibration result into an artefact.
#'
#' These three definitions are this package's own, and no published source
#' offers a ready-made ground truth to borrow. The closest is Shojania et al.
#' (2007), whose survival analysis counts a review as needing an update on
#' "changes in statistical significance or relative changes in effect magnitude
#' of at least 50\%", plus qualitative signals a program cannot infer.
#' [truth_conclusion()] is the definition closest to that quantitative half.
#' [truth_shift()] and [truth_surprise()] deliberately go elsewhere: they never
#' look at a significance threshold or an effect-size ratio, which is what
#' keeps them free of every detector's own logic. Treat them as three
#' operational choices, stated so they can be argued with, rather than as a
#' standard.
#'
#' @param theta_t Pooled effect at the cut point.
#' @param se_t Standard error at the cut point.
#' @param theta_final Pooled effect using all evidence.
#' @param se_final Standard error using all evidence.
#' @param p_t,p_final Two-sided p-values at the cut point and at the end.
#' @param threshold Number of standard errors for `truth_shift`.
#' @param alpha Significance level for `truth_conclusion`.
#'
#' @return `truth_shift()`, `truth_surprise()` and `truth_conclusion()` each
#'   return a single logical value: `TRUE` if the meta-analysis was out of date
#'   at the cut point under that definition, `FALSE` if it was not, and `NA`
#'   when the question cannot be answered from the values supplied.
#'
#'   `NA` is a deliberate third answer rather than a failure. `truth_shift()`
#'   and `truth_surprise()` divide by a standard error, so one that is zero,
#'   negative or non-finite makes the distance impossible to judge — not
#'   infinite.
#'   Returning `TRUE` there would score an unknowable truth as a certain event
#'   and bias every metric computed from it, which is why [calibration()] and
#'   [lead_time()] drop those rows instead.
#'
#'   `available_truths()` returns a character vector of the truth names the
#'   rest of the package accepts: `"shift"`, `"surprise"` and `"conclusion"`.
#'   These are the strings [backtest()], [calibration()] and [lead_time()] take
#'   in their `truth` arguments, and the suffixes of the three functions above.
#'
#' @section An estimate of exactly zero:
#' `truth_conclusion()` compares `sign(theta_t)` with `sign(theta_final)`, and
#' `sign(0)` is `0`, which differs from both `+1` and `-1`. An estimate
#' sitting exactly on the null therefore reads as a sign change against any
#' non-zero one. That is the intended reading — moving from "no effect at
#' all" to a definite direction is a change in the practical conclusion — and
#' on real data an exactly-zero pooled estimate is measure-zero in any case.
#' @examples
#' # Truth is measured against the FINAL body of evidence, which is why it
#' # can only ever be computed in hindsight -- and why no detector is
#' # allowed to see it.
#' theta_1970 <- -0.30; theta_final <- -0.75; se_final <- 0.12
#'
#' # Did the pooled effect end up more than 1.96 final SEs away?
#' truth_shift(theta_1970, theta_final, se_final)
#'
#' # Was that move a surprise given how precise the 1970 review thought it was?
#' truth_surprise(theta_1970, se_t = 0.30, theta_final)
#'
#' # Did the conclusion itself change -- sign or significance?
#' truth_conclusion(theta_t = -0.30, p_t = 0.21,
#'                  theta_final = -0.75, p_final = 0.001)
#'
#' available_truths()
#' @name truth
NULL

#' @rdname truth
#' @export
truth_shift <- function(theta_t, theta_final, se_final, threshold = 1.96) {
  check_scalar_input(theta_t, theta_final, se_final, arg = "truth_shift()")
  # A standard error of zero, negative or non-finite means the distance cannot
  # be judged, not that it is infinite. Dividing anyway returned TRUE, scoring
  # an unknowable truth as a certain event; NA is the honest answer, and the
  # metrics already drop NA rows for exactly this reason.
  if (!is.finite(se_final) || se_final <= 0) return(NA)
  abs(theta_final - theta_t) / se_final > threshold
}

#' @rdname truth
#' @export
truth_surprise <- function(theta_t, se_t, theta_final, threshold = 1.96) {
  check_scalar_input(theta_t, se_t, theta_final, arg = "truth_surprise()")
  if (!is.finite(se_t) || se_t <= 0) return(NA)
  abs(theta_final - theta_t) / se_t > threshold
}

#' @rdname truth
#' @export
truth_conclusion <- function(theta_t, p_t, theta_final, p_final, alpha = 0.05) {
  check_scalar_input(theta_t, p_t, theta_final, p_final,
                     arg = "truth_conclusion()")
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
#' @format A data frame with one row per circular detector-truth pair and two
#'   character columns:
#' \describe{
#'   \item{method}{Name of the detector, one of [available_methods()].}
#'   \item{truth}{Name of the truth definition whose logic it shares, one of
#'     [available_truths()].}
#' }
#'   It currently holds a single row: `ottawa` against `truth_conclusion`.
#'   The table is this package's own judgement about its own detectors, not
#'   data from an external source, so it is stated here to be argued with.
#'
#' @return A data frame; see **Format**.
#'
#' @examples
#' CONTAMINATED_PAIRS
#'
#' # Flagged, not dropped: the circular pair is still computed, and the
#' # warning travels next to the number rather than as a missing row.
#' # Filtering it out is one line for anyone who wants that instead.
#' pair <- CONTAMINATED_PAIRS[1, ]
#' paste(pair$method, "scored against truth_", pair$truth, sep = "")
#' @export
CONTAMINATED_PAIRS <- data.frame(
  method = "ottawa",
  truth  = "conclusion",
  stringsAsFactors = FALSE
)
