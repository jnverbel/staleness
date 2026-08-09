#' Operational evaluation targets for backtesting
#'
#' Three ways of asking "did this meta-analysis change by time t?", used to
#' score detector verdicts. They are **operational choices, not ground truth**,
#' and the distinction is not modesty: nothing here observes what a review
#' team actually did, whether a recommendation changed, or whether anyone was
#' harmed by acting on the old estimate. They observe the estimate moving.
#'
#' @section What each one measures, and what they share:
#' All three read the same underlying event -- the pooled estimate moved
#' between the cut and the target -- and differ in how they scale or classify
#' it:
#'
#' * [truth_shift()] and [truth_surprise()] use the **identical numerator**,
#'   `|theta_target - theta_t|`, and differ only in the denominator: the
#'   standard error at the target, or the standard error at the cut. They are
#'   one distance on two scales, not two independent facts, and a result that
#'   holds under one and not the other is a statement about precision, not
#'   about the detector.
#' * [truth_conclusion()] is categorical rather than metric: a flip in the sign
#'   of the effect, or across a significance threshold.
#'
#' Because the first two share a numerator, treat agreement between them as
#' expected and disagreement as informative -- the reverse of how independent
#' definitions would be read. Use them for sensitivity analysis across scalings
#' rather than as three votes.
#'
#' @section The 1.96 is a scaling, not a test:
#' `truth_shift()` divides by the standard error of the TARGET estimate, and
#' `theta_t` is nested inside it -- the cut's studies are part of the target's.
#' The variance of the difference between two nested estimates is not the
#' variance of either, so `|theta_target - theta_t| / se_target` is not a test
#' statistic and `1.96` is not a critical value. It is a threshold expressed in
#' units of final precision, chosen because those units are interpretable, and
#' it should be read that way.
#'
#' @section Why not a published definition:
#' No published source offers one to borrow. The closest is Shojania et al.
#' (2007), whose survival analysis counts a review as needing an update on
#' "changes in statistical significance or relative changes in effect magnitude
#' of at least 50%", plus qualitative signals a program cannot infer.
#' [truth_conclusion()] is the closest to that quantitative half.
#'
#' What all three do avoid is circularity: scoring a detector against a target
#' built from the same rule it uses makes it correct by construction.
#' [truth_shift()] and [truth_surprise()] never look at a significance
#' threshold or an effect-size ratio, so no detector's logic leaks into them;
#' `truth_conclusion` does share logic with [ottawa()], which is exactly why
#' that pair is named in [CONTAMINATED_PAIRS].
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
  check_positive_number(threshold, "threshold")
  # A standard error of zero, negative or non-finite means the distance cannot
  # be judged, not that it is infinite. Dividing anyway returned TRUE, scoring
  # an unknowable truth as a certain event; NA is the honest answer, and the
  # metrics already drop NA rows for exactly this reason.
  if (!is.finite(se_final) || se_final <= 0) return(NA)
  # And the same treatment for the effects, which did not get it. An infinite
  # theta made the distance Inf, which cleared any threshold and returned TRUE
  # -- a certain event conjured from a quantity that cannot be compared. That
  # is worse than a wrong answer here: a spurious TRUE is a free hit for every
  # detector that fires at that cut, so it inflates sensitivity and biases the
  # calibration the package exists to produce.
  if (!is.finite(theta_t) || !is.finite(theta_final)) return(NA)
  abs(theta_final - theta_t) / se_final > threshold
}

#' @rdname truth
#' @export
truth_surprise <- function(theta_t, se_t, theta_final, threshold = 1.96) {
  check_scalar_input(theta_t, se_t, theta_final, arg = "truth_surprise()")
  check_positive_number(threshold, "threshold")
  if (!is.finite(se_t) || se_t <= 0) return(NA)
  if (!is.finite(theta_t) || !is.finite(theta_final)) return(NA)
  abs(theta_final - theta_t) / se_t > threshold
}

#' @rdname truth
#' @export
truth_conclusion <- function(theta_t, p_t, theta_final, p_final, alpha = 0.05) {
  check_scalar_input(theta_t, p_t, theta_final, p_final,
                     arg = "truth_conclusion()")
  check_p_value(p_t, "p_t")
  check_p_value(p_final, "p_final")
  check_probability(alpha, "alpha")
  # sign(Inf) is 1, so an infinite effect read as a definite direction and any
  # comparison against it "changed sign". Same reasoning as the other two: a
  # quantity that cannot be located cannot have moved.
  if (!is.finite(theta_t) || !is.finite(theta_final)) return(NA)
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
