#' Build a time-ordered evidence stream from a meta-analysis
#'
#' The only part of the package that knows about time. Detectors receive
#' snapshots, never dates, which is what allows [backtest()] to reuse them
#' unmodified.
#'
#' @section What the stream carries from the model, and what it refuses:
#' Every snapshot is refitted from `yi` and `vi`, so anything the caller set
#' that changes the estimator has to travel with the stream or the snapshots
#' answer under a model nobody fitted. Carried: `method`, `test` (the
#' Knapp-Hartung adjustment moves p-values by orders of magnitude, and
#' [ottawa()] decides on p-values), `weighted`, and a `tau2` the caller fixed.
#' Verified against a direct refit for each.
#'
#' Refused, with an explanation rather than a silent drop: models with
#' moderators, whose `beta` is a vector of coefficients rather than a pooled
#' effect; models with custom per-study `weights`, which cannot follow a
#' subset in any defensible way; and fits that are not `rma.uni`.
#'
#' Not carried, because no detector reads them: `level`, `digits`, `slab`.
#'
#' @param ma An `rma.uni` object from `metafor`. Other fits (`rma.mh`,
#'   `rma.peto`) are refused with an explanatory error: every snapshot is
#'   refitted with [metafor::rma()], so accepting them would return
#'   inverse-variance estimates under a Mantel-Haenszel label, and refitting
#'   with Mantel-Haenszel is impossible from `yi` and `vi` alone. This matters
#'   in practice: Mantel-Haenszel is RevMan's default for binary outcomes, so
#'   it is what Cochrane reviews report. Both historical cases this package is
#'   validated against need it to reproduce their published intervals to the
#'   digit -- Lau et al. (1992) and the Cochrane review behind
#'   `metadat::dat.li2007`. Inverse variance gets the point estimate right in
#'   both and rounds an interval bound differently in both.
#' @param date Numeric vector of publication **years**, one entry per study,
#'   same order as the data used to fit `ma`. Missing values are an error; they
#'   are never imputed.
#'
#'   A `Date` is refused rather than converted. Every window in this package is
#'   denominated in years — `cuts = "yearly"` steps by one, and `horizon`,
#'   `window` and [lead_time()] all read as years — and a `Date` passed through
#'   `as.numeric()` becomes days since 1970, so `"yearly"` would cut once per
#'   *day* and both windows would silently become days. Convert explicitly with
#'   `as.numeric(format(date, "%Y"))`.
#' @param study_id Identifier of the study each estimate came from, one per
#'   row. **Required.** Every row used to be treated as an independent study
#'   and nothing in the stream could tell otherwise, so several outcomes, time
#'   points or arms from one trial entered as separate studies: their
#'   participants summed twice in [barrowman()], their weights counted twice in
#'   every pooled estimate, and the cumulative series moved on evidence that had
#'   not been collected twice.
#'
#'   It is an identifier rather than a promise of independence because a
#'   promise cannot be checked. If every row really is a distinct study, pass
#'   the study labels or `seq_len(k)` — the point is that the choice is made
#'   deliberately.
#' @param allow_dependence Set `TRUE` to build a stream whose `study_id` has
#'   duplicates. Refused by default: dependent estimates make one trial look
#'   like several. When allowed, the stream records it, and every metric
#'   computed from it should be read as optimistic.
#' @param ni Optional numeric vector of sample sizes per study, needed by
#'   [barrowman()]. Defaults to `ma$ni` when `metafor` recorded it.
#' @return An object of class `staleness_stream`.
#' @examples
#' library(metafor)
#' # The BCG vaccine trials, in the order metafor stores them: log risk
#' # ratios, their variances, the publication year and the total enrolled.
#' bcg <- data.frame(
#'   yi   = c(-0.89, -1.59, -1.35, -1.44, -0.22, -0.79, -1.62,
#'             0.01, -0.47, -1.37, -0.34,  0.45, -0.02),
#'   vi   = c(0.326, 0.195, 0.415, 0.020, 0.051, 0.007, 0.223,
#'            0.004, 0.056, 0.073, 0.012, 0.533, 0.071),
#'   year = c(1948, 1949, 1960, 1977, 1973, 1953, 1973,
#'            1980, 1968, 1961, 1974, 1969, 1976),
#'   ni   = c(262, 609, 451, 26465, 10877, 2992, 3174,
#'            176782, 14776, 3381, 77972, 4839, 34767)
#' )
#' ma <- rma(yi, vi, data = bcg, measure = "RR")
#'
#' # The stream sorts by date; the meta-analysis itself is order-blind.
#' stream <- evidence_stream(ma, date = bcg$year, study_id = seq_along(bcg$year), ni = bcg$ni)
#' stream
#' @export
evidence_stream <- function(ma, date, study_id, ni = NULL,
                            allow_dependence = FALSE) {
  if (!inherits(ma, "rma.uni")) {
    # Other metafor fits -- rma.mh, rma.peto -- carry every field read here,
    # but snapshot_at() refits each cut with rma.uni, and silently returning
    # inverse-variance snapshots for a Mantel-Haenszel input would be a
    # different analysis wearing the caller's label. Refitting with MH is not
    # possible either: it needs the 2x2 counts, and a stream holds only yi and
    # vi. So the limitation is real, and stated instead of papered over.
    if (inherits(ma, "rma")) {
      stop("`ma` is a ", class(ma)[1], " fit; this package refits every ",
           "snapshot with `rma()`, so it only accepts `rma.uni`. Convert the ",
           "data first -- `escalc()` then `rma()` -- and note that the pooled ",
           "estimates will be inverse-variance, not ", class(ma)[1], ".",
           call. = FALSE)
    }
    stop("`ma` must be an rma.uni object from metafor", call. = FALSE)
  }
  # A meta-regression's `beta` is a vector of coefficients, not a pooled
  # effect, and every snapshot below is refitted without moderators. Accepting
  # one would return a plain pooled analysis wearing the caller's model label.
  if (!is.null(ma$formula.mods) || NROW(ma$beta) > 1L) {
    stop("`ma` includes moderators. Snapshots are refitted without them, so ",
         "the backtest would answer a different question from the model you ",
         "fitted. Supply a model with no `mods` argument.", call. = FALSE)
  }
  # Per-study weights cannot follow a subset in any defensible way: dropping
  # studies changes what the supplied vector means. Refused rather than
  # silently ignored.
  if (!is.null(ma$weights)) {
    stop("`ma` was fitted with custom `weights`, which cannot be carried ",
         "through the subsetting each snapshot does. Fit without them, or ",
         "build the stream from `yi` and `vi` directly.", call. = FALSE)
  }
  k <- ma$k
  if (length(date) != k) {
    stop("`date` has length ", length(date), " but the meta-analysis has ",
         k, " studies", call. = FALSE)
  }
  if (anyNA(date)) {
    stop("`date` has missing values; dates are never imputed", call. = FALSE)
  }
  # Everything downstream is denominated in YEARS: cuts = "yearly" steps by
  # one, and `horizon`, `window` and lead_time() are all read as years. A Date
  # used to pass straight through as.numeric() into days since 1970, so
  # "yearly" cut once per day and the two windows silently became days. The
  # unit is refused here, at the single point every date enters, rather than
  # reinterpreted somewhere it cannot be seen.
  if (inherits(date, c("Date", "POSIXct", "POSIXlt"))) {
    stop("`date` must be the publication year as a number, not a Date: every ",
         "window in this package is measured in years. Convert with ",
         "as.numeric(format(date, \"%Y\"))", call. = FALSE)
  }
  if (!is.numeric(date)) {
    stop("`date` must be numeric, giving the publication year of each study; ",
         "got ", class(date)[1], call. = FALSE)
  }
  # anyNA() above catches NA and NaN but not the infinities, and an infinite
  # year used to travel all the way into backtest() before seq() rejected it
  # with "'to' must be a finite number" -- R's message about its own argument,
  # raised three functions away from the call that caused it. Worse, the stream
  # was usable in between: snapshot_at() happily returned a k.
  if (!all(is.finite(date))) {
    stop("`date` has infinite values; every study needs a real publication ",
         "year", call. = FALSE)
  }
  # Every row was treated as an independent study, and nothing in the stream
  # could tell otherwise: it carried yi, vi, a date and optionally ni, with no
  # identity attached. Several outcomes, time points or arms from one trial
  # therefore entered as separate studies -- their participants summed twice
  # in barrowman(), their weights counted twice in every pooled estimate, and
  # the cumulative series moved on evidence that had not been collected twice.
  #
  # The package cannot detect that without being told who is who, so it asks.
  # A required argument rather than an optional one, and rather than a boolean
  # promising independence: a promise cannot be checked, and an identifier can.
  if (missing(study_id)) {
    stop("`study_id` is required: one identifier per row, naming the study ",
         "each estimate came from. Without it, several outcomes or time ",
         "points from one trial enter as separate studies and are counted ",
         "twice. If every row really is a distinct study, pass the study ",
         "labels or seq_len(k).", call. = FALSE)
  }
  if (length(study_id) != k) {
    stop("`study_id` has length ", length(study_id), " but the meta-analysis ",
         "has ", k, " studies", call. = FALSE)
  }
  if (anyNA(study_id)) {
    stop("`study_id` has missing values; an estimate of unknown provenance ",
         "cannot be checked for independence", call. = FALSE)
  }
  dup <- duplicated(study_id)
  if (any(dup) && !allow_dependence) {
    n_dup <- length(unique(study_id[dup]))
    stop(n_dup, " stud", if (n_dup == 1) "y contributes" else "ies contribute",
         " more than one estimate (", paste(utils::head(unique(study_id[dup]), 3),
         collapse = ", "), if (n_dup > 3) ", ..." else "", "). Those estimates ",
         "are dependent: pooling them treats one trial as several, and ",
         "barrowman() would count its participants more than once. Select one ",
         "estimate per study before building the stream, or pass ",
         "`allow_dependence = TRUE` to proceed knowing the metrics will be ",
         "optimistic.", call. = FALSE)
  }
  dependent <- any(dup)

  ni_supplied <- !is.null(ni)
  if (is.null(ni)) ni <- ma$ni
  if (!is.null(ni) && length(ni) != k) {
    stop("`ni` has length ", length(ni), " but the meta-analysis has ", k,
         " studies", call. = FALSE)
  }
  # Only when the caller actually passed `ni`. An NA in a sample size someone
  # typed is a malformed argument and gets rejected next to the identical
  # check on `date`. An NA in the `ni` metafor derived on its own is a fact
  # about the dataset, and refusing to build the stream over it would take
  # down the four detectors that never look at `ni` -- for a value the caller
  # never asked to use. barrowman() is the only consumer, and it already
  # answers `not_applicable` naming the non-finite n, so nothing downstream
  # reports a missing argument that was never missing.
  if (ni_supplied && anyNA(ni)) {
    stop("`ni` has missing values; sample sizes are never imputed",
         call. = FALSE)
  }
  # Same side of the line, and for the same reason: an infinite or negative
  # sample size someone typed is a malformed argument. barrowman() sums these
  # over a snapshot, so a negative one silently shrinks the total it divides
  # by, and an infinite one makes it Inf. Checked only when supplied, so a
  # derived `ni` still cannot take down the four detectors that never read it.
  if (ni_supplied && (!all(is.finite(ni)) || any(ni <= 0))) {
    stop("`ni` must be finite and positive: it is a count of participants, ",
         "and barrowman() sums it across a snapshot", call. = FALSE)
  }

  ord <- order(date, seq_along(date))  # stable: ties keep input order
  structure(
    list(
      yi      = as.numeric(ma$yi)[ord],
      vi      = as.numeric(ma$vi)[ord],
      ni      = if (is.null(ni)) NULL else as.numeric(ni)[ord],
      date    = date[ord],
      study_id = study_id[ord],
      # TRUE only when the caller opted into dependent estimates. It travels
      # with the stream so that print(), backtest() and the metrics can say
      # so, rather than reporting rates that quietly assume independence.
      dependent = dependent,
      measure = ma$measure,
      method  = ma$method,
      # Carried so every snapshot is tested the way the caller asked. The
      # Knapp-Hartung adjustment moves p-values by orders of magnitude, and
      # ottawa() decides on p-values.
      test    = if (is.null(ma$test)) "z" else ma$test,
      # weighted = FALSE and a user-fixed tau2 both change beta, se and pval.
      # Anything that moves the estimator has to travel with the stream, or
      # the snapshots answer under a model the caller never fitted.
      weighted = if (is.null(ma$weighted)) TRUE else isTRUE(ma$weighted),
      tau2_fix = if (isTRUE(ma$tau2.fix)) ma$tau2 else NULL,
      k       = k
    ),
    class = "staleness_stream"
  )
}

#' Refit the meta-analysis as it stood at a point in time
#'
#' @param stream A `staleness_stream`.
#' @param cut Numeric or Date, the cut point. Studies with `date <= cut` are used.
#' @return An `rma.uni` object.
#' @examples
#' library(metafor)
#' bcg <- data.frame(
#'   yi   = c(-0.89, -1.59, -1.35, -1.44, -0.22, -0.79, -1.62,
#'             0.01, -0.47, -1.37, -0.34,  0.45, -0.02),
#'   vi   = c(0.326, 0.195, 0.415, 0.020, 0.051, 0.007, 0.223,
#'            0.004, 0.056, 0.073, 0.012, 0.533, 0.071),
#'   year = c(1948, 1949, 1960, 1977, 1973, 1953, 1973,
#'            1980, 1968, 1961, 1974, 1969, 1976)
#' )
#' stream <- evidence_stream(rma(yi, vi, data = bcg, measure = "RR"),
#'                           date = bcg$year, study_id = seq_along(bcg$year))
#'
#' # The review as it stood in 1970: 7 of the 13 trials had been published.
#' prev <- snapshot_at(stream, 1970)
#' prev$k
#'
#' # Refitting is what keeps a backtest honest -- the 1970 snapshot knows
#' # nothing about the six trials that came later, tau^2 included.
#' c(tau2_1970 = prev$tau2, tau2_all = snapshot_at(stream, 1980)$tau2)
#' @export
snapshot_at <- function(stream, cut) {
  check_class(stream, "staleness_stream", "stream", "evidence_stream()")
  check_cut_point(cut, "cut")
  keep <- stream$date <= cut
  if (sum(keep) < 2) {
    stop("a snapshot needs at least 2 studies; found ", sum(keep),
         " at cut ", format(cut), call. = FALSE)
  }
  metafor::rma(yi = stream$yi[keep], vi = stream$vi[keep],
               measure  = stream$measure, method = stream$method,
               test     = if (is.null(stream$test)) "z" else stream$test,
               weighted = if (is.null(stream$weighted)) TRUE else stream$weighted,
               tau2     = stream$tau2_fix)
}

#' Studies published in a half-open time window
#'
#' @param stream A `staleness_stream`.
#' @param from,to Bounds of the interval `(from, to]`.
#' @return A list with `yi`, `vi`, `ni` and `k`.
#' @examples
#' library(metafor)
#' bcg <- data.frame(
#'   yi   = c(-0.89, -1.59, -1.35, -1.44, -0.22, -0.79, -1.62,
#'             0.01, -0.47, -1.37, -0.34,  0.45, -0.02),
#'   vi   = c(0.326, 0.195, 0.415, 0.020, 0.051, 0.007, 0.223,
#'            0.004, 0.056, 0.073, 0.012, 0.533, 0.071),
#'   year = c(1948, 1949, 1960, 1977, 1973, 1953, 1973,
#'            1980, 1968, 1961, 1974, 1969, 1976),
#'   ni   = c(262, 609, 451, 26465, 10877, 2992, 3174,
#'            176782, 14776, 3381, 77972, 4839, 34767)
#' )
#' stream <- evidence_stream(rma(yi, vi, data = bcg, measure = "RR"),
#'                           date = bcg$year, study_id = seq_along(bcg$year), ni = bcg$ni)
#'
#' # The interval is half-open, (from, to], so a trial published exactly in
#' # 1970 would belong to the snapshot, never to the new evidence.
#' new <- window_between(stream, 1970, 1980)
#' new$k
#' sum(new$ni)
#' @export
window_between <- function(stream, from, to) {
  check_class(stream, "staleness_stream", "stream", "evidence_stream()")
  check_cut_point(from, "from")
  check_cut_point(to, "to")
  # A backwards interval is empty by construction. Returning k = 0 for it
  # reads as "we looked and found nothing", which is not what happened.
  if (from > to) {
    stop("`from` must be before `to`; got ", format(from), " and ", format(to),
         call. = FALSE)
  }
  keep <- stream$date > from & stream$date <= to
  list(
    yi = stream$yi[keep],
    vi = stream$vi[keep],
    ni = if (is.null(stream$ni)) NULL else stream$ni[keep],
    k  = sum(keep)
  )
}

#' @export
print.staleness_stream <- function(x, ...) {
  cat("<staleness_stream>\n")
  cat("  studies :", x$k, "\n")
  cat("  measure :", x$measure, "\n")
  cat("  span    :", format(min(x$date)), "to", format(max(x$date)), "\n")
  # Said out loud, every time. A stream built with allow_dependence = TRUE
  # produces metrics that assume independent studies over evidence that is
  # not, and the only place a reader might notice is here.
  if (isTRUE(x$dependent)) {
    n_studies <- length(unique(x$study_id))
    cat("  NOTE    :", x$k, "estimates from", n_studies, "studies;",
        "dependence was allowed.\n")
    cat("            Rates computed from this stream are optimistic:",
        "one trial\n            counts as several.\n")
  }
  invisible(x)
}
