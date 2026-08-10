#' Rows of a backtest that may be scored against a truth definition
#'
#' The one filter shared by [calibration()] and [lead_time()]. Three kinds of
#' row are dropped, none of them counted as either a hit or a miss: `censored`
#' cuts, `not_applicable` verdicts, and rows whose truth column is `NA`
#' because the standard error it divides by was degenerate.
#'
#' @param bt A `staleness_backtest`.
#' @param truth One of [available_truths()].
#' @return The eligible subset of `bt$results`.
#' @keywords internal
eligible_rows <- function(bt, truth) {
  col <- paste0("truth_", truth)
  res <- bt$results
  res[!res$censored &
      res$verdict != "not_applicable" &
      !is.na(res[[col]]), ]
}

#' Calibration of the detectors against an evaluation target
#'
#' Turns a backtest's raw results into sensitivity, specificity and false
#' alarm rate, one row per method, scored against a chosen operational target
#' (see [truth]). Those targets observe the pooled estimate moving; they do not
#' observe what a review team did, so these are rates against a stated
#' criterion rather than against an external outcome.
#'
#' Three kinds of rows are excluded from the count before any metric is
#' computed, none of them counted as either a hit or a miss:
#' \itemize{
#'   \item `censored` cuts, too close to the end of the series to be fairly
#'     evaluated (see [backtest()]).
#'   \item `not_applicable` verdicts, where the detector declined to answer.
#'   \item Rows whose truth column for the requested `truth` is `NA`. This
#'     happens when [truth_shift()] or [truth_surprise()] divide by a
#'     degenerate (zero) standard error and cannot say what should have
#'     happened; `backtest()` deliberately lets that `NA` propagate rather
#'     than guessing. Left unfiltered, it would flow into `sum(hit & ev)` and
#'     silently turn a whole cell into `NA` instead of merely omitting one
#'     row.
#' }
#'
#' `ottawa` signals on a change of significance, which is exactly what
#' `truth_conclusion` measures: scored against that truth it is correct by
#' construction. Rather than silently drop that pair or footnote it, every
#' row this function returns carries a `contaminated` flag looked up from
#' [CONTAMINATED_PAIRS], so no downstream reader can miss it.
#'
#' The rows are one per method **requested in the backtest** (`bt$methods`),
#' not one per method that happened to survive the filter. A detector that was
#' `not_applicable` at every cut — as `barrowman()` and `simulation()` are on
#' any consistently significant series — gets a row with `n = 0` and `NA`
#' metrics. "This detector never applied to this evidence" is itself a result
#' about the detector, and belongs in the table as a row rather than as an
#' absence the reader has to notice.
#'
#' @section These rates describe one series; they do not estimate a rate:
#' `n` counts **cuts**, and consecutive cuts of one review share almost every
#' study — the snapshot at 1970 and the one at 1971 differ by whatever appeared
#' in a single year. Nineteen cuts are not nineteen observations, so what comes
#' back is a description of this body of evidence and not an estimate of how a
#' detector behaves in general.
#'
#' That is also why no interval is returned. A binomial one computed on `n`
#' would treat overlapping snapshots as independent draws and come out far too
#' narrow: the honest answer on a single series is no interval, not a
#' confident-looking one.
#'
#' For a rate that generalises, use [pooled_calibration()] over several
#' independent reviews. Reviews share no studies, so resampling *them* gives an
#' interval that means something.
#'
#' A second caution, measured rather than theoretical: against
#' `truth_target = "final"` a review whose effect kept moving has **no true
#' negatives at all** — every earlier cut counts as out of date — so its
#' specificity is `NA` rather than perfect. That is the case in 6 of the 17
#' reviews swept in `inst/applicability/`.
#'
#' @param bt A `staleness_backtest`, see [backtest()].
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`, see
#'   [available_truths()].
#' @return A data frame, one row per method in `bt$methods`, with columns
#'   `method`, `truth`, `sensitivity`, `specificity`, `false_alarm`, `n` and
#'   `contaminated` and `dependent`. `dependent` is `TRUE` when the stream
#'   behind the backtest was built with `allow_dependence = TRUE`, in which
#'   case one trial contributed several estimates and every rate in the row is
#'   optimistic. Descriptive of this series; see the section above.
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
#' bt <- backtest(evidence_stream(rma(yi, vi, data = bcg, measure = "RR"),
#'                                date = bcg$year, study_id = seq_along(bcg$year), ni = bcg$ni))
#'
#' # Every requested method gets a row. barrowman and simulation never apply
#' # to this evidence -- the prior review is significant throughout -- so they
#' # come back with n = 0 and NA metrics rather than quietly vanishing.
#' calibration(bt)
#'
#' # The truth definition is a choice, not a constant, and the answer moves
#' # with it. `contaminated` marks the pairs where a detector is being scored
#' # against a truth built from its own rule.
#' calibration(bt, truth = "conclusion")
#' CONTAMINATED_PAIRS
#' @export
calibration <- function(bt, truth = "shift") {
  check_class(bt, "staleness_backtest", "bt", "backtest()")
  truth <- match.arg(truth, available_truths())
  col <- paste0("truth_", truth)
  res <- eligible_rows(bt, truth)

  out <- lapply(bt$methods, function(m) {
    d  <- res[res$method == m, ]
    hit <- d$verdict == "out_of_date"
    ev  <- d[[col]]
    tp <- sum(hit & ev); fn <- sum(!hit & ev)
    tn <- sum(!hit & !ev); fp <- sum(hit & !ev)
    data.frame(
      method       = m,
      truth        = truth,
      sensitivity  = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
      specificity  = if (tn + fp > 0) tn / (tn + fp) else NA_real_,
      false_alarm  = if (fp + tn > 0) fp / (fp + tn) else NA_real_,
      n            = nrow(d),
      contaminated = any(CONTAMINATED_PAIRS$method == m &
                         CONTAMINATED_PAIRS$truth  == truth),
      # Travels beside the rate it qualifies, for the same reason
      # `contaminated` does: a reader who receives this table and not the
      # stream has no other way to learn that one trial counted as several.
      dependent    = isTRUE(bt$dependent),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

#' How early a detector fired before the evidence actually changed
#'
#' The metric that decides practical usefulness. A detector that only fires
#' in the same period the evidence already moved looks perfect in a
#' contingency table and is useless in practice: `calibration()` cannot tell
#' the two apart, because it scores each cut in isolation.
#'
#' `lead_time()` looks at every true event (per the chosen `truth`
#' definition), not just the first. For each true event it finds the
#' detector's earliest `out_of_date` firing at or before that event's cut
#' and records `event_cut - firing_cut` as that event's lead; an event the
#' detector never flagged in time (no firing at or before it) contributes no
#' lead at all rather than a miss disguised as a number. `median_lead` is
#' the median of the leads that *are* defined, i.e. the median time-to-event
#' across the events the detector caught early or on time — the events it
#' missed entirely are excluded from that median, not folded into it as
#' zero. A `median_lead` of 0 is a real, meaningful value: it means the
#' detector, in the middle of its caught events, only ever fired in the same
#' period the evidence had already moved, not that it failed to catch
#' anything.
#'
#' @section Why the two sides of the comparison are filtered differently:
#' This is the only metric here that relates *different* rows: the event comes
#' from one cut and the firing that preceded it from another. [calibration()]
#' scores each row on its own, so it needs the verdict and the truth from the
#' same row and must drop the row when either is missing. That does not carry
#' over.
#'
#' Censored cuts and `not_applicable` verdicts are excluded from both sides,
#' as in [calibration()]. A truth of `NA` is excluded from the **events**
#' only: an event is defined by `truth == TRUE`, and unknown is not true. It
#' is *not* excluded from the **firings**. A firing is an observed fact about
#' the detector, and its role here is to precede a later event; whether the
#' truth of its own cut could be computed says nothing about that.
#'
#' Filtering the firings too — which earlier versions did, by reusing
#' `calibration()`'s eligibility filter unchanged — deletes real early
#' warnings. A detector that fired at 2000 where truth was indeterminate, for
#' an event at 2001, reported `median_lead = NA` instead of 1. The bias has a
#' direction: dropping a firing can only ever lengthen a lead or erase it,
#' never shorten it, so it always reads against the detector.
#'
#' A method with no true event under `truth` in the uncensored window gets
#' `n_events = 0` and `median_lead = NA`: there is nothing to lead. A method
#' that never fired `out_of_date` at or before any of its true events also
#' gets `median_lead = NA`, with `n_events` still reporting how many events
#' it missed (visible in detail via `calibration()`'s sensitivity). As in
#' [calibration()], every method in `bt$methods` gets a row, including one
#' that was `not_applicable` at every cut and therefore has no eligible rows
#' at all.
#'
#' @section How long a warning stays relevant:
#' A firing counts as advance warning for an event if it happened at or before
#' that event and no earlier than `within` units before it. The default,
#' `Inf`, keeps the definition above: the earliest firing at or before the
#' event, however distant.
#'
#' That default is generous, and worth seeing plainly. With events at 2005,
#' 2007 and 2009, a detector that fires once in 2000 and stays silent for nine
#' years scores a median lead of 7 -- the same as one that fires at every
#' single cut -- while its sensitivity is 0, because it flagged none of those
#' events at their own cut. A warning that old is not advance notice of
#' anything in particular.
#'
#' Set `within` to bound it. `within = bt$horizon` restricts a warning to the
#' span over which the backtest considers truth evaluable, and turns that
#' 7 into `NA`. Nothing here picks a value for you, because the right one
#' depends on how fast the field moves; the default is stated rather than
#' hidden, and `calibration()`'s sensitivity should always be read alongside
#' this column.
#'
#' @param bt A `staleness_backtest`, see [backtest()].
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`, see
#'   [available_truths()].
#' @param within Maximum age of a firing for it to count as advance warning of
#'   an event, in the same units as the cuts. Defaults to the backtest's own
#'   `horizon`, which is the window the backtest was built around.
#'
#'   The old default was `Inf`, and it flattered every detector: one firing in
#'   1955 counted as advance warning for an event in 1980, so a detector that
#'   fired once and then went quiet scored the same median lead as one that
#'   warned before each event in turn. `Inf` is still available and still
#'   meaningful — "did this detector ever fire before the event" is a real
#'   question — but as something a caller asks for, not something they get
#'   without choosing it.
#' @return A data frame, one row per method in `bt$methods`, with columns
#'   `method`, `median_lead`, `n_events` and `dependent` (see
#'   [calibration()] for what the last one means).
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
#' bt <- backtest(evidence_stream(rma(yi, vi, data = bcg, measure = "RR"),
#'                                date = bcg$year, study_id = seq_along(bcg$year), ni = bcg$ni))
#'
#' # Whether a detector eventually fires is the easy question. This is the
#' # one the methods papers never answer: how far ahead of the evidence.
#' # `median_lead` is NA when a detector caught no true event at all, and
#' # `n_events` still says how many it had the chance to catch.
#' lead_time(bt)
#' @export
lead_time <- function(bt, truth = "shift", within = NULL) {
  check_class(bt, "staleness_backtest", "bt", "backtest()")
  # The backtest's own horizon, not Inf: a warning older than the window the
  # backtest was built around is not advance warning of this event, it is a
  # firing that happened to precede it.
  if (is.null(within)) within <- bt$horizon
  truth <- match.arg(truth, available_truths())
  if (!is.numeric(within) || length(within) != 1L || is.na(within) ||
      within <= 0) {
    stop("`within` must be a single positive number, or Inf", call. = FALSE)
  }
  col <- paste0("truth_", truth)

  # This is the one metric that relates DIFFERENT rows -- an event in one row,
  # the firing that preceded it in another -- and the two sides do not carry
  # the same requirement.
  #
  # An event is defined by `truth == TRUE`, so a row whose truth is NA cannot
  # supply one: unknown is not true. A firing is an observed fact about the
  # detector, and its job here is to precede a later event. Whether the truth
  # of the firing's OWN cut could be computed says nothing about that.
  #
  # eligible_rows() applies `!is.na(truth)` to both sides, which is right for
  # calibration() -- that scores each row in isolation and needs the verdict
  # and the truth from the same row -- and wrong here. Applied to the firings
  # it deletes real early warnings and biases every lead in one direction,
  # since it can only ever remove a candidate, never add one.
  usable <- bt$results
  usable <- usable[!usable$censored & usable$verdict != "not_applicable", ]
  scored <- usable[!is.na(usable[[col]]), ]

  out <- lapply(bt$methods, function(m) {
    d <- scored[scored$method == m, ]
    d <- d[order(d$cut), ]
    events <- d$cut[d[[col]]]
    f <- usable[usable$method == m, ]
    fired  <- f$cut[f$verdict == "out_of_date"]

    # One lead per true event: the gap to that event's earliest on-time
    # firing, or NA when the detector never fired at or before it. A miss
    # must never contribute a numeric lead (e.g. 0), or it would be
    # indistinguishable from genuine same-period detection.
    leads <- vapply(events, function(event_cut) {
      # `within` bounds how long a firing stays relevant. With the default of
      # Inf a single early firing counts as advance warning for every later
      # event, however distant -- see the note in the documentation.
      on_time <- fired[fired <= event_cut & (event_cut - fired) <= within]
      if (length(on_time)) event_cut - min(on_time) else NA_real_
    }, numeric(1))

    median_lead <- if (any(!is.na(leads))) {
      stats::median(leads, na.rm = TRUE)
    } else {
      NA_real_
    }

    data.frame(method = m, median_lead = median_lead, n_events = length(events),
               # Same reason as in calibration(): the caveat travels with the
               # number, not with the object the number came from.
               dependent = isTRUE(bt$dependent),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Calibration against all three truth definitions, stacked
#'
#' Runs [calibration()] once per entry of [available_truths()] and stacks the
#' results, so a single call shows sensitivity, specificity and contamination
#' side by side across `shift`, `surprise` and `conclusion`.
#'
#' @param object A `staleness_backtest`, see [backtest()].
#' @param ... Unused; present for S3 consistency with [summary()].
#' @return A data frame, the row-bound output of [calibration()] for each
#'   truth definition.
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
#' bt <- backtest(evidence_stream(rma(yi, vi, data = bcg, measure = "RR"),
#'                                date = bcg$year, study_id = seq_along(bcg$year)))
#'
#' # All three truth definitions at once, stacked. Reading them side by side
#' # is the point: a detector that looks good under one and bad under another
#' # is telling you about the truth definition, not only about itself.
#' summary(bt)
#' @export
summary.staleness_backtest <- function(object, ...) {
  do.call(rbind, lapply(available_truths(), function(t) calibration(object, t)))
}

#' Calibration pooled across independent reviews, with intervals
#'
#' [calibration()] describes one body of evidence. Its `n` counts cuts, and
#' consecutive cuts of a single review share almost every study: the snapshot
#' at 1970 and the snapshot at 1971 differ by whatever appeared in one year.
#' Nineteen cuts are not nineteen observations, so a rate computed from them
#' is a **description of that series**, not an estimate of how a detector
#' behaves in general, and it carries no interval for the same reason — a
#' binomial one would be built on a denominator that is not what it looks
#' like, and would come out far too narrow.
#'
#' Independence is available at the level of **reviews**. Different reviews
#' share no studies, so this function pools across them and resamples whole
#' reviews to get an interval. What varies between bootstrap replicates is
#' which reviews were drawn, which is the variability that matters when asking
#' whether a result would hold on other evidence.
#'
#' @section Two estimands, and why the choice cannot be made for you:
#' "Pooled sensitivity" is not one quantity. Summing the 2x2 counts across
#' reviews and dividing once answers **"across all the cuts in this corpus,
#' what share of the truly out-of-date ones did the detector flag?"** — a
#' review with thirty cuts contributes six times what a review with five
#' contributes. Averaging each review's own rate answers **"for a typical
#' review, what share?"** — every review counts once, however long its series.
#'
#' Those are different questions with different answers, and neither is a
#' better version of the other. `weighting = "cut"` gives the first,
#' `weighting = "review"` the second, and the returned table carries a
#' `weighting` column so a reader can never be left guessing which one is in
#' front of them.
#'
#' The gap between them is a finding, not noise: they diverge exactly when
#' detector behaviour depends on series length, which is the case for
#' [barrowman()] and [simulation()], applicable only while the prior review is
#' non-significant and therefore mostly on short early series. Reporting both
#' is cheap.
#'
#' @section Why you have to declare that the reviews are independent:
#' The bootstrap draws whole reviews and the pooling adds their counts. Both
#' steps assume the reviews share no studies. Nothing here can check that: a
#' `study_id` is a label chosen by whoever built each stream, and two reviews
#' that each numbered their studies `1:14` would look identical and be
#' unrelated, while two reviews of overlapping literature might use different
#' labels for the same trial. An automatic check would be wrong in both
#' directions.
#'
#' So `reviews_independent` has no default and must be given as `TRUE`. It is
#' the same reasoning as `study_id` in [evidence_stream()] — a promise that
#' cannot be checked is asked for explicitly rather than assumed — with the
#' difference that here there is nothing to identify, so a declaration is all
#' that is left. Two backtests of the same evidence, or of two outcomes from
#' one review, are not two reviews.
#'
#' @param bts A list of `staleness_backtest` objects, one per independent
#'   review. Two backtests of the same evidence are not two reviews.
#' @param truth One of [available_truths()].
#' @param R Bootstrap replicates.
#' @param seed Integer seed, or `NULL`. Without one the interval is not
#'   reproducible; the caller's random stream is restored either way.
#' @param conf Interval coverage.
#' @param weighting `"cut"` weights every cut equally, so long reviews count
#'   more; `"review"` weights every review equally. See the section above:
#'   these estimate different things, and the answer is reported with the
#'   choice attached.
#' @param reviews_independent Must be given as `TRUE`. There is no default and
#'   no check; see the section above for why the declaration is asked for
#'   rather than inferred.
#' @param accept_dependence Set `TRUE` to receive bootstrap bounds even when a
#'   pooled review was built with `allow_dependence = TRUE`. Withheld by
#'   default, with a warning: the bootstrap resamples whole reviews, which is
#'   the right unit only when each review is one independent body of studies,
#'   and drawing a dependent review again cannot undo the double counting
#'   inside it. Point estimates are descriptive and always returned; this
#'   argument governs only the interval, which is the one output here that
#'   reads as an inferential claim.
#' @return A data frame with one row per method: `weighting`, `n_reviews`,
#'   `n_cuts`, the pooled `sensitivity` and `specificity`, percentile bounds
#'   for each, and `contaminated` and `dependent` flags.
#'   Bounds are `NA` when fewer than two reviews contribute a defined rate --
#'   an interval from one review would describe nothing but that review -- and
#'   when `dependent` is `TRUE` and `accept_dependence` is `FALSE`.
#' @examples
#' # Two independent reviews, pooled. The interval comes from resampling the
#' # reviews, so it says what would happen on other bodies of evidence -- not
#' # what would happen on other cuts of these two.
#' library(metafor)
#' mk <- function(seed) {
#'   set.seed(seed)
#'   yi <- cumsum(rnorm(14, -0.05, 0.15)); vi <- runif(14, 0.02, 0.08)
#'   evidence_stream(rma(yi, vi, method = "FE"), date = 1990:2003,
#'                   study_id = seq_len(14))
#' }
#' bts <- lapply(c(1, 2), function(s) {
#'   suppressWarnings(backtest(mk(s), cuts = "yearly", horizon = 2,
#'                             window = 3, min_k = 3, seed = 1,
#'                             methods = c("rcma", "ottawa")))
#' })
#' # `reviews_independent` has no default: these two streams were simulated
#' # separately and share no studies, which only the caller can know.
#' pooled_calibration(bts, "shift", R = 200, seed = 1,
#'                    reviews_independent = TRUE)
#'
#' # The same evidence under the other estimand. "cut" above asks what share
#' # of all out-of-date cuts were flagged; "review" asks what share a typical
#' # review has flagged, counting each review once however long its series.
#' pooled_calibration(bts, "shift", R = 200, seed = 1, weighting = "review",
#'                    reviews_independent = TRUE)
#' @export
pooled_calibration <- function(bts, truth = "shift", R = 2000, seed = NULL,
                               conf = 0.95, weighting = c("cut", "review"),
                               reviews_independent,
                               accept_dependence = FALSE) {
  if (!is.list(bts) || !length(bts)) {
    stop("`bts` must be a non-empty list of staleness_backtest objects",
         call. = FALSE)
  }
  for (i in seq_along(bts)) {
    check_class(bts[[i]], "staleness_backtest", paste0("bts[[", i, "]]"),
                "backtest()")
  }
  truth <- match.arg(truth, available_truths())
  weighting <- match.arg(weighting)
  # No default, and TRUE is the only accepted value. Pooling and resampling
  # both assume the reviews share no studies, and nothing here can check it:
  # study_id labels are chosen per stream, so two unrelated reviews numbering
  # their studies 1:14 look identical, and two overlapping ones can use
  # different labels for the same trial. An automatic check would be wrong in
  # both directions, so the promise is asked for instead of inferred.
  if (missing(reviews_independent)) {
    stop("`reviews_independent` must be supplied. Pooling these backtests and ",
         "resampling them both assume the reviews share no studies, and that ",
         "cannot be checked from `study_id` labels chosen separately per ",
         "stream. Pass `reviews_independent = TRUE` to declare it. Two ",
         "backtests of the same evidence, or of two outcomes from one review, ",
         "are not two reviews.", call. = FALSE)
  }
  if (!isTRUE(reviews_independent)) {
    stop("`reviews_independent` must be TRUE. Over reviews that share studies ",
         "there is no defined estimand here: the counts are added as if every ",
         "study appeared once, and the bootstrap draws a review as if drawing ",
         "an independent one. Pool only reviews that share no studies.",
         call. = FALSE)
  }
  check_count(R, "R")
  check_probability(conf, "conf")
  check_seed(seed)
  if (!is.logical(accept_dependence) || length(accept_dependence) != 1L ||
      is.na(accept_dependence)) {
    stop("`accept_dependence` must be TRUE or FALSE", call. = FALSE)
  }

  # The bootstrap resamples whole reviews, which is exactly the right unit
  # WHEN each review is one independent body of studies. If a review was built
  # with allow_dependence = TRUE, its own cuts count one trial several times,
  # and drawing that review again cannot undo it: the interval comes out
  # narrower than the evidence supports, in the one output of this package
  # that looks like an inferential statement. So the point estimates are still
  # returned -- they are descriptive and the caller asked for them -- and the
  # bounds are withheld unless the caller says, in the call, that they want
  # them anyway.
  dep <- vapply(bts, function(b) isTRUE(b$dependent), logical(1))
  dependent_any <- any(dep)
  if (dependent_any && !accept_dependence) {
    warning(sum(dep), " of ", length(bts), " backtests came from a stream ",
            "built with `allow_dependence = TRUE`. Pooled point estimates are ",
            "returned; the bootstrap bounds are not, because resampling ",
            "reviews cannot repair dependence within one and the interval ",
            "would be too narrow. Pass `accept_dependence = TRUE` to get them ",
            "anyway, reading them as optimistic.", call. = FALSE)
  }

  targets <- unique(vapply(bts, function(b) {
    if (is.null(b$truth_target)) "final" else b$truth_target
  }, character(1)))
  if (length(targets) > 1) {
    stop("the backtests use different `truth_target` values (",
         paste(targets, collapse = ", "), "); they answer different ",
         "questions and cannot be pooled", call. = FALSE)
  }

  col <- paste0("truth_", truth)
  # One 2x2 count per review per method, kept separate so that a bootstrap
  # replicate can draw whole reviews.
  counts <- lapply(bts, function(b) {
    res <- eligible_rows(b, truth)
    stats::setNames(lapply(b$methods, function(m) {
      d <- res[res$method == m, ]
      hit <- d$verdict == "out_of_date"; ev <- d[[col]]
      c(tp = sum(hit & ev), fn = sum(!hit & ev),
        tn = sum(!hit & !ev), fp = sum(hit & !ev), n = nrow(d))
    }), b$methods)
  })
  methods <- unique(unlist(lapply(bts, function(b) b$methods)))

  # "cut": add the 2x2 counts over the drawn reviews and divide once, so a
  # thirty-cut review contributes six times a five-cut one. "review": take
  # each review's own rate and average the ones that are defined, so every
  # review counts once. Different questions; see the estimand section above.
  rate <- function(idx, m, num, den) {
    if (weighting == "cut") {
      tot <- c(tp = 0, fn = 0, tn = 0, fp = 0)
      for (i in idx) {
        cm <- counts[[i]][[m]]
        if (!is.null(cm)) tot <- tot + cm[c("tp", "fn", "tn", "fp")]
      }
      d <- sum(tot[den])
      return(if (d > 0) sum(tot[num]) / d else NA_real_)
    }
    # A review whose denominator is zero has no rate to average -- not a rate
    # of zero. Dropping it is what makes this the mean over the reviews that
    # could answer, rather than a mean pulled towards zero by the ones that
    # could not.
    per <- vapply(idx, function(i) {
      cm <- counts[[i]][[m]]
      if (is.null(cm)) return(NA_real_)
      d <- sum(cm[den])
      if (d > 0) sum(cm[num]) / d else NA_real_
    }, numeric(1))
    if (any(!is.na(per))) mean(per, na.rm = TRUE) else NA_real_
  }

  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  boot_ci <- with_preserved_seed(seed = seed, {
    lapply(methods, function(m) {
      reps <- vapply(seq_len(R), function(r) {
        idx <- sample(seq_along(bts), replace = TRUE)
        c(rate(idx, m, "tp", c("tp", "fn")), rate(idx, m, "tn", c("tn", "fp")))
      }, numeric(2))
      list(sens = stats::quantile(reps[1, ], probs, na.rm = TRUE),
           spec = stats::quantile(reps[2, ], probs, na.rm = TRUE))
    })
  })
  names(boot_ci) <- methods

  all_idx <- seq_along(bts)
  out <- lapply(methods, function(m) {
    contributing <- vapply(all_idx, function(i) {
      cm <- counts[[i]][[m]]
      !is.null(cm) && cm[["n"]] > 0
    }, logical(1))
    # An interval drawn from one review describes that review, not a class of
    # them, so it is withheld rather than printed narrow.
    enough <- sum(contributing) >= 2 && (!dependent_any || accept_dependence)
    ci <- boot_ci[[m]]
    data.frame(
      method      = m,
      truth       = truth,
      # Named in the table, not only in the call, so a saved data frame still
      # says which of the two questions it answers.
      weighting   = weighting,
      n_reviews   = sum(contributing),
      n_cuts      = sum(vapply(all_idx, function(i) {
                      cm <- counts[[i]][[m]]; if (is.null(cm)) 0L else cm[["n"]]
                    }, numeric(1))),
      sensitivity = rate(all_idx, m, "tp", c("tp", "fn")),
      sens_lo     = if (enough) unname(ci$sens[1]) else NA_real_,
      sens_hi     = if (enough) unname(ci$sens[2]) else NA_real_,
      specificity = rate(all_idx, m, "tn", c("tn", "fp")),
      spec_lo     = if (enough) unname(ci$spec[1]) else NA_real_,
      spec_hi     = if (enough) unname(ci$spec[2]) else NA_real_,
      contaminated = any(CONTAMINATED_PAIRS$method == m &
                         CONTAMINATED_PAIRS$truth  == truth),
      # TRUE when at least one pooled review allowed dependent estimates. With
      # accept_dependence = FALSE the bounds beside it are NA; with TRUE they
      # are present and this column is what says how to read them.
      dependent    = dependent_any,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
