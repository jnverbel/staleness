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

#' Calibration of the detectors against a truth definition
#'
#' Turns a backtest's raw results into sensitivity, specificity and false
#' alarm rate, one row per method, scored against a chosen definition of
#' ground truth (see [truth]).
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
#' @param bt A `staleness_backtest`, see [backtest()].
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`, see
#'   [available_truths()].
#' @return A data frame, one row per method in `bt$methods`, with columns
#'   `method`, `truth`, `sensitivity`, `specificity`, `false_alarm`, `n` and
#'   `contaminated`.
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
#'                                date = bcg$year, ni = bcg$ni))
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
#' The same three exclusions as [calibration()] apply, through the same
#' internal filter: censored cuts, `not_applicable` verdicts, and rows whose
#' truth column is `NA` are dropped before either the events or the firings
#' are located. Both sides of the comparison are read from that already
#' filtered frame, so a degenerate standard error removes an event and the
#' verdict that would have been scored against it together, rather than
#' leaving one without the other.
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
#' @param bt A `staleness_backtest`, see [backtest()].
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`, see
#'   [available_truths()].
#' @return A data frame, one row per method in `bt$methods`, with columns
#'   `method`, `median_lead` and `n_events`.
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
#'                                date = bcg$year, ni = bcg$ni))
#'
#' # Whether a detector eventually fires is the easy question. This is the
#' # one the methods papers never answer: how far ahead of the evidence.
#' # `median_lead` is NA when a detector caught no true event at all, and
#' # `n_events` still says how many it had the chance to catch.
#' lead_time(bt)
#' @export
lead_time <- function(bt, truth = "shift") {
  truth <- match.arg(truth, available_truths())
  col <- paste0("truth_", truth)
  res <- eligible_rows(bt, truth)

  out <- lapply(bt$methods, function(m) {
    d <- res[res$method == m, ]
    d <- d[order(d$cut), ]
    events <- d$cut[d[[col]]]
    fired  <- d$cut[d$verdict == "out_of_date"]

    # One lead per true event: the gap to that event's earliest on-time
    # firing, or NA when the detector never fired at or before it. A miss
    # must never contribute a numeric lead (e.g. 0), or it would be
    # indistinguishable from genuine same-period detection.
    leads <- vapply(events, function(event_cut) {
      on_time <- fired[fired <= event_cut]
      if (length(on_time)) event_cut - min(on_time) else NA_real_
    }, numeric(1))

    median_lead <- if (any(!is.na(leads))) {
      stats::median(leads, na.rm = TRUE)
    } else {
      NA_real_
    }

    data.frame(method = m, median_lead = median_lead, n_events = length(events),
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
#'                                date = bcg$year))
#'
#' # All three truth definitions at once, stacked. Reading them side by side
#' # is the point: a detector that looks good under one and bad under another
#' # is telling you about the truth definition, not only about itself.
#' summary(bt)
#' @export
summary.staleness_backtest <- function(object, ...) {
  do.call(rbind, lapply(available_truths(), function(t) calibration(object, t)))
}
