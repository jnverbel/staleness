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
#' @param bt A `staleness_backtest`, see [backtest()].
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`, see
#'   [available_truths()].
#' @return A data frame, one row per method, with columns `method`, `truth`,
#'   `sensitivity`, `specificity`, `false_alarm`, `n` and `contaminated`.
#' @export
calibration <- function(bt, truth = "shift") {
  truth <- match.arg(truth, available_truths())
  col <- paste0("truth_", truth)
  res <- bt$results[!bt$results$censored &
                    bt$results$verdict != "not_applicable" &
                    !is.na(bt$results[[col]]), ]

  out <- lapply(unique(res$method), function(m) {
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
#' the two apart, because it scores each cut in isolation. `lead_time()`
#' reports, per method, how many time units before the first true event (per
#' the chosen `truth` definition) the detector first flagged `out_of_date`.
#'
#' The same three exclusions as [calibration()] apply before the event is
#' located: censored cuts, `not_applicable` verdicts, and rows whose truth
#' column is `NA` are dropped first, so a degenerate standard error can
#' neither manufacture nor hide an event.
#'
#' A method with no true event under `truth` in the uncensored window gets
#' `n_events = 0` and `median_lead = NA`: there is nothing to lead. A method
#' that never fired `out_of_date` at or before its first true event also gets
#' `median_lead = NA`: it never gave advance warning, so there is no lead
#' time to report, only a miss (already visible in `calibration()`).
#'
#' @param bt A `staleness_backtest`, see [backtest()].
#' @param truth One of `"shift"`, `"surprise"`, `"conclusion"`, see
#'   [available_truths()].
#' @return A data frame, one row per method, with columns `method`,
#'   `median_lead` and `n_events`.
#' @export
lead_time <- function(bt, truth = "shift") {
  truth <- match.arg(truth, available_truths())
  col <- paste0("truth_", truth)
  res <- bt$results[!bt$results$censored &
                    bt$results$verdict != "not_applicable" &
                    !is.na(bt$results[[col]]), ]

  out <- lapply(unique(res$method), function(m) {
    d <- res[res$method == m, ]
    d <- d[order(d$cut), ]
    events <- d$cut[d[[col]]]
    if (!length(events)) {
      return(data.frame(method = m, median_lead = NA_real_, n_events = 0L,
                        stringsAsFactors = FALSE))
    }
    first_event <- min(events)
    fired <- d$cut[d$verdict == "out_of_date" & d$cut <= first_event]
    lead <- if (length(fired)) first_event - min(fired) else NA_real_
    data.frame(method = m, median_lead = lead, n_events = length(events),
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
#' @export
summary.staleness_backtest <- function(object, ...) {
  do.call(rbind, lapply(available_truths(), function(t) calibration(object, t)))
}
