#' Backtest the detectors against historical evidence
#'
#' Walks the evidence stream, reconstructs what was known at each point in the
#' past, asks every detector the same question it would be asked today, and
#' compares the answer with what actually happened afterwards.
#'
#' At each cut `t`, the snapshot (`prev`, from [snapshot_at()]) and the new
#' evidence window (`new`, from [window_between()]) are built using only
#' studies with `date <= t` and `date` in `(t, t + window]` respectively.
#' Nothing past the cut, including tau-squared, feeds into `prev`: it is
#' refit from scratch on the restricted data by [snapshot_at()], never
#' derived from the full series. The `final` meta-analysis, fit once over the
#' complete stream, exists solely to score truth after the fact and is never
#' passed to a detector.
#'
#' When the stream carries sample sizes (`stream$ni`), they are summed on
#' both sides of the cut and passed to [check_currency()] as `n_prev`/`n_new`
#' so that [barrowman()] can be evaluated like every other detector. When the
#' stream has no `ni`, both are passed as `NULL` and `barrowman()` degrades to
#' `not_applicable` through [check_currency()]'s own handling, rather than
#' erroring.
#'
#' Cuts too close to the end of the stream cannot be fairly evaluated: they
#' are kept in `results` (marked `censored = TRUE`) so that a caller can see
#' them, but must be excluded from any accuracy metric computed later. "Too
#' close" means within `max(horizon, window)` of the last study date, not just
#' within `horizon`. Both bounds have to hold, and for different reasons: a
#' cut needs `horizon` units of future for its truth to have had time to
#' materialise, *and* it needs `window` units of future for the detector to be
#' handed the same amount of new evidence every other uncensored cut was
#' handed. Censoring on `horizon` alone would let a run with
#' `window > horizon` score late cuts as full observations while their
#' evidence window ran off the end of the data — a truncated window shows a
#' detector less change than really happened, biasing exactly those cuts
#' toward `current`, invisibly. A backtest needs at least 3 uncensored cuts to
#' say anything useful about a detector's calibration and is refused
#' otherwise.
#'
#' @param stream A `staleness_stream` from [evidence_stream()].
#' @param cuts `"yearly"`, or a numeric vector of explicit cut points. Sorted
#'   and de-duplicated: a repeated cut is the same point in time, and counting
#'   it twice would inflate the denominator of every metric.
#' @param methods Detector names, see [available_methods()].
#' @param horizon Units of future required for a cut's truth to be evaluated.
#' @param window Length of the window of new evidence assessed at each cut.
#'   Cuts with fewer than `max(horizon, window)` units of future are marked
#'   censored and excluded from metrics — see the note on censoring above.
#' @param min_k Minimum studies required in a snapshot. At least 2, since a
#'   meta-analysis cannot be fitted from fewer.
#' @param seed Integer seed for [simulation()].
#' @return An object of class `staleness_backtest` with elements `results` (a
#'   data.frame with columns `cut`, `method`, `verdict`, `signal`, `reason`,
#'   `truth_shift`, `truth_surprise`, `truth_conclusion`, `censored`),
#'   `stream`, `methods`, `horizon`, `window`, `n_cuts` and `n_censored`. The
#'   `truth_shift` and `truth_surprise` columns can be `NA` when the standard
#'   error they divide by (`se_final` or the snapshot's own `se`) is
#'   degenerate; that `NA` means "truth could not be determined here", not
#'   "no shift occurred", and is deliberately left unresolved rather than
#'   guessed at. Consumers must exclude those rows from any accuracy metric,
#'   the same treatment `censored == TRUE` rows get — never count an `NA`
#'   truth as a miss. [calibration()] and [lead_time()] already do this
#'   filtering for you.
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
#'                           date = bcg$year, ni = bcg$ni)
#'
#' bt <- backtest(stream, cuts = "yearly")
#' bt
#'
#' # Every cut refits from scratch, so no detector ever sees the future.
#' # Cuts without `max(horizon, window)` units of future left are censored
#' # rather than scored against a truth that cannot be known yet.
#' c(cuts = bt$n_cuts, censored = bt$n_censored)
#'
#' # One row per cut per method. `not_applicable` is recorded, not dropped.
#' head(bt$results[, c("cut", "method", "verdict", "censored")])
#' @export
backtest <- function(stream, cuts = "yearly", methods = available_methods(),
                     horizon = 5, window = 3, min_k = 3, seed = NULL) {
  if (!inherits(stream, "staleness_stream")) {
    stop("`stream` must come from evidence_stream()", call. = FALSE)
  }
  # Checked here rather than left to surface as an internal R error further
  # down. A caller who passes cuts containing NA used to get "missing value
  # where TRUE/FALSE needed", which says nothing about what they did wrong.
  if (!is.numeric(horizon) || length(horizon) != 1L || !is.finite(horizon) ||
      horizon < 0) {
    stop("`horizon` must be a single, finite, non-negative number; a negative ",
         "horizon would ask the truth to be evaluated before the cut",
         call. = FALSE)
  }
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) ||
      window <= 0) {
    stop("`window` must be a single, finite, positive number: a cut with no ",
         "window has no new evidence to assess", call. = FALSE)
  }
  if (!is.numeric(min_k) || length(min_k) != 1L || !is.finite(min_k) ||
      min_k < 2) {
    stop("`min_k` must be at least 2: a snapshot cannot be fitted from fewer ",
         "than two studies", call. = FALSE)
  }
  if (!length(methods)) {
    stop("`methods` is empty; name at least one of: ",
         paste(available_methods(), collapse = ", "), call. = FALSE)
  }

  dates <- as.numeric(stream$date)
  if (identical(cuts, "yearly")) {
    cuts <- seq(min(dates), max(dates), by = 1)
  }
  if (!is.numeric(cuts) || !length(cuts) || anyNA(cuts)) {
    stop("`cuts` must be \"yearly\" or a numeric vector with no missing values",
         call. = FALSE)
  }
  # A repeated cut is the same point in time. Left in, it entered the results
  # once per repetition and inflated the denominator of every metric.
  cuts <- sort(unique(cuts))

  # a cut is usable when the snapshot has enough studies AND new evidence exists
  usable <- vapply(cuts, function(cut) {
    sum(dates <= cut) >= min_k && sum(dates > cut & dates <= cut + window) >= 1
  }, logical(1))
  cuts <- cuts[usable]

  # A cut is uncensored only if BOTH its truth horizon and its evidence window
  # fit inside the data. Censoring on `horizon` alone silently truncates the
  # new-evidence window of late cuts whenever `window > horizon`, then scores
  # them as if they had seen a full window.
  last_uncensored <- max(dates) - max(horizon, window)
  if (sum(cuts <= last_uncensored) < 3) {
    stop("backtesting needs at least 3 uncensored cuts; found ",
         sum(cuts <= last_uncensored), call. = FALSE)
  }

  # final state, used only to evaluate truth, never fed back into a snapshot
  final <- metafor::rma(yi = stream$yi, vi = stream$vi,
                        measure = stream$measure, method = stream$method)
  theta_final <- as.numeric(final$beta)
  se_final    <- final$se
  p_final     <- final$pval

  rows <- list()
  for (cut in cuts) {
    prev <- snapshot_at(stream, cut)
    new  <- window_between(stream, cut, cut + window)

    # Correction to the brief: thread sample sizes through so barrowman() is
    # actually exercised when the stream has them, instead of being blind at
    # every cut. sum(new$ni) mirrors sum(stream$ni[...]) below: both are NULL
    # together, both are numeric together, because window_between() carries
    # `ni` (or its absence) forward from the same stream.
    if (is.null(stream$ni)) {
      n_prev <- NULL
      n_new  <- NULL
    } else {
      n_prev <- sum(stream$ni[dates <= cut])
      n_new  <- sum(new$ni)
    }

    chk <- check_currency(prev, new, methods = methods,
                          n_prev = n_prev, n_new = n_new, seed = seed)
    df  <- as.data.frame(chk)

    theta_t <- as.numeric(prev$beta)
    df$cut      <- cut
    df$censored <- cut > last_uncensored
    # truth_shift/truth_surprise divide by a standard error; if that SE is
    # ever exactly zero (a degenerate snapshot or final fit) the result is
    # NA, not TRUE/FALSE. That NA is left as-is rather than coerced: it means
    # "truth could not be determined here", which is a different fact from
    # "no shift occurred", and collapsing the two would quietly bias any
    # calibration metric computed downstream. Consumers of `results` should
    # treat an NA truth column the same way they treat `censored == TRUE`:
    # excluded from accuracy metrics, not counted as a miss.
    df$truth_shift      <- truth_shift(theta_t, theta_final, se_final)
    df$truth_surprise   <- truth_surprise(theta_t, prev$se, theta_final)
    df$truth_conclusion <- truth_conclusion(theta_t, prev$pval,
                                            theta_final, p_final)
    rows[[length(rows) + 1L]] <- df
  }

  results <- do.call(rbind, rows)
  results <- results[, c("cut", "method", "verdict", "signal", "reason",
                         "truth_shift", "truth_surprise", "truth_conclusion",
                         "censored")]

  structure(
    list(results = results, stream = stream, methods = methods,
         horizon = horizon, window = window,
         n_cuts = length(cuts),
         n_censored = sum(results$censored) / length(methods)),
    class = "staleness_backtest"
  )
}

#' @export
print.staleness_backtest <- function(x, ...) {
  cat("<staleness_backtest>\n")
  cat("  cuts     :", x$n_cuts, "(", x$n_censored, "censored )\n")
  cat("  methods  :", paste(x$methods, collapse = ", "), "\n")
  cat("  horizon  :", x$horizon, " window:", x$window, "\n")
  invisible(x)
}
