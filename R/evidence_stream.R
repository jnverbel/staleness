#' Build a time-ordered evidence stream from a meta-analysis
#'
#' The only part of the package that knows about time. Detectors receive
#' snapshots, never dates, which is what allows [backtest()] to reuse them
#' unmodified.
#'
#' @param ma An `rma.uni` object from `metafor`.
#' @param date Numeric or Date vector, one entry per study, same order as the
#'   data used to fit `ma`. Missing values are an error; they are never imputed.
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
#' stream <- evidence_stream(ma, date = bcg$year, ni = bcg$ni)
#' stream
#' @export
evidence_stream <- function(ma, date, ni = NULL) {
  if (!inherits(ma, "rma.uni")) {
    stop("`ma` must be an rma.uni object from metafor", call. = FALSE)
  }
  k <- ma$k
  if (length(date) != k) {
    stop("`date` has length ", length(date), " but the meta-analysis has ",
         k, " studies", call. = FALSE)
  }
  if (anyNA(date)) {
    stop("`date` has missing values; dates are never imputed", call. = FALSE)
  }
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

  ord <- order(date, seq_along(date))  # stable: ties keep input order
  structure(
    list(
      yi      = as.numeric(ma$yi)[ord],
      vi      = as.numeric(ma$vi)[ord],
      ni      = if (is.null(ni)) NULL else as.numeric(ni)[ord],
      date    = date[ord],
      measure = ma$measure,
      method  = ma$method,
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
#'                           date = bcg$year)
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
  keep <- stream$date <= cut
  if (sum(keep) < 2) {
    stop("a snapshot needs at least 2 studies; found ", sum(keep),
         " at cut ", format(cut), call. = FALSE)
  }
  metafor::rma(yi = stream$yi[keep], vi = stream$vi[keep],
               measure = stream$measure, method = stream$method)
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
#'                           date = bcg$year, ni = bcg$ni)
#'
#' # The interval is half-open, (from, to], so a trial published exactly in
#' # 1970 would belong to the snapshot, never to the new evidence.
#' new <- window_between(stream, 1970, 1980)
#' new$k
#' sum(new$ni)
#' @export
window_between <- function(stream, from, to) {
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
  invisible(x)
}
