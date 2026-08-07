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
  if (is.null(ni)) ni <- ma$ni
  if (!is.null(ni) && length(ni) != k) {
    stop("`ni` has length ", length(ni), " but the meta-analysis has ", k,
         " studies", call. = FALSE)
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
