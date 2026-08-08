#' Check whether a meta-analysis is still current
#'
#' Applies the selected detectors to a prior meta-analysis and the evidence
#' published since.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new A list with `yi`, `vi` and `k`, as returned by [window_between()].
#' @param methods Character vector of detector names, see [available_methods()].
#' @param n_prev,n_new Sample sizes, required by [barrowman()].
#' @param qualitative Character vector of qualitative signals, see [ottawa()].
#' @param seed Integer seed for [simulation()].
#' @return An object of class `staleness_check`.
#' @export
check_currency <- function(prev, new, methods = available_methods(),
                           n_prev = NULL, n_new = NULL,
                           qualitative = character(), seed = NULL) {
  unknown <- setdiff(methods, available_methods())
  if (length(unknown)) {
    stop("unknown method: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (is.null(new$k) || new$k < 1) {
    # Absence of new evidence is not evidence of currency. It gets its own
    # class so that no caller can mistake it for a verified "current".
    verdicts <- lapply(methods, function(m) {
      verdict_na(m, "no new evidence supplied; currency cannot be assessed")
    })
    names(verdicts) <- methods
    return(structure(
      list(prev = prev, updated = prev, verdicts = verdicts,
           disagreement = FALSE, i2 = prev$I2),
      class = c("staleness_no_evidence", "staleness_check")
    ))
  }

  updated <- metafor::rma(
    yi      = c(as.numeric(prev$yi), new$yi),
    vi      = c(as.numeric(prev$vi), new$vi),
    measure = prev$measure,
    method  = prev$method
  )

  verdicts <- lapply(methods, function(m) {
    switch(m,
      rcma        = rcma(prev, updated),
      ottawa      = ottawa(prev, updated, qualitative = qualitative),
      sufficiency = sufficiency(prev, updated),
      simulation  = simulation(prev, new, seed = seed),
      barrowman   = if (is.null(n_prev) || is.null(n_new)) {
                      verdict_na("barrowman",
                        "sample size not supplied for prior and new studies")
                    } else {
                      barrowman(prev, n_prev = n_prev, n_new = n_new)
                    }
    )
  })
  names(verdicts) <- methods

  calls <- vapply(verdicts, function(v) v$verdict, character(1))
  decided <- calls[calls != "not_applicable"]

  structure(
    list(prev = prev, updated = updated, verdicts = verdicts,
         disagreement = length(unique(decided)) > 1,
         i2 = updated$I2),
    class = "staleness_check"
  )
}

#' @export
as.data.frame.staleness_check <- function(x, ...) {
  data.frame(
    method  = vapply(x$verdicts, function(v) v$method,  character(1)),
    verdict = vapply(x$verdicts, function(v) v$verdict, character(1)),
    signal  = vapply(x$verdicts, function(v) as.numeric(v$signal), numeric(1)),
    reason  = vapply(x$verdicts, function(v) v$reason,  character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @export
print.staleness_check <- function(x, ...) {
  cat("<staleness_check>\n")
  cat("  studies:", x$prev$k, "prior +", x$updated$k - x$prev$k, "new\n")
  cat("  I2 (updated):", format(x$i2, digits = 3), "%\n")
  if (!is.na(x$i2) && x$i2 > 75) {
    cat("  note: heterogeneity above 75%; the pooled estimate is itself debatable\n")
  }
  cat("\n")
  for (v in x$verdicts) print(v)
  if (x$disagreement) cat("\n  detectors disagree\n")
  invisible(x)
}
