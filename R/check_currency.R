#' Check whether a meta-analysis is still current
#'
#' Applies the selected detectors to a prior meta-analysis and the evidence
#' published since.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new A list with `yi`, `vi` and `k`, as returned by
#'   [window_between()]. The three must agree: `yi` and `vi` of the same
#'   length, and `k` equal to that length. `k` decides whether there is new
#'   evidence at all and `yi`/`vi` are what the updated model is fitted on, so
#'   a mismatch is refused rather than resolved in favour of either. The values
#'   must also be usable: `yi` finite, `vi` finite and strictly positive. A
#'   missing or infinite variance is dropped or ignored by the model fit, which
#'   returns the prior estimate unchanged and so a verdict of `"current"`
#'   computed from evidence that carried nothing.
#' @param methods Character vector of detector names, see [available_methods()].
#' @param n_prev,n_new Sample sizes, required by [barrowman()].
#' @param qualitative Character vector of qualitative signals, see [ottawa()].
#' @param seed Integer seed for [simulation()].
#' @return An object of class `staleness_check`.
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
#' prev <- snapshot_at(stream, 1970)
#' new  <- window_between(stream, 1970, 1980)
#'
#' # `new` is the new evidence on its own, never the updated model: this
#' # function refits it internally. Handing it a fitted rma is an error, not
#' # a quietly doubled dataset.
#' check_currency(prev, new, methods = c("rcma", "ottawa", "sufficiency"))
#'
#' # barrowman() needs the participant counts, and says so when they are absent.
#' check_currency(prev, new, methods = "barrowman")
#' check_currency(prev, new, methods = "barrowman",
#'                n_prev = sum(bcg$ni[bcg$year <= 1970]), n_new = sum(new$ni))
#' @export
check_currency <- function(prev, new, methods = available_methods(),
                           n_prev = NULL, n_new = NULL,
                           qualitative = character(), seed = NULL) {
  # backtest() has raised on an empty `methods` since its argument sweep;
  # this function quietly built a staleness_check holding zero verdicts -- an
  # object that answers no question, for the same mistake, in the same package.
  if (!length(methods)) {
    stop("`methods` is empty; name at least one of: ",
         paste(available_methods(), collapse = ", "), call. = FALSE)
  }
  unknown <- setdiff(methods, available_methods())
  if (length(unknown)) {
    stop("unknown method: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  # Naming a detector twice is naming it once: `verdicts` is keyed by name,
  # and duplicates would produce duplicate rows downstream.
  methods <- unique(methods)
  # `new` is the new evidence alone, but rcma() and ottawa() take the updated
  # model, so passing the updated rma here is the natural mistake. An rma.uni
  # carries $yi, $vi and $k, which is everything the guards below look for, so
  # duck typing would let it through and pool the prior studies twice -- a
  # quietly wrong answer rather than an error.
  if (inherits(new, "rma")) {
    stop("`new` must be the new evidence on its own -- a list with `yi`, `vi` ",
         "and `k`, as returned by `window_between()` -- not a fitted ",
         "meta-analysis. Passing the updated model here would count every ",
         "prior study twice.", call. = FALSE)
  }
  check_rma_uni(prev, "prev")
  check_new_evidence(new)
  check_seed(seed)
  if (new$k < 1) {
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

  # Refitted with the same options as `prev`, not just the same method. A
  # caller who fitted with test = "knha" must not have the updated model --
  # the one ottawa() reads p-values off -- scored under the default z test.
  updated <- metafor::rma(
    yi       = c(as.numeric(prev$yi), new$yi),
    vi       = c(as.numeric(prev$vi), new$vi),
    measure  = prev$measure,
    method   = prev$method,
    test     = if (is.null(prev$test)) "z" else prev$test,
    weighted = if (is.null(prev$weighted)) TRUE else isTRUE(prev$weighted),
    tau2     = if (isTRUE(prev$tau2.fix)) prev$tau2 else NULL
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
