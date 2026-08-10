#' What series is actually being evaluated
#'
#' Every rate this package reports is a rate over one body of evidence, and
#' the reader of a calibration table cannot see what that body was. Eight
#' facts decide how the numbers should be read, and each of them has silently
#' changed an answer at some point in this package's history: how many
#' estimates there are, how many distinct studies they came from, how many
#' estimates a study contributed, over what span, on what measure, under what
#' model, whether dependence was allowed, and whether sample sizes are there
#' at all.
#'
#' Nothing here is a new requirement. It reports what the stream already knows,
#' in one row, so that a saved result can be audited later against the series
#' it came from.
#'
#' @section Why these eight and not others:
#' * `k` against `n_studies` is the dependence question in its raw form. Equal
#'   means one estimate per study; unequal means some trial is being counted
#'   more than once, and every rate is optimistic by an amount that
#'   `max_per_study` bounds.
#' * `from` and `to` say what "yearly cuts" was actually cutting. A stream
#'   spanning six years supports a different claim from one spanning forty.
#' * `measure` decides which detectors can answer at all: [ottawa()]'s effect
#'   criterion is defined on comparative ratio measures and returns
#'   `not_applicable` elsewhere, and [rcma()] is undefined for difference
#'   measures with a prior effect near zero.
#' * `model`, `test`, `weighted` and `tau2_fixed` are the fit every snapshot is
#'   refitted under. Fixed effects and REML disagree on the same evidence by
#'   factors that change verdicts, and the Knapp-Hartung adjustment moves
#'   p-values by orders of magnitude, which [ottawa()] reads directly.
#' * `ni_available` decides whether [barrowman()] can answer. Without sample
#'   sizes it returns `not_applicable` at every cut, which looks in a results
#'   table exactly like a detector that never fired.
#'
#' @param x A `staleness_stream` or a `staleness_backtest`.
#' @return A one-row data frame with columns `k`, `n_studies`,
#'   `max_per_study`, `from`, `to`, `measure`, `model`, `test`, `weighted`,
#'   `tau2_fixed`, `dependent` and `ni_available`.
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
#'                           date = bcg$year,
#'                           study_id = seq_along(bcg$year), ni = bcg$ni)
#' eligibility(stream)
#'
#' # A backtest answers for the stream it ran on, so a saved result can be
#' # audited against the series it came from without keeping the stream.
#' eligibility(backtest(stream))
#' @export
eligibility <- function(x) {
  stream <- if (inherits(x, "staleness_backtest")) x$stream else x
  check_class(stream, "staleness_stream", "x",
              "evidence_stream() (or backtest())")

  per <- table(stream$study_id)
  data.frame(
    k             = stream$k,
    n_studies     = length(per),
    max_per_study = max(as.integer(per)),
    from          = min(stream$date),
    to            = max(stream$date),
    measure       = if (is.null(stream$measure)) NA_character_
                    else as.character(stream$measure),
    model         = if (is.null(stream$method)) NA_character_
                    else as.character(stream$method),
    test          = stream$test,
    weighted      = isTRUE(stream$weighted),
    # NA rather than FALSE: the model either fixes tau2 at a value the caller
    # chose, in which case that value is the fact worth reporting, or it
    # estimates it, in which case there is no number here at all.
    tau2_fixed    = if (is.null(stream$tau2_fix)) NA_real_
                    else as.numeric(stream$tau2_fix),
    dependent     = isTRUE(stream$dependent),
    # Counted, not merely present. metafor can derive an `ni` that is NA for
    # some studies, and barrowman() sums over a snapshot, so "sample sizes are
    # available" is a matter of how many.
    ni_available  = if (is.null(stream$ni)) 0L else sum(is.finite(stream$ni)),
    stringsAsFactors = FALSE
  )
}
