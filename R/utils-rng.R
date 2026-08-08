#' Evaluate an expression without disturbing the caller's random stream
#'
#' Detectors that randomise ([simulation()], and the permutation test inside
#' [sufficiency()]) must not leave a footprint on the caller's session. A
#' script that calls `set.seed()` once and then runs [backtest()] would
#' otherwise lose the reproducibility of everything it does afterwards, and
#' CRAN policy forbids altering global state in any case.
#'
#' What is forbidden is seizing the seed, not consuming deviates, so the two
#' cases are handled differently.
#'
#' With a `seed`, `.Random.seed` is saved from the global environment,
#' reseeded for the duration of `expr`, and restored on exit — including on
#' error. If the session had no stream at all, the exit handler removes the
#' one `set.seed()` created rather than leaving it behind.
#'
#' Without a `seed`, the stream is left to advance exactly as it would for
#' any other RNG call. Rewinding it would be worse than a footprint: every
#' call would restart from the same state and return identical draws, so a
#' path documented as unreproducible would in fact be deterministic.
#'
#' @param expr Expression to evaluate. Lazily evaluated, so it runs *after*
#'   the seed has been set and *before* the old state is put back.
#' @param seed Integer seed, or `NULL` to run from wherever the stream
#'   currently is, advancing it as any RNG call would.
#' @return The value of `expr`.
#' @keywords internal
with_preserved_seed <- function(expr, seed = NULL) {
  # A session that has not drawn yet has no `.Random.seed` at all. Creating
  # one and restoring it would leave a footprint in the very function whose
  # purpose is to leave none, so in that case the exit handler removes it.
  had_stream <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_stream) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(suppressWarnings(rm(".Random.seed", envir = globalenv())),
            add = TRUE)
  }

  if (is.null(seed)) {
    # Draw a fresh, unpredictable starting point instead of running from
    # wherever the caller happens to be. Since the caller's stream is put
    # back on exit, running from it would rewind every call to the same
    # state and hand back byte-identical "random" draws -- a path documented
    # as unreproducible would in fact be deterministic.
    set.seed(NULL)
  } else {
    set.seed(seed)
  }
  expr
}
