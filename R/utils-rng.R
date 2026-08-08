#' Evaluate an expression without disturbing the caller's random stream
#'
#' Detectors that randomise ([simulation()], and the permutation test inside
#' [sufficiency()]) must not leave a footprint on the caller's session. A
#' script that calls `set.seed()` once and then runs [backtest()] would
#' otherwise lose the reproducibility of everything it does afterwards, and
#' CRAN policy forbids altering global state in any case.
#'
#' `.Random.seed` lives in the global environment, so it is saved from there,
#' optionally reseeded for the duration of `expr`, and restored on exit —
#' including on error. When no stream exists yet, one is initialised with
#' `set.seed(NULL)` first; that consumes no draws, and leaves the caller in
#' exactly the state it would have reached by generating its own first random
#' number.
#'
#' @param expr Expression to evaluate. Lazily evaluated, so it runs *after*
#'   the seed has been set and *before* the old state is put back.
#' @param seed Integer seed, or `NULL` to run from wherever the stream
#'   currently is (unreproducible, but still non-destructive).
#' @return The value of `expr`.
#' @keywords internal
with_preserved_seed <- function(expr, seed = NULL) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    set.seed(NULL)
  }
  old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  expr
}
