# Argument checks shared by the five detectors.
#
# Each helper raises with a single message that names what the argument must
# be, rather than reporting the first condition that happened to fail. The
# reader of an error wants the contract, not the branch.
#
# These guard the CALL. A datum that is present but unusable -- an NA sample
# size, a non-finite p-value -- is a fact about the evidence and belongs in a
# "not_applicable" verdict with a reason, which is why the sample-size check
# below deliberately lets NA and Inf through to the detector.

check_scalar_numeric <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L) {
    stop("`", arg, "` must be a single number", call. = FALSE)
  }
  invisible(NULL)
}

check_probability <- function(x, arg, closed = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    if (closed) x >= 0 && x <= 1 else x > 0 && x < 1
  if (!ok) {
    stop("`", arg, "` must be a single number ",
         if (closed) "between 0 and 1 inclusive" else "strictly between 0 and 1",
         call. = FALSE)
  }
  invisible(NULL)
}

check_count <- function(x, arg) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    x >= 1 && x == round(x)
  if (!ok) {
    stop("`", arg, "` must be a single whole number of at least 1",
         call. = FALSE)
  }
  invisible(NULL)
}

# `NULL` is a valid seed and means "run from wherever the stream is", which is
# documented behaviour. Anything else must be a whole number: set.seed()
# truncates towards zero, so seed = 1.5 used to be accepted in silence and
# produce the identical stream to seed = 1 -- two values a reader would record
# as different runs, giving the same result.
check_seed <- function(x, arg = "seed") {
  if (is.null(x)) return(invisible(NULL))
  # The bound is R's integer range, and its lower end is -.Machine$integer.max,
  # not the int32 minimum: R reserves -2147483648 for NA_integer_, so
  # set.seed() rejects it. Checked by running both ends, not read off the type.
  lim <- .Machine$integer.max
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x) ||
      abs(x) > lim) {
    stop("`", arg, "` must be NULL or a single whole number in [-", lim, ", ",
         lim, "]; `set.seed()` truncates anything else, so two different ",
         "values can silently give the same stream, and anything past that ",
         "range it rejects outright", call. = FALSE)
  }
  invisible(NULL)
}

# The new-evidence object carries two things that must agree: `k`, which
# decides whether there is anything to assess at all, and `yi`/`vi`, which are
# what actually gets fitted. Nothing tied them together, so the object could
# say one thing and carry another -- and the failure was silent in both
# directions. k = 1 with no studies produced a "current" verdict from a model
# refitted on the prior evidence alone, walking past the guard written to stop
# exactly that; k = 0 with real studies threw them away.
#
# window_between() builds `k` as sum(keep) over the same subset it takes yi
# and vi from, so this contract is what the canonical source already produces.
check_new_evidence <- function(new, arg = "new") {
  if (!is.list(new)) {
    stop("`", arg, "` must be a list with `yi`, `vi` and `k`, as returned by ",
         "`window_between()`", call. = FALSE)
  }
  if (!is.numeric(new$yi) || !is.numeric(new$vi)) {
    stop("`", arg, "$yi` and `", arg, "$vi` must both be numeric vectors",
         call. = FALSE)
  }
  if (length(new$yi) != length(new$vi)) {
    stop("`", arg, "$yi` and `", arg, "$vi` must be the same length; got ",
         length(new$yi), " and ", length(new$vi), call. = FALSE)
  }
  if (!is.numeric(new$k) || length(new$k) != 1L || !is.finite(new$k) ||
      new$k < 0 || new$k != round(new$k)) {
    stop("`", arg, "$k` must be a single whole number of 0 or more",
         call. = FALSE)
  }
  if (new$k != length(new$yi)) {
    stop("`", arg, "$k` says ", new$k, " but `", arg, "$yi` holds ",
         length(new$yi), ". `k` decides whether there is new evidence and ",
         "`yi`/`vi` are what gets fitted, so a disagreement between them is ",
         "answered from one and reported from the other", call. = FALSE)
  }
  invisible(NULL)
}

check_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive number", call. = FALSE)
  }
  invisible(NULL)
}
