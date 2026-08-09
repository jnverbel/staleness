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
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != round(x)) {
    stop("`", arg, "` must be NULL or a single whole number; `set.seed()` ",
         "truncates anything else, so two different values can silently give ",
         "the same stream", call. = FALSE)
  }
  invisible(NULL)
}

check_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive number", call. = FALSE)
  }
  invisible(NULL)
}
