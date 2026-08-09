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

check_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive number", call. = FALSE)
  }
  invisible(NULL)
}
