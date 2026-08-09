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
  # Shape is not usability. An NA or infinite variance makes metafor drop the
  # study or ignore it, so the "updated" model comes back identical to the
  # prior one, every ratio is 1, and the answer is a confident "current" from
  # evidence that carried nothing -- the same failure the k/yi mismatch above
  # produced, reached by a different route. Zero-length vectors are vacuously
  # finite and positive, so the empty case still passes.
  if (!all(is.finite(new$yi))) {
    stop("`", arg, "$yi` must be finite; a missing or infinite effect is ",
         "dropped by the model fit and returns a verdict computed without it",
         call. = FALSE)
  }
  if (!all(is.finite(new$vi)) || any(new$vi <= 0)) {
    stop("`", arg, "$vi` must be finite and strictly positive; a missing, ",
         "infinite or non-positive variance gives the study no weight, or ",
         "none that the model can use", call. = FALSE)
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

# Both detectors read $measure off `prev` to choose their branch, and neither
# checked that the updated model was on the same scale. A risk ratio against a
# mean difference returned a signal of 2.98 and a verdict, which is a number
# with no meaning wearing the same clothes as every other verdict.
check_same_measure <- function(prev, new_ma) {
  a <- if (is.null(prev$measure)) NA_character_ else as.character(prev$measure)
  b <- if (is.null(new_ma$measure)) NA_character_ else as.character(new_ma$measure)
  if (!identical(a, b)) {
    stop("`prev` and `new_ma` must be on the same scale; got ",
         if (is.na(a)) "none" else a, " and ", if (is.na(b)) "none" else b,
         ". Comparing effects across measures does not define a ratio",
         call. = FALSE)
  }
  invisible(NULL)
}

# The truth functions each document "a single logical value". They took
# vectors: truth_shift() returned one per element, truth_conclusion() died
# inside R's own coercion. Documentation stronger than code, again.
check_scalar_input <- function(..., arg) {
  vals <- list(...)
  bad <- vapply(vals, function(x) !is.numeric(x) || length(x) != 1L,
                logical(1))
  if (any(bad)) {
    stop("`", arg, "` takes a single number for each argument, not a vector; ",
         "apply it one row at a time", call. = FALSE)
  }
  invisible(NULL)
}

# A cut point is one number on the date scale. `from`, `to` and `cut` went
# straight into a comparison against the whole date vector, so NA made every
# element NA, a character compared lexically, and a length-two vector recycled
# against thirteen dates with only a warning.
check_cut_point <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", arg, "` must be a single finite number on the same scale as the ",
         "stream's dates", call. = FALSE)
  }
  invisible(NULL)
}

# Every exported detector documents an `rma.uni` and none of them checked.
# Four died with R's own "argument is of length zero"; sufficiency() was worse
# and returned a verdict of not_applicable from an empty list.
#
# The subclass case is the one that matters in practice. evidence_stream() has
# always refused rma.mh with a reasoned message -- it refits each snapshot with
# rma(), which needs yi and vi, and a Mantel-Haenszel fit cannot be reproduced
# from those -- while the detectors took the same object and answered
# "current". Reproducing a Cochrane review to the digit requires
# Mantel-Haenszel, so rma.mh is exactly what a user arrives with.
check_rma_uni <- function(x, arg) {
  if (inherits(x, "rma.uni")) return(invisible(NULL))
  if (inherits(x, "rma")) {
    stop("`", arg, "` is a ", class(x)[1], " fit; the detectors are defined on ",
         "the pooled effect and its standard error as `rma.uni` reports them. ",
         "Convert the data first -- `escalc()` then `rma()` -- and note that ",
         "the pooled estimates will be inverse-variance, not ", class(x)[1],
         ".", call. = FALSE)
  }
  stop("`", arg, "` must be an rma.uni object from metafor; got ",
       class(x)[1], call. = FALSE)
}

check_positive_number <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", arg, "` must be a single positive number", call. = FALSE)
  }
  invisible(NULL)
}
