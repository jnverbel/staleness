VALID_VERDICTS <- c("out_of_date", "current", "not_applicable")

#' Construct a detector verdict
#'
#' @param method Character, detector name.
#' @param verdict One of `"out_of_date"`, `"current"`, `"not_applicable"`.
#' @param signal Numeric, the detector's continuous statistic, or `NA_real_`.
#' @param reason Character, mandatory when `verdict` is `"not_applicable"`.
#' @param detail List, method-specific extras.
#' @return An object of class `staleness_verdict`.
#' @keywords internal
new_verdict <- function(method, verdict, signal = NA_real_, reason = "",
                        detail = list()) {
  # Checked here rather than left to print.staleness_verdict, where nzchar()
  # on a vector dies with "the condition has length > 1" -- an error about an
  # if(), raised far from the call that caused it, naming nothing the caller
  # can act on.
  if (!is.character(reason) || length(reason) != 1L || is.na(reason)) {
    stop("`reason` must be a single, non-NA string (a scalar)", call. = FALSE)
  }
  if (!verdict %in% VALID_VERDICTS) {
    stop("`verdict` must be one of: ", paste(VALID_VERDICTS, collapse = ", "),
         call. = FALSE)
  }
  if (verdict == "not_applicable" && !nzchar(reason)) {
    stop("a not_applicable verdict requires a non-empty `reason`", call. = FALSE)
  }
  structure(
    list(method = method, verdict = verdict, signal = signal,
         reason = reason, detail = detail),
    class = "staleness_verdict"
  )
}

#' @rdname new_verdict
#' @keywords internal
verdict_na <- function(method, reason) {
  new_verdict(method, "not_applicable", signal = NA_real_, reason = reason)
}

#' @export
print.staleness_verdict <- function(x, ...) {
  label <- switch(x$verdict,
    out_of_date    = "OUT OF DATE",
    current        = "current",
    not_applicable = "not applicable"
  )
  cat(format(x$method, width = 12), label, "\n")
  if (nzchar(x$reason)) cat("  reason:", x$reason, "\n")
  if (!is.na(x$signal)) cat("  signal:", format(x$signal, digits = 3), "\n")
  invisible(x)
}
