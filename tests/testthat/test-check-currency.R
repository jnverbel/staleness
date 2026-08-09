fit_rr <- function(yi, vi) metafor::rma(yi = yi, vi = vi, measure = "RR")

test_that("available_methods lists the five detectors", {
  expect_setequal(available_methods(),
                  c("rcma", "ottawa", "barrowman", "sufficiency_changepoint", "simulation"))
})

test_that("check_currency runs the requested detectors and returns one row each", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  new  <- list(yi = rep(log(0.20), 4), vi = rep(0.02, 4), k = 4)
  res <- check_currency(prev, new, methods = c("rcma", "ottawa"))
  expect_s3_class(res, "staleness_check")
  df <- as.data.frame(res)
  expect_equal(nrow(df), 2)
  expect_setequal(df$method, c("rcma", "ottawa"))
  expect_true(all(df$verdict %in% c("out_of_date", "current", "not_applicable")))
})

test_that("check_currency refits the updated meta-analysis from prior plus new", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  new  <- list(yi = rep(log(0.20), 4), vi = rep(0.02, 4), k = 4)
  res <- check_currency(prev, new, methods = "rcma")
  expect_equal(res$updated$k, 10)
})

test_that("barrowman is skipped with a reason when sample sizes are absent", {
  prev <- fit_rr(c(0.02, -0.01, 0.03, -0.02), rep(0.05, 4))
  new  <- list(yi = c(0.01, 0.02), vi = c(0.05, 0.05), k = 2)
  df <- as.data.frame(check_currency(prev, new, methods = "barrowman"))
  expect_equal(df$verdict, "not_applicable")
  expect_match(df$reason, "sample size")
})

test_that("an unknown method name is an error", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  new  <- list(yi = rep(log(0.20), 4), vi = rep(0.02, 4), k = 4)
  expect_error(check_currency(prev, new, methods = "nope"), "unknown method")
})

test_that("no new evidence is its own class, not a clean bill of health", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  none <- list(yi = numeric(0), vi = numeric(0), k = 0)
  res <- check_currency(prev, none)
  expect_s3_class(res, "staleness_no_evidence")
  df <- as.data.frame(res)
  expect_true(all(df$verdict == "not_applicable"))
  expect_true(all(nzchar(df$reason)))
  # crucially: nothing reads as "current"
  expect_false(any(df$verdict == "current"))
})

test_that("disagreement between detectors is reported", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  new  <- list(yi = rep(log(0.49), 4), vi = rep(0.02, 4), k = 4)
  res <- check_currency(prev, new, methods = c("rcma", "ottawa"))
  expect_false(res$disagreement)   # both say current
})

test_that("an rma.uni passed as `new` is refused, not silently double-counted", {
  prev <- fit_rr(rep(log(0.50), 4), rep(0.05, 4))
  # The updated meta-analysis: the four prior studies plus three new ones.
  updated <- fit_rr(c(rep(log(0.50), 4), log(0.90), log(1.05), log(0.98)),
                    c(rep(0.05, 4), 0.03, 0.03, 0.02))
  # `new` wants the new evidence alone, but rcma()/ottawa() want the updated
  # model, so confusing the two is the natural mistake. An rma.uni carries
  # $yi, $vi and $k, so duck typing lets it through every guard: the prior
  # studies would be concatenated onto themselves and pooled twice.
  expect_error(check_currency(prev, updated, methods = "rcma"),
               "new evidence")
})

test_that("the guard does not reject a legitimate list", {
  prev <- fit_rr(rep(log(0.50), 4), rep(0.05, 4))
  new  <- list(yi = c(log(0.90), log(1.05)), vi = c(0.03, 0.03), k = 2)
  res  <- check_currency(prev, new, methods = "rcma")
  expect_s3_class(res, "staleness_check")
  expect_equal(res$updated$k, 6)   # 4 prior + 2 new, each counted once
})

# --- The `new` contract: k decides, yi/vi are fitted, nothing tied them -----
#
# check_currency() reads `new$k` to decide whether any new evidence exists and
# `new$yi`/`new$vi` to fit the updated model. Nothing required the two to
# agree, so the object could say one thing and carry another.

test_that("check_currency refuses evidence whose k contradicts its data", {
  prev <- metafor::rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4),
                       measure = "RR")

  # The worst case. This one returned a verified-looking "current" from no new
  # evidence at all: the updated model was refitted on the prior studies
  # alone, so every ratio was 1. The guard above it exists precisely to stop
  # that -- "absence of new evidence is not evidence of currency" -- and
  # setting k = 1 walked straight past it.
  expect_error(
    check_currency(prev, list(yi = numeric(), vi = numeric(), k = 1)),
    "k")

  # Real evidence discarded in silence, because k said there was none.
  expect_error(
    check_currency(prev, list(yi = log(0.9), vi = 0.03, k = 0)),
    "k")

  # Declared 99, carried 1. The verdict was unaffected but detail$k_new
  # reported 99, so the output stated a study count that was never used.
  expect_error(
    check_currency(prev, list(yi = log(0.9), vi = 0.03, k = 99)),
    "k")

  # Mismatched lengths used to fail inside metafor with "Length of 'yi' and
  # 'vi' (or 'sei') are not the same" -- a message about metafor's arguments,
  # not about this package's.
  expect_error(
    check_currency(prev, list(yi = c(log(0.9), log(0.8)), vi = 0.03, k = 2)),
    "same length")

  # What window_between() actually produces still passes untouched.
  expect_silent(check_currency(prev, list(yi = rep(log(0.9), 3),
                                          vi = rep(0.03, 3), k = 3),
                               methods = "rcma"))
  # And k = 0 with nothing in it is the honest empty case, which keeps its
  # own class rather than becoming an error.
  empty <- check_currency(prev, list(yi = numeric(), vi = numeric(), k = 0),
                          methods = "rcma")
  expect_s3_class(empty, "staleness_no_evidence")
})

test_that("check_currency refuses an empty method list, as backtest does", {
  # backtest() raises "`methods` is empty; name at least one of: ..." while
  # check_currency() built a staleness_check holding zero verdicts -- an
  # object that answers no question, from the same package, for the same
  # mistake.
  prev <- metafor::rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4),
                       measure = "RR")
  new <- list(yi = rep(log(0.9), 3), vi = rep(0.03, 3), k = 3)
  expect_error(check_currency(prev, new, methods = character()), "empty")
  expect_error(check_currency(prev, new, methods = NULL), "empty")
})

test_that("new evidence must be statistically usable, not just well shaped", {
  # The structural contract (k == length(yi), matching lengths) said nothing
  # about the values. Every degenerate case below returned "current" from both
  # rcma and ottawa -- the same failure mode as k = 1 with an empty yi, which
  # the previous fix closed: an NA or infinite variance makes metafor drop or
  # ignore the study, the "updated" model comes back identical to the prior
  # one, every ratio is 1, and the verdict looks verified.
  prev <- metafor::rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4),
                       measure = "RR")
  expect_error(check_currency(prev, list(yi = NA_real_, vi = 0.03, k = 1)),
               "finite")
  expect_error(check_currency(prev, list(yi = Inf, vi = 0.03, k = 1)),
               "finite")
  expect_error(check_currency(prev, list(yi = log(0.9), vi = NA_real_, k = 1)),
               "positive")
  expect_error(check_currency(prev, list(yi = log(0.9), vi = Inf, k = 1)),
               "positive")
  expect_error(check_currency(prev, list(yi = log(0.9), vi = 0, k = 1)),
               "positive")
  expect_error(check_currency(prev, list(yi = log(0.9), vi = -1, k = 1)),
               "positive")

  # simulation() takes the same object directly. Four of these threw R's own
  # "missing value where TRUE/FALSE needed"; yi = NA returned a verdict of
  # current with a signal, because simulation reads vi and k but never yi.
  ns <- metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
                     vi = c(0.16, 0.20, 0.18, 0.15, 0.22), measure = "RR")
  expect_error(simulation(ns, list(yi = log(0.9), vi = NA_real_, k = 1), B = 5),
               "positive")
  expect_error(simulation(ns, list(yi = log(0.9), vi = 0, k = 1), B = 5),
               "positive")
  expect_error(simulation(ns, list(yi = NA_real_, vi = 0.03, k = 1), B = 5),
               "finite")

  # The empty case is still empty, not invalid: zero-length yi and vi are
  # vacuously finite and positive, and k = 0 keeps its own class.
  expect_s3_class(
    check_currency(prev, list(yi = numeric(), vi = numeric(), k = 0),
                   methods = "rcma"),
    "staleness_no_evidence")
})

# The registry and the dispatch are two lists of the same names in two files,
# and nothing tied them together. Renaming sufficiency() to
# sufficiency_changepoint() moved available_methods() and left the switch()
# label behind -- a label is a NAME, not a string, so no search for the old
# string found it. switch() answered NULL, check_currency() carried the NULL
# to a vapply four frames later, and the error named a length, not a method.
test_that("every name in available_methods() actually dispatches to its own detector", {
  prev <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.05, 6), measure = "RR")
  new  <- list(yi = rep(log(1.20), 6), vi = rep(0.02, 6), k = 6)

  res <- check_currency(prev, new, methods = available_methods(),
                        n_prev = 600, n_new = 600)

  # One verdict per requested method, each carrying the name it was asked for.
  expect_named(res$verdicts, available_methods())
  for (m in available_methods()) {
    v <- res$verdicts[[m]]
    expect_s3_class(v, "staleness_verdict")
    expect_identical(v$method, m)
  }
})

test_that("a method with no dispatch branch fails by name, not by length", {
  prev <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.05, 6), measure = "RR")
  new  <- list(yi = rep(log(1.20), 6), vi = rep(0.02, 6), k = 6)

  # Force the desync the test above exists to prevent: a registry that offers a
  # name the switch() has never heard of.
  local_mocked_bindings(available_methods = function() c("rcma", "ghost"))
  expect_error(check_currency(prev, new, methods = "ghost"),
               "no dispatch branch for method \"ghost\"")
})
