fit_rr <- function(yi, vi) metafor::rma(yi = yi, vi = vi, measure = "RR")

test_that("available_methods lists the five detectors", {
  expect_setequal(available_methods(),
                  c("rcma", "ottawa", "barrowman", "sufficiency", "simulation"))
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
