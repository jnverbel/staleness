fit_rr <- function(yi, vi) metafor::rma(yi = yi, vi = vi, measure = "RR")

test_that("rcma signals when the pooled effect halves on the natural scale", {
  prev <- fit_rr(rep(log(0.50), 4), rep(0.04, 4))
  upd  <- fit_rr(rep(log(0.25), 4), rep(0.04, 4))
  v <- rcma(prev, upd)
  expect_equal(v$verdict, "out_of_date")
  expect_equal(v$signal, 0.5, tolerance = 1e-6)
})

test_that("rcma reports current when the effect barely moves", {
  prev <- fit_rr(rep(log(0.50), 4), rep(0.04, 4))
  upd  <- fit_rr(rep(log(0.48), 4), rep(0.04, 4))
  v <- rcma(prev, upd)
  expect_equal(v$verdict, "current")
})

test_that("rcma signals on the upper threshold too", {
  prev <- fit_rr(rep(log(0.50), 4), rep(0.04, 4))
  upd  <- fit_rr(rep(log(0.80), 4), rep(0.04, 4))
  v <- rcma(prev, upd)
  expect_equal(v$signal, 1.6, tolerance = 1e-6)
  expect_equal(v$verdict, "out_of_date")
})

test_that("rcma refuses a near-null difference effect instead of guessing", {
  prev <- metafor::rma(yi = rep(0.02, 4), vi = rep(0.25, 4), measure = "MD")
  upd  <- metafor::rma(yi = rep(0.40, 4), vi = rep(0.25, 4), measure = "MD")
  v <- rcma(prev, upd)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "indistinguishable from zero")
})

test_that("rcma is computed on the natural scale, not the log scale", {
  # On the log scale log(0.25)/log(0.50) = 2.0, which would read as a signal
  # in the wrong direction. On the natural scale the ratio is 0.5.
  prev <- fit_rr(rep(log(0.50), 4), rep(0.04, 4))
  upd  <- fit_rr(rep(log(0.25), 4), rep(0.04, 4))
  expect_equal(rcma(prev, upd)$signal, 0.5, tolerance = 1e-6)
})
