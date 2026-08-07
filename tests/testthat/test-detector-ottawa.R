fit_rr <- function(yi, vi) metafor::rma(yi = yi, vi = vi, measure = "RR")

test_that("ottawa signals when a null result becomes significant", {
  prev <- fit_rr(c(0.02, -0.01, 0.03, -0.02), rep(0.05, 4))   # p >> 0.04
  upd  <- fit_rr(rep(log(0.40), 8), rep(0.01, 8))             # p << 0.04
  v <- ottawa(prev, upd)
  expect_equal(v$verdict, "out_of_date")
  expect_true(v$detail$signal_significance)
})

test_that("ottawa uses the published 0.04 threshold, not 0.05", {
  prev <- fit_rr(c(0.02, -0.01, 0.03, -0.02), rep(0.05, 4))
  expect_equal(formals(ottawa)$alpha, 0.04)
})

test_that("ottawa signals on a 50 percent effect change alone", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))   # already significant
  upd  <- fit_rr(rep(log(0.20), 6), rep(0.02, 6))   # still significant, big shift
  v <- ottawa(prev, upd)
  expect_equal(v$verdict, "out_of_date")
  expect_false(v$detail$signal_significance)
  expect_true(v$detail$signal_effect)
})

test_that("ottawa reports current when neither signal fires", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  upd  <- fit_rr(rep(log(0.49), 6), rep(0.02, 6))
  v <- ottawa(prev, upd)
  expect_equal(v$verdict, "current")
})

test_that("sig_change = 'gain' ignores loss of significance, 'any' catches it", {
  prev <- fit_rr(rep(log(0.50), 8), rep(0.01, 8))   # significant
  upd  <- fit_rr(c(rep(log(0.50), 8), rep(log(1.6), 12)), c(rep(0.01, 8), rep(0.01, 12)))
  expect_false(ottawa(prev, upd, sig_change = "gain")$detail$signal_significance)
  expect_true(ottawa(prev, upd, sig_change = "any")$detail$signal_significance)
})

test_that("qualitative signals are carried but never inferred", {
  prev <- fit_rr(rep(log(0.50), 6), rep(0.02, 6))
  upd  <- fit_rr(rep(log(0.49), 6), rep(0.02, 6))
  v <- ottawa(prev, upd, qualitative = "substantial harm reported")
  expect_equal(v$verdict, "out_of_date")
  expect_equal(v$detail$qualitative, "substantial harm reported")
  # with no declared qualitative signal the same data reads as current
  expect_equal(ottawa(prev, upd)$verdict, "current")
})
