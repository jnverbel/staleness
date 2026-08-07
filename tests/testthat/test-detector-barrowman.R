fit_md <- function(yi, vi) metafor::rma(yi = yi, vi = vi, measure = "MD")

test_that("barrowman computes q = m / n with n = Zc^2 * N / Z^2", {
  # prior non-significant: zval small
  prev <- fit_md(c(0.10, -0.05, 0.08, -0.02), rep(0.25, 4))
  N <- 400; m <- 900
  z <- prev$zval
  n_expected <- (1.96^2 * N) / z^2
  v <- barrowman(prev, n_prev = N, n_new = m)
  expect_equal(v$detail$n_required, n_expected, tolerance = 1e-8)
  expect_equal(v$signal, m / n_expected, tolerance = 1e-8)
})

test_that("barrowman signals when q exceeds 1", {
  prev <- fit_md(c(0.10, -0.05, 0.08, -0.02), rep(0.25, 4))
  n_req <- (1.96^2 * 400) / prev$zval^2
  expect_equal(barrowman(prev, 400, ceiling(n_req * 1.5))$verdict, "out_of_date")
  expect_equal(barrowman(prev, 400, floor(n_req * 0.5))$verdict, "current")
})

test_that("barrowman does not apply to an already significant meta-analysis", {
  prev <- metafor::rma(yi = rep(log(0.4), 8), vi = rep(0.01, 8), measure = "RR")
  v <- barrowman(prev, n_prev = 400, n_new = 900)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "already significant")
  expect_true(is.na(v$signal))
})

test_that("barrowman refuses a prior effect at exactly zero", {
  prev <- fit_md(c(0.5, -0.5, 0.5, -0.5), rep(0.25, 4))  # zval ~ 0
  v <- barrowman(prev, n_prev = 400, n_new = 900)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "too close to zero")
})

test_that("barrowman requires sample sizes", {
  prev <- fit_md(c(0.10, -0.05, 0.08, -0.02), rep(0.25, 4))
  expect_error(barrowman(prev, n_prev = NULL, n_new = 900), "sample size")
})
