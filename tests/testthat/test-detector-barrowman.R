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

# --- a sample size that is present but unusable is a verdict, not an error --
#
# barrowman() stopped on any non-finite size, which turned one NA anywhere in a
# stream into a hard failure of the whole backtest. Omitting an argument is a
# malformed call and stays an error; a sample size that is present but is NA,
# NaN or infinite is a fact about the data, and gets the same "not_applicable
# with a reason" treatment as every other question the method cannot answer.

test_that("a non-finite sample size yields not_applicable, not an error", {
  prev <- metafor::rma(yi = c(0.10, -0.05, 0.08, -0.02), vi = rep(0.05, 4),
                       measure = "MD")
  for (bad in list(NA_real_, NaN, Inf)) {
    v <- expect_no_error(barrowman(prev, n_prev = bad, n_new = 500))
    expect_equal(v$verdict, "not_applicable")
    expect_match(v$reason, "not a finite number")
    w <- expect_no_error(barrowman(prev, n_prev = 500, n_new = bad))
    expect_equal(w$verdict, "not_applicable")
  }
})

test_that("omitting a sample size altogether is still an error", {
  prev <- metafor::rma(yi = c(0.10, -0.05, 0.08, -0.02), vi = rep(0.05, 4),
                       measure = "MD")
  expect_error(barrowman(prev, n_prev = NULL, n_new = 500), "needs the sample size")
  expect_error(barrowman(prev, n_prev = 500, n_new = NULL), "needs the sample size")
})

test_that("a non-finite prior p-value is declined, not fatal", {
  # `if (prev$pval < alpha)` died with "missing value where TRUE/FALSE needed"
  # -- an error about an if(), naming nothing the caller can act on.
  prev <- metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
                       vi = c(0.16, 0.20, 0.18, 0.15, 0.22))
  prev$pval <- NA_real_
  res <- barrowman(prev, n_prev = 555, n_new = 2265)
  expect_equal(res$verdict, "not_applicable")
  expect_match(res$reason, "p-value")
})
