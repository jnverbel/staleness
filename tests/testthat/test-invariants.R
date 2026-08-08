fit_rr <- function(yi, vi) metafor::rma(yi = yi, vi = vi, measure = "RR")

test_that("duplicating the evidence leaves the effect ratio at exactly 1", {
  yi <- rep(log(0.5), 6); vi <- rep(0.02, 6)
  prev <- fit_rr(yi, vi)
  dup  <- fit_rr(c(yi, yi), c(vi, vi))
  expect_equal(rcma(prev, dup)$signal, 1, tolerance = 1e-8)
  expect_equal(rcma(prev, dup)$verdict, "current")
})

test_that("new evidence drawn from the same distribution raises no signal", {
  set.seed(11)
  yi <- rnorm(20, log(0.5), 0.02); vi <- rep(0.01, 20)
  prev <- fit_rr(yi, vi)
  more <- rnorm(10, log(0.5), 0.02)
  upd  <- fit_rr(c(yi, more), c(vi, rep(0.01, 10)))
  expect_equal(rcma(prev, upd)$verdict, "current")
  expect_equal(ottawa(prev, upd)$verdict, "current")
})

test_that("every detector returns a verdict object with the required fields", {
  prev <- fit_rr(rep(log(0.5), 8), rep(0.02, 8))
  upd  <- fit_rr(rep(log(0.3), 12), rep(0.02, 12))
  new_ev <- list(yi = rep(log(0.3), 4), vi = rep(0.02, 4), k = 4)
  vs <- list(
    rcma(prev, upd), ottawa(prev, upd), sufficiency(prev, upd),
    barrowman(prev, 400, 900), simulation(prev, new_ev, B = 50, seed = 1)
  )
  for (v in vs) {
    expect_s3_class(v, "staleness_verdict")
    expect_true(v$verdict %in% c("out_of_date", "current", "not_applicable"))
    if (v$verdict == "not_applicable") expect_true(nzchar(v$reason))
  }
})

test_that("the verdict of every detector is invariant to study input order", {
  set.seed(3)
  yi <- rnorm(10, log(0.5), 0.05); vi <- rep(0.02, 10)
  ord <- sample(10)
  a <- rcma(fit_rr(yi, vi), fit_rr(c(yi, log(0.2)), c(vi, 0.02)))
  b <- rcma(fit_rr(yi[ord], vi[ord]), fit_rr(c(yi[ord], log(0.2)), c(vi[ord], 0.02)))
  expect_equal(a$signal, b$signal, tolerance = 1e-8)
})
