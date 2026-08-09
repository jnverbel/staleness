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

test_that("an infinite ratio is not a signal", {
  # exp(theta_new - theta_prev) overflows when the prior effect underflows to
  # zero on the natural scale. Inf passed is.na(), so the detector reported
  # out_of_date with a signal of Inf -- a number no reader can act on.
  prev <- metafor::rma(yi = rep(-500, 4), vi = rep(0.05, 4), measure = "RR")
  new  <- metafor::rma(yi = rep(500, 4), vi = rep(0.05, 4), measure = "RR")
  res <- rcma(prev, new)
  expect_equal(res$verdict, "not_applicable")
  expect_true(is.na(res$signal))
  expect_match(res$reason, "finite")
})

test_that("rcma is not contained in ottawa, in both directions", {
  # The documentation asserted containment for several versions, and the
  # assertion survived the correction in two places because nobody had pinned
  # it to a run. A claim about two detectors' relationship belongs in a test:
  # prose cannot go red.
  #
  # Direction one: rcma fires and ottawa does not. The pooled effect moves
  # from RR 0.20 to 0.33 -- a ratio of 1.63 -- while the risk reduction it
  # implies goes from 80% to 67%, a ratio of 0.84, well inside ottawa's band.
  strong <- metafor::rma(yi = rep(log(0.20), 4), vi = rep(0.05, 4),
                         measure = "RR")
  weaker <- metafor::rma(yi = c(rep(log(0.20), 4), rep(log(0.50), 4)),
                         vi = c(rep(0.05, 4), rep(0.02, 4)), measure = "RR")
  expect_equal(rcma(strong, weaker)$verdict, "out_of_date")
  expect_equal(ottawa(strong, weaker)$verdict, "current")

  # Direction two is already covered by the ten published reviews in
  # test-external-validation.R, where all ten fire on ottawa's RRR ratio and
  # not one fires on rcma's ratio of effects.

  # And the signals are different numbers even when the verdicts agree, which
  # is the cheapest way to see that one rule is not the other.
  prev <- metafor::rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4),
                       measure = "RR")
  upd  <- metafor::rma(yi = c(rep(log(0.50), 4), rep(log(1.10), 4)),
                       vi = c(rep(0.05, 4), rep(0.02, 4)), measure = "RR")
  expect_equal(rcma(prev, upd)$verdict, ottawa(prev, upd)$verdict)
  expect_false(isTRUE(all.equal(rcma(prev, upd)$signal,
                                ottawa(prev, upd)$signal)))
})
