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

test_that("an unstable RRR ratio is flagged without altering the verdict", {
  # The Ottawa effect criterion divides by 1 - RR_prev, which approaches zero
  # exactly on the null reviews the method targets. The instability is a
  # property of the published method and is NOT corrected here -- correcting
  # it would mean implementing a different method, and the point of this
  # package is to find out how the published ones behave.
  #
  # What the package must not do is hand back a ratio of -19 with no sign
  # that it came from dividing by -0.005. Flagged, like CONTAMINATED_PAIRS
  # flags a circular pair: visible in the result, not silently removed.
  prev <- metafor::rma(yi = c(0.02, -0.01, 0.03, -0.02), vi = rep(0.05, 4),
                       measure = "RR")
  new  <- metafor::rma(yi = c(0.02, -0.01, 0.03, -0.02, -0.30, -0.25),
                       vi = c(rep(0.05, 4), 0.04, 0.04), measure = "RR")
  res <- ottawa(prev, new)

  expect_true(res$detail$effect_unstable)
  expect_true(abs(res$detail$rrr_prev) < 0.05)
  # The verdict and the signal are untouched: fidelity to the method.
  expect_equal(res$verdict, "out_of_date")
  expect_true(is.finite(res$signal))

  # A prior effect that is clearly away from no-effect is not flagged.
  prev2 <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.01, 6),
                        measure = "RR")
  new2  <- metafor::rma(yi = c(rep(log(0.50), 6), rep(log(0.80), 4)),
                        vi = c(rep(0.01, 6), rep(0.01, 4)), measure = "RR")
  expect_false(ottawa(prev2, new2)$detail$effect_unstable)
})

test_that("difference measures report the flag as FALSE, not missing", {
  # On difference measures the effect half defers to the rcma rule, which has
  # its own guard, so there is no RRR to be unstable. The field must still be
  # present and answerable.
  prev <- metafor::rma(yi = c(0.40, 0.35, 0.45), vi = rep(0.01, 3),
                       measure = "MD")
  new  <- metafor::rma(yi = c(0.40, 0.35, 0.45, 0.10), vi = rep(0.01, 4),
                       measure = "MD")
  res <- ottawa(prev, new)
  expect_false(res$detail$effect_unstable)
  expect_true(is.na(res$detail$rrr_prev))
})
