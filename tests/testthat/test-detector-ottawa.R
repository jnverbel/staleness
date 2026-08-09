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

test_that("ottawa says when it stops reproducing the published arithmetic", {
  # Shojania et al. (2007), AHRQ Technical Review 16, Methods: "we performed
  # the updated meta-analyses by combining the original pooled result with the
  # individual results of eligible new trials". The prior pooled estimate goes
  # in as ONE point. check_currency() instead refits over all studies, which
  # is the better estimate but not automatically the same arithmetic.
  #
  # Read from the open-access report rather than from a summary of it, and
  # then checked numerically, because the report also claims the two coincide
  # under fixed effects and that claim is testable.
  set.seed(4)
  yi_old <- rnorm(8, -0.4, 0.2); vi_old <- runif(8, 0.02, 0.10)
  yi_new <- rnorm(4, -0.1, 0.2); vi_new <- runif(4, 0.02, 0.08)

  same <- function(m) {
    prev  <- metafor::rma(yi_old, vi_old, method = m)
    all_s <- metafor::rma(c(yi_old, yi_new), c(vi_old, vi_new), method = m)
    point <- metafor::rma(c(as.numeric(prev$beta), yi_new),
                          c(prev$se^2, vi_new), method = m)
    isTRUE(all.equal(as.numeric(all_s$beta), as.numeric(point$beta),
                     tolerance = 1e-8))
  }
  # Under fixed effects the report is right and there is nothing to fix.
  expect_true(same("FE"))
  # Under REML it is not: refitting all studies is a different computation.
  expect_false(same("REML"))

  # So the verdict carries which of the two situations it is in.
  prev_fe <- metafor::rma(yi_old, vi_old, measure = "RR", method = "FE")
  upd_fe  <- metafor::rma(c(yi_old, yi_new), c(vi_old, vi_new),
                          measure = "RR", method = "FE")
  expect_true(ottawa(prev_fe, upd_fe)$detail$reproduces_published)

  prev_re <- metafor::rma(yi_old, vi_old, measure = "RR", method = "REML")
  upd_re  <- metafor::rma(c(yi_old, yi_new), c(vi_old, vi_new),
                          measure = "RR", method = "REML")
  expect_false(ottawa(prev_re, upd_re)$detail$reproduces_published)
})
