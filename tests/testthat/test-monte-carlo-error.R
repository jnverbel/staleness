# Two detectors report a number produced by a finite simulation: simulation()
# a power over B replicates, sufficiency_changepoint() a p-value over n_perm
# permutations. Both printed that number to three decimals with nothing saying
# it would have been a different number under a different seed -- and the
# variability matters most exactly where the estimate sits beside the
# threshold that decides the verdict.
#
# Nothing here may change a verdict. That is the first test.

sim_prev <- function() {
  metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
               vi = c(0.16, 0.20, 0.18, 0.15, 0.22), measure = "RR")
}
sim_new <- function(k = 6, yi = log(0.55)) {
  list(yi = rep(yi, k), vi = rep(0.02, k), k = k)
}

test_that("Monte Carlo reporting leaves every verdict exactly as it was", {
  # The verdict is the published rule applied to the point estimate. Compared
  # against the rule recomputed here, not against a recorded string, so this
  # stays true if the rule itself is ever changed on purpose.
  prev <- sim_prev()
  for (k in c(4, 6, 10)) {
    for (thr in c(0.2, 0.5, 0.8, 0.95)) {
      v <- simulation(prev, sim_new(k), B = 300, power_threshold = thr,
                      seed = 7)
      expect_equal(v$verdict,
                   if (v$signal > thr) "out_of_date" else "current",
                   info = paste(k, thr))
    }
  }

  st <- metafor::rma(yi = c(rep(log(0.50), 6), rep(log(1.30), 6)),
                     vi = rep(0.02, 12), measure = "RR")
  pr <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.02, 6), measure = "RR")
  for (a in c(0.01, 0.05, 0.2)) {
    v <- sufficiency_changepoint(pr, st, alpha_stability = a, n_perm = 199,
                                 seed = 3)
    stable <- v$detail$p_stability >= a
    expect_equal(v$verdict,
                 if (v$detail$sufficient && !stable) "out_of_date" else "current",
                 info = a)
  }
})

test_that("simulation() reports the Monte Carlo error of its power", {
  v <- simulation(sim_prev(), sim_new(), B = 400, seed = 7)
  d <- v$detail
  n <- d$n_evaluable
  p <- v$signal

  expect_equal(d$mc_se, sqrt(p * (1 - p) / n))
  # The divisor is `n_evaluable`, not `B`: a replicate whose model would not
  # refit contributes to neither numerator nor denominator, and dividing by B
  # would shrink the standard error using replicates that produced nothing.
  #
  # Declared limit: this assertion cannot currently tell the two apart,
  # because no input found makes them differ. A mutation swapping `evaluable`
  # for `B` passes every test in this file. Since simulation() was changed to
  # hold tau2 at the prior's value, refits do not fail: 36 configurations
  # (k_prev 3-5 x vi 0.001-2 x k_new 1-3, B = 200) returned n_nonconverged = 0
  # in every one. The code uses `evaluable` because that is correct if a
  # failure ever occurs; the test says so rather than implying coverage.
  expect_equal(d$n_nonconverged, 0)
  expect_equal(n, d$B - d$n_nonconverged)
  # The interval brackets the estimate and stays inside [0, 1].
  expect_lte(d$mc_lo, p)
  expect_gte(d$mc_hi, p)
  expect_gte(d$mc_lo, 0)
  expect_lte(d$mc_hi, 1)

  # More replicates, a narrower interval. Quadrupling B should roughly halve
  # the standard error; checked loosely, since the estimate moves too.
  a <- simulation(sim_prev(), sim_new(), B = 200, seed = 7)$detail
  b <- simulation(sim_prev(), sim_new(), B = 3200, seed = 7)$detail
  expect_lt(b$mc_hi - b$mc_lo, a$mc_hi - a$mc_lo)
})

test_that("near_threshold marks the powers a rerun could have flipped", {
  # A threshold set to the estimate itself is inside any interval around it,
  # so it must be flagged; one far outside must not be. Both verdicts stand.
  v <- simulation(sim_prev(), sim_new(), B = 400, seed = 7)
  p <- v$signal
  skip_if(p <= 0.02 || p >= 0.98, "power too extreme to place thresholds")

  at <- simulation(sim_prev(), sim_new(), B = 400, seed = 7,
                   power_threshold = p)
  expect_true(at$detail$near_threshold)

  far_lo <- simulation(sim_prev(), sim_new(), B = 400, seed = 7,
                       power_threshold = 0)
  far_hi <- simulation(sim_prev(), sim_new(), B = 400, seed = 7,
                       power_threshold = 1)
  expect_false(far_lo$detail$near_threshold)
  expect_false(far_hi$detail$near_threshold)
})

test_that("the permutation p-value carries its own Monte Carlo interval", {
  pr <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.02, 6), measure = "RR")
  st <- metafor::rma(yi = c(rep(log(0.50), 6), rep(log(1.30), 6)),
                     vi = rep(0.02, 12), measure = "RR")
  v <- sufficiency_changepoint(pr, st, n_perm = 999, seed = 3)
  d <- v$detail

  expect_true(is.finite(d$mc_se))
  expect_lte(d$mc_lo, d$p_stability)
  expect_gte(d$mc_hi, d$p_stability)
  # The bounds live on the scale of the reported p-value, which is never below
  # 1 / (n_perm + 1) by construction, so neither is the lower bound.
  expect_gte(d$mc_lo, 1 / (999 + 1))
  expect_lte(d$mc_hi, 1)

  # near_threshold tracks alpha_stability, and the verdict does not move with
  # it: the same p is compared against the same alpha either way.
  at <- sufficiency_changepoint(pr, st, n_perm = 999, seed = 3,
                                alpha_stability = d$p_stability)
  expect_true(at$detail$near_threshold)
  expect_equal(at$detail$p_stability, d$p_stability)
})

test_that("a p-value that was never simulated reports no Monte Carlo error", {
  # A cumulative series constant to rounding short-circuits to p = 1 without
  # drawing a single permutation. That 1 is determined, not estimated, and
  # reporting an interval around it would invent a simulation that never ran.
  flat <- metafor::rma(yi = rep(log(0.50), 12), vi = rep(0.02, 12),
                       measure = "RR")
  pr   <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.02, 6),
                       measure = "RR")
  v <- sufficiency_changepoint(pr, flat, n_perm = 199, seed = 3)
  expect_equal(v$detail$p_stability, 1)
  expect_true(is.na(v$detail$mc_se))
  expect_true(is.na(v$detail$mc_lo))
  expect_true(is.na(v$detail$mc_hi))
  expect_false(v$detail$near_threshold)
})

# The two helpers, tested where the choice between them actually shows. An
# earlier version of this file reached the p = 1 case only through a
# simulation that had to happen to return 1, and it did not: a mutation
# swapping Wilson for Wald passed every assertion here.
test_that("the interval is Wilson, which keeps a width at 0 and at 1", {
  for (n in c(50, 200, 999)) {
    # Wald has width exactly zero at both ends -- it would report a simulated
    # certainty from a finite number of draws. Wilson must not.
    for (x in c(0, n)) {
      ci <- mc_interval(x, n)
      expect_gt(ci[2] - ci[1], 0)
      expect_gte(ci[1], 0)
      expect_lte(ci[2], 1)
    }
    # At 0 the lower bound is 0 and the upper is not; at n, the mirror image.
    expect_equal(mc_interval(0, n)[1], 0)
    expect_gt(mc_interval(0, n)[2], 0)
    expect_equal(mc_interval(n, n)[2], 1)
    expect_lt(mc_interval(n, n)[1], 1)
  }

  # And it is the Wilson formula, not merely something with a width: checked
  # against the closed form written out independently.
  n <- 400; x <- 137; z <- stats::qnorm(0.975); ph <- x / n
  d <- 1 + z^2 / n
  lo <- (ph + z^2 / (2 * n) - z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / d
  hi <- (ph + z^2 / (2 * n) + z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2))) / d
  expect_equal(mc_interval(x, n), c(lo, hi))
  # Which is not the Wald interval, on this same input.
  expect_false(isTRUE(all.equal(
    mc_interval(x, n), c(ph - z * sqrt(ph * (1 - ph) / n),
                         ph + z * sqrt(ph * (1 - ph) / n)))))

  # Nothing to interval: no draws, or a count that is not a number.
  expect_true(all(is.na(mc_interval(3, 0))))
  expect_true(all(is.na(mc_interval(NA_integer_, 100))))
  expect_true(is.na(mc_se(NA_integer_, 100)))
  expect_true(is.na(mc_se(3, 0)))
  expect_false(mc_near_threshold(c(NA_real_, NA_real_), 0.5))
})

test_that("the permutation SE is scaled by the (1 + x) / (n + 1) map", {
  # p_stability is not count/n_perm: it is (1 + count) / (n_perm + 1), an
  # affine map with slope n / (n + 1). Reporting the unscaled standard error
  # would overstate it, slightly and systematically. Without this, a mutation
  # dropping the factor passed the whole file.
  pr <- metafor::rma(yi = rep(log(0.50), 6), vi = rep(0.02, 6), measure = "RR")
  st <- metafor::rma(yi = c(rep(log(0.50), 6), rep(log(1.30), 6)),
                     vi = rep(0.02, 12), measure = "RR")
  for (n_perm in c(199, 999)) {
    v <- sufficiency_changepoint(pr, st, n_perm = n_perm, seed = 3)
    d <- v$detail
    # Recover the exceedance count from the reported p-value, then rebuild the
    # standard error from it without going through the package.
    count <- d$p_stability * (n_perm + 1) - 1
    ph <- count / n_perm
    expect_equal(d$mc_se,
                 sqrt(ph * (1 - ph) / n_perm) * n_perm / (n_perm + 1))
    # And the bounds are the raw Wilson bounds put through the same map.
    raw <- mc_interval(count, n_perm)
    expect_equal(c(d$mc_lo, d$mc_hi), (1 + raw * n_perm) / (n_perm + 1))
  }
})
