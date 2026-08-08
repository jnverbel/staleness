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

# --- The invariant, across many seeds and ALL FIVE detectors ---------------
#
# The single-seed version above is what let the sufficiency defect ship. Design
# section 7 states the invariant as "statistically identical new evidence => no
# detector signals", but the test only ever checked two of the five, at one
# seed. sufficiency() fired on 209 of 300 such samples: its stability test was
# an OLS t-test on a cumulative mean, a series that is autocorrelated by
# construction and converges by the law of large numbers, so the test detected
# convergence and reported instability. Replaced by a permutation test over
# study order (see ?sufficiency); the invariant is now checked the way it
# should have been from the start.
#
# The right assertion is a RATE, not "never fires": rcma and ottawa are
# deterministic threshold rules and must never fire on unchanging evidence,
# while sufficiency is now an honest alpha-level test and is entitled to fire
# on about alpha of samples. Demanding zero firings from it would only be
# satisfiable by a broken test.

#' Run all five detectors over `n` samples of evidence that does not change.
#'
#' `effect` and `sd_study` pick which regime is exercised: a strong, precise
#' effect makes `sufficiency` applicable and *sufficient* (so its stability
#' half is genuinely under test) while making `barrowman` and `simulation`
#' inapplicable; a null effect does the reverse. Both regimes are needed to
#' put all five detectors under the invariant, which is why there are two
#' tests below and not one.
no_change_run <- function(seed0, n, effect, vi_study, n_prev = 2000, n_new = 1000) {
  verdicts <- matrix(NA_character_, nrow = n, ncol = 5,
                     dimnames = list(NULL, available_methods()))
  for (i in seq_len(n)) {
    set.seed(seed0 + i)
    yi   <- rnorm(20, effect, 0.05)
    vi   <- rep(vi_study, 20)
    more <- rnorm(10, effect, 0.05)          # same distribution: nothing changed
    prev <- fit_rr(yi, vi)
    upd  <- fit_rr(c(yi, more), c(vi, rep(vi_study, 10)))
    new  <- list(yi = more, vi = rep(vi_study, 10), k = 10)
    vs <- list(
      rcma        = rcma(prev, upd),
      ottawa      = ottawa(prev, upd),
      barrowman   = barrowman(prev, n_prev = n_prev, n_new = n_new),
      sufficiency = sufficiency(prev, upd),
      simulation  = simulation(prev, new, B = 200, seed = 1)
    )
    for (m in available_methods()) verdicts[i, m] <- vs[[m]]$verdict
  }
  verdicts
}

test_that("unchanging evidence never signals: null-effect regime, all five detectors apply", {
  # A null effect keeps the prior meta-analysis non-significant, which is the
  # only regime where barrowman() and simulation() have anything to say. All
  # five detectors answer at every seed, and none of them may fire.
  v <- no_change_run(seed0 = 1000, n = 60, effect = 0, vi_study = 0.05)
  for (m in available_methods()) {
    expect_equal(sum(v[, m] == "not_applicable"), 0, info = m)
    expect_equal(sum(v[, m] == "out_of_date"), 0, info = m)
  }
})

test_that("unchanging evidence signals no more often than chance: strong-effect regime", {
  # A strong, precise effect is what makes sufficiency() *sufficient*
  # (index ~ 100), so `sufficient && !stable` turns entirely on the stability
  # test -- the regime the defect lived in. It is also, necessarily, a regime
  # where barrowman() and simulation() decline: both require a prior that was
  # not yet significant. Declining is the correct answer and is asserted, so
  # neither can collect an unearned "current".
  v <- no_change_run(seed0 = 0, n = 60, effect = log(0.5), vi_study = 0.01)

  expect_equal(sum(v[, "rcma"]   == "out_of_date"), 0)
  expect_equal(sum(v[, "ottawa"] == "out_of_date"), 0)
  expect_equal(sum(v[, "barrowman"]  != "not_applicable"), 0)
  expect_equal(sum(v[, "simulation"] != "not_applicable"), 0)

  # sufficiency was applicable and sufficient throughout, so every "current"
  # here is the stability test declining to fire.
  expect_equal(sum(v[, "sufficiency"] == "not_applicable"), 0)

  # It is now an honest alpha = 0.05 permutation test, so it is *entitled* to
  # fire on about 5% of samples; demanding zero would only be satisfiable by a
  # broken test. Over 60 samples the upper 99.9% binomial bound is 8. The
  # pre-fix OLS implementation fired on roughly 70% of such samples (209 of
  # 300 in the measured experiment), so this bound has real teeth.
  expect_lte(sum(v[, "sufficiency"] == "out_of_date"), 8)
})

# --- The invariant under a REALISTIC variance schedule ---------------------
#
# no_change_run() above gives every study the same variance, and that is what
# let the second sufficiency defect ship. A permutation test over study order
# is exact only if the studies are exchangeable, and in a real evidence stream
# they are not: trials get larger over time, so the variance schedule is
# informative about position even when nothing about the effect is changing.
# Under that schedule the cumulative series wanders early -- when only small
# trials have accumulated -- and then settles as the large ones arrive, which
# looks exactly like drift. Permuting orders puts the large trials early half
# the time, so the null is flat and the real series looks extreme against it.
#
# Measured on the slope statistic that shipped before: 83/300 = 28% false
# alarms with variance falling over time, and 101/300 = 34% on a step schedule
# of twenty small trials followed by ten large ones. Five to seven times
# nominal, on evidence with no drift in it whatsoever. The fix was to make the
# statistic pivotal so it carries no imprint of the schedule (see
# ?sufficiency); this is the test that holds it to that.

# Counts how often the STABILITY half fires, not how often the verdict does:
# a verdict-level count would be diluted by the occasional sample whose `prev`
# came out insufficient, which says nothing about the property under test.
hetero_false_alarms <- function(vi, n = 200, effect = log(0.5), seed0 = 9000) {
  k <- length(vi); k_prev <- 20
  unstable <- 0L
  for (i in seq_len(n)) {
    set.seed(seed0 + i)
    yi <- rnorm(k, effect, sqrt(vi))     # no drift: one mean throughout
    prev <- metafor::rma(yi = yi[seq_len(k_prev)], vi = vi[seq_len(k_prev)],
                         measure = "RR", method = "FE")
    upd  <- metafor::rma(yi = yi, vi = vi, measure = "RR", method = "FE")
    unstable <- unstable + !sufficiency(prev, upd)$detail$stable
  }
  unstable
}

test_that("time-correlated heteroscedasticity does not manufacture false alarms", {
  # 200 samples, no drift, variance falling steadily over time: early small
  # trials, later large ones. At the nominal 5% the expected count is 10 and
  # the upper 99.9% binomial bound is 20; the shipped slope statistic produced
  # about 56 here, so the bound has teeth.
  expect_lte(hetero_false_alarms(seq(0.16, 0.01, length.out = 30)), 20)
})

test_that("a step change in study size does not manufacture false alarms", {
  # The harder version of the same shape, and the one the slope statistic was
  # worst on (34%): twenty small trials, then ten large ones, no drift.
  expect_lte(hetero_false_alarms(c(rep(0.25, 20), rep(0.005, 10))), 20)
})

test_that("the reverse schedule does not destroy the test either", {
  # Variance *rising* over time is the mirror image, and under the slope
  # statistic it was worse than a false alarm: the rate fell to 0/200, i.e. the
  # test had no power left at all in that regime -- a detector that can never
  # fire trivially passes an upper bound. So this one is bounded on BOTH sides.
  # At a true 5% rate the chance of seeing zero firings in 200 samples is
  # 3.5e-5, so the lower bound is safe and it is the half with the teeth here.
  fired <- hetero_false_alarms(seq(0.01, 0.16, length.out = 30))
  expect_lte(fired, 20)
  expect_gte(fired, 1)
})

test_that("sufficiency does not fire on a cumulative series that is only rounding noise", {
  # 12 byte-identical studies. The cumulative effect spans 2.2e-16 -- pure
  # floating-point noise around a constant -- and the old OLS test fitted that
  # ULP, returning slope -7.4e-17 with p = 0.014 and a verdict of out_of_date.
  # A constant series is maximally stable, by definition.
  prev <- fit_rr(rep(log(0.5), 12), rep(0.02, 12))
  cum <- cumulative_effect(as.numeric(prev$yi), as.numeric(prev$vi))
  expect_lt(diff(range(cum)), 1e-15)     # the series really is degenerate
  v <- sufficiency(prev, prev)
  expect_true(v$detail$sufficient)       # sufficiency is not the thing failing
  expect_true(v$detail$stable)
  expect_equal(v$detail$slope, 0)
  expect_equal(v$detail$p_stability, 1)
  expect_equal(v$verdict, "current")
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
