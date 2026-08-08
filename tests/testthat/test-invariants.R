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
  expect_equal(v$detail$p_slope, 1)
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
