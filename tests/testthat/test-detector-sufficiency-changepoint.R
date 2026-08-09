test_that("failsafe_n follows Rosenthal's formula", {
  yi <- rep(log(0.5), 5); vi <- rep(0.04, 5)
  z  <- yi / sqrt(vi)
  expected <- (sum(z)^2) / (1.645^2) - length(z)
  expect_equal(failsafe_n(yi, vi), expected, tolerance = 1e-8)
})

test_that("sufficiency_changepoint refuses fewer than five studies", {
  prev <- metafor::rma(yi = rep(log(0.5), 4), vi = rep(0.04, 4), measure = "RR")
  v <- sufficiency_changepoint(prev, prev)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "at least 5")
})

test_that("a large, stable body of evidence reads as current", {
  # Per the primary source (Pattanittum et al. 2012, Table 1), an out-of-date
  # verdict requires sufficiency AND instability together; sufficient-and-stable
  # is squarely "current". seed = 8 was originally picked to dodge the old OLS
  # test's serially correlated errors, which made a no-trend sample look
  # significant for many seeds. Under the permutation test the seed no longer
  # has to be hand-picked -- the multi-seed invariant in test-invariants.R
  # shows the property holds at the nominal rate across 60 samples -- but it is
  # kept so the assertion below stays a fixed, inspectable case.
  set.seed(8)
  yi <- rnorm(30, log(0.5), 0.02); vi <- rep(0.01, 30)
  prev <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  v <- sufficiency_changepoint(prev, prev)
  expect_equal(v$verdict, "current")
  expect_true(v$detail$sufficient)
  expect_true(v$detail$stable)
})

test_that("a drifting cumulative effect reads as out of date", {
  # cumulative effect trends steadily: unstable by construction. It is also
  # sufficient (index > 1), which is required by the source's criterion for
  # an out-of-date verdict (sufficient AND unstable).
  yi <- seq(log(0.9), log(0.2), length.out = 20); vi <- rep(0.01, 20)
  prev <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  v <- sufficiency_changepoint(prev, prev)
  expect_true(v$detail$sufficient)
  expect_false(v$detail$stable)
  expect_equal(v$verdict, "out_of_date")
})

test_that("the sufficiency index is Nfs / (5k + 10)", {
  yi <- rep(log(0.5), 10); vi <- rep(0.04, 10)
  prev <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  v <- sufficiency_changepoint(prev, prev)
  expect_equal(v$detail$index, failsafe_n(yi, vi) / (5 * 10 + 10),
               tolerance = 1e-8)
})

test_that("insufficient evidence never reads as out of date, regardless of stability", {
  # Confirmed empirically in the primary source: all 80 reviews in Pattanittum
  # et al. (2012) had a failsafe ratio below 1, and precisely because of that
  # none were flagged out-of-date by this method. Insufficiency alone must
  # never trigger "out_of_date".
  yi <- rep(log(0.99), 6); vi <- rep(0.04, 6)  # tiny effect: low index
  prev <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  v <- sufficiency_changepoint(prev, prev)
  expect_false(v$detail$sufficient)
  expect_equal(v$verdict, "current")
})

test_that("sufficiency_changepoint is read off prev, not new_ma; stability off new_ma", {
  # Per Table 1 of the primary source (Pattanittum et al. 2012), Nfs and its
  # k in "5k+10" are "of the previous meta-analysis"; the stability slope is
  # "of the updated meta-analysis". This is the case where that distinction is
  # load-bearing: prev alone is weak (few studies, small effect) and reads as
  # insufficient on its own, but folding in a large, strongly-effective,
  # trending batch of new studies would look sufficient AND unstable if
  # (wrongly) computed from new_ma throughout -- which would read as
  # "out_of_date". Reading sufficiency from prev alone, as the source
  # specifies, must keep this "current": insufficiency at the time of
  # publication blocks the out-of-date verdict regardless of what the new
  # studies later show.
  prev_yi <- rep(log(0.85), 4); prev_vi <- rep(0.05, 4)
  prev <- metafor::rma(yi = prev_yi, vi = prev_vi, measure = "RR")

  new_yi <- seq(log(0.5), log(0.05), length.out = 16)  # strong, drifting
  new_vi <- rep(0.02, 16)
  new_ma <- metafor::rma(yi = c(prev_yi, new_yi), vi = c(prev_vi, new_vi),
                          measure = "RR")

  # Sanity check: computed the wrong way (off new_ma throughout), this
  # fixture would be sufficient and unstable -- i.e. "out_of_date" under the
  # bug this test guards against.
  wrong_index <- failsafe_n(as.numeric(new_ma$yi), as.numeric(new_ma$vi)) /
    (5 * new_ma$k + 10)
  expect_gt(wrong_index, 1)

  v <- sufficiency_changepoint(prev, new_ma)
  expect_false(v$detail$sufficient)     # prev alone: index <= 1
  expect_equal(v$detail$k, prev$k)      # sufficiency's k is prev's k
  expect_equal(v$detail$k_new, new_ma$k) # stability's k is new_ma's k
  expect_equal(v$verdict, "current")    # insufficiency blocks out_of_date
})

# --- The stability test ----------------------------------------------------
#
# Three implementations, each replaced after measurement; see ?sufficiency_changepoint for
# the full numbers and vignette("methods") for the declaration.
#
#   1. summary.lm() t-test on the slope of the cumulative effect against study
#      index. No valid null distribution -- a cumulative mean is autocorrelated
#      by construction and converges by the law of large numbers, so the test
#      detects convergence and reports instability. 209 firings on 300 samples
#      of unchanging evidence, against 0 for rcma and 0 for ottawa.
#   2. Permutation over study order, same slope statistic. Calibrated on
#      average (16/300) but blind to late change -- power 1/200 against 10 new
#      studies at RR 0.30, and 0 as the shift got larger -- and its null is not
#      valid when study variances change over calendar time (83/300 = 28% false
#      alarms with no drift at all).
#   3. Permutation over study order, statistic replaced with the largest
#      standardised movement left in the cumulative series. Nominal everywhere
#      measured and full power against late change.
#
# The fixtures below pin down each of those failures so none can come back.

test_that("p_stability is a permutation p-value, and the permutation is deterministic", {
  yi <- seq(log(0.9), log(0.2), length.out = 20); vi <- rep(0.01, 20)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  a <- sufficiency_changepoint(ma, ma)
  b <- sufficiency_changepoint(ma, ma)
  expect_equal(a$detail$p_stability, b$detail$p_stability)  # same answer every call
  # A permutation p-value with 999 draws lives on the grid (1 + j) / 1000.
  expect_equal(a$detail$p_stability * 1000, round(a$detail$p_stability * 1000))
  expect_gte(a$detail$p_stability, 1 / 1000)            # never exactly zero
  expect_lte(a$detail$p_stability, 1)
  # A monotone drift is the case the test must catch.
  expect_false(a$detail$stable)
  expect_equal(a$verdict, "out_of_date")
})

test_that("alpha_stability is the permutation p-value cutoff", {
  set.seed(8)
  yi <- rnorm(30, log(0.5), 0.02); vi <- rep(0.01, 30)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  p <- sufficiency_changepoint(ma, ma)$detail$p_stability
  expect_true(sufficiency_changepoint(ma, ma, alpha_stability = p - 1e-9)$detail$stable)
  expect_false(sufficiency_changepoint(ma, ma, alpha_stability = p + 1e-9)$detail$stable)
})

test_that("sufficiency_changepoint leaves the caller's random stream untouched", {
  # The permutation draw must not be visible to the caller: a script that seeds
  # once and then backtests would otherwise lose reproducibility downstream.
  ma <- metafor::rma(yi = seq(log(0.9), log(0.2), length.out = 20),
                     vi = rep(0.01, 20), measure = "RR")
  set.seed(123); expected <- runif(3)
  set.seed(123); invisible(sufficiency_changepoint(ma, ma)); got <- runif(3)
  expect_equal(got, expected)
})

test_that("detail still carries every documented field", {
  ma <- metafor::rma(yi = seq(log(0.9), log(0.2), length.out = 20),
                     vi = rep(0.01, 20), measure = "RR")
  d <- sufficiency_changepoint(ma, ma)$detail
  expect_setequal(names(d), c("index", "sufficient", "stable", "slope",
                              "z_shift", "split", "p_stability", "k", "k_new"))
  expect_true(is.numeric(d$slope) && is.finite(d$slope))
  expect_true(is.numeric(d$z_shift) && is.finite(d$z_shift) && d$z_shift >= 0)
  expect_true(d$split %in% seq_len(ma$k - 1))
})

# --- The late-drift regime, which is the whole point of the package --------
#
# A mature review updated with a smaller batch of new studies is this package's
# canonical case, and it is exactly the case implementation 2 could not see:
# its statistic, the slope of the cumulative series, is dominated by the first
# few points where almost no information has accumulated and the series swings
# hardest, so a change confined to the tail never cleared the permutation null.
# The suite up to this point only ever exercised a fully monotone trend, which
# is the regime that statistic handles best -- which is why the blind spot
# shipped. These fixtures are the ones that would have caught it.

test_that("a shift confined to the last studies reads as out of date", {
  # 15 studies at RR 0.5, then 5 at RR 0.05. Under implementation 2 this
  # returned "current": observed |slope| 0.0247 against a null whose median was
  # 0.022, because a permuted order that happens to start on one of the RR 0.05
  # studies swings the cumulative series far harder at low index than a genuine
  # late cliff ever moves it at high index.
  yi <- c(rep(log(0.5), 15), rep(log(0.05), 5)); vi <- rep(0.01, 20)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  v <- sufficiency_changepoint(ma, ma)
  expect_true(v$detail$sufficient)
  expect_false(v$detail$stable)
  expect_equal(v$verdict, "out_of_date")
  # And it says where: the estimate was still moving after study 15.
  expect_equal(v$detail$split, 15L)
})

test_that("the shift is caught wherever it falls in the series, early or late", {
  # The shift-position scan from the review, as a fixture. Implementation 2
  # fired at 2, 5, 8 and 10 and was silent at 12, 15 and 18 -- the later the
  # change, the more relevant it is to a staleness question and the less this
  # detector could see it.
  fired <- vapply(c(2, 5, 8, 10, 12, 15, 18), function(s) {
    yi <- c(rep(log(0.5), s), rep(log(0.05), 20 - s)); vi <- rep(0.01, 20)
    ma <- metafor::rma(yi = yi, vi = vi, measure = "RR")
    sufficiency_changepoint(ma, ma)$verdict == "out_of_date"
  }, logical(1))
  expect_true(all(fired))
})

test_that("power against a smaller late batch is not merely non-zero", {
  # 20 prior studies at RR 0.5 plus 10 new at RR 0.30. Implementation 2 caught
  # 1 of 200 such samples, and 0 of 200 once the new studies moved to RR 0.15
  # or beyond -- power that *fell* as the change grew. Ten samples here, all of
  # which must be caught; the full 200-sample measurement is in ?sufficiency_changepoint.
  caught <- vapply(1:10, function(i) {
    set.seed(5000 + i)
    yi <- rnorm(20, log(0.5), 0.05); vi <- rep(0.01, 20)
    more <- rnorm(10, log(0.30), 0.05)
    prev <- metafor::rma(yi = yi, vi = vi, measure = "RR")
    upd  <- metafor::rma(yi = c(yi, more), vi = c(vi, rep(0.01, 10)),
                         measure = "RR")
    sufficiency_changepoint(prev, upd)$verdict == "out_of_date"
  }, logical(1))
  expect_true(all(caught))
})

test_that("the reported slope is fitted against accumulated information", {
  # Pattanittum et al. (2012, Table 1) specifies the regression as being
  # "versus information increment". An earlier round of this package noticed
  # that the code used the study index instead and resolved it by rewriting the
  # documentation to match the code. This test resolves it the other way round.
  # It needs unequal variances to have any teeth: when every vi is the same,
  # cumsum(1 / vi) is an affine function of the index and the two slopes differ
  # only by a constant factor.
  yi <- seq(log(0.9), log(0.2), length.out = 12)
  vi <- seq(0.20, 0.01, length.out = 12)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR", method = "FE")
  cum <- cumulative_effect(yi, vi)
  y <- cum[-1]
  x_info  <- cumsum(1 / vi)[-1]
  x_index <- seq_along(y)
  ols <- function(x, y) {
    xc <- x - mean(x); sum(xc * (y - mean(y))) / sum(xc^2)
  }
  expect_equal(sufficiency_changepoint(ma, ma)$detail$slope, ols(x_info, y),
               tolerance = 1e-10)
  # The two really are different numbers here, so the assertion above is not
  # satisfied by the index version.
  expect_false(isTRUE(all.equal(ols(x_info, y), ols(x_index, y))))
})

test_that("the stability statistic is pivotal in the noise scale", {
  # This is the property that makes the permutation null valid when study sizes
  # change systematically over calendar time: the statistic is a ratio of the
  # movement to its own standard error, so inflating the noise and the declared
  # variances together -- which is exactly what a smaller study is -- leaves it
  # untouched. The published slope, by contrast, carries the units and moves by
  # a^3 under the same transformation, which is how the variance schedule got
  # imprinted on it and why implementation 2 false-alarmed at 28%.
  set.seed(21)
  vi <- seq(0.15, 0.01, length.out = 15)
  th <- log(0.5)
  yi <- th + rnorm(15, 0, sqrt(vi))
  a  <- 3
  yi_a <- th + a * (yi - th); vi_a <- vi * a^2

  expect_equal(stability_shift_z(yi, vi), stability_shift_z(yi_a, vi_a),
               tolerance = 1e-10)
  s1 <- cum_drift_slope(cumulative_effect(yi, vi), cumsum(1 / vi))
  s2 <- cum_drift_slope(cumulative_effect(yi_a, vi_a), cumsum(1 / vi_a))
  expect_false(isTRUE(all.equal(s1, s2)))
  expect_equal(s2 / s1, a^3, tolerance = 1e-8)
})

# --- C2: an exactly constant cumulative series must not crash the caller ----
#
# When cum_theta is constant in floating point, summary.lm() returns
# Pr(>|t|) = NaN, so `stable <- NaN >= 0.05` is NA and `if (sufficient &&
# !stable)` became `if (NA)` -> "missing value where TRUE/FALSE needed". That
# error propagated out of backtest() and killed the whole run. Both
# reproductions below come straight from the review.

test_that("a byte-identical body of evidence returns a verdict instead of erroring", {
  # Reproduction 2 from the review, at exactly min_k.
  prev <- metafor::rma(yi = rep(log(0.4), 5), vi = rep(0.01, 5), measure = "RR")
  v <- expect_no_error(sufficiency_changepoint(prev, prev))
  expect_true(v$verdict %in% c("current", "out_of_date"))
  expect_true(v$detail$stable)            # a constant series is maximally stable
  expect_equal(v$verdict, "current")
})

test_that("a constant cumulative series does not take down a whole backtest", {
  # Reproduction 1 from the review: this exact call died before the fix, and
  # succeeded as soon as "sufficiency_changepoint" was dropped from `methods`.
  stream <- evidence_stream(
    metafor::rma(yi = rep(0.5, 20), vi = rep(0.02, 20), measure = "MD"),
    date = 2000:2019, study_id = seq_along(2000:2019))
  bt <- expect_no_error(
    backtest(stream, methods = c("rcma", "ottawa", "sufficiency_changepoint"),
             horizon = 3, window = 3))
  expect_s3_class(bt, "staleness_backtest")
  suff <- bt$results[bt$results$method == "sufficiency_changepoint", ]
  expect_true(all(suff$verdict %in% c("current", "not_applicable")))
})

test_that("a split that carries no information does not swallow the statistic", {
  # A study with infinite variance contributes no weight. If it sits last, the
  # "after" block of the final split holds no information and its Z is 0/0 --
  # one NaN that would otherwise take the maximum down with it and leave
  # `stable <- NaN >= alpha` to resolve. The undefined split is dropped and the
  # remaining ones still answer.
  yi <- c(rep(log(0.5), 6), rep(log(0.05), 3), log(0.5))
  vi <- c(rep(0.01, 9), Inf)
  expect_true(any(!is.finite(split_z(yi, vi))))   # the fixture really is degenerate
  z <- stability_shift_z(yi, vi)
  expect_true(is.finite(z) && z > 0)
  expect_equal(stability_shift_at(yi, vi), 6L)

  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR", method = "FE")
  v <- expect_no_error(sufficiency_changepoint(ma, ma))
  expect_true(v$verdict %in% c("current", "out_of_date"))
  expect_true(is.logical(v$detail$stable) && !is.na(v$detail$stable))
})

test_that("evidence with no usable split reads as stable rather than unstable", {
  # Every study after the first weightless: no split carries information, so
  # there is nothing to test. The statistic must say NA, and NA must resolve to
  # "stable" -- if it were left to the permutation p-value, `abs(perm) >= NA` is
  # all NA, `na.rm = TRUE` counts zero exceedances, and p = 1/1000 would report
  # "unstable" from an absence of evidence.
  yi <- c(log(0.5), rep(log(0.05), 5)); vi <- c(0.01, rep(Inf, 5))
  expect_true(all(!is.finite(split_z(yi, vi))))
  expect_true(is.na(stability_shift_z(yi, vi)))
  expect_true(is.na(stability_shift_at(yi, vi)))
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR", method = "FE")
  v <- expect_no_error(sufficiency_changepoint(ma, ma))
  expect_true(v$detail$stable)
  expect_equal(v$verdict, "current")
})

test_that("a non-finite cumulative effect is declined, not guessed at", {
  # A zero study variance makes the weights infinite and the cumulative series
  # NaN. Whatever else happens, `if (NA)` must never be reached.
  prev <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6), measure = "RR")
  broken <- prev
  broken$vi <- c(0, rep(0.02, 5))
  v <- expect_no_error(sufficiency_changepoint(prev, broken))
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "not finite")
})

test_that("a flat cumulative tail does not hide a split at the first study", {
  # cum_theta[-1] is exactly constant here, which is what the OLD statistic
  # was fitted to -- but the change-point statistic reads splits, and the one
  # isolating study 1 is enormous. Keying the degenerate short circuit to the
  # cumulative series would report perfect stability over a 10-SE split.
  yi <- c(10, -10, rep(0, 28))
  vi <- rep(1, 30)
  expect_lt(diff(range(cumulative_effect(yi, vi)[-1])), 1e-12)  # the trap
  expect_gt(stability_shift_z(yi, vi), 10)                      # the truth

  prev <- metafor::rma(yi = rep(0, 5), vi = rep(1, 5))
  new  <- metafor::rma(yi = yi, vi = vi)
  res  <- sufficiency_changepoint(prev, new)
  expect_equal(res$detail$z_shift, stability_shift_z(yi, vi), tolerance = 1e-8)
  expect_false(is.na(res$detail$split))
  expect_lt(res$detail$p_stability, 1)
})

test_that("byte-identical studies are still maximally stable, with no p-value", {
  # The case the short circuit exists for: rounding noise must never be
  # tested. Its spread in yi is zero, so the statistic has nothing to read.
  yi <- rep(log(0.5), 12)
  vi <- rep(0.02, 12)
  prev <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6))
  new  <- metafor::rma(yi = yi, vi = vi)
  res  <- sufficiency_changepoint(prev, new)
  expect_true(res$detail$stable)
  expect_equal(res$detail$z_shift, 0)
  expect_equal(res$detail$p_stability, 1)
})

test_that("a single new study resolves cleanly instead of warning three times", {
  prev <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6))
  new  <- metafor::rma(yi = c(log(0.5), log(0.6)), vi = c(0.02, 0.02))
  # min_k = 1 is user-settable, and k_new = 1 used to reach max(numeric(0)).
  one <- metafor::rma(yi = rep(log(0.5), 2), vi = rep(0.02, 2))
  expect_no_warning(sufficiency_changepoint(prev, one, min_k = 1))
})

test_that("scale pivotality is about yi and vi TOGETHER, not vi alone", {
  yi <- c(1, 2, 3)
  base <- stability_shift_z(yi, rep(1, 3))

  # The property the statistic actually has: rescale the deviations and the
  # variances together and nothing moves.
  a <- 3; theta <- 2
  expect_equal(stability_shift_z(theta + a * (yi - theta), rep(1, 3) * a^2),
               base)

  # The property it does NOT have, asserted here so the documentation cannot
  # drift back to claiming it: vi alone divides every Z_m by sqrt(c).
  expect_equal(stability_shift_z(yi, rep(4, 3)), base / 2)
  expect_false(isTRUE(all.equal(stability_shift_z(yi, rep(4, 3)), base)))
})

test_that("failsafe_n on empty input is undefined, not zero", {
  # (sum(z)^2)/z_crit^2 - length(z) evaluates to 0 on empty input, and 0 is a
  # meaningful fail-safe N: it flows into index = 0/(5*0+10) = 0 and reports
  # sufficient = FALSE. A quiet, wrong "not sufficient" beats no answer only
  # if you never look.
  expect_true(is.na(failsafe_n(numeric(0), numeric(0))))
  expect_false(is.na(failsafe_n(c(0.5, 0.4), c(0.05, 0.05))))
})

test_that("min_k gates the updated evidence only, and that is load-bearing", {
  # The two halves read different objects: sufficiency from `prev`, stability
  # from `new_ma`. min_k guards the stability permutation, which needs enough
  # studies to have a null at all, so it is checked against new_ma. Nothing
  # guards `prev`, which means a two-study prior review can drive the
  # sufficiency half. Pinned so the asymmetry stays a documented decision.
  prev2 <- metafor::rma(yi = c(0.5, 0.4), vi = c(0.05, 0.05))
  new8  <- metafor::rma(yi = rep(0.5, 8), vi = rep(0.05, 8))
  res <- sufficiency_changepoint(prev2, new8)
  expect_equal(res$detail$k, 2)       # the prior review really is that small
  expect_false(is.na(res$detail$index))
  expect_equal(res$verdict, "current")

  # And the gate does fire on the updated side.
  new3 <- metafor::rma(yi = rep(0.5, 3), vi = rep(0.05, 3))
  expect_equal(sufficiency_changepoint(prev2, new3)$verdict, "not_applicable")
})

test_that("an undeterminable index is resolved, not left to break the if()", {
  # The code already carries this lesson in a comment, for p_stability: a
  # non-finite value "must be resolved deliberately, never left to turn
  # `if (sufficient && !stable)` into `if (NA)`". The same reasoning was never
  # applied to `sufficient` itself, one line above. With a degenerate prior
  # the fail-safe N is NaN, so sufficient is NA; NA && TRUE is NA, and the
  # detector died with "missing value where TRUE/FALSE needed" -- but only
  # when the evidence also happened to be unstable, which is why it survived.
  prev <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6))
  prev$yi <- rep(0, 6)
  prev$vi <- rep(0, 6)                       # fail-safe N -> NaN
  unstable <- metafor::rma(yi = c(rep(log(0.5), 6), rep(log(1.4), 6)),
                           vi = rep(0.02, 12))

  res <- sufficiency_changepoint(prev, unstable)
  expect_equal(res$verdict, "not_applicable")
  expect_match(res$reason, "sufficiency index is not finite")

  # And a stable one takes the same path, rather than reaching "current" by
  # the accident of NA && FALSE being FALSE.
  stable <- metafor::rma(yi = rep(log(0.5), 12), vi = rep(0.02, 12))
  expect_equal(sufficiency_changepoint(prev, stable)$verdict, "not_applicable")
})
