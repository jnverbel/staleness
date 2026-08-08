test_that("failsafe_n follows Rosenthal's formula", {
  yi <- rep(log(0.5), 5); vi <- rep(0.04, 5)
  z  <- yi / sqrt(vi)
  expected <- (sum(z)^2) / (1.645^2) - length(z)
  expect_equal(failsafe_n(yi, vi), expected, tolerance = 1e-8)
})

test_that("sufficiency refuses fewer than five studies", {
  prev <- metafor::rma(yi = rep(log(0.5), 4), vi = rep(0.04, 4), measure = "RR")
  v <- sufficiency(prev, prev)
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
  v <- sufficiency(prev, prev)
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
  v <- sufficiency(prev, prev)
  expect_true(v$detail$sufficient)
  expect_false(v$detail$stable)
  expect_equal(v$verdict, "out_of_date")
})

test_that("the sufficiency index is Nfs / (5k + 10)", {
  yi <- rep(log(0.5), 10); vi <- rep(0.04, 10)
  prev <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  v <- sufficiency(prev, prev)
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
  v <- sufficiency(prev, prev)
  expect_false(v$detail$sufficient)
  expect_equal(v$verdict, "current")
})

test_that("sufficiency is read off prev, not new_ma; stability off new_ma", {
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

  v <- sufficiency(prev, new_ma)
  expect_false(v$detail$sufficient)     # prev alone: index <= 1
  expect_equal(v$detail$k, prev$k)      # sufficiency's k is prev's k
  expect_equal(v$detail$k_new, new_ma$k) # stability's k is new_ma's k
  expect_equal(v$verdict, "current")    # insufficiency blocks out_of_date
})

# --- The stability test: permutation over study order ----------------------
#
# The first implementation tested the slope of the cumulative pooled effect
# with the t-test from summary.lm(). A cumulative mean is near-perfectly
# autocorrelated by construction and drifts toward the pooled effect by the law
# of large numbers alone, so that test has no valid null distribution: it
# detects convergence and reports instability. Measured over 300 samples of
# genuinely unchanging evidence it returned out_of_date 209 times, against 0
# for rcma and 0 for ottawa. Replaced by a permutation test over study order,
# which assumes nothing about independence and asks the question the method
# actually means: is this drift larger than the same studies in a random order
# would produce? See ?sufficiency and vignette("methods") for the declaration.

test_that("p_slope is a permutation p-value, and the permutation is deterministic", {
  yi <- seq(log(0.9), log(0.2), length.out = 20); vi <- rep(0.01, 20)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  a <- sufficiency(ma, ma)
  b <- sufficiency(ma, ma)
  expect_equal(a$detail$p_slope, b$detail$p_slope)  # same answer every call
  # A permutation p-value with 999 draws lives on the grid (1 + j) / 1000.
  expect_equal(a$detail$p_slope * 1000, round(a$detail$p_slope * 1000))
  expect_gte(a$detail$p_slope, 1 / 1000)            # never exactly zero
  expect_lte(a$detail$p_slope, 1)
  # A monotone drift is the case the test must catch.
  expect_false(a$detail$stable)
  expect_equal(a$verdict, "out_of_date")
})

test_that("alpha_slope is the permutation p-value cutoff", {
  set.seed(8)
  yi <- rnorm(30, log(0.5), 0.02); vi <- rep(0.01, 30)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "RR")
  p <- sufficiency(ma, ma)$detail$p_slope
  expect_true(sufficiency(ma, ma, alpha_slope = p - 1e-9)$detail$stable)
  expect_false(sufficiency(ma, ma, alpha_slope = p + 1e-9)$detail$stable)
})

test_that("sufficiency leaves the caller's random stream untouched", {
  # The permutation draw must not be visible to the caller: a script that seeds
  # once and then backtests would otherwise lose reproducibility downstream.
  ma <- metafor::rma(yi = seq(log(0.9), log(0.2), length.out = 20),
                     vi = rep(0.01, 20), measure = "RR")
  set.seed(123); expected <- runif(3)
  set.seed(123); invisible(sufficiency(ma, ma)); got <- runif(3)
  expect_equal(got, expected)
})

test_that("detail still carries every documented field", {
  ma <- metafor::rma(yi = seq(log(0.9), log(0.2), length.out = 20),
                     vi = rep(0.01, 20), measure = "RR")
  d <- sufficiency(ma, ma)$detail
  expect_setequal(names(d), c("index", "sufficient", "stable", "slope",
                              "p_slope", "k", "k_new"))
  expect_true(is.numeric(d$slope) && is.finite(d$slope))
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
  v <- expect_no_error(sufficiency(prev, prev))
  expect_true(v$verdict %in% c("current", "out_of_date"))
  expect_true(v$detail$stable)            # a constant series is maximally stable
  expect_equal(v$verdict, "current")
})

test_that("a constant cumulative series does not take down a whole backtest", {
  # Reproduction 1 from the review: this exact call died before the fix, and
  # succeeded as soon as "sufficiency" was dropped from `methods`.
  stream <- evidence_stream(
    metafor::rma(yi = rep(0.5, 20), vi = rep(0.02, 20), measure = "MD"),
    date = 2000:2019)
  bt <- expect_no_error(
    backtest(stream, methods = c("rcma", "ottawa", "sufficiency"),
             horizon = 3, window = 3))
  expect_s3_class(bt, "staleness_backtest")
  suff <- bt$results[bt$results$method == "sufficiency", ]
  expect_true(all(suff$verdict %in% c("current", "not_applicable")))
})

test_that("a non-finite cumulative effect is declined, not guessed at", {
  # A zero study variance makes the weights infinite and the cumulative series
  # NaN. Whatever else happens, `if (NA)` must never be reached.
  prev <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6), measure = "RR")
  broken <- prev
  broken$vi <- c(0, rep(0.02, 5))
  v <- expect_no_error(sufficiency(prev, broken))
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "not finite")
})
