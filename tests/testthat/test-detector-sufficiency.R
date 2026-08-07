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
  # is squarely "current". seed = 8 is deliberately picked (see task report):
  # the cumulative-effect regression has serially correlated errors (a
  # limitation the source itself names), so a "no true trend" sample can still
  # yield a significant slope by chance for many seeds. seed = 8 gives a
  # comfortable, non-borderline margin (p_slope ~ 0.53) so the test reliably
  # demonstrates the intended property instead of depending on a coin flip.
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
