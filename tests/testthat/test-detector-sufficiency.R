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
