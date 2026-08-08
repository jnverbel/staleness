# Anti-look-ahead guarantee.
#
# The entire scientific value of backtest() rests on one property: at a cut
# in the past, every detector must see only what was known at that time. If
# information from the future leaks in (the classic mistake is estimating
# tau-squared on the full series), every calibration number the package
# produces is worthless. This file is the adversarial test for that property.

test_that("no detector sees a step change before it happens", {
  # 12 studies with effect 0.0, then 12 with effect 0.8, from 2010 onwards.
  # Any signal at a cut before 2010 is look-ahead leakage.
  set.seed(99)
  yi <- c(rnorm(12, 0.0, 0.02), rnorm(12, 0.8, 0.02))
  vi <- rep(0.01, 24)
  dates <- c(1998:2009, 2010:2021)

  ma <- metafor::rma(yi = yi, vi = vi, measure = "MD")
  s  <- evidence_stream(ma, date = dates)

  bt <- backtest(s, methods = c("rcma", "ottawa", "sufficiency"),
                 horizon = 3, window = 3, seed = 1)
  before <- bt$results[bt$results$cut < 2007, ]

  # cuts before 2007 have a window ending in 2009 at the latest: pure pre-change
  # evidence. No detector may call out_of_date there.
  expect_false(any(before$verdict == "out_of_date"),
               info = paste("leakage at cuts:",
                            paste(unique(before$cut[before$verdict == "out_of_date"]),
                                  collapse = ", ")))
})

test_that("the same detectors do fire after the change, proving the test has teeth", {
  set.seed(99)
  yi <- c(rnorm(12, 0.0, 0.02), rnorm(12, 0.8, 0.02))
  vi <- rep(0.01, 24)
  dates <- c(1998:2009, 2010:2021)

  ma <- metafor::rma(yi = yi, vi = vi, measure = "MD")
  s  <- evidence_stream(ma, date = dates)

  bt <- backtest(s, methods = c("rcma", "ottawa"), horizon = 3, window = 3, seed = 1)
  after <- bt$results[bt$results$cut >= 2010, ]
  expect_true(any(after$verdict == "out_of_date"))
})

test_that("snapshot_at never uses a study published after the cut", {
  set.seed(7)
  yi <- rnorm(20, 0.3, 0.05); vi <- rep(0.01, 20)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "MD")
  s  <- evidence_stream(ma, date = 2000:2019)

  for (cut in c(2005, 2010, 2015)) {
    snap <- snapshot_at(s, cut)
    expect_equal(snap$k, sum(s$date <= cut))
    expect_setequal(round(as.numeric(snap$yi), 10),
                    round(s$yi[s$date <= cut], 10))
  }
})
