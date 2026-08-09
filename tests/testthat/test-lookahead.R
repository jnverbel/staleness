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
  s  <- evidence_stream(ma, date = dates, study_id = seq_along(dates))

  bt <- backtest(s, methods = c("rcma", "ottawa", "sufficiency_changepoint"),
                 horizon = 3, window = 3, seed = 1)
  before <- bt$results[bt$results$cut < 2007, ]

  # Guard against passing vacuously: expect_false(any(character(0) == ...))
  # is also FALSE, so a future change to backtest()'s cut-selection (the
  # `usable` filter, or the min_k default) that silently emptied `before`
  # would leave the assertion below green while testing nothing. Anchor it.
  expect_true(nrow(before) > 0)

  # cuts before 2007 have a window ending in 2009 at the latest: pure pre-change
  # evidence. No detector may call out_of_date there.
  expect_false(any(before$verdict == "out_of_date"),
               info = paste("leakage at cuts:",
                            paste(unique(before$cut[before$verdict == "out_of_date"]),
                                  collapse = ", ")))
})

test_that("ottawa and sufficiency_changepoint do fire after the change, proving test 1 has teeth", {
  set.seed(99)
  yi <- c(rnorm(12, 0.0, 0.02), rnorm(12, 0.8, 0.02))
  vi <- rep(0.01, 24)
  dates <- c(1998:2009, 2010:2021)

  ma <- metafor::rma(yi = yi, vi = vi, measure = "MD")
  s  <- evidence_stream(ma, date = dates, study_id = seq_along(dates))

  bt <- backtest(s, methods = c("rcma", "ottawa", "sufficiency_changepoint"),
                 horizon = 3, window = 3, seed = 1)
  after <- bt$results[bt$results$cut >= 2010, ]
  expect_true(any(after$verdict == "out_of_date"))

  # Every method test 1's silence relies on must be shown capable of firing
  # here — otherwise a regression that made a detector always answer
  # "current" would sail through test 1 trivially and never be caught. ottawa
  # and sufficiency both genuinely reach out_of_date after the change; rcma's
  # own math (effect_ratio() in R/utils-measure.R) never lets its ratio cross
  # the 0.5/1.5 threshold for this particular synthetic series in this
  # seeded run, so it is not asserted here.
  expect_true(any(after$verdict[after$method == "ottawa"] == "out_of_date"))
  expect_true(any(after$verdict[after$method == "sufficiency_changepoint"] == "out_of_date"))
})

test_that("snapshot_at never uses a study published after the cut", {
  set.seed(7)
  yi <- rnorm(20, 0.3, 0.05); vi <- rep(0.01, 20)
  ma <- metafor::rma(yi = yi, vi = vi, measure = "MD")
  s  <- evidence_stream(ma, date = 2000:2019, study_id = seq_along(2000:2019))

  for (cut in c(2005, 2010, 2015)) {
    snap <- snapshot_at(s, cut)
    expect_equal(snap$k, sum(s$date <= cut))
    expect_equal(round(as.numeric(snap$yi), 10),
                 round(s$yi[s$date <= cut], 10))
  }
})
