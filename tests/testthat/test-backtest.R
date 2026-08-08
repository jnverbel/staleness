# Tests for backtest(): the retro-validation engine.
#
# Walks a historical evidence stream, reconstructs what was known at each cut,
# asks every detector the same question it would be asked today, and compares
# the answer against what actually happened afterwards.

make_bcg_stream <- function() {
  dat <- metadat::dat.bcg
  es  <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat)
  ma  <- metafor::rma(yi, vi, data = es)
  evidence_stream(ma, date = es$year)
}

test_that("backtest returns one row per usable cut and method", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  bt <- backtest(s, methods = c("rcma", "ottawa"), horizon = 5, seed = 1)
  expect_s3_class(bt, "staleness_backtest")
  expect_true(all(c("cut", "method", "verdict",
                    "truth_shift", "truth_surprise", "truth_conclusion",
                    "censored") %in% names(bt$results)))
  expect_setequal(unique(bt$results$method), c("rcma", "ottawa"))
})

test_that("cuts too close to the end of the series are censored, not counted", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  bt <- backtest(s, methods = "rcma", horizon = 5, seed = 1)
  last_ok <- max(s$date) - 5
  expect_true(all(bt$results$censored[bt$results$cut > last_ok]))
  expect_true(all(!bt$results$censored[bt$results$cut <= last_ok]))
})

test_that("backtest refuses a series with fewer than three usable cuts", {
  ma <- metafor::rma(yi = rep(log(0.5), 3), vi = rep(0.04, 3), measure = "RR")
  s  <- evidence_stream(ma, date = c(2000, 2001, 2002))
  expect_error(backtest(s, horizon = 5), "at least 3")
})

test_that("truth columns are logical and computed against the full series", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  bt <- backtest(s, methods = "rcma", horizon = 5, seed = 1)
  expect_type(bt$results$truth_shift, "logical")
  expect_type(bt$results$truth_surprise, "logical")
  expect_type(bt$results$truth_conclusion, "logical")
})

# --- ni correction (binding correction to the brief) ------------------------
#
# The brief's backtest() calls check_currency() without n_prev/n_new, so
# barrowman() would return not_applicable at every single cut regardless of
# whether sample sizes were ever available. The stream already carries `ni`,
# and window_between() already returns it, so backtest() must sum it on both
# sides of the cut and thread it through. When the stream has no `ni`,
# check_currency() must still degrade barrowman() to not_applicable without
# erroring.

make_ni_stream <- function(with_ni = TRUE) {
  yi   <- c(0.05, 0.08, -0.02, 0.10, 0.04, 0.12, 0.15, 0.20, 0.18, 0.25,
            0.30, 0.35, 0.40, 0.45, 0.50)
  vi   <- rep(0.05, length(yi))
  date <- 2000:2014
  ma   <- metafor::rma(yi = yi, vi = vi, method = "REML")
  if (with_ni) {
    # Small early studies, large later studies: n_new eventually exceeds the
    # sample size barrowman requires, producing real out_of_date/current
    # verdicts rather than a wall of not_applicable.
    ni <- c(rep(50, 8), rep(2000, 7))
    evidence_stream(ma, date = date, ni = ni)
  } else {
    evidence_stream(ma, date = date)
  }
}

test_that("backtest threads n_prev/n_new through so barrowman produces real verdicts", {
  s  <- make_ni_stream(with_ni = TRUE)
  bt <- backtest(s, methods = "barrowman", horizon = 5, seed = 1)
  decided <- bt$results$verdict[bt$results$verdict != "not_applicable"]
  expect_true(length(decided) > 0)
  expect_false(any(grepl("sample size not supplied", bt$results$reason)))
})

test_that("backtest degrades barrowman to not_applicable, without erroring, when the stream has no ni", {
  s <- make_ni_stream(with_ni = FALSE)
  expect_null(s$ni)
  bt <- backtest(s, methods = "barrowman", horizon = 5, seed = 1)
  expect_true(all(bt$results$verdict == "not_applicable"))
  expect_true(all(grepl("sample size not supplied", bt$results$reason)))
})

# --- window > horizon must censor too --------------------------------------
#
# Censoring depended only on `cut > max(dates) - horizon`, with nothing
# checking that the evidence window fitted inside the data. On dat.bcg with
# horizon = 3 and window = 5 -- the settings the backtesting vignette itself
# uses -- cuts 1976 and 1977 were marked UNCENSORED while seeing only 2 and 1
# new studies, because their five-year windows ran past 1980. They were then
# scored as full observations alongside cuts that got a complete window. A
# truncated window shows a detector less change than really happened, so those
# cuts are biased toward "current", invisibly. Fixed by censoring on
# max(horizon, window).

test_that("cuts whose evidence window runs past the end of the data are censored", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  bt <- backtest(s, methods = "rcma", horizon = 3, window = 5, seed = 1)
  last_date <- max(s$date)                       # 1980

  # The two cuts the review named. Under horizon alone they were uncensored.
  affected <- bt$results[bt$results$cut %in% c(1976, 1977), ]
  expect_equal(nrow(affected), 2)
  expect_true(all(affected$censored))

  # And they really do have a truncated window, which is the reason: their
  # five-year window ends after the last study in the stream, so they see less
  # new evidence than an earlier cut with the same nominal window.
  for (cut in c(1976, 1977)) {
    expect_gt(cut + 5, last_date)
    expect_lt(sum(s$date > cut & s$date <= cut + 5),
              sum(s$date > 1970 & s$date <= 1975))
  }

  # The general rule, not just those two cuts.
  expect_true(all(bt$results$censored[bt$results$cut >  last_date - 5]))
  expect_true(all(!bt$results$censored[bt$results$cut <= last_date - 5]))

  # Every uncensored cut sees a full window's worth of future.
  uncens <- unique(bt$results$cut[!bt$results$censored])
  expect_true(all(uncens + 5 <= last_date))
})

test_that("window <= horizon leaves the censoring boundary at horizon", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  bt <- backtest(s, methods = "rcma", horizon = 5, window = 3, seed = 1)
  expect_true(all(bt$results$censored[bt$results$cut >  max(s$date) - 5]))
  expect_true(all(!bt$results$censored[bt$results$cut <= max(s$date) - 5]))
})
