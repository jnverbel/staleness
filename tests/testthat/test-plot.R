test_that("plotting a backtest runs and returns the object invisibly", {
  skip_if_not_installed("metadat")
  dat <- metadat::dat.bcg
  es  <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat)
  ma  <- metafor::rma(yi, vi, data = es)
  s   <- evidence_stream(ma, date = es$year)
  bt  <- backtest(s, methods = c("rcma", "ottawa"), horizon = 5, seed = 1)

  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  out <- plot(bt)
  grDevices::dev.off()

  expect_true(file.exists(f))
  expect_s3_class(out, "staleness_backtest")
  # `plot()` must hand back the same backtest object it was given, unchanged
  # (sensitivity/specificity values themselves are covered by
  # test-metrics.R's calibration() tests, not re-verified here).
  expect_identical(out, bt)
})

test_that("plotting a backtest restores the caller's graphics parameters", {
  # A minimal hand-built backtest, so this check does not depend on metadat
  # or on running a real backtest: it only needs calibration() to succeed.
  bt <- structure(
    list(
      results = data.frame(
        cut     = c(2000, 2001, 2000, 2001),
        method  = rep(c("rcma", "ottawa"), each = 2),
        verdict = c("out_of_date", "current", "current", "current"),
        signal  = NA_real_,
        reason  = "",
        truth_shift      = c(TRUE, FALSE, TRUE, FALSE),
        truth_surprise   = rep(FALSE, 4),
        truth_conclusion = rep(TRUE, 4),
        censored = rep(FALSE, 4),
        stringsAsFactors = FALSE
      ),
      methods = c("rcma", "ottawa"), horizon = 5, window = 3,
      n_cuts = 2, n_censored = 0
    ),
    class = "staleness_backtest"
  )

  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit(grDevices::dev.off(), add = TRUE)

  before <- graphics::par(no.readonly = TRUE)
  plot(bt)
  after <- graphics::par(no.readonly = TRUE)

  expect_identical(before, after)
})

test_that("plotting a backtest with no eligible rows fails with a clear message", {
  # Every row is ineligible for calibration(), by one of its three exclusion
  # rules: not_applicable, censored, or an undeterminable (NA) truth value.
  # calibration() then returns NULL (do.call(rbind, list()) on zero groups),
  # which must not reach graphics::barplot() as-is -- that produces the
  # opaque "'height' must be a vector or a matrix" error instead of saying
  # what actually went wrong.
  bt <- structure(
    list(
      results = data.frame(
        cut     = c(2000, 2001, 2002),
        method  = c("rcma", "ottawa", "sufficiency"),
        verdict = c("not_applicable", "current", "current"),
        signal  = NA_real_,
        reason  = c("insufficient studies", "", ""),
        truth_shift      = c(TRUE, NA, TRUE),
        truth_surprise   = c(FALSE, FALSE, FALSE),
        truth_conclusion = c(TRUE, TRUE, TRUE),
        censored = c(FALSE, FALSE, TRUE),
        stringsAsFactors = FALSE
      ),
      methods = c("rcma", "ottawa", "sufficiency"), horizon = 5, window = 3,
      n_cuts = 3, n_censored = 1
    ),
    class = "staleness_backtest"
  )

  expect_error(
    plot(bt, truth = "shift"),
    "nothing to plot for truth = \"shift\""
  )
})
