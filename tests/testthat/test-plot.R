test_that("plotting a backtest runs and returns the object invisibly", {
  skip_if_not_installed("metadat")
    es <- bcg_es()
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
  # calibration() now returns one row per requested method with NA metrics
  # rather than NULL (see ?calibration), so the all-NA matrix must not reach
  # graphics::barplot() either: it would draw an empty pair of axes that looks
  # like a result. Either way the caller is told what actually went wrong.
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

  # A null device, so a plot() that got as far as drawing cannot leave an
  # Rplots.pdf behind in the test directory.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(
    plot(bt, truth = "shift"),
    "nothing to plot for truth = \"shift\""
  )
})

test_that("a method with no eligible rows is drawn as a gap, not dropped or fatal", {
  # calibration() emits an NA row for a never-applicable detector; barplot()
  # must tolerate it (an empty slot in the chart) rather than error, so that a
  # backtest where four detectors answered and one never did still plots.
  bt <- structure(
    list(
      results = data.frame(
        cut     = c(2000, 2001, 2000, 2001),
        method  = rep(c("rcma", "barrowman"), each = 2),
        verdict = c("out_of_date", "current", "not_applicable", "not_applicable"),
        signal  = NA_real_,
        reason  = "",
        truth_shift      = c(TRUE, FALSE, TRUE, FALSE),
        truth_surprise   = rep(FALSE, 4),
        truth_conclusion = rep(TRUE, 4),
        censored = rep(FALSE, 4),
        stringsAsFactors = FALSE
      ),
      methods = c("rcma", "barrowman"), horizon = 5, window = 3,
      n_cuts = 2, n_censored = 0
    ),
    class = "staleness_backtest"
  )

  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(bt, truth = "shift"))
})

test_that("graphical arguments the caller supplies override the defaults", {
  # ?plot.staleness_backtest says `...` is passed to barplot(), but main, ylab
  # and ylim were also set inside the call, so supplying any of them failed
  # with "formal argument matched by multiple actual arguments" -- a message
  # about R's argument matching, not about anything the caller did wrong.
  st <- evidence_stream(
    suppressWarnings(metafor::rma(yi = seq(-1, 1, length.out = 12),
                                  vi = rep(0.05, 12), measure = "RR")),
    date = 2000:2011, ni = rep(100, 12))
  bt <- backtest(st, cuts = 2003:2008, horizon = 2, window = 2, seed = 1)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_silent(plot(bt, main = "my title"))
  expect_silent(plot(bt, ylab = "my label"))
  expect_silent(plot(bt, ylim = c(0, 2)))
  expect_silent(plot(bt, col = "red"))
  # And the defaults still apply when nothing is supplied.
  expect_silent(plot(bt))
})
