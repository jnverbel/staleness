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
})
