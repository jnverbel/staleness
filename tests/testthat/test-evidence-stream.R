make_bcg_stream <- function() {
  dat <- metadat::dat.bcg
  es  <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat)
  ma  <- metafor::rma(yi, vi, data = es)
  evidence_stream(ma, date = es$year)
}

test_that("a stream is built from an rma object and sorted by date", {
  skip_if_not_installed("metadat")
  s <- make_bcg_stream()
  expect_s3_class(s, "staleness_stream")
  expect_equal(s$k, 13)
  expect_equal(s$measure, "RR")
  expect_false(is.unsorted(s$date))
})

test_that("missing dates are an explicit error, never imputed", {
  skip_if_not_installed("metadat")
  dat <- metadat::dat.bcg
  es  <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat)
  ma  <- metafor::rma(yi, vi, data = es)
  bad <- es$year; bad[3] <- NA
  expect_error(evidence_stream(ma, date = bad), "missing")
})

test_that("date length must match the number of studies", {
  skip_if_not_installed("metadat")
  dat <- metadat::dat.bcg
  es  <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg,
                         ci = cpos, di = cneg, data = dat)
  ma  <- metafor::rma(yi, vi, data = es)
  expect_error(evidence_stream(ma, date = c(1950, 1960)), "length")
})

test_that("snapshot_at refits using only studies up to the cut", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  m1 <- snapshot_at(s, 1960)
  expect_true(m1$k < s$k)
  expect_equal(m1$k, sum(s$date <= 1960))
  expect_equal(m1$measure, "RR")
})

test_that("window_between returns studies in the half-open interval", {
  skip_if_not_installed("metadat")
  s <- make_bcg_stream()
  w <- window_between(s, 1960, 1970)
  expect_equal(w$k, sum(s$date > 1960 & s$date <= 1970))
})

test_that("a snapshot with fewer than 2 studies is refused", {
  skip_if_not_installed("metadat")
  s <- make_bcg_stream()
  expect_error(snapshot_at(s, min(s$date) - 1), "at least 2")
})
