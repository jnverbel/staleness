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

# --- ni is validated for NA, like date -------------------------------------
#
# evidence_stream() checked `ni`'s length but not for missing values, two lines
# after checking `date` for exactly that. A single NA survived into
# backtest()'s sum(stream$ni[...]) -- sum() over an NA is NA -- and surfaced,
# many frames later, as barrowman() stopping with "needs the sample size of
# both the prior meta-analysis and the new studies". That message describes a
# call that supplied no sample sizes at all, which is not what happened.

test_that("evidence_stream refuses missing sample sizes, with a message that says so", {
  ma <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6), measure = "RR")
  expect_error(
    evidence_stream(ma, date = 2000:2005, ni = c(100, 100, NA, 100, 100, 100)),
    "missing values"
  )
  # And the same stream with a complete `ni` is accepted, so the check is not
  # rejecting the shape of the input.
  expect_s3_class(evidence_stream(ma, date = 2000:2005, ni = rep(100, 6)),
                  "staleness_stream")
})
