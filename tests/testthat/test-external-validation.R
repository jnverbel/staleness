# Cross-validation against independent implementations.
#
# Everything else in this suite checks the package against itself. These
# check it against code written by other people for the same quantities:
# metafor's fsn() and cumul(). An arithmetic mistake that is consistent
# throughout our own code would survive every other test in this directory
# and die here.

skip_if_not_installed("metafor")

bcg <- function() {
  skip_if_not_installed("metadat")
  d <- metadat::dat.bcg
  metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                  data = d)
}

test_that("failsafe_n agrees with metafor::fsn on Rosenthal's estimator", {
  es <- bcg()
  z05 <- stats::qnorm(0.05, lower.tail = FALSE)

  # metafor reports the fail-safe N as a whole number of studies; we keep the
  # continuous value, because it feeds the index N_fs / (5k + 10) where
  # rounding first would move the ratio. Ceiling ours and they must agree.
  ours <- failsafe_n(es$yi, es$vi, z_crit = z05)
  theirs <- metafor::fsn(yi, vi, data = es, type = "Rosenthal")$fsnum
  expect_equal(ceiling(ours), theirs)

  # Not a one-off: same agreement across randomly drawn evidence bodies.
  set.seed(20260808)
  for (i in 1:5) {
    k  <- sample(6:20, 1)
    yi <- stats::rnorm(k, 0.3, 0.2)
    vi <- stats::runif(k, 0.01, 0.2)
    expect_equal(max(0, ceiling(failsafe_n(yi, vi, z_crit = z05))),
                 metafor::fsn(yi, vi, type = "Rosenthal")$fsnum,
                 info = paste("k =", k))
  }
})

test_that("cumulative_effect agrees with metafor::cumul step for step", {
  es  <- bcg()
  ord <- order(es$year, seq_len(nrow(es)))
  esf <- es[ord, ]

  # cumul() on a fixed-effect fit is the same running pooled estimate our
  # cumulative_effect() computes in closed form.
  theirs <- as.numeric(
    metafor::cumul(metafor::rma(yi, vi, data = esf, method = "FE"))$estimate)
  ours <- cumulative_effect(esf$yi, esf$vi)

  expect_equal(length(ours), nrow(esf))
  expect_equal(ours, theirs, tolerance = 1e-12)
})

test_that("the snapshot at a cut equals metafor refitted on the same studies", {
  # No-look-ahead rests on this: a snapshot must be exactly what you would
  # have got by fitting only the studies published by then, tau^2 included.
  es <- bcg()
  st <- evidence_stream(metafor::rma(yi, vi, data = es), date = es$year)
  for (cut in c(1960, 1970, 1975)) {
    keep <- es$year <= cut
    direct <- metafor::rma(yi = es$yi[keep], vi = es$vi[keep])
    snap   <- snapshot_at(st, cut)
    expect_equal(as.numeric(snap$beta), as.numeric(direct$beta),
                 tolerance = 1e-10, info = paste("cut", cut))
    expect_equal(snap$tau2, direct$tau2, tolerance = 1e-10)
    expect_equal(snap$k, sum(keep))
  }
})
