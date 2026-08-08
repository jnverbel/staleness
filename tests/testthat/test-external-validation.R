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

# --- ottawa's effect criterion, against the published application ----------
#
# Pattanittum et al. (2012), Table 1, defines the Ottawa effect signal for
# ratio measures as the ratio of RELATIVE RISK REDUCTIONS -- (1 - RR_updated)
# over (1 - RR_previous) -- not the ratio of the risk ratios themselves. For
# mean differences it falls back to the rCMA ratio. Appendix S3, Table S6
# publishes the ten reviews with the largest indicator, with both effects and
# the RRR ratio the authors computed, so the criterion can be checked against
# real numbers rather than against our reading of the sentence.

OTTAWA_S6 <- data.frame(
  review    = c("Alfirevic 2010", "Dare 2006", "Hofmeyr 2008", "Tooher 2010",
                "Anotayanonth 2004", "Porter 2006", "Boulvain 2005",
                "Empson 2005", "Brown 2008", "Hopkins 1999"),
  rr_prev   = c(0.995, 0.998, 1.070, 1.057, 0.953, 1.095, 0.967, 1.021, 0.940, 0.929),
  rr_new    = c(0.848, 0.944, 1.367, 1.244, 0.803, 1.332, 0.903, 1.054, 0.854, 0.872),
  rrr_ratio = c(33.07, 22.80, 5.22, 4.30, 4.14, 3.51, 2.93, 2.54, 2.42, 1.79),
  stringsAsFactors = FALSE
)

test_that("ottawa reproduces the published RRR ratios and verdicts", {
  # The published effects are rounded to three decimals, and RRR = 1 - RR is
  # tiny when RR sits near 1, so the ratio is violently sensitive to that
  # rounding: at RR_prev = 0.998 the reduction is 0.002, and half a unit in
  # the last published decimal moves the ratio by a quarter. Comparing
  # against a fixed tolerance would be arbitrary. Instead propagate the
  # rounding interval and require the published figure to lie inside the band
  # it implies -- which is the strongest statement the published precision
  # can support.
  band <- function(rr_prev, rr_new, h = 0.0005) {
    lo_p <- rr_prev - h; hi_p <- rr_prev + h
    lo_n <- rr_new  - h; hi_n <- rr_new  + h
    vals <- c()
    for (p in c(lo_p, hi_p)) for (n in c(lo_n, hi_n))
      vals <- c(vals, (1 - n) / (1 - p))
    range(vals)
  }

  for (i in seq_len(nrow(OTTAWA_S6))) {
    row <- OTTAWA_S6[i, ]
    prev <- metafor::rma(yi = rep(log(row$rr_prev), 2), vi = c(0.001, 0.001),
                         measure = "RR")
    new  <- metafor::rma(yi = rep(log(row$rr_new), 2), vi = c(0.001, 0.001),
                         measure = "RR")
    res <- ottawa(prev, new)
    b <- band(row$rr_prev, row$rr_new)

    expect_gte(row$rrr_ratio, b[1] * 0.999)
    expect_lte(row$rrr_ratio, b[2] * 1.001)
    expect_gte(res$signal, b[1] * 0.999)
    expect_lte(res$signal, b[2] * 1.001)
    # All ten were published as out-of-date, and all ten must fire on effect.
    expect_equal(res$verdict, "out_of_date", info = row$review)
    expect_true(res$detail$signal_effect, info = row$review)
  }
})

test_that("the ratio of risk ratios would have missed all ten", {
  # The criterion the package used before: none of the ten fires, because the
  # risk ratios themselves barely move. Kept as a regression guard so the
  # implementation cannot drift back to it silently.
  fired_on_rr <- vapply(seq_len(nrow(OTTAWA_S6)), function(i) {
    r <- OTTAWA_S6$rr_new[i] / OTTAWA_S6$rr_prev[i]
    r <= 0.5 || r >= 1.5
  }, logical(1))
  expect_equal(sum(fired_on_rr), 0)
})
