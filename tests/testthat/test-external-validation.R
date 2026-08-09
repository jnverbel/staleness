# Cross-validation against independent implementations.
#
# Everything else in this suite checks the package against itself. These
# check it against code written by other people for the same quantities:
# metafor's fsn() and cumul(). An arithmetic mistake that is consistent
# throughout our own code would survive every other test in this directory
# and die here.

skip_if_not_installed("metafor")

# Wraps the shared fixture so this file's tests skip cleanly when metadat is
# absent, which every one of them needs.
bcg <- function() {
  skip_if_not_installed("metadat")
  bcg_es()
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
  st <- evidence_stream(metafor::rma(yi, vi, data = es), date = es$year, study_id = seq_along(es$year))
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

# --- The remaining three detectors, against the same published appendix ----

test_that("sufficiency's benchmark is 5k + 10 on the PRIOR review, as published", {
  # Appendix S3, Table S4 publishes Nfs, the benchmark and the failsafe ratio
  # for the ten reviews with the largest indicator. The benchmark is 5k + 10
  # computed on the number of studies in the PREVIOUS meta-analysis, which is
  # independent confirmation of the argument order this package had to correct
  # during development: sufficiency reads `prev`, stability reads `new_ma`.
  s4 <- data.frame(
    review  = c("Boulvain 2008", "Kelly 2009", "Alfirevic 2006", "Abalos 2007",
                "Alfirevic 2009", "Kenyon 2010", "Boulvain 2005",
                "Alfirevic 2010", "Haas 2008", "French 2001"),
    k_prev  = c(24, 22, 22, 20, 21, 16, 15, 15, 14, 13),
    nfs     = c(23.7, 21.7, 21.3, 20.0, 19.7, 15.8, 15.0, 14.4, 13.9, 12.9),
    bench   = c(130, 120, 120, 110, 115, 90, 85, 85, 80, 75),
    ratio   = c(0.18, 0.18, 0.18, 0.18, 0.17, 0.18, 0.18, 0.17, 0.17, 0.17)
  )
  expect_equal(5 * s4$k_prev + 10, s4$bench)
  expect_equal(s4$nfs / s4$bench, s4$ratio, tolerance = 0.02)

  # Every one of them is well below 1, i.e. insufficient -- which is why the
  # published study flagged none of its eighty reviews by this method, and why
  # `sufficient && !stable` is the correct reading of the rule rather than its
  # inverse.
  expect_true(all(s4$ratio < 1))

  # Our index is built the same way.
  prev <- metafor::rma(yi = rep(log(0.9), 24), vi = rep(0.05, 24))
  new  <- metafor::rma(yi = rep(log(0.9), 30), vi = rep(0.05, 30))
  res  <- sufficiency(prev, new)
  expect_equal(res$detail$k, 24)
  expect_equal(res$detail$index,
               failsafe_n(as.numeric(prev$yi), as.numeric(prev$vi)) / (5 * 24 + 10))
})

test_that("barrowman's participant ratio matches the published quotient", {
  # Table S5 publishes n_actual, n_expected and the ratio for ten reviews.
  s5 <- data.frame(
    n_act = c(1734, 126, 377, 4725, 494, 6241, 300, 116, 85, 375),
    n_exp = c(50, 12, 52, 863, 249, 3173, 185, 143, 119, 681),
    ratio = c(34.92, 10.22, 7.25, 5.48, 1.98, 1.97, 1.62, 0.81, 0.71, 0.55)
  )
  expect_equal(s5$n_act / s5$n_exp, s5$ratio, tolerance = 0.03)
  # Seven of the ten exceed 1 and were published as out-of-date; three did not.
  expect_equal(sum(s5$ratio > 1), 7)

  # Our q is that same quotient, with n_expected from the published formula
  # n = (Z_crit^2 * N) / Z^2.
  prev <- metafor::rma(yi = c(0.10, -0.05, 0.08, -0.02), vi = rep(0.05, 4),
                       measure = "MD")
  n_req <- (1.96^2 * 400) / prev$zval^2
  expect_equal(barrowman(prev, n_prev = 400, n_new = 2 * n_req)$signal, 2,
               tolerance = 1e-8)
})

test_that("simulation's 80% threshold is what left the published cohort empty", {
  # Table S7: the ten reviews with the highest simulated power top out at
  # 63.4%, below the method's 80% threshold -- which is why the study reports
  # zero out-of-date reviews by this method. Our default threshold is the same.
  powers <- c(63.4, 62.0, 60.7, 60.3, 49.5, 44.7, 37.9, 29.0, 27.5, 27.0) / 100
  expect_true(all(powers < 0.80))
  expect_equal(formals(simulation)$power_threshold, 0.80)
})

test_that("simulation follows the published procedure, step for step", {
  # Pattanittum et al. (2012), Appendix S1, sets out the simulation-based
  # power method: (a) draw the new study's effect from a t-distribution with
  # the parameters of the previous meta-analysis; the new study carries the
  # combined size of the recent studies; (b-d) re-pool and test at 5%;
  # (e) repeat 10,000 times, and power ABOVE 80% means out-of-date.
  expect_equal(formals(simulation)$B, 10000)
  expect_equal(formals(simulation)$alpha, 0.05)
  expect_equal(formals(simulation)$power_threshold, 0.80)

  prev <- metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
                       vi = c(0.16, 0.20, 0.18, 0.15, 0.22))
  new  <- list(yi = c(-0.45, -0.38, -0.52, -0.29, -0.41),
               vi = c(0.05, 0.04, 0.06, 0.05, 0.04), k = 5)
  res <- simulation(prev, new, B = 2000, seed = 1)

  # One simulated study, carrying the summed precision of the recent ones.
  expect_equal(res$detail$k_simulated, 1)
  expect_equal(res$detail$vi_new, 1 / sum(1 / new$vi), tolerance = 1e-12)
  expect_equal(res$detail$df, length(prev$yi) - 1)

  # The threshold is strict: power exactly at the threshold is not a signal.
  expect_equal(simulation(prev, new, B = 100, seed = 1,
                          power_threshold = 0)$verdict, "out_of_date")
  at <- simulation(prev, new, B = 2000, seed = 1)
  expect_equal(simulation(prev, new, B = 2000, seed = 1,
                          power_threshold = at$signal)$verdict, "current")
})

test_that("a second source's worked examples pick the same reading", {
  # Kuhnisch et al. (2013, PMC3881834) apply the Ottawa method citing Shojania
  # and work two examples, both declared to meet the 50% criterion. Only one
  # reading of "relative effect size" makes both statements true.
  ex <- data.frame(prev = c(2.10, 2.61), new = c(1.51, 1.66))
  fires <- function(x) x <= 0.5 || x >= 1.5

  ratio_of_effects <- ex$new / ex$prev             # 0.719, 0.636
  pct_change       <- abs(ex$new - ex$prev) / ex$prev  # 0.281, 0.364
  ratio_of_rrr     <- (1 - ex$new) / (1 - ex$prev) # 0.464, 0.410

  expect_false(any(vapply(ratio_of_effects, fires, logical(1))))
  expect_false(any(pct_change > 0.5))
  expect_true(all(vapply(ratio_of_rrr, fires, logical(1))))

  # And the detector agrees with the source on both.
  for (i in 1:2) {
    prev <- metafor::rma(yi = rep(log(ex$prev[i]), 2), vi = c(0.001, 0.001),
                         measure = "RR")
    new  <- metafor::rma(yi = rep(log(ex$new[i]), 2), vi = c(0.001, 0.001),
                         measure = "RR")
    expect_true(ottawa(prev, new)$detail$signal_effect)
  }
})

# --- The historical case, end to end -------------------------------------
#
# Lau et al. (1992), NEJM 327:248-254, introduced cumulative meta-analysis on
# 33 streptokinase trials and reported two things precisely enough to test
# against: a statistically significant mortality reduction was reached in 1973
# after eight trials (OR 0.74, 95% CI 0.59 to 0.92), and the 25 later trials
# "had little or no effect on the odds ratio". The dataset is public, as
# metadat::dat.lau1992, so the whole pipeline can be run against a documented
# historical fact rather than against itself.

test_that("the package reproduces Lau's 1973 cut point", {
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.lau1992)
  # Mantel-Haenszel, which is what Lau used, and what reproduces all three
  # published figures exactly. Inverse-variance fixed effects gives 0.744
  # (0.595, 0.929) -- the same conclusion, but it rounds to a different
  # interval, which is why the estimator has to match to claim reproduction.
  mh73 <- metafor::rma.mh(measure = "OR", ai = ai, n1i = n1i, ci = ci,
                          n2i = n2i,
                          data = metadat::dat.lau1992[
                            metadat::dat.lau1992$year <= 1973, ])
  expect_equal(mh73$k, 8)
  expect_equal(round(exp(as.numeric(mh73$beta)), 2), 0.74)
  expect_equal(round(exp(mh73$ci.lb), 2), 0.59)
  expect_equal(round(exp(mh73$ci.ub), 2), 0.92)

  # The package's own pipeline reaches the same conclusion on the same cut.
  fit73 <- metafor::rma(yi, vi, data = es[es$year <= 1973, ], method = "FE")
  expect_equal(fit73$k, 8)
  expect_lt(fit73$pval, 0.05)

  # And the cut before it is not yet significant, so 1973 really is the year.
  fit71 <- metafor::rma(yi, vi, data = es[es$year <= 1971, ], method = "FE")
  expect_gt(fit71$pval, 0.05)
})

test_that("the backtest puts the turning point in 1973 and stays quiet after", {
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.lau1992)
  ma <- metafor::rma(yi, vi, data = es, measure = "OR", method = "FE")
  st <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = es$n1i + es$n2i)
  bt <- backtest(st, cuts = "yearly", horizon = 3, window = 5, min_k = 3,
                 seed = 42)
  r <- bt$results[!bt$results$censored, ]

  # The conclusion still had to change from every cut before 1973, and from
  # none after it.
  tc <- unique(r[, c("cut", "truth_conclusion")])
  expect_true(all(tc$truth_conclusion[tc$cut < 1973]))
  expect_false(any(tc$truth_conclusion[tc$cut >= 1973]))

  # rcma watches the size of the effect, and Lau reports the effect did not
  # move after 1973. It must not fire.
  rc <- r[r$method == "rcma", ]
  expect_equal(sum(rc$verdict == "out_of_date"), 0)

  # ottawa watches significance, and fires on exactly the cuts before 1973.
  ot <- r[r$method == "ottawa", ]
  expect_setequal(ot$cut[ot$verdict == "out_of_date"], c(1969, 1970, 1971, 1972))

  # It warns ahead of the event rather than alongside it. Pinned to the figure
  # the validation report cites, not to a bound: `> 0` would still pass if the
  # lead collapsed from 1.5 years to a single month, and a validation whose
  # numbers can drift without going red is not validating anything.
  lt <- lead_time(bt, "conclusion")
  expect_equal(lt$median_lead[lt$method == "ottawa"], 1.5)

  # And the circular pair is flagged without being asked.
  cal <- calibration(bt, "conclusion")
  expect_true(cal$contaminated[cal$method == "ottawa"])
})

test_that("the estimator choice moves the answer by years, not decimals", {
  # Worth pinning: Lau's 1973 conclusion depends on fixed effects. Under REML
  # the same evidence is not significant until 1979. staleness inherits the
  # method from the rma it is handed, so this is the analyst's decision, and a
  # substantive one.
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.lau1992)
  first_sig <- function(method) {
    for (y in sort(unique(es$year))) {
      d <- es[es$year <= y, ]
      if (nrow(d) < 3) next
      f <- suppressWarnings(metafor::rma(yi, vi, data = d, method = method))
      if (f$pval < 0.05) return(y)
    }
    NA_integer_
  }
  expect_equal(first_sig("FE"), 1973)
  # 1979, not merely "later": the gap is six years, and that size is the point
  # of the test.
  expect_equal(first_sig("REML"), 1979)
})

test_that("a Mantel-Haenszel fit is refused with a usable explanation", {
  skip_if_not_installed("metadat")
  mh <- metafor::rma.mh(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.lau1992)
  expect_error(evidence_stream(mh, date = metadat::dat.lau1992$year, study_id = seq_along(metadat::dat.lau1992$year)),
               "rma.mh")
  expect_error(evidence_stream(mh, date = metadat::dat.lau1992$year, study_id = seq_along(metadat::dat.lau1992$year)),
               "inverse-variance")
})

# --- The mirror case: a review that WAS dangerously out of date -----------
#
# dat.li2007 is intravenous magnesium for acute myocardial infarction, the
# textbook counterexample to Lau. Small trials suggested a large benefit; the
# 58,050-patient ISIS-4 trial in 1995 removed it. The Cochrane review these
# data come from reports both pooled results, so the reproduction is checkable.
#
# The pair matters more than either case alone. In Lau the correct behaviour
# for rcma is silence -- the effect never moved. Here it is to speak, loudly
# and early. A detector that only ever fires, or only ever stays quiet, would
# look right in one of the two and wrong in the other.

test_that("the package reproduces the Cochrane review's two pooled results", {
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.li2007)
  # Mantel-Haenszel for the fixed-effect result: RevMan's default for binary
  # outcomes, and the second case in this file where inverse variance gets the
  # point estimate right but rounds the interval differently (1.05 vs 1.04).
  # Reproducing a Cochrane review to the digit means matching its estimator.
  fe <- metafor::rma.mh(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.li2007)
  re <- metafor::rma(yi, vi, data = es, method = "DL")   # Cochrane's default

  # "a fixed-effect meta-analysis showed no difference (OR 0.99, 0.94 to 1.04)"
  expect_equal(round(exp(as.numeric(fe$beta)), 2), 0.99)
  expect_equal(round(exp(fe$ci.lb), 2), 0.94)
  expect_equal(round(exp(fe$ci.ub), 2), 1.04)

  # "a random-effects meta-analysis showed a significant reduction
  #  (OR 0.66, 0.53 to 0.82)"
  expect_equal(round(exp(as.numeric(re$beta)), 2), 0.66)
  expect_equal(round(exp(re$ci.lb), 2), 0.53)
  expect_equal(round(exp(re$ci.ub), 2), 0.82)
})

test_that("ISIS-4 removes the effect, and rcma saw it coming", {
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.li2007)

  # Before ISIS-4: a significant 35% reduction. After it: nothing.
  before <- metafor::rma(yi, vi, data = es[es$year <= 1994, ], method = "FE")
  after  <- metafor::rma(yi, vi, data = es[es$year <= 1995, ], method = "FE")
  expect_equal(round(exp(as.numeric(before$beta)), 3), 0.650)
  expect_equal(round(before$pval, 5), 0.00047)
  expect_equal(round(exp(as.numeric(after$beta)), 3), 1.021)
  expect_equal(round(after$pval, 4), 0.4948)

  ma <- metafor::rma(yi, vi, data = es, measure = "OR", method = "FE")
  st <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = es$n1i + es$n2i)
  bt <- backtest(st, cuts = "yearly", horizon = 3, window = 5, min_k = 3,
                 seed = 42)
  r <- bt$results[!bt$results$censored, ]

  # rcma watches effect size, and here the effect moved from 0.65 to 1.02. It
  # fires on the cuts before ISIS-4 and on none after -- the opposite of its
  # behaviour on dat.lau1992, where the effect never moved and it never fired.
  rc <- r[r$method == "rcma", ]
  fired <- rc$cut[rc$verdict == "out_of_date"]
  # The exact set, not a count with a floor: nine consecutive cuts, every one
  # of them before ISIS-4.
  expect_equal(fired, c(1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994))

  # Scored against a truth that shares no logic with it: perfect, and early.
  cal <- calibration(bt, "shift")
  expect_equal(cal$sensitivity[cal$method == "rcma"], 1)
  expect_equal(cal$specificity[cal$method == "rcma"], 1)
  expect_false(cal$contaminated[cal$method == "rcma"])
  expect_equal(lead_time(bt, "shift")$median_lead[
    lead_time(bt, "shift")$method == "rcma"], 4)
})

# --- A third case, for the branch the other two never reach ---------------
#
# dat.lau1992 and dat.li2007 are both odds ratios in cardiology. Every ratio
# measure takes one branch of effect_ratio(); the difference-measure branch,
# with its guard against a prior effect indistinguishable from zero, had never
# been exercised on real data.
#
# dat.bangertdrowns2004 is 48 writing-to-learn studies, 1926-1998, on the
# standardized mean difference scale, with a small pooled effect (0.22). Note
# what this case is NOT for: metadat states that these values are
# bias-corrected and so "differ slightly from the values reported in the
# article", and that the variances assume equal group sizes. Reproducing the
# published pooled estimate is therefore impossible by construction, and not
# attempted. This case is here for branch coverage on real data.

test_that("difference measures take the other branch, and the guard fires", {
  skip_if_not_installed("metadat")
  d <- metadat::dat.bangertdrowns2004
  d <- d[!is.na(d$yi) & !is.na(d$vi) & !is.na(d$year), ]
  ma <- metafor::rma(yi, vi, data = d, measure = "SMD")
  st <- evidence_stream(ma, date = d$year, study_id = seq_along(d$year), ni = d$ni)

  expect_false(is_comparative_ratio(st$measure))

  bt <- backtest(st, cuts = "yearly", horizon = 3, window = 5, min_k = 5,
                 seed = 42)
  r <- bt$results[!bt$results$censored, ]

  # The guard documented in effect_ratio() -- a prior effect too close to zero
  # for a ratio to mean anything -- actually fires here, on a real body of
  # evidence, and says why.
  rc <- r[r$method == "rcma", ]
  na <- rc[rc$verdict == "not_applicable", ]
  expect_gt(nrow(na), 0)
  expect_true(all(grepl("indistinguishable from zero", na$reason)))

  # A ratio it refuses to compute must never be dressed up as a signal.
  expect_equal(sum(rc$verdict == "out_of_date"), 0)
})

test_that("ottawa still works when its effect half cannot be computed", {
  # On difference measures the Ottawa effect criterion defers to the rcma
  # rule, so when the guard blanks the ratio, only the significance half is
  # left. The detector has to keep answering on that half rather than falling
  # over or going silent.
  skip_if_not_installed("metadat")
  d <- metadat::dat.bangertdrowns2004
  d <- d[!is.na(d$yi) & !is.na(d$vi) & !is.na(d$year), ]
  ma <- metafor::rma(yi, vi, data = d, measure = "SMD")
  st <- evidence_stream(ma, date = d$year, study_id = seq_along(d$year), ni = d$ni)

  prev <- snapshot_at(st, 1985)
  new  <- snapshot_at(st, 1990)
  res  <- ottawa(prev, new)

  expect_true(is.na(res$signal))              # the ratio was refused
  expect_false(res$detail$signal_effect)      # so the effect half is silent
  expect_true(res$detail$signal_significance) # and the other half carried it
  expect_equal(res$verdict, "out_of_date")
})

# --- A fourth case, found by searching for one the others could not give ---
#
# barrowman() and simulation() only speak when the prior meta-analysis was not
# significant, and in the three cases above that condition rarely held, so both
# spent most cuts returning not_applicable. Sweeping metadat for datasets with
# several non-significant cuts turned up dat.laopaiboon2015: azithromycin for
# lower respiratory tract infection, 15 trials 1991-2002, pooled OR 1.12
# (0.61-2.04), p = 0.72. A null meta-analysis from end to end.

test_that("barrowman and simulation answer on every cut of a null review", {
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.laopaiboon2015)
  es <- es[!is.na(es$yi) & !is.na(es$vi) & es$vi > 0, ]
  ma <- metafor::rma(yi, vi, data = es, measure = "OR")
  expect_gt(ma$pval, 0.05)                      # the precondition both need

  st <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = es$n1i + es$n2i)
  bt <- backtest(st, cuts = "yearly", horizon = 2, window = 3, min_k = 3,
                 seed = 42)
  r <- bt$results[!bt$results$censored, ]

  for (m in c("barrowman", "simulation")) {
    d <- r[r$method == m, ]
    expect_equal(sum(d$verdict == "not_applicable"), 0, info = m)
    expect_equal(nrow(d), 8, info = m)
  }
  # Neither raises a false alarm on evidence that never moved.
  cal <- calibration(bt, "shift")
  expect_equal(cal$specificity[cal$method == "barrowman"], 1)
  expect_equal(cal$specificity[cal$method == "simulation"], 1)
})

test_that("ottawa's instability on null reviews shows up in real data too", {
  # Measured by simulation elsewhere in this suite: under a null effect the
  # RRR criterion fires on most samples containing no change, because
  # 1 - RR_prev approaches zero. Here is the same failure on a real review,
  # which is the kind of confirmation a simulated result cannot give itself.
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.laopaiboon2015)
  es <- es[!is.na(es$yi) & !is.na(es$vi) & es$vi > 0, ]
  ma <- metafor::rma(yi, vi, data = es, measure = "OR")
  st <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = es$n1i + es$n2i)
  bt <- backtest(st, cuts = "yearly", horizon = 2, window = 3, min_k = 3,
                 seed = 42)
  cal <- calibration(bt, "shift")

  # Specificity collapses for ottawa and holds for everything else.
  # 1/7 exactly -- one correct silence in seven chances. Pinned as the figure
  # the report cites rather than as "below 0.3".
  expect_equal(cal$specificity[cal$method == "ottawa"], 1 / 7)
  expect_equal(cal$specificity[cal$method == "rcma"], 1)
  expect_equal(cal$specificity[cal$method == "sufficiency"], 1)
})

# --- Every refit in the package must honour the same options ---------------

test_that("all three refit sites propagate the caller's model options", {
  # snapshot_at() was fixed to carry test/weighted/tau2. Two other places
  # refit: backtest() builds the `final` model that every truth column is
  # scored against, and check_currency() builds the `updated` model handed to
  # rcma, ottawa and sufficiency. Fixing one and not the others leaves the
  # verdicts and the truth computed under different models.
  skip_if_not_installed("metadat")
  es <- bcg_es()
  kn <- metafor::rma(yi, vi, data = es, test = "knha")
  ni <- es$tpos + es$tneg + es$cpos + es$cneg
  st <- evidence_stream(kn, date = es$year, study_id = seq_along(es$year), ni = ni)

  # 1. snapshot_at
  expect_equal(snapshot_at(st, 1970)$test, "knha")

  # 2. check_currency's updated model, reached through its verdicts: with the
  #    Knapp-Hartung adjustment the updated p-value is much larger, and
  #    ottawa() decides on p-values.
  prev <- snapshot_at(st, 1970)
  new  <- window_between(st, 1970, 1980)
  chk  <- check_currency(prev, new, methods = "ottawa")
  direct <- metafor::rma(yi = c(as.numeric(prev$yi), new$yi),
                         vi = c(as.numeric(prev$vi), new$vi),
                         measure = prev$measure, method = prev$method,
                         test = prev$test)
  expect_equal(chk$updated$test, "knha")
  expect_equal(chk$updated$pval, direct$pval, tolerance = 1e-10)

  # 3. backtest's final model, which every truth column is scored against.
  bt <- suppressWarnings(backtest(st, cuts = "yearly", horizon = 3,
                                  window = 3, seed = 1))
  fin <- metafor::rma(yi = st$yi, vi = st$vi, measure = st$measure,
                      method = st$method, test = st$test)
  # truth_shift divides by se_final; recompute it and require a match.
  prev70 <- snapshot_at(st, 1970)
  expected <- truth_shift(as.numeric(prev70$beta), as.numeric(fin$beta), fin$se)
  got <- unique(bt$results$truth_shift[bt$results$cut == 1970])
  expect_equal(got, expected)
})

# --- Two more real cases, with the criterion matched to the detector -------
#
# The four cases above all validate against what happened to a review's
# CONCLUSION. That works when magnitude and conclusion move together, as they
# do in dat.li2007. When they come apart, a conclusion-level criterion cannot
# judge rcma(), which measures the size of the pooled effect and nothing else.
#
# These two are a pair on purpose: in the first the correct behaviour is
# silence, in the second it is to fire. Neither case alone separates a working
# detector from one that is simply mute, or simply noisy.
#
# Both criteria were read from the raw abstracts of the successive Cochrane
# versions via Europe PMC, not from a summarising fetch.

test_that("damico2009: the detectors stay silent on a review that held up", {
  # metadat::dat.damico2009 is the evidence behind Cochrane CD000022.pub3
  # (2009): topical plus systemic antibiotics to prevent respiratory tract
  # infections in intensive care. The review was updated as .pub4 in 2021 and
  # its conclusion did not change -- combined prophylaxis reduces both RTIs
  # and mortality, topical alone reduces RTIs but not mortality, in both
  # versions. Twelve years and a full update later, staying quiet was right.
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = xt, n1i = nt, ci = xc, n2i = nc,
                        data = metadat::dat.damico2009)
  keep <- is.finite(es$yi) & is.finite(es$vi) & es$vi > 0
  es <- es[keep, ]
  yr <- metadat::dat.damico2009$year[keep]
  ni <- (metadat::dat.damico2009$nt + metadat::dat.damico2009$nc)[keep]

  ma <- metafor::rma(yi, vi, data = es, measure = "OR", method = "FE")
  st <- evidence_stream(ma, date = yr, study_id = seq_along(yr), ni = ni)
  bt <- suppressWarnings(backtest(st, cuts = "yearly", horizon = 3,
                                  window = 5, min_k = 3, seed = 42))
  r <- bt$results[!bt$results$censored, ]

  for (m in c("rcma", "ottawa", "sufficiency")) {
    d <- r[r$method == m, ]
    expect_equal(nrow(d), 12, info = m)
    expect_equal(sum(d$verdict == "out_of_date"), 0, info = m)
  }
  # The effect barely moves, which is why: a 7% shift over the whole series.
  last <- r[r$method == "rcma" & r$cut == max(r$cut), ]
  expect_equal(round(last$signal, 3), 1.069)
})

test_that("lee2004: rcma fires on a real change of magnitude, and is right", {
  # metadat::dat.lee2004 is P6 acupoint stimulation for postoperative nausea,
  # the evidence behind Cochrane CD003281.pub2 (2004).
  #
  # The conclusion never changed: P6 still works in .pub3 (2009) and .pub5
  # (2025). The MAGNITUDE did, inside the window this dataset covers -- the
  # pooled odds ratio goes from 0.38 up to 0.60, a 55% move, past rcma's 50%
  # threshold. rcma measures magnitude, so it fires, and it is correct to.
  #
  # A conclusion-level criterion would have scored this as a false positive.
  # It is not one: the two Cochrane versions either side report RR 0.72 and
  # RR 0.71, so the effect settled AFTER the period this dataset ends in,
  # which is exactly the story the firing tells.
  skip_if_not_installed("metadat")
  es <- metafor::escalc(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.lee2004)
  keep <- is.finite(es$yi) & is.finite(es$vi) & es$vi > 0
  es <- es[keep, ]
  yr <- metadat::dat.lee2004$year[keep]

  early <- metafor::rma(yi, vi, data = es[yr <= 1996, ], method = "FE")
  whole <- metafor::rma(yi, vi, data = es, method = "FE")
  expect_equal(round(exp(as.numeric(early$beta)), 3), 0.385)
  expect_equal(round(exp(as.numeric(whole$beta)), 3), 0.598)
  # The move that rcma is built to see, computed here independently of it.
  expect_gt(exp(as.numeric(whole$beta)) / exp(as.numeric(early$beta)), 1.5)

  ma <- metafor::rma(yi, vi, data = es, measure = "OR", method = "FE")
  st <- evidence_stream(ma, date = yr, study_id = seq_along(yr),
                        ni = (metadat::dat.lee2004$n1i +
                              metadat::dat.lee2004$n2i)[keep])
  bt <- suppressWarnings(backtest(st, cuts = "yearly", horizon = 3,
                                  window = 5, min_k = 3, seed = 42))
  r <- bt$results[!bt$results$censored, ]

  rc <- r[r$method == "rcma", ]
  expect_equal(rc$cut[rc$verdict == "out_of_date"], c(1995, 1996, 1997))

  # And ottawa stays quiet on the same evidence, because significance never
  # changed. The two disagree here and both are right -- which is the whole
  # reason a criterion has to match the detector it judges.
  ot <- r[r$method == "ottawa", ]
  expect_equal(sum(ot$verdict == "out_of_date"), 0)
})
