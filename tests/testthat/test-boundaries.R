# Exact-boundary tests.
#
# Every detector here decides on a threshold, and each threshold is either
# inclusive or strict. Testing 1.6 tells you nothing about 1.5. These pin the
# boundary itself, so that flipping < to <= anywhere below turns a test red
# instead of quietly moving every verdict that lands exactly on the line.

fit <- function(yi, vi, measure = "RR") {
  metafor::rma(yi = yi, vi = vi, measure = measure)
}

test_that("rcma's thresholds are inclusive at exactly the boundary", {
  # Floating point never lands ON a threshold: building data for a ratio of
  # 1.5 yields 1.5000000000000004, which sits above the line and so cannot
  # tell >= from >. The boundary is reachable from the other side -- both the
  # ratio and the threshold are exact values, so move the threshold onto the
  # observed ratio and the comparison becomes an equality.
  prev <- fit(rep(log(0.40), 4), rep(0.01, 4))
  new  <- fit(rep(log(0.40 * 1.5), 4), rep(0.01, 4))
  r <- rcma(prev, new)$signal
  expect_gt(r, 1.4)                       # a sane ratio, not the point

  # ratio == upper exactly: inclusive, so it signals.
  expect_equal(rcma(prev, new, upper = r)$verdict, "out_of_date")
  # A hair above it does not.
  expect_equal(rcma(prev, new, upper = r * (1 + 1e-12))$verdict, "current")

  # Same on the lower side, where the comparison is ratio <= lower.
  down <- fit(rep(log(0.40 * 0.5), 4), rep(0.01, 4))
  rl <- rcma(prev, down)$signal
  expect_equal(rcma(prev, down, lower = rl)$verdict, "out_of_date")
  expect_equal(rcma(prev, down, lower = rl * (1 - 1e-12))$verdict, "current")
})

test_that("ottawa's significance test is strict at exactly alpha", {
  # Same trick: alpha is a parameter, so set it to the p-value itself and the
  # comparison `pval < alpha` becomes a strict test of equality. Under `<` a
  # p sitting exactly on alpha is NOT significant; under `<=` it would be,
  # and the verdict would move.
  prev <- fit(c(0.05, -0.03, 0.04, -0.02), rep(0.20, 4))  # not significant
  new  <- fit(rep(log(0.20), 6), rep(0.02, 6))            # significant
  expect_gt(prev$pval, 0.04)

  # alpha set exactly to the updated p-value: `pval < alpha` is FALSE, so no
  # gain of significance is recorded and only the effect-size signal remains.
  at_alpha <- ottawa(prev, new, alpha = new$pval)
  expect_false(isTRUE(at_alpha$detail$signal_significance))

  # A hair above, and the very same data does count as newly significant.
  above <- ottawa(prev, new, alpha = new$pval * (1 + 1e-12))
  expect_true(isTRUE(above$detail$signal_significance))

  # The other side of the same comparison. `was_sig` reads the PRIOR p-value,
  # so its boundary is only reachable with alpha set onto that one. Under `<`
  # the prior review is not significant, the updated one is, and "gain" holds;
  # under `<=` the prior would count as significant too and the gain vanishes.
  on_prev <- ottawa(prev, new, alpha = prev$pval)
  expect_lt(new$pval, prev$pval)
  expect_true(isTRUE(on_prev$detail$signal_significance))
})

test_that("barrowman is strict at exactly q == 1", {
  prev <- fit(c(0.10, -0.05, 0.08, -0.02), rep(0.05, 4), measure = "MD")
  n_req <- (1.96^2 * 400) / prev$zval^2
  at    <- barrowman(prev, n_prev = 400, n_new = n_req)
  just  <- barrowman(prev, n_prev = 400, n_new = n_req * 1.0001)
  expect_equal(at$signal, 1, tolerance = 1e-9)
  expect_equal(at$verdict, "current")       # q == 1 is not "more than needed"
  expect_equal(just$verdict, "out_of_date")
})

test_that("studies sharing a date keep their input order", {
  # order(date, seq_along(date)) is what makes this stable. Without the tie
  # breaker the order would depend on the sort algorithm, and every snapshot
  # boundary with a tie could shift between R versions.
  ma <- fit(c(0.10, 0.20, 0.30, 0.40), rep(0.05, 4))
  st <- evidence_stream(ma, date = c(2000, 2000, 2000, 2001), study_id = seq_along(c(2000, 2000, 2000, 2001)))
  expect_equal(st$yi[1:3], c(0.10, 0.20, 0.30))
  expect_equal(st$date, c(2000, 2000, 2000, 2001))
})

test_that("effect_ratio refuses a prior effect it cannot tell from zero", {
  # A difference measure whose prior estimate is within min_z standard errors
  # of zero has no meaningful ratio; the guard needs se_prev to be usable.
  expect_true(is.na(effect_ratio(0.40, 0.02, "MD", se_prev = 0.5)$ratio))
  expect_true(is.na(effect_ratio(0.40, 0.02, "MD", se_prev = NA_real_)$ratio))
  expect_true(is.na(effect_ratio(0.40, 0.02, "MD", se_prev = Inf)$ratio))
  expect_true(is.na(effect_ratio(0.40, 0.02, "MD", se_prev = 0)$ratio))
  # Far enough from zero, it computes.
  expect_false(is.na(effect_ratio(0.40, 0.20, "MD", se_prev = 0.01)$ratio))
})

test_that("disagreement is reported when detectors actually disagree", {
  prev <- fit(rep(log(0.50), 6), rep(0.02, 6))
  new  <- list(yi = rep(log(1.20), 6), vi = rep(0.02, 6), k = 6)
  res  <- check_currency(prev, new, methods = c("rcma", "sufficiency"))
  calls <- vapply(res$verdicts, function(v) v$verdict, character(1))
  decided <- calls[calls != "not_applicable"]
  expect_equal(res$disagreement, length(unique(decided)) > 1)
  if (res$disagreement) expect_output(print(res), "disagree")
})

test_that("print methods cover their conditional branches", {
  prev <- fit(rep(log(0.50), 4), rep(0.30, 4))   # high heterogeneity path
  new  <- list(yi = c(log(0.20), log(1.60)), vi = c(0.02, 0.02), k = 2)
  expect_output(print(check_currency(prev, new, methods = "rcma")), "I2")

  st <- evidence_stream(fit(seq(0.10, 0.60, length.out = 12), rep(0.05, 12)),
                        date = 2000:2011, study_id = seq_along(2000:2011), ni = rep(100, 12))
  bt <- backtest(st, cuts = 2003:2008, horizon = 2, window = 2, seed = 3)
  expect_output(print(bt), "staleness_backtest")
  expect_output(print(st), "staleness_stream")
})

test_that("the plot marks contaminated methods on the axis", {
  # calibration(truth = "conclusion") flags ottawa, and the plot is supposed
  # to say so under the bars. Without this the warning exists only in the
  # data frame, where a reader looking at the picture never sees it.
  st <- evidence_stream(fit(seq(0.10, 0.70, length.out = 14), rep(0.05, 14)),
                        date = 2000:2013, study_id = seq_along(2000:2013), ni = rep(100, 14))
  bt <- backtest(st, cuts = 2003:2010, horizon = 2, window = 2, seed = 5)
  expect_true(any(calibration(bt, "conclusion")$contaminated))

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp); on.exit(unlink(tmp), add = TRUE)
  expect_silent(plot(bt, truth = "conclusion"))
  grDevices::dev.off()
})

test_that("a missing measure is not a ratio measure", {
  # rma objects built from yi/vi alone carry measure = "GEN" or nothing at
  # all, and effect_ratio() has to decide without one.
  expect_false(is_comparative_ratio(NULL))
  expect_false(is_comparative_ratio(""))
  expect_true(is_comparative_ratio("RR"))
})
