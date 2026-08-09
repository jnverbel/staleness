test_that("comparative ratio measures are identified", {
  expect_true(is_comparative_ratio("RR"))
  expect_true(is_comparative_ratio("OR"))
  expect_true(is_comparative_ratio("IRR"))
  expect_true(is_comparative_ratio("HR"))
  expect_true(is_comparative_ratio("ROM"))
  expect_true(is_comparative_ratio("PETO"))
  expect_false(is_comparative_ratio("MD"))
  expect_false(is_comparative_ratio("SMD"))
  expect_false(is_comparative_ratio("RD"))
})

test_that("to_natural exponentiates only ratio measures", {
  # Empirically verified with metadat::dat.bcg: beta = -0.7145323 is log(RR)
  expect_equal(to_natural(-0.7145323, "RR"), 0.4894209, tolerance = 1e-6)
  expect_equal(to_natural(-0.7145323, "MD"), -0.7145323)
})

test_that("effect_ratio works on the natural scale for ratio measures", {
  # Prior RR 0.50, new RR 0.25 -> ratio 0.5, signal of rCMA
  r <- effect_ratio(log(0.25), log(0.50), "RR", se_prev = 0.1)
  expect_equal(r$ratio, 0.5, tolerance = 1e-8)
  expect_equal(r$reason, "")
})

test_that("effect_ratio refuses to divide by a near-null difference effect", {
  # Prior MD 0.05 with se 0.10: indistinguishable from zero, the ratio is unstable
  r <- effect_ratio(0.40, 0.05, "MD", se_prev = 0.10)
  expect_true(is.na(r$ratio))
  expect_match(r$reason, "indistinguishable from zero")
})

test_that("effect_ratio computes difference ratios when the prior effect is solid", {
  r <- effect_ratio(0.20, 0.80, "MD", se_prev = 0.05)
  expect_equal(r$ratio, 0.25, tolerance = 1e-8)
  expect_equal(r$reason, "")
})

test_that("single-group summaries are not comparative effects", {
  # The old RATIO_MEASURES answered two questions with one list: how a measure
  # is STORED, and whether it is a comparative effect at all. It was documented
  # as the first and used as the second, and it contained PLO and IRLN.
  #
  # PLO is a proportion on the logit scale, IRLN an incidence rate on the log
  # scale. Neither compares a treated arm against a control, so the Ottawa
  # criterion's `1 - exp(theta)` is not a risk reduction of anything: on a
  # pooled proportion of 0.24 it returned 0.68, and a verdict was issued on it.
  for (m in c("PLO", "IRLN", "PR", "IR", "MN", "COR", "ZCOR")) {
    expect_false(is_comparative_ratio(m), info = m)
    expect_true(is_single_group(m), info = m)
  }
  for (m in c("RR", "OR", "PETO", "IRR", "HR", "ROM")) {
    expect_true(is_comparative_ratio(m), info = m)
    expect_false(is_single_group(m), info = m)
  }
})

test_that("the criteria refuse a measure with no treatment effect to compare", {
  # Refused before either branch: running a single-group summary through the
  # difference branch would be as meaningless as through the ratio branch.
  for (m in c("PLO", "IRLN", "PR", "COR")) {
    r <- effect_ratio(log(0.30), log(0.25), m, se_prev = 0.1)
    expect_true(is.na(r$ratio), info = m)
    expect_match(r$reason, "one group", info = m)

    q <- rrr_ratio(log(0.30), log(0.25), m, se_prev = 0.1)
    expect_true(is.na(q$ratio), info = m)
    expect_match(q$reason, "one group", info = m)
  }
  # And a detector handed one answers not_applicable rather than a number.
  prev <- metafor::rma(yi = rep(qlogis(0.20), 4), vi = rep(0.05, 4),
                       measure = "PLO")
  upd  <- metafor::rma(yi = c(rep(qlogis(0.20), 4), rep(qlogis(0.30), 3)),
                       vi = c(rep(0.05, 4), rep(0.02, 3)), measure = "PLO")
  expect_equal(rcma(prev, upd)$verdict, "not_applicable")
  expect_match(rcma(prev, upd)$reason, "one group")

  # ottawa is different, and deliberately so: it carries three signals, and
  # the significance half is still computable here. It answers, but its EFFECT
  # criterion does not fire and its reason says why -- the same behaviour
  # already validated on metadat::dat.bangertdrowns2004, where the ratio comes
  # back NA and ottawa keeps answering on significance alone.
  o <- ottawa(prev, upd)
  expect_true(is.na(o$signal))
  expect_false(o$detail$signal_effect)
  expect_match(o$reason, "one group")
})

test_that("a logit is inverted by plogis, not by exp", {
  # PLO was wrong on the storage question too. metafor stores logit(p), and
  # exp() on it returns the ODDS: for a pooled proportion of 0.2434 the old
  # code reported 0.3216.
  p <- 0.2434
  theta <- qlogis(p)
  expect_equal(to_natural(theta, "PLO"), p, tolerance = 1e-9)
  expect_false(isTRUE(all.equal(to_natural(theta, "PLO"), exp(theta))))
  # A log-scale measure still exponentiates.
  expect_equal(to_natural(log(0.5), "RR"), 0.5, tolerance = 1e-9)
  expect_equal(to_natural(log(0.05), "IRLN"), 0.05, tolerance = 1e-9)
  # And a difference measure is left alone.
  expect_equal(to_natural(-0.71, "SMD"), -0.71)
})
