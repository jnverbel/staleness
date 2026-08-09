test_that("truth_shift measures movement relative to final precision", {
  expect_true(truth_shift(theta_t = 0.0, theta_final = 0.5, se_final = 0.1))
  expect_false(truth_shift(theta_t = 0.48, theta_final = 0.5, se_final = 0.1))
})

test_that("truth_surprise uses the precision of the time, not the final one", {
  # final estimate far outside the old confidence interval
  expect_true(truth_surprise(theta_t = 0.0, se_t = 0.05, theta_final = 0.5))
  # same final estimate, but back then the interval was wide enough to contain it
  expect_false(truth_surprise(theta_t = 0.0, se_t = 0.40, theta_final = 0.5))
})

test_that("truth_conclusion catches sign and significance changes", {
  expect_true(truth_conclusion(0.5, 0.01, -0.5, 0.01))   # sign flip
  expect_true(truth_conclusion(0.5, 0.20, 0.5, 0.01))    # became significant
  expect_false(truth_conclusion(0.5, 0.01, 0.6, 0.01))   # neither
})

test_that("the contaminated pairs are declared in data, not in a footnote", {
  expect_true(any(CONTAMINATED_PAIRS$method == "ottawa" &
                  CONTAMINATED_PAIRS$truth  == "conclusion"))
})

test_that("shift and surprise are not the same test", {
  # identical inputs, opposite answers: se_t is wide, se_final is narrow
  expect_true(truth_shift(0.0, 0.5, se_final = 0.10))
  expect_false(truth_surprise(0.0, se_t = 0.40, theta_final = 0.5))
})

test_that("an effect of exactly zero counts as a sign change, deliberately", {
  # sign(0) is 0, so it differs from both +1 and -1 and truth_conclusion()
  # reads a flip. Pinned here so the behaviour is a decision rather than an
  # accident of sign(): an estimate sitting exactly on the null moving to a
  # definite direction IS a change in the practical conclusion, and an exact
  # zero is measure-zero on real data anyway.
  expect_true(truth_conclusion(0, 0.9, -0.5, 0.01))
  expect_true(truth_conclusion(-0.5, 0.01, 0, 0.9))
  # Two zeros do not flip: sign(0) == sign(0).
  expect_false(truth_conclusion(0, 0.9, 0, 0.9))
})

test_that("a degenerate standard error makes truth unknown, not certain", {
  # ?backtest promises these columns are NA when the standard error they
  # divide by is degenerate, and calibration()/lead_time() drop NA rows for
  # exactly that reason. Division by zero was returning TRUE instead: an
  # unknowable truth scored as a certain event, silently inflating every
  # detector's apparent miss rate at that cut.
  expect_true(is.na(truth_shift(0.5, 1.0, 0)))
  expect_true(is.na(truth_shift(0.5, 1.0, -1)))
  expect_true(is.na(truth_shift(0.5, 1.0, NA_real_)))
  expect_true(is.na(truth_shift(0.5, 1.0, Inf)))
  expect_true(is.na(truth_surprise(0.5, 0, 1.0)))
  expect_true(is.na(truth_surprise(0.5, NA_real_, 1.0)))

  # A usable standard error still answers.
  expect_true(truth_shift(0.5, 1.0, 0.1))
  expect_false(truth_shift(0.5, 0.51, 0.1))
})

test_that("an effect that cannot be located cannot have moved", {
  # The three truths already returned NA for a degenerate standard error --
  # "the distance cannot be judged, not that it is infinite" -- and then
  # returned TRUE for a degenerate EFFECT. Same class of degeneracy, opposite
  # treatment, inside the same function:
  #
  #   truth_shift(0, 1, se_final = 0)   -> NA    (SE degenerate)
  #   truth_shift(Inf, 0, se_final = 1) -> TRUE  (effect degenerate)
  #
  # A spurious TRUE is worse than a wrong answer here: it is a free hit for
  # every detector that fires at that cut, so it inflates sensitivity and
  # biases the calibration this package exists to produce. rcma() had always
  # treated a non-finite effect as not_applicable; the truths now agree.
  for (bad in c(Inf, -Inf)) {
    expect_true(is.na(truth_shift(bad, 0, 1)))
    expect_true(is.na(truth_shift(0, bad, 1)))
    expect_true(is.na(truth_surprise(bad, 1, 0)))
    expect_true(is.na(truth_surprise(0, 1, bad)))
    expect_true(is.na(truth_conclusion(bad, 0.1, 0, 0.1)))
    expect_true(is.na(truth_conclusion(0, 0.1, bad, 0.1)))
  }
  # NaN and NA already propagated to NA; they still do.
  expect_true(is.na(truth_shift(NaN, 0, 1)))
  expect_true(is.na(truth_shift(NA_real_, 0, 1)))
  expect_true(is.na(truth_conclusion(NaN, 0.1, 0, 0.1)))

  # Finite effects are untouched, including the boundary the sign test uses.
  expect_true(truth_shift(-0.30, -0.75, 0.12))
  expect_false(truth_shift(-0.30, -0.31, 0.12))
  expect_true(truth_conclusion(-0.30, 0.21, -0.75, 0.001))
  expect_true(truth_conclusion(0, 0.21, -0.75, 0.001))   # sign(0) is its own
})
