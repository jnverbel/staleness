test_that("ratio measures are identified", {
  expect_true(is_ratio_measure("RR"))
  expect_true(is_ratio_measure("OR"))
  expect_true(is_ratio_measure("IRR"))
  expect_true(is_ratio_measure("HR"))
  expect_true(is_ratio_measure("ROM"))
  expect_true(is_ratio_measure("PETO"))
  expect_false(is_ratio_measure("MD"))
  expect_false(is_ratio_measure("SMD"))
  expect_false(is_ratio_measure("RD"))
})

test_that("to_natural exponentiates only ratio measures", {
  # verificado empiricamente con metadat::dat.bcg: beta = -0.7145323 es log(RR)
  expect_equal(to_natural(-0.7145323, "RR"), 0.4894209, tolerance = 1e-6)
  expect_equal(to_natural(-0.7145323, "MD"), -0.7145323)
})

test_that("effect_ratio works on the natural scale for ratio measures", {
  # RR previo 0.50, RR nuevo 0.25 -> ratio 0.5, senal de rCMA
  r <- effect_ratio(log(0.25), log(0.50), "RR", se_prev = 0.1)
  expect_equal(r$ratio, 0.5, tolerance = 1e-8)
  expect_equal(r$reason, "")
})

test_that("effect_ratio refuses to divide by a near-null difference effect", {
  # MD previo 0.05 con se 0.10: indistinguible de cero, el ratio es inestable
  r <- effect_ratio(0.40, 0.05, "MD", se_prev = 0.10)
  expect_true(is.na(r$ratio))
  expect_match(r$reason, "indistinguishable from zero")
})

test_that("effect_ratio computes difference ratios when the prior effect is solid", {
  r <- effect_ratio(0.20, 0.80, "MD", se_prev = 0.05)
  expect_equal(r$ratio, 0.25, tolerance = 1e-8)
  expect_equal(r$reason, "")
})
