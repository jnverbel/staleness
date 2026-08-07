test_that("a verdict carries method, verdict, signal and reason", {
  v <- new_verdict("ottawa", "out_of_date", signal = 0.41, reason = "")
  expect_s3_class(v, "staleness_verdict")
  expect_equal(v$method, "ottawa")
  expect_equal(v$verdict, "out_of_date")
  expect_equal(v$signal, 0.41)
})

test_that("invalid verdict values are rejected", {
  expect_error(new_verdict("ottawa", "maybe"), "must be one of")
})

test_that("a not_applicable verdict requires a reason", {
  expect_error(verdict_na("barrowman", ""), "reason")
  v <- verdict_na("barrowman", "prior meta-analysis was already significant")
  expect_equal(v$verdict, "not_applicable")
  expect_true(is.na(v$signal))
})
