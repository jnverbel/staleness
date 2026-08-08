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

test_that("print method returns invisible and formats verdict correctly", {
  # Test out_of_date verdict with signal
  v_ood <- new_verdict("ottawa", "out_of_date", signal = 0.41, reason = "")
  expect_invisible(print(v_ood))
  expect_output(print(v_ood), "ottawa.*OUT OF DATE")
  expect_output(print(v_ood), "signal.*0\\.41")

  # Test current verdict without signal (NA signal should not print signal line)
  v_current <- new_verdict("barrowman", "current", signal = NA_real_, reason = "")
  expect_invisible(print(v_current))
  expect_output(print(v_current), "barrowman.*current")
  out_current <- capture.output(print(v_current))
  expect_false(any(grepl("signal", out_current)))

  # Test not_applicable verdict with reason (NA signal should not print signal line)
  v_na <- verdict_na("seo", "sample size too small")
  expect_invisible(print(v_na))
  expect_output(print(v_na), "seo.*not applicable")
  expect_output(print(v_na), "reason.*sample size too small")
  out_na <- capture.output(print(v_na))
  expect_false(any(grepl("signal", out_na)))
})

test_that("a non-scalar reason is refused where it is built, not where it prints", {
  # nzchar() on a length-2 vector makes print.staleness_verdict die with
  # "the condition has length > 1" -- an error about an if(), raised far from
  # the call that caused it, naming nothing the caller can act on.
  expect_error(new_verdict("rcma", "current", reason = c("a", "b")), "scalar")
  expect_error(new_verdict("rcma", "current", reason = character(0)), "scalar")
  expect_error(new_verdict("rcma", "current", reason = NA), "scalar")
  # The legitimate cases still work, and still print.
  expect_s3_class(new_verdict("rcma", "current"), "staleness_verdict")
  expect_output(print(new_verdict("rcma", "current", reason = "because")),
                "because")
})
