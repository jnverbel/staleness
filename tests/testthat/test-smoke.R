test_that("the package loads and metafor is available", {
  expect_true(requireNamespace("metafor", quietly = TRUE))
})
