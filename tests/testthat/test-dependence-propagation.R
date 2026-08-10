# allow_dependence = TRUE used to be a fact about the stream and nothing else.
# The analyst who receives a calibration table receives numbers computed as if
# every estimate were a separate trial, with no column saying otherwise -- and
# the rates are optimistic, because the denominator counts one trial several
# times. The flag now travels: stream -> backtest -> calibration/lead_time,
# and pooled_calibration(, reviews_independent = TRUE) withholds its interval rather than printing one
# that resampling reviews cannot justify.

dep_stream <- function(n = 14, seed = 11, dependent = TRUE) {
  set.seed(seed)
  yi <- cumsum(stats::rnorm(n, -0.05, 0.15))
  vi <- stats::runif(n, 0.02, 0.08)
  ma <- metafor::rma(yi, vi, method = "FE")
  ids <- if (dependent) c(1L, 1L, seq.int(2L, n - 1L)) else seq_len(n)
  suppressWarnings(
    evidence_stream(ma, date = seq.int(1990, length.out = n), study_id = ids,
                    allow_dependence = dependent)
  )
}

dep_backtest <- function(...) {
  suppressWarnings(backtest(dep_stream(...), cuts = "yearly", horizon = 2,
                            window = 3, min_k = 3, seed = 1,
                            methods = c("rcma", "ottawa")))
}

test_that("the dependence flag reaches the backtest and both metrics tables", {
  bt <- dep_backtest()
  expect_true(bt$dependent)

  cal <- calibration(bt)
  expect_true("dependent" %in% names(cal))
  expect_true(all(cal$dependent))

  lt <- lead_time(bt)
  expect_true("dependent" %in% names(lt))
  expect_true(all(lt$dependent))

  # summary() stacks calibration(), so it inherits the column.
  expect_true(all(summary(bt)$dependent))
})

test_that("an independent stream reports FALSE, so the column discriminates", {
  bt <- dep_backtest(dependent = FALSE)
  expect_false(bt$dependent)
  expect_false(any(calibration(bt)$dependent))
  expect_false(any(lead_time(bt)$dependent))
})

test_that("pooled_calibration(, reviews_independent = TRUE) withholds its interval over dependent reviews", {
  bts <- list(dep_backtest(seed = 11), dep_backtest(seed = 12))

  expect_warning(res <- pooled_calibration(bts, "shift", R = 100, seed = 1, reviews_independent = TRUE),
                 "accept_dependence")
  # Descriptive figures are still computed -- the caller asked for them and
  # they describe what happened. It is the interval that is withheld.
  expect_true(all(res$dependent))
  expect_true(any(!is.na(res$sensitivity) | !is.na(res$specificity)))
  expect_true(all(is.na(res$sens_lo)))
  expect_true(all(is.na(res$sens_hi)))
  expect_true(all(is.na(res$spec_lo)))
  expect_true(all(is.na(res$spec_hi)))

  # Expert use is not blocked: say so in the call and the bounds come back.
  ok <- pooled_calibration(bts, "shift", R = 100, seed = 1,
                           accept_dependence = TRUE, reviews_independent = TRUE)
  expect_true(all(ok$dependent))
  expect_true(any(!is.na(ok$sens_lo)))
})

test_that("independent reviews keep their interval and their FALSE flag", {
  bts <- list(dep_backtest(seed = 11, dependent = FALSE),
              dep_backtest(seed = 12, dependent = FALSE))
  res <- expect_no_warning(
    pooled_calibration(bts, "shift", R = 100, seed = 1,
                       reviews_independent = TRUE))
  expect_false(any(res$dependent))
  expect_true(any(!is.na(res$sens_lo)))
})

test_that("accept_dependence refuses anything that is not one TRUE or FALSE", {
  bts <- list(dep_backtest(seed = 11, dependent = FALSE))
  for (bad in list(NA, "yes", 1, c(TRUE, TRUE), NULL)) {
    expect_error(pooled_calibration(bts, "shift", R = 20, seed = 1,
                                    accept_dependence = bad,
                                    reviews_independent = TRUE),
                 "must be TRUE or FALSE")
  }
})
