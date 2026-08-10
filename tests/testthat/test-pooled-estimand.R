# "Pooled sensitivity" is not one quantity, and the function used to compute
# one of the two without saying which. Summing the 2x2 counts across reviews
# answers "across all cuts in this corpus, what share of the out-of-date ones
# were flagged?", where a thirty-cut review counts six times a five-cut one.
# Averaging each review's own rate answers "for a typical review, what share?".
#
# The tests below hold three things: that the two are computed differently,
# that a corpus of unequal review lengths actually separates them, and that
# the choice is recorded in the table rather than only in the call.

mk_bt <- function(seed, n, methods = c("rcma", "ottawa")) {
  set.seed(seed)
  yi <- cumsum(stats::rnorm(n, -0.05, 0.15))
  vi <- stats::runif(n, 0.02, 0.08)
  st <- evidence_stream(metafor::rma(yi, vi, method = "FE"),
                        date = seq.int(1990, length.out = n),
                        study_id = seq_len(n))
  suppressWarnings(backtest(st, cuts = "yearly", horizon = 2, window = 3,
                            min_k = 3, seed = 1, methods = methods))
}

test_that("the weighting is reported in the table, not only chosen in the call", {
  bts <- list(mk_bt(11, 14), mk_bt(12, 14))
  for (w in c("cut", "review")) {
    res <- pooled_calibration(bts, "shift", R = 100, seed = 1, weighting = w,
                              reviews_independent = TRUE)
    expect_true("weighting" %in% names(res))
    expect_true(all(res$weighting == w))
  }
  expect_error(
    pooled_calibration(bts, "shift", R = 50, seed = 1, weighting = "study",
                       reviews_independent = TRUE),
    "'arg' should be one of")
})

test_that("cut and review weighting give different answers on unequal reviews", {
  # A long review and a short one, chosen so that both contribute a defined
  # rate: 20 and 23 eligible cuts against 5 and 7, with per-review
  # sensitivities of 0 and 0.188 against 1 and 1. Under "cut" the long review
  # dominates; under "review" they count the same. If these two ever agree to
  # the digit on a corpus this lopsided, one of them is not being computed.
  bts <- list(mk_bt(1, 28), mk_bt(14, 12))

  cut <- pooled_calibration(bts, "shift", R = 100, seed = 1, weighting = "cut",
                            reviews_independent = TRUE)
  rev <- pooled_calibration(bts, "shift", R = 100, seed = 1,
                            weighting = "review", reviews_independent = TRUE)

  expect_equal(cut$method, rev$method)
  both <- !is.na(cut$sensitivity) & !is.na(rev$sensitivity)
  skip_if_not(any(both), "no method has a defined sensitivity in this corpus")
  expect_false(isTRUE(all.equal(cut$sensitivity[both], rev$sensitivity[both])))

  # And the review-weighted figure is the plain mean of the per-review rates,
  # computed here from calibration() rather than from the function under test.
  per <- lapply(bts, function(b) calibration(b, "shift"))
  for (i in which(both)) {
    m <- rev$method[i]
    rates <- vapply(per, function(p) p$sensitivity[p$method == m], numeric(1))
    expect_equal(rev$sensitivity[i], mean(rates, na.rm = TRUE), info = m)
  }
})

test_that("a review that cannot answer is dropped, not averaged in as zero", {
  # barrowman() is not applicable on a consistently significant series, so it
  # has no rate there. Averaging that in as 0 would be an answer invented from
  # a review that gave none.
  bts <- list(mk_bt(31, 20, methods = c("rcma", "barrowman")),
              mk_bt(32, 20, methods = c("rcma", "barrowman")))
  rev <- pooled_calibration(bts, "shift", R = 50, seed = 1,
                            weighting = "review", reviews_independent = TRUE)
  b <- rev[rev$method == "barrowman", ]
  per <- vapply(bts, function(x) {
    p <- calibration(x, "shift"); p$sensitivity[p$method == "barrowman"]
  }, numeric(1))
  if (all(is.na(per))) {
    expect_true(is.na(b$sensitivity))
  } else {
    expect_equal(b$sensitivity, mean(per, na.rm = TRUE))
  }
})

test_that("reviews_independent has to be declared, and TRUE is the only value", {
  bts <- list(mk_bt(11, 14), mk_bt(12, 14))
  expect_error(pooled_calibration(bts, "shift", R = 50, seed = 1),
               "`reviews_independent` must be supplied")
  for (bad in list(FALSE, NA, "yes", 1, c(TRUE, TRUE))) {
    expect_error(pooled_calibration(bts, "shift", R = 50, seed = 1,
                                    reviews_independent = bad),
                 "must be TRUE")
  }
})
