# Every rate this package reports is a rate over one body of evidence, and the
# reader of a calibration table cannot see what that body was. Each of the
# eight facts below has changed an answer at some point in this package's
# history, and none of them was visible in an output a user could save.

mk <- function(n = 13, ids = seq_len(n), ni = NULL, method = "REML",
               dependent = FALSE, measure = "RR") {
  set.seed(4)
  yi <- stats::rnorm(n, -0.4, 0.3); vi <- stats::runif(n, 0.02, 0.2)
  ma <- metafor::rma(yi, vi, method = method, measure = measure)
  suppressWarnings(
    evidence_stream(ma, date = seq.int(1970, length.out = n), study_id = ids,
                    ni = ni, allow_dependence = dependent))
}

test_that("eligibility() reports the eight facts, one row, from the stream", {
  st <- mk(ni = rep(500, 13))
  e  <- eligibility(st)

  expect_s3_class(e, "data.frame")
  expect_equal(nrow(e), 1L)
  expect_named(e, c("k", "n_studies", "max_per_study", "from", "to",
                    "measure", "model", "test", "weighted", "tau2_fixed",
                    "dependent", "ni_available"))
  expect_equal(e$k, 13L)
  expect_equal(e$n_studies, 13L)
  expect_equal(e$max_per_study, 1L)
  expect_equal(e$from, 1970)
  expect_equal(e$to, 1982)
  expect_equal(e$measure, "RR")
  expect_equal(e$model, "REML")
  expect_equal(e$test, "z")
  expect_true(e$weighted)
  expect_true(is.na(e$tau2_fixed))
  expect_false(e$dependent)
  expect_equal(e$ni_available, 13L)
})

test_that("k against n_studies is the dependence question in its raw form", {
  # Three estimates from one trial. k stays 13, the study count drops, and
  # max_per_study says by how much any single trial is being over-counted.
  ids <- c("A", "A", "A", as.character(4:13))
  st  <- mk(ids = ids, dependent = TRUE)
  e   <- eligibility(st)
  expect_equal(e$k, 13L)
  expect_equal(e$n_studies, 11L)
  expect_equal(e$max_per_study, 3L)
  expect_true(e$dependent)
})

test_that("the model the snapshots are refitted under is reported, not assumed", {
  # FE and REML disagree on the same evidence by factors that change verdicts,
  # and Knapp-Hartung moves the p-values ottawa() reads directly.
  expect_equal(eligibility(mk(method = "FE"))$model, "FE")

  set.seed(4)
  yi <- stats::rnorm(13, -0.4, 0.3); vi <- stats::runif(13, 0.02, 0.2)
  kn <- metafor::rma(yi, vi, test = "knha", measure = "RR")
  expect_equal(eligibility(evidence_stream(kn, date = 1970:1982,
                                           study_id = 1:13))$test, "knha")

  fx <- metafor::rma(yi, vi, tau2 = 0.05, measure = "RR")
  expect_equal(eligibility(evidence_stream(fx, date = 1970:1982,
                                           study_id = 1:13))$tau2_fixed, 0.05)

  uw <- metafor::rma(yi, vi, weighted = FALSE, measure = "RR")
  expect_false(eligibility(evidence_stream(uw, date = 1970:1982,
                                           study_id = 1:13))$weighted)
})

test_that("ni_available counts sizes, because barrowman sums them", {
  # Absent sample sizes make barrowman() answer not_applicable at every cut,
  # which in a results table looks exactly like a detector that never fired.
  expect_equal(eligibility(mk())$ni_available, 0L)
  expect_equal(eligibility(mk(ni = rep(500, 13)))$ni_available, 13L)

  # The case the count exists for, and the one a presence check cannot see: an
  # `ni` metafor derived on its own, with holes in it. A supplied `ni` with an
  # NA is refused as a malformed argument; a derived one is a fact about the
  # dataset and is let through, so "sample sizes are available" is a matter of
  # how many. Without this, reporting length(ni) instead of the finite count
  # passed every other test in this file.
  set.seed(4)
  yi <- stats::rnorm(13, -0.4, 0.3); vi <- stats::runif(13, 0.02, 0.2)
  ma <- metafor::rma(yi, vi, measure = "RR")
  ma$ni <- c(rep(500, 10), NA, NA, NA)
  st <- evidence_stream(ma, date = 1970:1982, study_id = 1:13)
  expect_equal(eligibility(st)$ni_available, 10L)
  expect_equal(eligibility(st)$k, 13L)
  # And print() reports the shortfall rather than implying full coverage.
  expect_true(any(grepl("10 of 13 studies", utils::capture.output(print(st)))))
})

test_that("a backtest answers for the stream it ran on", {
  st <- mk(ni = rep(500, 13))
  bt <- suppressWarnings(backtest(st, cuts = "yearly", horizon = 2,
                                  window = 3, min_k = 3, seed = 1,
                                  methods = c("rcma", "ottawa")))
  expect_equal(eligibility(bt), eligibility(st))
  expect_error(eligibility(bt$results), "staleness_stream")
})

test_that("print() shows the series, and says when it is dependent", {
  out <- utils::capture.output(print(mk(ni = rep(500, 13))))
  expect_true(any(grepl("13 from 13 studies", out)))
  expect_true(any(grepl("1970 to 1982", out)))
  expect_true(any(grepl("REML", out)))
  expect_true(any(grepl("13 of 13 studies", out)))

  none <- utils::capture.output(print(mk()))
  expect_true(any(grepl("barrowman cannot answer", none)))

  st <- mk(ids = c("A", "A", as.character(3:13)), dependent = TRUE)
  bt <- suppressWarnings(backtest(st, cuts = "yearly", horizon = 2,
                                  window = 3, min_k = 3, seed = 1,
                                  methods = c("rcma", "ottawa")))
  ob <- utils::capture.output(print(bt))
  expect_true(any(grepl("13 estimates from 12 studies", ob)))
  expect_true(any(grepl("dependence allowed", ob)))
  # The claim the note makes is bounded: it says the rates cannot be read as
  # coming from independent studies, not that they are all biased upwards.
  # An earlier version said "every rate below is optimistic", which asserted
  # a direction dependence does not guarantee.
  expect_true(any(grepl("independent studies", ob)))
  expect_false(any(grepl("optimistic", ob)))
  # And the target, which decides what the rates below it even mean.
  expect_true(any(grepl("target", ob)))
})
