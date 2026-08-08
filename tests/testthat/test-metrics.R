# Tests for calibration(), lead_time() and summary.staleness_backtest(): the
# metrics that turn a backtest's raw results into an answer to "does this
# method actually work?".
#
# calibration() must surface contamination (ottawa/conclusion is correct by
# construction, per CONTAMINATED_PAIRS in R/truth.R) rather than let it hide
# inside an otherwise-clean sensitivity number. lead_time() is the metric
# nobody has reported for these methods: a detector that only fires in the
# same period the evidence already moved looks perfect in a contingency
# table and is useless in practice.

fake_bt <- function() {
  res <- data.frame(
    cut     = c(2000, 2001, 2002, 2003, 2000, 2001, 2002, 2003),
    method  = rep(c("rcma", "ottawa"), each = 4),
    verdict = c("out_of_date", "out_of_date", "current", "current",
                "current", "current", "current", "current"),
    signal  = NA_real_,
    reason  = "",
    truth_shift      = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
    truth_surprise   = rep(FALSE, 8),
    truth_conclusion = rep(TRUE, 8),
    censored = rep(FALSE, 8),
    stringsAsFactors = FALSE
  )
  structure(list(results = res, methods = c("rcma", "ottawa"),
                 horizon = 5, window = 3, n_cuts = 4, n_censored = 0),
            class = "staleness_backtest")
}

test_that("calibration computes sensitivity and specificity per method", {
  cal <- calibration(fake_bt(), truth = "shift")
  rcma_row <- cal[cal$method == "rcma", ]
  expect_equal(rcma_row$sensitivity, 1)    # 2 of 2 true events caught
  expect_equal(rcma_row$specificity, 1)    # 2 of 2 non-events correct
  ott_row <- cal[cal$method == "ottawa", ]
  expect_equal(ott_row$sensitivity, 0)     # missed both
})

test_that("censored cuts are excluded from the metrics", {
  bt <- fake_bt()
  bt$results$censored[bt$results$cut == 2003] <- TRUE
  cal <- calibration(bt, truth = "shift")
  expect_equal(cal$n[cal$method == "rcma"], 3)
})

test_that("contaminated detector-truth pairs are flagged", {
  cal <- calibration(fake_bt(), truth = "conclusion")
  expect_true(cal$contaminated[cal$method == "ottawa"])
  expect_false(cal$contaminated[cal$method == "rcma"])
})

test_that("not_applicable verdicts are excluded, never counted as current", {
  bt <- fake_bt()
  bt$results$verdict[1] <- "not_applicable"
  cal <- calibration(bt, truth = "shift")
  expect_equal(cal$n[cal$method == "rcma"], 3)
})

test_that("lead_time reports how early a detector fired before the event", {
  bt <- fake_bt()
  lt <- lead_time(bt, truth = "shift")
  expect_true("median_lead" %in% names(lt))
  expect_equal(nrow(lt), 2)
})

# --- NA correction (binding correction to the brief) ------------------------
#
# truth_shift() and truth_surprise() return NA, not TRUE/FALSE, when the
# standard error they divide by is degenerate (see R/backtest.R), and
# backtest() deliberately lets that NA propagate into `results` rather than
# guessing. The brief's filter (`!censored & verdict != "not_applicable"`)
# says nothing about this case: an NA truth value flowing into
# `sum(hit & ev)` silently turns the whole cell into NA, corrupting
# sensitivity/specificity instead of merely omitting one row. A row whose
# truth column is NA must be excluded from that truth's metrics exactly like
# a censored row: never counted as a hit or a miss.

test_that("calibration excludes NA-truth rows instead of propagating NA into the metrics", {
  bt <- fake_bt()
  bt$results$truth_shift[bt$results$method == "rcma" & bt$results$cut == 2000] <- NA
  cal <- calibration(bt, truth = "shift")
  rcma_row <- cal[cal$method == "rcma", ]
  expect_equal(rcma_row$n, 3)
  expect_true(is.finite(rcma_row$sensitivity))
  expect_true(is.finite(rcma_row$specificity))
})

test_that("lead_time excludes NA-truth rows, n_events reflects the exclusion", {
  bt <- fake_bt()
  bt$results$truth_shift[bt$results$method == "rcma" & bt$results$cut == 2000] <- NA
  lt <- lead_time(bt, truth = "shift")
  rcma_row <- lt[lt$method == "rcma", ]
  expect_equal(rcma_row$n_events, 1)  # only cut 2001 remains a true event
  expect_true(is.finite(rcma_row$median_lead))
})
