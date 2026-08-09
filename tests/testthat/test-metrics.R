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
  # rcma: events at 2000 and 2001, fired out_of_date at both. Event 2000's
  # earliest on-time firing is 2000 itself (lead 0); event 2001's earliest
  # on-time firing is 2000 (lead 1). median(c(0, 1)) = 0.5.
  rcma_row <- lt[lt$method == "rcma", ]
  expect_equal(rcma_row$median_lead, 0.5)
  expect_equal(rcma_row$n_events, 2)
})

# --- Fix: lead_time() must look at every true event, not just the first ----
#
# Coordinator review found the first implementation collapsed to
# first_event <- min(events) and only ever checked firings at or before that
# single earliest event. With true events at 2000, 2001 and 2002 and the
# detector firing out_of_date only at 2002, that version returned
# n_events = 3 but median_lead = NA: a firing that genuinely anticipated a
# later event was invisible, and the column promised a median that was
# never computed. Fixed to compute one lead per true event (NA when that
# event was never caught in time) and take the median of the defined leads.

test_that("lead_time takes the median across all true events, not just the first", {
  # The exact scenario from review: true events at cuts 2000, 2001 and 2002;
  # the detector only ever fires out_of_date at 2002. The two earlier events
  # were missed outright (no on-time firing, contribute no lead at all); the
  # 2002 event was caught in the very period it happened, contributing a
  # lead of 0 -- a real value, not a stand-in for "missed".
  res <- data.frame(
    cut     = c(2000, 2001, 2002),
    method  = "rcma",
    verdict = c("current", "current", "out_of_date"),
    signal  = NA_real_,
    reason  = "",
    truth_shift      = c(TRUE, TRUE, TRUE),
    truth_surprise   = FALSE,
    truth_conclusion = FALSE,
    censored = FALSE,
    stringsAsFactors = FALSE
  )
  bt <- structure(list(results = res, methods = "rcma",
                        horizon = 5, window = 3, n_cuts = 3, n_censored = 0),
                   class = "staleness_backtest")
  lt <- lead_time(bt, truth = "shift")
  expect_equal(lt$median_lead, 0)
  expect_equal(lt$n_events, 3)
})

test_that("lead_time reports the actual lead when a detector fires well before the event", {
  res <- data.frame(
    cut     = c(2000, 2001, 2002, 2003, 2004),
    method  = "rcma",
    verdict = c("out_of_date", "current", "current", "current", "current"),
    signal  = NA_real_,
    reason  = "",
    truth_shift      = c(FALSE, FALSE, FALSE, FALSE, TRUE),
    truth_surprise   = FALSE,
    truth_conclusion = FALSE,
    censored = FALSE,
    stringsAsFactors = FALSE
  )
  bt <- structure(list(results = res, methods = "rcma",
                        horizon = 5, window = 3, n_cuts = 5, n_censored = 0),
                   class = "staleness_backtest")
  lt <- lead_time(bt, truth = "shift")
  expect_equal(lt$median_lead, 4)  # fired at 2000, the single event is at 2004
  expect_equal(lt$n_events, 1)
})

test_that("lead_time returns NA when the detector never fires in time; n_events still counts the misses", {
  bt <- fake_bt()
  lt <- lead_time(bt, truth = "shift")
  ott_row <- lt[lt$method == "ottawa", ]
  expect_true(is.na(ott_row$median_lead))
  expect_equal(ott_row$n_events, 2)  # both true events counted, neither caught in time
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
  # rcma at cut 2000: was a true positive (out_of_date, truth TRUE); making
  # its truth NA removes that row entirely, leaving one true positive (2001)
  # and two true negatives (2002, 2003) rather than two and two.
  bt$results$truth_shift[bt$results$method == "rcma" & bt$results$cut == 2000] <- NA
  cal <- calibration(bt, truth = "shift")
  rcma_row <- cal[cal$method == "rcma", ]
  expect_equal(rcma_row$n, 3)
  expect_true(is.finite(rcma_row$sensitivity))
  expect_true(is.finite(rcma_row$specificity))
  expect_equal(rcma_row$sensitivity, 1)  # the one remaining true event, caught
  expect_equal(rcma_row$specificity, 1)  # the two remaining non-events, correct
})

test_that("an NA truth removes the event but not the firing at that cut", {
  bt <- fake_bt()
  bt$results$truth_shift[bt$results$method == "rcma" & bt$results$cut == 2000] <- NA
  lt <- lead_time(bt, truth = "shift")
  rcma_row <- lt[lt$method == "rcma", ]

  # 2000 stops being a true event, because unknown is not true.
  expect_equal(rcma_row$n_events, 1)

  # But rcma DID fire at 2000, and the surviving event is at 2001, so the
  # lead is 1. This assertion used to read 0, which is what the shared
  # eligibility filter produced by deleting the 2000 firing along with its
  # unknown truth -- and 0 does not mean "we lost the firing", it means "the
  # detector only spoke once the evidence had already moved". The two are
  # opposite verdicts on the detector, so the old expectation stated the
  # defect rather than the intent.
  expect_equal(rcma_row$median_lead, 1)
})

# --- Every requested method gets a row, even one that never applied ---------
#
# calibration() and lead_time() iterated unique(res$method) -- the POST-filter
# set -- while promising "one row per method". A detector that was
# not_applicable at every cut simply vanished from the table. On the shipped
# BCG backtest that dropped 2 of the 5 detectors from summary(bt), and the
# backtesting vignette had to explain their absence in prose. "This detector
# never applied to this evidence" is itself a result about the detector: it
# belongs in the table as a row with n = 0, not as an absence the reader has to
# notice. Both functions now iterate bt$methods, which the backtest object
# already carries.

na_everywhere_bt <- function() {
  res <- data.frame(
    cut     = c(2000, 2001, 2002, 2000, 2001, 2002),
    method  = rep(c("rcma", "barrowman"), each = 3),
    verdict = c("out_of_date", "current", "current",
                "not_applicable", "not_applicable", "not_applicable"),
    signal  = NA_real_,
    reason  = "",
    truth_shift      = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
    truth_surprise   = FALSE,
    truth_conclusion = FALSE,
    censored = FALSE,
    stringsAsFactors = FALSE
  )
  structure(list(results = res, methods = c("rcma", "barrowman"),
                 horizon = 5, window = 3, n_cuts = 3, n_censored = 0),
            class = "staleness_backtest")
}

test_that("calibration keeps a row for a method that was never applicable", {
  cal <- calibration(na_everywhere_bt(), truth = "shift")
  expect_setequal(cal$method, c("rcma", "barrowman"))
  b <- cal[cal$method == "barrowman", ]
  expect_equal(nrow(b), 1)
  expect_equal(b$n, 0)
  expect_true(is.na(b$sensitivity))
  expect_true(is.na(b$specificity))
  expect_true(is.na(b$false_alarm))
  # The method that did answer is unaffected.
  expect_equal(cal$sensitivity[cal$method == "rcma"], 1)
})

test_that("lead_time keeps a row for a method that was never applicable", {
  lt <- lead_time(na_everywhere_bt(), truth = "shift")
  expect_setequal(lt$method, c("rcma", "barrowman"))
  b <- lt[lt$method == "barrowman", ]
  expect_equal(b$n_events, 0)
  expect_true(is.na(b$median_lead))
})

test_that("summary() reports all five methods across all three truths", {
  bt <- na_everywhere_bt()
  s <- summary(bt)
  expect_equal(nrow(s), 2 * length(available_truths()))
  for (t in available_truths()) {
    expect_setequal(s$method[s$truth == t], bt$methods)
  }
})

test_that("rows come out in the order the backtest requested the methods", {
  bt <- na_everywhere_bt()
  expect_equal(calibration(bt, "shift")$method, bt$methods)
  expect_equal(lead_time(bt, "shift")$method, bt$methods)
})

test_that("a firing still counts as early warning when its own truth is NA", {
  # lead_time() is the only metric that relates DIFFERENT rows: an event in
  # one row, the firing that preceded it in another. eligible_rows() was
  # built for calibration(), which scores each row in isolation and therefore
  # must drop a row whose truth is unknown. Reused here it dropped the
  # firings too, and a firing is an observed fact about the detector -- its
  # role is to precede a later event, not to be scored against its own cut.
  mk <- function(res) structure(
    list(results = res, methods = "rcma", horizon = 3, window = 3,
         n_cuts = nrow(res), n_censored = 0),
    class = "staleness_backtest")

  res <- data.frame(
    cut = c(2000, 2001), method = "rcma",
    verdict = c("out_of_date", "current"), signal = c(2.1, 1.0), reason = "",
    truth_shift = c(NA, TRUE), truth_surprise = c(NA, TRUE),
    truth_conclusion = c(NA, TRUE), censored = FALSE,
    stringsAsFactors = FALSE)

  lt <- lead_time(mk(res))
  expect_equal(lt$n_events, 1)
  expect_equal(lt$median_lead, 1)   # was NA: the 2000 firing had been erased

  # The identical history with 2000's truth merely FALSE must agree. The two
  # differ in what is known about 2000, not in what the detector did.
  res_known <- res; res_known$truth_shift[1] <- FALSE
  expect_equal(lead_time(mk(res_known))$median_lead, lt$median_lead)
})

test_that("an unknown truth still cannot manufacture an event", {
  # The other half of the asymmetry: NA must never be read as TRUE. Events
  # require a known truth; only the firings are exempt.
  mk <- function(res) structure(
    list(results = res, methods = "rcma", horizon = 3, window = 3,
         n_cuts = nrow(res), n_censored = 0),
    class = "staleness_backtest")
  res <- data.frame(
    cut = c(2000, 2001), method = "rcma",
    verdict = c("out_of_date", "current"), signal = c(2.1, 1.0), reason = "",
    truth_shift = c(NA, NA), truth_surprise = c(NA, NA),
    truth_conclusion = c(NA, NA), censored = FALSE,
    stringsAsFactors = FALSE)
  lt <- lead_time(mk(res))
  expect_equal(lt$n_events, 0)
  expect_true(is.na(lt$median_lead))
})

test_that("an ancient firing counts as warning only when `within` allows it", {
  # This fixture used to demonstrate the defect and call it a definition: with
  # `within = Inf`, a single firing nine years before an event counted as
  # advance warning of it, so a detector that fired once in 2000 and went
  # quiet scored the SAME median lead as one that warned before every event.
  # The test recorded that and the default kept it. The default is now the
  # backtest's own horizon, and `Inf` has to be asked for.
  mk <- function(v, t) structure(list(
    results = data.frame(cut = 2000:2009, method = "rcma", verdict = v,
                         signal = NA_real_, reason = "", truth_shift = t,
                         truth_surprise = t, truth_conclusion = t,
                         censored = FALSE, stringsAsFactors = FALSE),
    methods = "rcma", horizon = 3, window = 1, n_cuts = 10, n_censored = 0),
    class = "staleness_backtest")
  ev <- c(rep(FALSE, 5), TRUE, FALSE, TRUE, FALSE, TRUE)  # 2005, 2007, 2009

  once  <- mk(c("out_of_date", rep("current", 9)), ev)   # fires only in 2000
  always <- mk(rep("out_of_date", 10), ev)

  # Asked for explicitly, the old behaviour is unchanged and still available:
  # "did this detector ever fire before the event" is a real question.
  expect_equal(lead_time(once, within = Inf)$median_lead, 7)
  expect_equal(lead_time(always, within = Inf)$median_lead, 7)

  # By default the two are no longer indistinguishable. The lone firing is too
  # old to be warning of anything, and drops out; the detector that warns each
  # time keeps its lead.
  expect_true(is.na(lead_time(once)$median_lead))
  expect_equal(lead_time(always)$median_lead, 3)
  expect_equal(lead_time(once)$n_events, lead_time(always)$n_events)
  # ...while its sensitivity is zero, which is why the two must be read together.
  expect_equal(calibration(once)$sensitivity, 0)

  # Bounded, the stale warning stops counting.
  expect_true(is.na(lead_time(once, within = 3)$median_lead))
  expect_equal(lead_time(always, within = 3)$median_lead, 3)

  expect_error(lead_time(once, within = 0), "within")
  expect_error(lead_time(once, within = -1), "within")
})

# --- Rates from one series describe it; they do not estimate anything -------

mk_stream <- function(seed) {
  set.seed(seed)
  yi <- cumsum(stats::rnorm(14, -0.05, 0.15))
  vi <- stats::runif(14, 0.02, 0.08)
  evidence_stream(metafor::rma(yi, vi, method = "FE"), date = 1990:2003,
                  study_id = seq_len(14))
}
mk_bt <- function(seed) {
  suppressWarnings(backtest(mk_stream(seed), cuts = "yearly", horizon = 2,
                            window = 3, min_k = 3, seed = 1,
                            methods = c("rcma", "ottawa")))
}

test_that("pooled_calibration resamples reviews, not cuts", {
  # calibration()'s n counts cuts, and consecutive cuts of one review share
  # almost every study, so a binomial interval on that denominator would be
  # far too narrow. Independence lives at the level of reviews: different
  # reviews share no studies, so those are what get resampled.
  bts <- lapply(1:4, mk_bt)
  p <- pooled_calibration(bts, "shift", R = 200, seed = 1)

  expect_true(all(c("n_reviews", "sens_lo", "sens_hi",
                    "spec_lo", "spec_hi") %in% names(p)))
  expect_lte(max(p$n_reviews), length(bts))
  # The interval brackets the point estimate it belongs to.
  ok <- !is.na(p$specificity) & !is.na(p$spec_lo)
  expect_true(all(p$spec_lo[ok] <= p$specificity[ok] + 1e-9))
  expect_true(all(p$spec_hi[ok] >= p$specificity[ok] - 1e-9))

  # Reproducible under a seed, and the caller's stream survives.
  expect_equal(pooled_calibration(bts, "shift", R = 200, seed = 1),
               pooled_calibration(bts, "shift", R = 200, seed = 1))
  set.seed(5); before <- runif(2)
  set.seed(5); invisible(pooled_calibration(bts, "shift", R = 100, seed = 2))
  expect_equal(runif(2), before)
})

test_that("one review gets a rate but no interval", {
  # An interval drawn from a single review describes that review and nothing
  # else, so it is withheld rather than printed narrow. The rate itself is
  # still reported, and agrees with calibration() on the same series.
  bts <- lapply(1:4, mk_bt)
  one <- pooled_calibration(bts[2], "shift", R = 100, seed = 1)
  expect_true(all(one$n_reviews == 1))
  expect_true(all(is.na(one$spec_lo)))
  expect_true(all(is.na(one$sens_lo)))

  solo <- calibration(bts[[2]], "shift")
  expect_equal(one$specificity[one$method == "rcma"],
               solo$specificity[solo$method == "rcma"])
  expect_equal(one$specificity[one$method == "ottawa"],
               solo$specificity[solo$method == "ottawa"])
})

test_that("backtests answering different questions are not pooled", {
  # "final" and "horizon" score a cut against different targets. Pooling them
  # would average two different questions into one number.
  a <- mk_bt(1)
  b <- suppressWarnings(backtest(mk_stream(2), cuts = "yearly", horizon = 2,
                                 window = 3, min_k = 3, seed = 1,
                                 methods = c("rcma", "ottawa"),
                                 truth_target = "horizon"))
  expect_error(pooled_calibration(list(a, b), "shift", R = 50, seed = 1),
               "different")
  expect_error(pooled_calibration(list(), "shift"), "non-empty")
  expect_error(pooled_calibration(list(a, "not a backtest"), "shift"),
               "staleness_backtest")
})
