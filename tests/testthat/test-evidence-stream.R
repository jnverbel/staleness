make_bcg_stream <- function() {
    es <- bcg_es()
  ma  <- metafor::rma(yi, vi, data = es)
  evidence_stream(ma, date = es$year, study_id = seq_along(es$year))
}

test_that("a stream is built from an rma object and sorted by date", {
  skip_if_not_installed("metadat")
  s <- make_bcg_stream()
  expect_s3_class(s, "staleness_stream")
  expect_equal(s$k, 13)
  expect_equal(s$measure, "RR")
  expect_false(is.unsorted(s$date))
})

test_that("missing dates are an explicit error, never imputed", {
  skip_if_not_installed("metadat")
    es <- bcg_es()
  ma  <- metafor::rma(yi, vi, data = es)
  bad <- es$year; bad[3] <- NA
  expect_error(evidence_stream(ma, date = bad, study_id = seq_along(bad)), "missing")
})

test_that("date length must match the number of studies", {
  skip_if_not_installed("metadat")
    es <- bcg_es()
  ma  <- metafor::rma(yi, vi, data = es)
  expect_error(evidence_stream(ma, date = c(1950, 1960), study_id = seq_along(c(1950, 1960))), "length")
})

test_that("snapshot_at refits using only studies up to the cut", {
  skip_if_not_installed("metadat")
  s  <- make_bcg_stream()
  m1 <- snapshot_at(s, 1960)
  expect_true(m1$k < s$k)
  expect_equal(m1$k, sum(s$date <= 1960))
  expect_equal(m1$measure, "RR")
})

test_that("window_between returns studies in the half-open interval", {
  skip_if_not_installed("metadat")
  s <- make_bcg_stream()
  w <- window_between(s, 1960, 1970)
  expect_equal(w$k, sum(s$date > 1960 & s$date <= 1970))
})

test_that("a snapshot with fewer than 2 studies is refused", {
  skip_if_not_installed("metadat")
  s <- make_bcg_stream()
  expect_error(snapshot_at(s, min(s$date) - 1), "at least 2")
})

# --- ni is validated for NA, like date -------------------------------------
#
# evidence_stream() checked `ni`'s length but not for missing values, two lines
# after checking `date` for exactly that. A single NA survived into
# backtest()'s sum(stream$ni[...]) -- sum() over an NA is NA -- and surfaced,
# many frames later, as barrowman() stopping with "needs the sample size of
# both the prior meta-analysis and the new studies". That message describes a
# call that supplied no sample sizes at all, which is not what happened.

test_that("evidence_stream refuses missing sample sizes, with a message that says so", {
  ma <- metafor::rma(yi = rep(log(0.5), 6), vi = rep(0.02, 6), measure = "RR")
  expect_error(
    evidence_stream(ma, date = 2000:2005, study_id = seq_along(2000:2005), ni = c(100, 100, NA, 100, 100, 100)),
    "missing values"
  )
  # And the same stream with a complete `ni` is accepted, so the check is not
  # rejecting the shape of the input.
  expect_s3_class(evidence_stream(ma, date = 2000:2005, study_id = seq_along(2000:2005), ni = rep(100, 6)),
                  "staleness_stream")
})

test_that("an NA in an explicitly supplied ni is still refused", {
  ma <- metafor::rma(yi = c(0.1, 0.2, 0.15, 0.05), vi = rep(0.05, 4))
  expect_error(
    evidence_stream(ma, date = 2001:2004, study_id = seq_along(2001:2004), ni = c(100, NA, 120, 130)),
    "never imputed"
  )
})

test_that("an NA in metafor's own ni does not block the whole stream", {
  # ni is filled in from ma$ni when the caller does not supply it, so a
  # dataset where one study never reported its sample size used to kill
  # construction outright -- taking the four detectors that never look at ni
  # down with it. barrowman() already answers "not_applicable" with the right
  # reason for a non-finite n, so the hard error was both redundant and more
  # destructive than the thing it prevented.
  d <- data.frame(ai = c(4, 5, 6, 7), bi = c(100, 110, 120, 130),
                  ci = c(11, 12, 13, 14), di = c(90, 95, 100, 105))
  es <- metafor::escalc(measure = "RR", ai = ai, bi = bi, ci = ci, di = di,
                        data = d)
  ma <- metafor::rma(yi, vi, data = es)
  ma$ni[2] <- NA                      # one study never reported its n

  expect_no_error(evidence_stream(ma, date = 2001:2004, study_id = seq_along(2001:2004)))
  stream <- evidence_stream(ma, date = 2001:2004, study_id = seq_along(2001:2004))
  expect_s3_class(stream, "staleness_stream")
  expect_true(anyNA(stream$ni))

  # And the detector that does need it says so precisely, rather than
  # claiming no sample sizes were supplied at all.
  prev <- snapshot_at(stream, 2002)
  res  <- barrowman(prev, n_prev = sum(stream$ni[1:2]), n_new = 250)
  expect_equal(res$verdict, "not_applicable")
  expect_match(res$reason, "finite")
})

test_that("a meta-regression is refused rather than silently flattened", {
  # rma(mods = ~ x) fits a meta-regression: beta is a vector of coefficients,
  # not a pooled effect. Every snapshot here is refitted WITHOUT moderators,
  # so accepting one would hand back a plain pooled analysis under the label
  # of the model the caller supplied. Different question, same-looking answer.
  es <- bcg_es()
  mr <- metafor::rma(yi, vi, mods = ~ ablat, data = es)
  expect_gt(length(mr$beta), 1)
  expect_error(evidence_stream(mr, date = es$year, study_id = seq_along(es$year)), "moderator")
})

test_that("the test statistic the caller chose survives into every snapshot", {
  # test = "knha" changes the p-value by orders of magnitude, and ottawa()
  # decides on p-values. Refitting snapshots with the default z test would
  # score the caller's evidence under a test they did not ask for.
  es <- bcg_es()
  kn <- metafor::rma(yi, vi, data = es, test = "knha")
  st <- evidence_stream(kn, date = es$year, study_id = seq_along(es$year))
  expect_equal(st$test, "knha")

  snap <- snapshot_at(st, 1975)
  expect_equal(snap$test, "knha")

  direct <- metafor::rma(yi, vi, data = es[es$year <= 1975, ], test = "knha")
  expect_equal(snap$pval, direct$pval, tolerance = 1e-10)

  # And the default is untouched for callers who never asked for anything.
  plain <- evidence_stream(metafor::rma(yi, vi, data = es), date = es$year, study_id = seq_along(es$year))
  expect_equal(snapshot_at(plain, 1975)$test, "z")
})

test_that("options that change the estimator survive into the snapshots", {
  # method and test were not the only ones. weighted = FALSE and a fixed tau2
  # both change beta, se and pval, and a snapshot refitted without them scores
  # the caller's evidence under a model they did not fit.
  es <- bcg_es()
  keep <- es$year <= 1975

  unw <- metafor::rma(yi, vi, data = es, weighted = FALSE)
  s1  <- snapshot_at(evidence_stream(unw, date = es$year, study_id = seq_along(es$year)), 1975)
  d1  <- metafor::rma(yi, vi, data = es[keep, ], weighted = FALSE)
  expect_false(s1$weighted)
  expect_equal(as.numeric(s1$beta), as.numeric(d1$beta), tolerance = 1e-10)

  fx <- metafor::rma(yi, vi, data = es, tau2 = 0.1)
  s2 <- snapshot_at(evidence_stream(fx, date = es$year, study_id = seq_along(es$year)), 1975)
  d2 <- metafor::rma(yi, vi, data = es[keep, ], tau2 = 0.1)
  expect_equal(s2$tau2, 0.1)
  expect_equal(as.numeric(s2$beta), as.numeric(d2$beta), tolerance = 1e-10)

  # Custom per-study weights cannot be carried through a subset sensibly, so
  # they are refused rather than dropped.
  wt <- metafor::rma(yi, vi, data = es, weights = rep(1, nrow(es)))
  expect_error(evidence_stream(wt, date = es$year, study_id = seq_along(es$year)), "weights")
})

test_that("evidence_stream requires an identifier and refuses dependence", {
  # Every row was treated as an independent study, and nothing in the stream
  # could tell otherwise: it held yi, vi, a date and optionally ni, with no
  # identity attached. Several outcomes, time points or arms from one trial
  # entered as separate studies -- their participants summed twice in
  # barrowman(), their weights counted twice in every pooled estimate.
  #
  # This was known before it was reported: the metadat sweep in
  # inst/applicability/ had to exclude seven datasets with several effects per
  # study BY HAND. Having to do that by hand was the defect.
  es <- bcg_es()
  ma <- metafor::rma(yi, vi, data = es, measure = "RR", method = "FE")

  expect_error(evidence_stream(ma, date = es$year), "study_id")
  expect_error(evidence_stream(ma, date = es$year,
                               study_id = seq_len(3)), "length")
  expect_error(evidence_stream(ma, date = es$year,
                               study_id = c(NA, seq_len(nrow(es) - 1))),
               "missing")

  # Duplicates are refused by default, and the message says which studies and
  # what it would have cost.
  dup <- c("trial A", "trial A", seq_len(nrow(es) - 2))
  expect_error(evidence_stream(ma, date = es$year, study_id = dup),
               "more than one estimate")
  expect_error(evidence_stream(ma, date = es$year, study_id = dup),
               "barrowman")

  # An identifier rather than a boolean promising independence, because a
  # promise cannot be checked and an identifier can. Distinct ids pass.
  st <- evidence_stream(ma, date = es$year, study_id = seq_len(nrow(es)))
  expect_s3_class(st, "staleness_stream")
  expect_false(st$dependent)
  expect_equal(length(st$study_id), nrow(es))
})

test_that("allowed dependence is recorded and announced, not forgotten", {
  es <- bcg_es()
  ma <- metafor::rma(yi, vi, data = es, measure = "RR", method = "FE")
  dup <- c("trial A", "trial A", seq_len(nrow(es) - 2))

  # Announced at the moment the decision is taken. Before this, opting in was
  # silent and the only trace was in print(stream) -- which the analyst who
  # receives a results table never calls.
  expect_warning(
    st <- evidence_stream(ma, date = es$year, study_id = dup,
                          allow_dependence = TRUE),
    "1 study contributes more than one estimate")
  # What the warning may claim: standard errors too small and intervals too
  # narrow -- that direction follows from an effective k below the nominal
  # one. What it may NOT claim is that every point rate moves upwards.
  expect_warning(
    evidence_stream(ma, date = es$year, study_id = dup,
                    allow_dependence = TRUE),
    "too narrow")
  expect_warning(
    evidence_stream(ma, date = es$year, study_id = dup,
                    allow_dependence = TRUE),
    "biased in either direction")
  expect_true(st$dependent)

  # And an independent stream says nothing at all.
  expect_no_warning(
    evidence_stream(ma, date = es$year, study_id = seq_len(nrow(es))))

  # print() is the one place a reader might notice, so it says so there.
  out <- utils::capture.output(print(st))
  expect_true(any(grepl("dependence was allowed", out)))
  expect_true(any(grepl("intervals too narrow", out)))
  expect_true(any(grepl("either direction", out)))
  # barrowman() is the one place the direction IS known, because it sums
  # participants across a snapshot, so the note is allowed to say so.
  expect_true(any(grepl("barrowman", out)))
  expect_false(any(grepl("optimistic", out)))

  # And the ids travel in the stream, sorted with everything else, so a caller
  # can check what got pooled with what.
  expect_equal(length(st$study_id), nrow(es))
  expect_equal(length(unique(st$study_id)), nrow(es) - 1)
})
