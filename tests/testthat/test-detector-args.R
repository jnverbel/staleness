# The engine's arguments were validated in an earlier pass; the detectors'
# were not. Sweeping all five turned up the same defect in every one, and in
# two of them it is silent: ottawa(alpha = NA) returned a verdict of
# out_of_date and sufficiency(alpha_stability = NA) returned current, both
# from an argument that means nothing. A wrong verdict that looks valid is
# worse than an error, which is why these are errors.
#
# The distinction each detector must keep: a malformed CALL is an error, while
# a datum that is present but unusable -- an NA sample size, a non-finite
# p-value -- is a fact about the evidence and yields "not_applicable" with a
# reason. barrowman() below tests both sides of that line.

skip_if_not_installed("metafor")

fit_prev <- function() {
  metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
               vi = c(0.16, 0.20, 0.18, 0.15, 0.22), measure = "RR")
}
fit_upd <- function() {
  metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10, -0.9, -0.8),
               vi = c(0.16, 0.20, 0.18, 0.15, 0.22, 0.05, 0.05), measure = "RR")
}
new_ev <- function() list(k = 3, yi = c(-0.9, -0.8, -0.85),
                          vi = c(0.05, 0.05, 0.05))

test_that("barrowman refuses impossible sample sizes", {
  prev <- fit_prev()
  # n_prev = 0 made n_required 0, so the ratio was Inf and the verdict came
  # back out_of_date with a signal no reader can act on. A negative size
  # produced a negative ratio, which is never > 1, so it always read current:
  # the failure was silent and biased in one direction.
  expect_error(barrowman(prev, n_prev = 0, n_new = 30), "positive")
  expect_error(barrowman(prev, n_prev = -555, n_new = 2265), "positive")
  expect_error(barrowman(prev, n_prev = 555, n_new = -100), "negative")
  expect_error(barrowman(prev, n_prev = c(1, 2), n_new = 30), "single")
  expect_error(barrowman(prev, n_prev = "555", n_new = 30), "single")
})

test_that("barrowman still answers when zero participants are new", {
  # Zero new participants is not malformed: it is the honest statement that
  # nothing has been added, and the correct answer to it is "current".
  prev <- fit_prev()
  v <- barrowman(prev, n_prev = 555, n_new = 0)
  expect_equal(v$verdict, "current")
  expect_equal(v$signal, 0)
})

test_that("barrowman keeps unusable data separate from a malformed call", {
  # This is the line the new validation must not blur. NA is a fact about the
  # evidence, so it is not_applicable; zero is a fact about the call.
  prev <- fit_prev()
  v <- barrowman(prev, n_prev = NA_real_, n_new = 2265)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "finite")
  expect_error(barrowman(prev, n_prev = 555, n_new = 2265, alpha = NA), "0 and 1")
  expect_error(barrowman(prev, n_prev = 555, n_new = 2265, alpha = 1.5), "0 and 1")
  expect_error(barrowman(prev, n_prev = 555, n_new = 2265, z_crit = -1), "positive")
})

test_that("simulation refuses replicate counts it cannot run", {
  prev <- fit_prev()
  # B = 0 ran no replicates, so power was 0/0 and the comparison threw R's
  # own "missing value where TRUE/FALSE needed" -- an internal error, not an
  # explanation. B = c(10, 20) silently used the first element.
  expect_error(simulation(prev, new_ev(), B = 0), "at least 1")
  expect_error(simulation(prev, new_ev(), B = -10), "at least 1")
  expect_error(simulation(prev, new_ev(), B = 2.5), "whole number")
  expect_error(simulation(prev, new_ev(), B = c(10, 20)), "single")
  expect_error(simulation(prev, new_ev(), B = NA), "single")
})

test_that("simulation refuses probabilities that are not probabilities", {
  prev <- fit_prev()
  expect_error(simulation(prev, new_ev(), B = 10, alpha = NA), "0 and 1")
  expect_error(simulation(prev, new_ev(), B = 10, alpha = 0), "0 and 1")
  expect_error(simulation(prev, new_ev(), B = 10, alpha = 1), "0 and 1")
  # power_threshold may sit at either end: 0 means every power fires, 1 means
  # none does. Both are coherent, so the interval is closed here and open for
  # alpha.
  expect_error(simulation(prev, new_ev(), B = 10, power_threshold = NA), "0 and 1")
  expect_error(simulation(prev, new_ev(), B = 10, power_threshold = 1.5), "0 and 1")
  expect_silent(simulation(prev, new_ev(), B = 10, power_threshold = 0, seed = 1))
  expect_silent(simulation(prev, new_ev(), B = 10, power_threshold = 1, seed = 1))
})

test_that("ottawa refuses an alpha that is not a probability", {
  # The worst of the five: ottawa(alpha = NA) returned out_of_date. Both
  # significance comparisons became NA, and the code downstream read them as
  # a change. Nothing warned.
  prev <- fit_prev(); upd <- fit_upd()
  expect_error(ottawa(prev, upd, alpha = NA), "0 and 1")
  expect_error(ottawa(prev, upd, alpha = 2), "0 and 1")
  expect_error(ottawa(prev, upd, alpha = c(0.04, 0.05)), "single")
})

test_that("sufficiency refuses arguments it cannot honour", {
  # sufficiency(alpha_stability = NA) returned current, equally silently.
  prev <- fit_prev(); upd <- fit_upd()
  expect_error(sufficiency(prev, upd, alpha_stability = NA), "0 and 1")
  expect_error(sufficiency(prev, upd, min_k = 0), "at least 1")
  expect_error(sufficiency(prev, upd, min_k = 2.5), "whole number")
  expect_error(sufficiency(prev, upd, n_perm = 0), "at least 1")
  expect_error(sufficiency(prev, upd, n_perm = NA), "single")
})

test_that("rcma refuses thresholds that cannot bracket a ratio", {
  prev <- fit_prev(); upd <- fit_upd()
  expect_error(rcma(prev, upd, upper = NA), "single")
  expect_error(rcma(prev, upd, lower = -1), "positive")
  expect_error(rcma(prev, upd, lower = 1.5, upper = 0.5), "below")
})

test_that("seed is validated wherever it enters", {
  # NEWS claimed all five detectors validated their arguments, and `seed` was
  # the exception that made the claim false. set.seed() truncates towards
  # zero, so seed = 1.5 was accepted in silence and produced the *same* stream
  # as seed = 1: two values a reader would record as different runs, giving
  # byte-identical results. A vector was accepted too, using its first element.
  prev <- fit_prev(); upd <- fit_upd()
  expect_error(simulation(prev, new_ev(), B = 10, seed = 1.5), "whole number")
  expect_error(simulation(prev, new_ev(), B = 10, seed = c(1, 2)), "whole number")
  expect_error(sufficiency(prev, upd, seed = 1.5), "whole number")
  expect_error(backtest(structure(list(), class = "staleness_stream"),
                        seed = 1.5), "whole number")

  # NULL stays valid: it is the documented way to ask for an unseeded run.
  expect_silent(simulation(prev, new_ev(), B = 10, seed = NULL))

  # And the check must fire even when the detector never reaches its RNG. Both
  # of these return not_applicable long before with_preserved_seed() is
  # called, so validating only at the point of use would miss them.
  significant <- metafor::rma(yi = rep(log(0.50), 4), vi = rep(0.05, 4),
                              measure = "RR")
  expect_error(simulation(significant, new_ev(), B = 10, seed = 1.5),
               "whole number")
  expect_error(sufficiency(prev, upd, min_k = 500, seed = 1.5), "whole number")
})

test_that("seed is bounded by the range set.seed() actually accepts", {
  # check_seed() promised "a whole number" and set.seed() accepts less than
  # that: past R's integer range it throws "supplied seed is not a valid
  # integer" with a coercion warning. The promise was wider than the code, so
  # the check handed the failure downstream instead of explaining it.
  #
  # The bound is +/- .Machine$integer.max, NOT the textbook int32 range: R
  # reserves -2147483648 for NA_integer_, so set.seed(-2147483648) fails too.
  # Verified by running both ends rather than assumed from the type.
  #
  # Unlike the power threshold, this boundary CAN be reached from the datum:
  # 2147483647 and 2147483648 are both exactly representable as doubles, so
  # the pair makes the operator visible.
  lim <- .Machine$integer.max
  expect_silent(check_seed(lim))
  expect_silent(check_seed(-lim))
  expect_error(check_seed(lim + 1), "whole number")
  expect_error(check_seed(-lim - 1), "whole number")

  prev <- fit_prev()
  expect_silent(simulation(prev, new_ev(), B = 5, seed = lim))
  expect_error(simulation(prev, new_ev(), B = 5, seed = lim + 1),
               "whole number")
  # And the error must come from us, not from set.seed() downstream.
  expect_error(simulation(prev, new_ev(), B = 5, seed = lim + 1),
               "set\\.seed", fixed = FALSE)
})

test_that("evidence_stream requires dates it can cut yearly", {
  # `date` was checked for length and for NAs, never for type. A Date vector
  # went through as.numeric() into DAYS SINCE 1970, so cuts = "yearly" stepped
  # one DAY at a time: ten studies three days apart produced 21 cuts numbered
  # 11329, 11330, ... and `horizon` and `window` silently became days too.
  # Everything here is denominated in years -- cuts, horizon, window,
  # lead_time -- so a Date is refused with the conversion spelled out rather
  # than reinterpreted.
  ma <- metafor::rma(yi = rnorm(6, -0.3, 0.1), vi = runif(6, 0.02, 0.05),
                     measure = "RR")
  d <- as.Date("2001-01-01") + seq(0, by = 3, length.out = 6)
  # The generic "must be numeric" guard would already reject a Date, so these
  # require the SPECIFIC message: one that names the unit and spells out the
  # conversion. Mutation showed the difference -- dropping the Date branch left
  # the tests green while the user lost the only sentence that tells them what
  # to do about it.
  expect_error(evidence_stream(ma, date = d, study_id = seq_along(d)), "as\\.numeric\\(format\\(")
  expect_error(evidence_stream(ma, date = as.POSIXct(d), study_id = seq_along(as.POSIXct(d))),
               "as\\.numeric\\(format\\(")
  expect_error(evidence_stream(ma, date = factor(2001:2006), study_id = seq_along(factor(2001:2006))), "numeric")
  expect_error(evidence_stream(ma, date = as.character(2001:2006), study_id = seq_along(as.character(2001:2006))), "numeric")
  # Years still work, including fractional ones.
  expect_silent(evidence_stream(ma, date = 2001:2006, study_id = seq_along(2001:2006)))
  expect_silent(evidence_stream(ma, date = seq(2001, 2003.5, length.out = 6), study_id = seq_along(seq(2001, 2003.5, length.out = 6))))
})

test_that("rcma and ottawa refuse models on different scales", {
  # Both read $measure off `prev` to decide which branch to take, and neither
  # checked that `new_ma` was on the same one. A risk ratio compared against a
  # mean difference returned a signal of 2.98 and a verdict of out_of_date --
  # a number with no meaning, presented like every other verdict.
  prev <- metafor::rma(yi = rep(log(0.5), 4), vi = rep(0.05, 4), measure = "RR")
  upd  <- metafor::rma(yi = rep(0.40, 4), vi = rep(0.05, 4), measure = "MD")
  expect_error(rcma(prev, upd), "same")
  expect_error(ottawa(prev, upd), "same")
  same <- metafor::rma(yi = rep(log(0.9), 4), vi = rep(0.05, 4), measure = "RR")
  expect_silent(rcma(prev, same))
})

test_that("ottawa refuses an unknown qualitative signal", {
  # nzchar(NA_character_) is TRUE, so an analyst recording "I don't know" for a
  # qualitative signal got out_of_date. Unknown is not the same as present,
  # and this is an argument rather than a datum, so it is refused.
  prev <- metafor::rma(yi = rep(log(0.5), 4), vi = rep(0.05, 4), measure = "RR")
  upd  <- metafor::rma(yi = c(rep(log(0.5), 4), rep(log(0.52), 3)),
                       vi = c(rep(0.05, 4), rep(0.02, 3)), measure = "RR")
  expect_equal(ottawa(prev, upd)$verdict, "current")
  expect_error(ottawa(prev, upd, qualitative = NA_character_), "missing")
  expect_error(ottawa(prev, upd, qualitative = c("harm", NA)), "missing")
  # An empty string is not a signal, and still is not.
  expect_equal(ottawa(prev, upd, qualitative = "")$verdict, "current")
  expect_equal(ottawa(prev, upd, qualitative = "substantial harm")$verdict,
               "out_of_date")
})

test_that("the truth functions return one logical, as they document", {
  # ?truth says each returns "a single logical value". truth_shift() happily
  # returned a vector of two, and truth_conclusion() died with R's own
  # "'length = 2' in coercion to 'logical(1)'". The documentation was stronger
  # than the code -- again.
  expect_error(truth_shift(c(-0.3, -0.2), -0.75, 0.12), "single")
  expect_error(truth_surprise(c(-0.3, -0.2), 0.3, -0.75), "single")
  expect_error(truth_conclusion(c(-0.3, -0.2), c(0.2, 0.3), -0.75, 0.001),
               "single")
  expect_error(truth_shift(-0.3, -0.75, c(0.12, 0.2)), "single")

  # Scalars still behave exactly as before.
  expect_true(truth_shift(-0.30, -0.75, 0.12))
  expect_true(truth_conclusion(-0.30, 0.21, -0.75, 0.001))
  expect_true(is.na(truth_shift(-0.30, -0.75, 0)))
})

# --- Contract sweep -------------------------------------------------------
#
# Four outside reviews in a row found the same shape of defect: a contract
# stated in the documentation and not enforced by the code. None of them can
# surface in a green suite, because the suite exercises the package correctly.
#
# So the contract was swept in bulk instead of one round at a time: every
# exported function, every argument, eight malformed inputs each, 384
# combinations, classified as "our error", "R's internal error" or "no
# complaint at all". These are what the sweep found.

test_that("ottawa's qualitative signals must be text", {
  # nzchar() coerces, so ANY non-empty value counted as a declared signal:
  # ottawa(qualitative = 0) fired, and so did a list. The argument carries an
  # analyst's judgement in words; a number is not one.
  prev <- metafor::rma(yi = rep(log(0.5), 4), vi = rep(0.05, 4), measure = "RR")
  upd  <- metafor::rma(yi = c(rep(log(0.5), 4), rep(log(0.52), 3)),
                       vi = c(rep(0.05, 4), rep(0.02, 3)), measure = "RR")
  expect_error(ottawa(prev, upd, qualitative = 0), "character")
  expect_error(ottawa(prev, upd, qualitative = c(1, 2)), "character")
  expect_error(ottawa(prev, upd, qualitative = list(1)), "character")
  expect_error(ottawa(prev, upd, qualitative = Inf), "character")
  expect_equal(ottawa(prev, upd, qualitative = "substantial harm")$verdict,
               "out_of_date")
})

test_that("the truth functions validate their thresholds", {
  # threshold and alpha were never checked. A negative threshold makes every
  # comparison TRUE, an NA makes every one NA, and truth_conclusion(alpha = NA)
  # still returned a verdict -- an evaluation target computed from a cutoff that is
  # not one.
  expect_error(truth_shift(-0.3, -0.75, 0.12, threshold = -1), "positive")
  expect_error(truth_shift(-0.3, -0.75, 0.12, threshold = NA), "positive")
  expect_error(truth_shift(-0.3, -0.75, 0.12, threshold = c(1, 2)), "positive")
  expect_error(truth_surprise(-0.3, 0.3, -0.75, threshold = Inf), "positive")
  expect_error(truth_conclusion(-0.3, 0.21, -0.75, 0.001, alpha = NA), "0 and 1")
  expect_error(truth_conclusion(-0.3, 0.21, -0.75, 0.001, alpha = 2), "0 and 1")
  expect_error(truth_conclusion(-0.3, 0.21, -0.75, 0.001, alpha = c(0.05, 0.1)),
               "0 and 1")
  # Defaults and sane values unchanged.
  expect_true(truth_shift(-0.30, -0.75, 0.12))
  expect_true(truth_conclusion(-0.30, 0.21, -0.75, 0.001, alpha = 0.05))
})

test_that("window_between and snapshot_at validate their cut points", {
  es <- bcg_es()
  ma <- metafor::rma(yi, vi, data = es, measure = "RR", method = "FE")
  st <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year))

  # from/to went straight into a comparison against the date vector. NA made
  # every element NA and the subset came back malformed; a character compared
  # lexically; a length-two vector recycled silently against 13 dates.
  expect_error(window_between(st, NA, 1970), "single")
  expect_error(window_between(st, 1960, NA), "single")
  expect_error(window_between(st, c(1960, 1965), 1970), "single")
  expect_error(window_between(st, "1960", 1970), "single")
  expect_error(window_between(st, Inf, 1970), "single")
  # An interval that runs backwards is empty by construction, and saying so
  # beats returning k = 0 as though the evidence had been examined.
  expect_error(window_between(st, 1975, 1970), "before")

  expect_error(snapshot_at(st, NA), "single")
  expect_error(snapshot_at(st, c(1970, 1975)), "single")
  expect_error(snapshot_at(st, "1970"), "single")
  expect_error(snapshot_at(st, Inf), "single")

  # Valid calls still work.
  expect_equal(window_between(st, 1960, 1970)$k,
               sum(es$year > 1960 & es$year <= 1970))
  expect_equal(snapshot_at(st, 1970)$k, sum(es$year <= 1970))
})

test_that("evidence_stream refuses infinite years at the entry point", {
  # anyNA() catches NA and NaN but not the infinities. An infinite year used to
  # reach backtest() before seq() rejected it with "'to' must be a finite
  # number" -- R's complaint about its own argument, three functions from the
  # call that caused it. In between, the stream looked usable:
  # snapshot_at(st, 2003) returned a k as though nothing were wrong.
  #
  # This one escaped the contract sweep because the sweep mutated arguments
  # with SCALARS, and for a vector argument the length check fires first and
  # hides the content. Mutating in place, keeping the length, is what finds it.
  ma <- metafor::rma(yi = rnorm(6, -0.3, 0.1), vi = runif(6, 0.02, 0.05),
                     measure = "RR")
  expect_error(evidence_stream(ma, date = c(2001:2005, Inf), study_id = seq_along(c(2001:2005, Inf))), "infinite")
  expect_error(evidence_stream(ma, date = c(-Inf, 2002:2006), study_id = seq_along(c(-Inf, 2002:2006))), "infinite")
  # NA and NaN keep their own message, which says something different.
  expect_error(evidence_stream(ma, date = c(NA, 2002:2006), study_id = seq_along(c(NA, 2002:2006))), "missing")
  expect_error(evidence_stream(ma, date = c(NaN, 2002:2006), study_id = seq_along(c(NaN, 2002:2006))), "missing")
  expect_silent(evidence_stream(ma, date = 2001:2006, study_id = seq_along(2001:2006)))
})

test_that("supplied sample sizes and cut points must be real numbers", {
  # Found by mutating vector arguments IN PLACE, keeping their length -- the
  # blind spot the scalar sweep had, since for a vector argument the length
  # check fires first and hides the content.
  es <- bcg_es()
  ma <- metafor::rma(yi, vi, data = es, measure = "RR", method = "FE")
  ni <- with(es, tpos + tneg + cpos + cneg)

  # barrowman() sums ni across a snapshot, so a negative element silently
  # shrinks the total it divides by and an infinite one makes it Inf.
  expect_error(evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = replace(ni, 1, Inf)),
               "finite")
  expect_error(evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = replace(ni, 1, -5)),
               "positive")
  expect_error(evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = replace(ni, 1, 0)),
               "positive")

  # But a DERIVED ni is a fact about the dataset, not about the call, and must
  # still not take down the four detectors that never read it. That line was
  # drawn deliberately for NA and holds here too.
  st <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year))
  expect_s3_class(st, "staleness_stream")

  # An infinite cut used to surface as "needs at least 3 uncensored cuts":
  # the error named a consequence and never the Inf that caused it.
  st2 <- evidence_stream(ma, date = es$year, study_id = seq_along(es$year), ni = ni)
  expect_error(backtest(st2, cuts = c(1960, 1965, Inf), horizon = 3,
                        window = 5, min_k = 3), "infinite")
  expect_silent(suppressWarnings(
    backtest(st2, cuts = c(1960, 1965, 1970), horizon = 3, window = 5,
             min_k = 3, seed = 1)))
})

test_that("the exported detectors require the model class they document", {
  # Every one of them documents an rma.uni and none checked for it. Four died
  # with R's own "argument is of length zero" or "missing value where
  # TRUE/FALSE needed", and sufficiency() was worse: it returned a verdict of
  # not_applicable from an empty list, as though it had examined something.
  for (f in list(rcma, ottawa, sufficiency)) {
    expect_error(f(list(), list()), "rma.uni")
  }
  expect_error(barrowman(list(), 100, 50), "rma.uni")
  expect_error(simulation(list(), list(k = 1, yi = 0.1, vi = 0.1), B = 5),
               "rma.uni")
  expect_error(check_currency(list(), list(yi = 0.1, vi = 0.1, k = 1),
                              methods = "rcma"), "rma.uni")
})

test_that("metafor subclasses are refused by the detectors, as by the stream", {
  # evidence_stream() has always refused rma.mh with a reasoned message: it
  # refits each snapshot with rma(), which needs yi and vi, and a
  # Mantel-Haenszel fit cannot be reproduced from those. The detectors accepted
  # the same object and returned "current".
  #
  # That inconsistency matters here specifically. Reproducing a Cochrane review
  # to the digit requires Mantel-Haenszel, so rma.mh is what a user arrives
  # with -- and they got a reasoned refusal from the stream and a silent
  # verdict from the detectors.
  skip_if_not_installed("metadat")
  mh <- metafor::rma.mh(measure = "OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i,
                        data = metadat::dat.lau1992)
  expect_false(inherits(mh, "rma.uni"))
  expect_error(rcma(mh, mh), "rma.mh")
  expect_error(ottawa(mh, mh), "rma.mh")
  expect_error(sufficiency(mh, mh), "rma.mh")
  # The stream's message is the one being made consistent, not replaced.
  expect_error(evidence_stream(mh, date = metadat::dat.lau1992$year, study_id = seq_along(metadat::dat.lau1992$year)), "rma.mh")

  # And what the package builds internally is still rma.uni, so nothing in the
  # engine trips over this.
  es <- bcg_es()
  st <- evidence_stream(metafor::rma(yi, vi, data = es, measure = "RR",
                                     method = "FE"), date = es$year, study_id = seq_along(es$year))
  expect_s3_class(snapshot_at(st, 1970), "rma.uni")
})

test_that("a factor cannot run one detector and be labelled another", {
  # The worst of the sweep. `methods` reaches switch(), and switch() on a
  # factor uses its INTEGER CODE, not its label. factor("ottawa") has one
  # level, so its code is 1 and switch() took the first branch:
  #
  #   check_currency(methods = factor("ottawa"))
  #     name in the list : ottawa
  #     verdict$method   : rcma
  #     signal           : 1.6272   (ottawa on the same data gives 0.8432)
  #
  # Two identities in one object, and plausible numbers under the wrong label.
  # In backtest() the results came back holding rcma rows while ottawa had
  # been asked for, with nothing to say so.
  prev <- metafor::rma(yi = rep(log(0.20), 4), vi = rep(0.05, 4), measure = "RR")
  nev  <- list(yi = rep(log(0.5), 4), vi = rep(0.02, 4), k = 4)
  expect_error(check_currency(prev, nev, methods = factor("ottawa")),
               "character")
  expect_error(check_currency(prev, nev, methods = c("rcma", NA)), "missing")

  es <- bcg_es()
  st <- evidence_stream(metafor::rma(yi, vi, data = es, measure = "RR",
                                     method = "FE"), date = es$year, study_id = seq_along(es$year))
  expect_error(backtest(st, methods = factor("ottawa")), "character")
  expect_error(backtest(st, methods = c("rcma", NA)), "missing")

  # The label and the detector must agree when it is asked properly.
  r <- check_currency(prev, nev, methods = "ottawa")
  expect_equal(names(r$verdicts), "ottawa")
  expect_equal(r$verdicts[[1]]$method, "ottawa")
})

test_that("truth_conclusion refuses p-values that cannot exist", {
  # p_t = -1 read as significant and returned TRUE: an evaluation target
  # manufactured from a number outside [0, 1]. Impossible is refused; NA stays
  # a datum and yields NA, as in the other two truths.
  expect_error(truth_conclusion(0.1, -1, 0.1, 0.10), "p-value")
  expect_error(truth_conclusion(0.1, 2, 0.1, 0.10), "p-value")
  expect_error(truth_conclusion(0.1, Inf, 0.1, 0.10), "p-value")
  expect_error(truth_conclusion(0.1, 0.21, 0.1, -1), "p-value")
  # The boundaries themselves are legitimate p-values.
  expect_false(is.na(truth_conclusion(0.1, 0, 0.1, 1)))
  expect_true(truth_conclusion(-0.30, 0.21, -0.75, 0.001))
})

test_that("the stream and backtest consumers check what they were handed", {
  # window_between(list(), 1, 2) returned k = 0 in silence, which reads as
  # "we looked and found no new evidence" rather than "that was not a stream".
  # snapshot_at(list(), 1) was worse: it reported "found 0 at cut 1", a
  # statement about evidence that had never been examined.
  expect_error(window_between(list(), 1, 2), "staleness_stream")
  expect_error(snapshot_at(list(), 1), "staleness_stream")
  # calibration() and lead_time() died with R's "invalid argument type".
  expect_error(calibration(list(), "shift"), "staleness_backtest")
  expect_error(lead_time(list(), "shift"), "staleness_backtest")
})

test_that("min_k counts studies, so it is a whole number", {
  # min_k = 2.1, 2.5 and 2.9 were all accepted and all behaved exactly like 3
  # -- 23 cuts where min_k = 2 gives 27 -- so two callers writing different
  # numbers got the same backtest and neither could tell which from the object.
  es <- bcg_es()
  st <- evidence_stream(metafor::rma(yi, vi, data = es, measure = "RR",
                                     method = "FE"), date = es$year, study_id = seq_along(es$year))
  for (bad in c(2.1, 2.5, 2.9)) {
    expect_error(backtest(st, min_k = bad), "whole number", info = bad)
  }
  expect_error(backtest(st, min_k = 1), "whole number")
  expect_silent(suppressWarnings(backtest(st, cuts = "yearly", horizon = 3,
                                          window = 5, min_k = 3, seed = 1)))
})
