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
