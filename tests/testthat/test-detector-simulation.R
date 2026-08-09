test_that("simulation does not apply to an already significant meta-analysis", {
  prev <- metafor::rma(yi = rep(log(0.4), 8), vi = rep(0.01, 8), measure = "RR")
  new_ev <- list(yi = rep(log(0.4), 3), vi = rep(0.01, 3), k = 3)
  v <- simulation(prev, new_ev, B = 50, seed = 1)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "already significant")
})

test_that("simulation needs at least one recent study", {
  prev <- metafor::rma(yi = c(0.02, -0.01, 0.03, -0.02), vi = rep(0.05, 4),
                       measure = "MD")
  v <- simulation(prev, list(yi = numeric(0), vi = numeric(0), k = 0),
                  B = 50, seed = 1)
  expect_equal(v$verdict, "not_applicable")
  expect_match(v$reason, "no recent studies")
})

test_that("simulation is reproducible under a fixed seed", {
  prev <- metafor::rma(yi = c(0.20, 0.10, 0.25, 0.05), vi = rep(0.05, 4),
                       measure = "MD")
  new_ev <- list(yi = c(0.22, 0.18), vi = c(0.02, 0.02), k = 2)
  a <- simulation(prev, new_ev, B = 200, seed = 42)
  b <- simulation(prev, new_ev, B = 200, seed = 42)
  expect_equal(a$signal, b$signal)
})

test_that("simulation reports power as its signal, bounded in [0, 1]", {
  prev <- metafor::rma(yi = c(0.20, 0.10, 0.25, 0.05), vi = rep(0.05, 4),
                       measure = "MD")
  new_ev <- list(yi = c(0.22, 0.18), vi = c(0.02, 0.02), k = 2)
  v <- simulation(prev, new_ev, B = 200, seed = 7)
  expect_gte(v$signal, 0)
  expect_lte(v$signal, 1)
  # Pinned to the verdict this fixture actually produces. The assertion used
  # to read `v$verdict %in% c("out_of_date", "current")`, which passes
  # whichever way the detector answers and so tests nothing about it: the
  # power here is 0.57, comfortably under the 0.80 threshold, and stays
  # `current` across seeds 1 to 5.
  expect_equal(v$verdict, "current")
  expect_lt(v$signal, 0.80)
})

test_that("high power crosses the threshold and signals", {
  # A solid prior trend plus precise new studies should reach significance
  # often. vi = 0.16 (not the original 0.04) keeps the prior meta-analysis
  # NOT significant (p = 0.137), which the method requires as a precondition
  # -- at vi = 0.04 the prior itself is already significant (p = 0.0029),
  # which would make this fixture hit the "not_applicable" branch instead of
  # exercising high power. True power under this fixture is 1.0 across
  # several seeds tried during design, comfortably above the 0.8 threshold.
  prev <- metafor::rma(yi = c(0.30, 0.28, 0.32, 0.29), vi = rep(0.16, 4),
                       measure = "MD")
  new_ev <- list(yi = c(0.30, 0.31, 0.29), vi = rep(0.005, 3), k = 3)
  v <- simulation(prev, new_ev, B = 500, seed = 3)
  expect_gt(v$signal, 0.8)
  expect_equal(v$verdict, "out_of_date")
})

# --- The caller's random stream is global state and must be left alone ------
#
# simulation() called set.seed() directly and then drew B * k_new normals, both
# of which write to .Random.seed in the global environment. Verified before the
# fix: `set.seed(123); runif(1)` gave 0.2875775, while
# `set.seed(123); simulation(...); runif(1)` gave 0.9074913. Any script that
# seeds once and then calls check_currency() or backtest() lost the
# reproducibility of everything downstream, and CRAN policy forbids altering
# global state. The seed is now applied inside with_preserved_seed(), which
# saves .Random.seed and restores it on exit.

test_that("simulation leaves the caller's random stream untouched, with a seed", {
  prev <- metafor::rma(yi = c(0.20, 0.10, 0.25, 0.05), vi = rep(0.05, 4),
                       measure = "MD")
  new_ev <- list(yi = c(0.22, 0.18), vi = c(0.02, 0.02), k = 2)
  set.seed(123); expected <- runif(3)
  set.seed(123); invisible(simulation(prev, new_ev, B = 200, seed = 42))
  got <- runif(3)
  expect_equal(got, expected)
})

test_that("simulation leaves the caller's random stream untouched, without a seed", {
  # seed = NULL consumed B * k_new draws from whatever stream the caller was
  # on. Unreproducible results are the documented price of seed = NULL;
  # silently advancing someone else's stream is not.
  prev <- metafor::rma(yi = c(0.20, 0.10, 0.25, 0.05), vi = rep(0.05, 4),
                       measure = "MD")
  new_ev <- list(yi = c(0.22, 0.18), vi = c(0.02, 0.02), k = 2)
  set.seed(123); expected <- runif(3)
  set.seed(123); invisible(simulation(prev, new_ev, B = 200, seed = NULL))
  got <- runif(3)
  expect_equal(got, expected)
})

test_that("a whole backtest leaves the caller's random stream untouched", {
  # The case that actually bites: seed once, backtest, and expect everything
  # afterwards to still be reproducible.
  yi <- c(0.05, 0.08, -0.02, 0.10, 0.04, 0.12, 0.15, 0.20, 0.18, 0.25,
          0.30, 0.35, 0.40, 0.45, 0.50)
  s <- evidence_stream(metafor::rma(yi = yi, vi = rep(0.05, length(yi))),
                       date = 2000:2014, ni = rep(100, length(yi)))
  set.seed(123); expected <- runif(3)
  set.seed(123); invisible(backtest(s, horizon = 5, seed = 4))
  got <- runif(3)
  expect_equal(got, expected)
})

test_that("without a seed, consecutive calls actually differ", {
  # Restoring the stream unconditionally made every seedless call start from
  # the same state, so simulation() returned byte-identical draws while its
  # documentation promised the opposite. In backtest(seed = NULL) that made
  # the Monte Carlo error perfectly correlated across cuts instead of
  # independent.
  prev <- metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
                       vi = c(0.16, 0.20, 0.18, 0.15, 0.22))
  new  <- list(yi = c(-0.45, -0.38, -0.52, -0.29, -0.41),
               vi = c(0.05, 0.04, 0.06, 0.05, 0.04), k = 5)
  set.seed(11)
  got <- replicate(6, simulation(prev, new, B = 200, seed = NULL)$signal)
  expect_gt(length(unique(got)), 1)
})

test_that("with a seed, the result repeats and the caller's stream survives", {
  prev <- metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
                       vi = c(0.16, 0.20, 0.18, 0.15, 0.22))
  new  <- list(yi = c(-0.45, -0.38, -0.52, -0.29, -0.41),
               vi = c(0.05, 0.04, 0.06, 0.05, 0.04), k = 5)
  a <- simulation(prev, new, B = 200, seed = 7)$signal
  b <- simulation(prev, new, B = 200, seed = 7)$signal
  expect_identical(a, b)

  set.seed(99); expected <- runif(1)
  set.seed(99); invisible(simulation(prev, new, B = 200, seed = 7))
  expect_identical(runif(1), expected)
})

test_that("a session with no random stream is left with none", {
  # The function whose stated purpose is to leave no footprint was creating
  # .Random.seed in globalenv() where there had been none.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(with_preserved_seed(stats::runif(1), seed = 3))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the power threshold is strict, as the source states", {
  # The vignette said "at or above 80%" while the code compares with `>`,
  # following the source's "Power >80%". The code was right and the prose was
  # wrong, but nothing could tell them apart, because no test touched the
  # boundary.
  #
  # And the boundary cannot be reached from the data: simulated power is
  # hits/B, so building evidence that lands on exactly 0.80 is guesswork, and
  # a value that misses by any amount makes `>` and `>=` behave identically.
  # Move the THRESHOLD to the datum instead -- it is a parameter, and a
  # parameter is exact.
  prev <- metafor::rma(yi = c(-0.20, -0.35, 0.05, -0.30, -0.10),
                       vi = c(0.16, 0.20, 0.18, 0.15, 0.22), measure = "RR")
  new_ev <- list(k = 3, yi = c(-0.9, -0.8, -0.85), vi = rep(0.02, 3))

  observed <- simulation(prev, new_ev, B = 200, seed = 7)$signal
  expect_equal(observed, 0.7)

  # Threshold exactly equal to the power: strict `>` says current.
  at <- simulation(prev, new_ev, B = 200, seed = 7,
                   power_threshold = observed)
  expect_equal(at$verdict, "current")

  # One replicate below it, the same power fires. The pair is what makes the
  # operator visible; either case alone passes under both `>` and `>=`.
  below <- simulation(prev, new_ev, B = 200, seed = 7,
                      power_threshold = observed - 1 / 200)
  expect_equal(below$verdict, "out_of_date")
})
