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
  expect_true(v$verdict %in% c("out_of_date", "current"))
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
