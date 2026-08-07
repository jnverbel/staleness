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
