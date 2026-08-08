# Reproduces the calibration figures quoted for sufficiency()'s stability test.
#
# The package departs from the published sufficiency method in one place: the
# stability half is tested with a change-point statistic under a permutation
# null instead of the OLS slope of the cumulative series that the source
# describes. That departure is justified by measurement, and this script is
# the measurement. Every figure quoted in ?sufficiency, vignette("methods")
# and the JOSS paper is produced here.
#
# Run it with:
#
#   Rscript inst/calibration/calibration.R
#
# It takes about a minute and needs only the package itself. Results are
# deterministic: every experiment fixes its seeds, and the seeds are the ones
# the original measurements used.

suppressMessages(library(staleness))

ALPHA <- 0.05  # the stability cutoff sufficiency() uses by default

# Does the stability test fire -- i.e. call this evidence unstable?
#
# method = "FE" throughout. sufficiency() reads only $yi and $vi off these
# objects -- the fail-safe N and the change-point statistic are both computed
# from the study-level values, never from the fitted model -- so the estimator
# cannot affect any figure below. Fixing it to FE keeps the script from dying
# where REML's Fisher scoring fails to converge, which it does on the strongly
# heteroscedastic series of E4.
fires <- function(yi, vi, prev_idx = NULL) {
  prev_idx <- if (is.null(prev_idx)) seq_along(yi) else prev_idx
  prev <- metafor::rma(yi = yi[prev_idx], vi = vi[prev_idx], method = "FE")
  new  <- metafor::rma(yi = yi, vi = vi, method = "FE")
  isTRUE(!sufficiency(prev, new)$detail$stable)
}

report <- function(label, got, published) {
  agree <- identical(as.character(got), as.character(published))
  cat(sprintf("  %-46s measured: %-12s published: %-12s %s\n",
              label, got, published, if (agree) "MATCH" else "<-- DIFFERS"))
  invisible(agree)
}

# --- E1: no drift at all. False alarms must sit at the nominal 5%. ----------
# yi ~ N(log 0.5, 0.05^2), vi = 0.01, prev = first 20 of 30, seeds 1001-1300.
e1 <- function(seeds = 1001:1300) {
  sum(vapply(seeds, function(s) {
    set.seed(s)
    yi <- stats::rnorm(30, log(0.5), 0.05)
    fires(yi, rep(0.01, 30), prev_idx = 1:20)
  }, logical(1)))
}

# --- E2: power against a shift confined to the last 10 studies. ------------
# The canonical regime for this package: a mature review plus a smaller batch
# of new trials. The previous statistic had power 0 here, and inverted.
e2 <- function(rr, seeds = 5001:5200) {
  sum(vapply(seeds, function(s) {
    set.seed(s)
    yi <- c(stats::rnorm(20, log(0.5), 0.05), stats::rnorm(10, log(rr), 0.05))
    fires(yi, rep(0.01, 30), prev_idx = 1:20)
  }, logical(1)))
}

# --- E3: deterministic scan of WHERE the shift sits. -----------------------
# s studies at log 0.5, then 20 - s at log 0.05. The old statistic went silent
# once the shift moved past study 10.
e3 <- function(s) {
  yi <- c(rep(log(0.50), s), rep(log(0.05), 20 - s))
  fires(yi, rep(0.01, 20))
}

# --- E4: no drift, but the variance schedule is correlated with time. ------
# The failure mode a permutation null on the raw series cannot see.
e4 <- function(vi, seeds = 9001:9300) {
  sum(vapply(seeds, function(s) {
    set.seed(s)
    yi <- stats::rnorm(length(vi), log(0.5), sqrt(vi))
    fires(yi, vi, prev_idx = 1:20)
  }, logical(1)))
}

# --- The statistic that was replaced, so the comparison is checkable too. ---
# The permutation test as it stood before the change: same 999 draws and the
# same (1 + count) / (n + 1) estimator, but the statistic is the OLS slope of
# the cumulative series against accumulated information. It is reconstructed
# here from the package's own cum_drift_slope(), which survives because
# sufficiency() still reports that slope in detail$slope. Without this the
# "before" column of the comparison would be a number nobody can re-derive.
old_fires <- function(yi, vi, n_perm = 999, seed = 20260807) {
  k <- length(yi)
  # Against the study INDEX: that is what the replaced test regressed on.
  # Refitting against cumsum(1/vi) was the hypothesis that got tested and
  # refuted, not the original behaviour.
  slope_of <- function(o) abs(staleness:::cum_drift_slope(
    staleness:::cumulative_effect(yi[o], vi[o]), seq_len(k)))
  obs <- slope_of(seq_len(k))
  if (!is.finite(obs)) return(FALSE)
  set.seed(seed)
  perm <- vapply(seq_len(n_perm), function(i) slope_of(sample(k)), numeric(1))
  (1 + sum(perm >= obs, na.rm = TRUE)) / (n_perm + 1) < ALPHA
}

# The statistic before THAT one: the plain OLS t-test on the cumulative
# series, which is what the source describes literally. It is the figure the
# JOSS paper leads with, so it has to be checkable too.
ols_fires <- function(yi, vi) {
  cum <- staleness:::cumulative_effect(yi, vi)
  y <- cum[-1]
  x <- seq_along(y)
  p <- tryCatch(summary(stats::lm(y ~ x))$coefficients[2, 4],
                error = function(e) NA_real_)
  isTRUE(p < ALPHA)
}

e1_ols <- function(seeds = 1001:1300) {
  sum(vapply(seeds, function(s) {
    set.seed(s)
    ols_fires(stats::rnorm(30, log(0.5), 0.05), rep(0.01, 30))
  }, logical(1)))
}

e4_old <- function(vi, seeds = 9001:9300) {
  sum(vapply(seeds, function(s) {
    set.seed(s)
    yi <- stats::rnorm(length(vi), log(0.5), sqrt(vi))
    old_fires(yi, vi)
  }, logical(1)))
}

# The step schedule runs 50:1, the ratio vignette("methods") records for it.
# Reconstructing it at 16:1 reproduced nothing in that regime; at 50:1 the
# replaced statistic lands on its published 127/300 exactly, which is what
# identifies the generator.
STEP <- c(rep(0.50, 20), rep(0.01, 10))

ok <- c()
cat("\nE1 - no drift, 20 prior + 10 new, equal variances, 300 seeds\n")
ok <- c(ok, report("false alarms", paste0(e1(), "/300"), "15/300"))

cat("\nE2 - power against a shift in the last 10 studies, 200 seeds each\n")
for (rr in c(0.40, 0.30, 0.15, 0.02))
  ok <- c(ok, report(sprintf("power at RR %.2f", rr),
                     paste0(e2(rr), "/200"), "200/200"))

cat("\nE3 - deterministic scan: shift after study s, of 20\n")
for (s in c(2, 5, 8, 10, 12, 15, 18))
  ok <- c(ok, report(sprintf("shift after study %2d", s),
                     if (e3(s)) "fires" else "silent", "fires"))

cat("\nE4 - no drift, variance correlated with time, 300 seeds each\n")
ok <- c(ok, report("variance falling (0.16 -> 0.01)",
                   paste0(e4(seq(0.16, 0.01, length.out = 30)), "/300"), "16/300"))
ok <- c(ok, report("variance rising  (0.01 -> 0.16)",
                   paste0(e4(seq(0.01, 0.16, length.out = 30)), "/300"), "20/300"))
ok <- c(ok, report("step: 20 small then 10 large trials",
                   paste0(e4(STEP), "/300"), "20/300"))

cat("\nThe two statistics this one replaced, so the comparison is checkable\n")
ok <- c(ok, report("OLS t-test (as published), E1 no drift",
                   paste0(e1_ols(), "/300"), "209/300"))
ok <- c(ok, report("OLD slope test, variance falling",
                   paste0(e4_old(seq(0.16, 0.01, length.out = 30)), "/300"), "83/300"))
ok <- c(ok, report("OLD slope test, step schedule",
                   paste0(e4_old(STEP), "/300"), "127/300"))

cat(sprintf("\n%d of %d figures reproduced.\n", sum(ok), length(ok)))
if (!all(ok)) {
  cat("A figure that no longer reproduces is a finding, not a nuisance:\n",
      "either the statistic changed or the published number was wrong.\n")
  quit(status = 1)
}
