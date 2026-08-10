# Does requiring the shift to PERSIST separate drift from a single outlier?
#
# The change-point statistic is max_m |Z_m| over split points, and the split
# isolating one study is among them. Measured in ?sufficiency_changepoint, a
# five-SE discordant study at either end forces `unstable` in 88.5% of samples
# at k = 25 and 99.5% at k = 40. So from moderate k the detector is
# substantially a single-outlier detector, and "the evidence has drifted" and
# "one study disagrees loudly" are not the same claim.
#
# The obvious variant: refuse splits whose smaller side holds fewer than r
# studies, so a shift has to persist across r studies to count at all. r = 1 is
# the current statistic.
#
# This asks whether that buys anything, on four regimes where the right answer
# is known by construction:
#
#   null       no change at all           -> should not fire
#   drift      a real late shift           -> should fire
#   outlier    one discordant study, last  -> should NOT fire, and does
#   outlier_mid  the same study in the middle
#
# Run:  Rscript inst/persistence/persistence.R
# About two minutes. Deterministic.

suppressMessages(library(metafor))
if (!requireNamespace("staleness", quietly = TRUE)) {
  pkgload::load_all(quiet = TRUE)
}

REPS   <- 400
ALPHA  <- 0.05
NPERM  <- 499
RS     <- c(1, 2, 3, 5)
KS     <- c(20, 30, 40)
SE_STUDY <- 0.2
VI     <- SE_STUDY^2

# max |Z_m| over splits where BOTH sides hold at least r studies. r = 1 is the
# statistic the package ships.
shift_z_r <- function(yi, vi, r = 1L) {
  k <- length(yi)
  if (k < 2 * r) return(NA_real_)
  w <- 1 / vi; sw <- cumsum(w); swy <- cumsum(w * yi)
  m <- seq.int(r, k - r)
  before <- sw[m]; after <- sw[k] - before
  z <- (swy[m] / before - (swy[k] - swy[m]) / after) / sqrt(1 / before + 1 / after)
  z <- z[is.finite(z)]
  if (!length(z)) NA_real_ else max(abs(z))
}

# Same permutation null, same estimator, so only r differs between arms.
p_perm <- function(yi, vi, r, n_perm = NPERM) {
  obs <- shift_z_r(yi, vi, r)
  if (!is.finite(obs)) return(NA_real_)
  k <- length(yi)
  null <- vapply(seq_len(n_perm), function(i) {
    o <- sample.int(k); shift_z_r(yi[o], vi[o], r)
  }, numeric(1))
  (1 + sum(null >= obs, na.rm = TRUE)) / (n_perm + 1)
}

make <- function(regime, k, seed) {
  set.seed(seed)
  yi <- rnorm(k, 0, SE_STUDY)
  if (regime == "drift") {
    # A real late shift: the last third moves by 1.5 SE and STAYS moved.
    j <- seq.int(floor(2 * k / 3) + 1, k)
    yi[j] <- yi[j] + 1.5 * SE_STUDY
  } else if (regime == "outlier") {
    yi[k] <- yi[k] + 5 * SE_STUDY            # one discordant study, last
  } else if (regime == "outlier_mid") {
    yi[k %/% 2] <- yi[k %/% 2] + 5 * SE_STUDY
  }
  list(yi = yi, vi = rep(VI, k))
}

regimes <- c("null", "drift", "outlier", "outlier_mid")
rows <- list()
for (k in KS) {
  for (reg in regimes) {
    for (r in RS) {
      fired <- vapply(seq_len(REPS), function(i) {
        d <- make(reg, k, seed = 1000 * k + 97 * match(reg, regimes) + i)
        p <- p_perm(d$yi, d$vi, r)
        isTRUE(p < ALPHA)
      }, logical(1))
      rows[[length(rows) + 1]] <- data.frame(
        k = k, regime = reg, r = r, rate = mean(fired))
      cat(sprintf("  k=%2d %-12s r=%d  %.3f\n", k, reg, r, mean(fired)))
      flush.console()
    }
  }
}
res <- do.call(rbind, rows)

cat("\n=== FIRING RATE BY PERSISTENCE REQUIREMENT ===\n")
cat("null and the two outlier rows should be LOW; drift should be HIGH.\n\n")
for (k in KS) {
  cat(sprintf("k = %d\n", k))
  cat(sprintf("  %-12s %s\n", "regime", paste(sprintf("r=%-6d", RS), collapse = "")))
  for (reg in regimes) {
    v <- res[res$k == k & res$regime == reg, ]
    v <- v[match(RS, v$r), ]
    cat(sprintf("  %-12s %s\n", reg, paste(sprintf("%-8.3f", v$rate), collapse = "")))
  }
  cat("\n")
}

cat("=== WHAT IT COSTS AND WHAT IT BUYS ===\n")
for (r in RS) {
  d <- res[res$r == r, ]
  cat(sprintf("  r=%d   drift %.3f   outlier(last) %.3f   outlier(mid) %.3f   null %.3f\n",
              r,
              mean(d$rate[d$regime == "drift"]),
              mean(d$rate[d$regime == "outlier"]),
              mean(d$rate[d$regime == "outlier_mid"]),
              mean(d$rate[d$regime == "null"])))
}
saveRDS(res, "inst/persistence/persistence.rds")
