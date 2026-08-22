#!/usr/bin/env Rscript
#
# POST HOC. Run 2026-08-22, AFTER the confirmatory results were seen.
#
# PROTOCOL.md says it plainly: "Anything run afterwards is post hoc and will be
# labelled so in the manuscript." This is that. It is not a second confirmatory
# run, it does not replace corpus/data/confirmatory.json, and no number it
# produces may be reported as pre-specified.
#
# Why it exists. The confirmatory run fixed its outcome threshold at the score
# reproducing a prevalence of 52%. That figure estimates the prevalence of the
# population the coded sample was drawn from -- 1,825 pairs parsing an effect
# at both ends -- and the detectors run on a different set: the 560 pairs whose
# two effects are COMPARABLE. Requiring comparability selects pairs whose
# abstracts repeat the same outcome, and those score low: 74% of the 560 sit in
# the bottom stratum against 33% of the sampled population. Rebuilt from its
# own stratum composition, the prevalence of the analysed set is 23%, not 52%.
#
# So the confirmatory labelling marked about 291 pairs as events where roughly
# 131 are expected, and half of those sat in the stratum where 8% of pairs are
# real events. This script asks the only question worth asking afterwards:
# does the negative result survive the correction, or was it an artefact of a
# threshold set too low?
#
#   Rscript corpus/06-posthoc-prevalence.R

suppressMessages({
  library(metafor)
  library(jsonlite)
})

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "corpus"
source(file.path(HERE, "05-confirmatory.R"), local = FALSE)

DATA <- file.path(HERE, "data")
OUT  <- file.path(DATA, "posthoc-prevalence.json")

CORRECTED_PREVALENCE <- 0.23   # rebuilt in 04-analyse-coding.py from the strata

main_posthoc <- function() {
  pairs_all <- read_pairs(file.path(DATA, "screened.jsonl"))
  pairs <- eligible_set(pairs_all)$pairs
  scores <- vapply(pairs, function(p) p$screen$score, numeric(1))
  contaminated <- CONTAMINATED_PAIRS$method

  threshold <- outcome_threshold(scores, CORRECTED_PREVALENCE)
  res <- run_all(pairs, threshold, contaminated)

  # The quantity that does not depend on the labelling at all. A detector
  # firing on 1.5% of the pairs it answers cannot reach a sensitivity worth
  # having against ANY prevalence in this range, which is the whole point.
  fire_rate <- lapply(names(res$detectors), function(d) {
    st <- vapply(res$rows, function(r) r[[d]]$state, character(1))
    fired <- vapply(res$rows, function(r) isTRUE(r[[d]]$fired), logical(1))
    list(answered = sum(st == "answered"),
         fired = sum(fired),
         rate = if (sum(st == "answered")) sum(fired) / sum(st == "answered") else NA)
  })
  names(fire_rate) <- names(res$detectors)

  res$rows <- NULL
  out <- list(
    label = "POST HOC -- run after the confirmatory results were seen",
    supersedes_nothing = "corpus/data/confirmatory.json remains the pre-specified run",
    reason = paste("confirmatory threshold used a prevalence estimated for the",
                   "1825-pair sampled population; the analysed set of 560 has a",
                   "reweighted prevalence of 23%"),
    prevalence = CORRECTED_PREVALENCE,
    threshold = threshold,
    result = res,
    firing_rate = fire_rate
  )
  write(toJSON(out, auto_unbox = TRUE, digits = 6, pretty = TRUE), OUT)

  cat("\n=== POST HOC: prevalence 0.23 instead of 0.52 ===\n")
  cat(sprintf("pairs %d | threshold %.3f | events %d (%.1f%%)\n",
              res$n_pairs, threshold, res$n_events,
              100 * res$n_events / res$n_pairs))
  for (d in names(res$detectors)) {
    m <- res$detectors[[d]]
    f <- fire_rate[[d]]
    cat(sprintf("\n%-10s %s%s\n", d, toupper(m$result),
                if (m$contaminated) "  [CONTAMINATED]" else ""))
    cat(sprintf("  sensitivity %.3f [%.3f-%.3f] on %d events\n",
                m$sensitivity, m$ci$sensitivity["lo"], m$ci$sensitivity["hi"], m$n_events))
    cat(sprintf("  specificity %.3f on %d non-events\n", m$specificity, m$n_non_events))
    cat(sprintf("  fires on %d of %d answered = %.1f%%  (label-independent)\n",
                f$fired, f$answered, 100 * f$rate))
  }
  cat("\nwritten:", OUT, "\n")
  invisible(TRUE)
}

main_posthoc()
