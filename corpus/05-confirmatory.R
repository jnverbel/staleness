#!/usr/bin/env Rscript
#
# Confirmatory evaluation of the detectors on Arm B.
#
# Implements corpus/PROTOCOL.md, frozen 2026-08-10 with SHA-256
# c55593ec5ca8eb27..., plus corpus/AMENDMENT-01.md, written 2026-08-21 before
# execution. Read both before reading this. Where they disagree with the code,
# they are right and the code is a bug.
#
# ONE RUN. §8: output is written once to corpus/data/confirmatory.json and not
# regenerated with altered options. The script refuses to overwrite an existing
# output file, because the way a held-out set gets spent twice is by someone
# re-running "just to check" after seeing a number they did not like.
#
#   Rscript corpus/05-confirmatory.R --selftest   # mechanics, on invented data
#   Rscript corpus/05-confirmatory.R              # the run
#
# --selftest touches nothing under corpus/data/. It exists so the machinery can
# be checked without opening the evidence: the alternative is a first execution
# that is simultaneously the debugging session, which is exactly what a
# pre-specified protocol is for.

suppressMessages({
  library(metafor)
  library(jsonlite)
})

# ---------------------------------------------------------------- constants
# Every number here is fixed by the protocol or the amendment. None is tuned.

TARGET_PREVALENCE <- 0.52    # AMENDMENT-01 §A: reweighted, published in the preprint
BOOT_R            <- 2000    # PROTOCOL §1
BOOT_LEVEL        <- 0.95    # PROTOCOL §1: percentile intervals
ERROR_LIMIT       <- 0.02    # PROTOCOL §3: above this a detector is reported failed
RATIO_MEASURES    <- c("RR", "OR", "HR", "IRR")  # PROTOCOL §2: ottawa's domain
LOG_SCALE         <- RATIO_MEASURES              # ratio measures are pooled in log

# PROTOCOL §7, fixed before execution so neither direction is the expected one.
POSITIVE <- list(sens = 0.60, spec = 0.60, lower_bound = 0.50, eligibility = 0.50)

# AMENDMENT-01 §A: stratum rates and their Wilson intervals, from the 120 coded
# pairs. S2 reweights with the extremes of these.
STRATA <- data.frame(
  stratum = c("alto", "medio", "bajo"),
  rate    = c(0.85, 0.65, 0.08),
  lo      = c(0.71, 0.50, 0.03),
  hi      = c(0.93, 0.78, 0.20),
  stringsAsFactors = FALSE
)

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || !nzchar(HERE)) HERE <- "corpus"
DATA <- file.path(HERE, "data")
OUT  <- file.path(DATA, "confirmatory.json")

# ------------------------------------------------------------ reconstruction

#' A pooled estimate, rebuilt as a one-study fixed-effect meta-analysis.
#'
#' The corpus carries what the abstract reported: a point estimate and its
#' interval. The detectors take `rma.uni` objects. Rebuilding one from the
#' interval is exact for the quantities they read -- the estimate and its
#' standard error -- and silent about everything it cannot know, which is the
#' per-study series. That silence is why sufficiency_changepoint and simulation
#' are excluded in §2 rather than approximated here.
as_pooled <- function(eff) {
  measure <- toupper(eff$measure)
  z <- qnorm(0.975)
  if (measure %in% LOG_SCALE) {
    if (!is.finite(eff$lo) || !is.finite(eff$hi) || eff$lo <= 0 || eff$est <= 0) {
      return(NULL)
    }
    yi  <- log(eff$est)
    sei <- (log(eff$hi) - log(eff$lo)) / (2 * z)
  } else {
    yi  <- eff$est
    sei <- (eff$hi - eff$lo) / (2 * z)
  }
  if (!is.finite(yi) || !is.finite(sei) || sei <= 0) return(NULL)
  suppressWarnings(
    tryCatch(rma(yi = yi, sei = sei, method = "FE", measure = measure),
             error = function(e) tryCatch(rma(yi = yi, sei = sei, method = "FE"),
                                          error = function(e2) NULL))
  )
}

# ------------------------------------------------------------------- verdict
#
# PROTOCOL §3. Three states, never merged:
#   answered       -- eligible, the detector returned out_of_date or current
#   not_applicable -- eligible, the detector was asked and declined
#   not_evaluable  -- the pair lacks the input the detector requires
# and errors, which are not verdicts: logged, excluded, counted.

ask <- function(fn, evaluable) {
  if (!isTRUE(evaluable)) return(list(state = "not_evaluable", fired = NA))
  res <- tryCatch(fn(), error = function(e) e)
  if (inherits(res, "error")) {
    return(list(state = "error", fired = NA, message = conditionMessage(res)))
  }
  v <- res$verdict
  if (is.null(v) || is.na(v)) return(list(state = "error", fired = NA,
                                          message = "detector returned no verdict"))
  if (grepl("not.?applicable", v)) return(list(state = "not_applicable", fired = NA))
  list(state = "answered", fired = identical(v, "out_of_date"))
}

# --------------------------------------------------------------------- data

read_pairs <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lapply(lines, fromJSON, simplifyVector = TRUE)
}

#' PROTOCOL §2 and §5, applied in that order and counted as we go.
eligible_set <- function(pairs) {
  drop <- c(protocol_or_missing = 0, same_abstract = 0,
            non_positive_interval = 0, not_comparable = 0)
  keep <- list()
  for (p in pairs) {
    if (is.null(p$conclusions_from) || is.null(p$conclusions_to) ||
        !nzchar(p$conclusions_from %||% "") || !nzchar(p$conclusions_to %||% "")) {
      drop["protocol_or_missing"] <- drop["protocol_or_missing"] + 1; next
    }
    if (isTRUE(p$same_abstract)) {
      drop["same_abstract"] <- drop["same_abstract"] + 1; next
    }
    if (!isTRUE(p$comparable_effects)) {
      drop["not_comparable"] <- drop["not_comparable"] + 1; next
    }
    from_date <- as.Date(p$from_date %||% NA)
    to_date   <- as.Date(p$to_date %||% NA)
    if (is.na(from_date) || is.na(to_date) || !(to_date > from_date)) {
      drop["non_positive_interval"] <- drop["non_positive_interval"] + 1; next
    }
    keep[[length(keep) + 1L]] <- p
  }
  list(pairs = keep, dropped = drop)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------------ outcome

#' AMENDMENT-01 §A. The threshold that reproduces the reweighted prevalence.
#'
#' Not the cut that maximises agreement with anything: this one never looks at
#' a detector, and the prevalence it reproduces was estimated and published
#' before the rule was chosen.
outcome_threshold <- function(scores, prevalence) {
  target <- floor(length(scores) * prevalence)   # nearest from below, per the amendment
  sorted <- sort(scores, decreasing = TRUE)
  sorted[target]
}

# ---------------------------------------------------------------- detectors

evaluate <- function(pairs, threshold) {
  rows <- vector("list", length(pairs))
  for (i in seq_along(pairs)) {
    p <- pairs[[i]]
    ef <- p$effect_from; et <- p$effect_to
    prev <- as_pooled(ef); new <- as_pooled(et)
    measure <- toupper(ef$measure)

    n_prev <- p$counts_from$n_participants
    n_to   <- p$counts_to$n_participants
    n_prev <- if (is.null(n_prev) || !length(n_prev)) NA_real_ else as.numeric(n_prev)
    n_to   <- if (is.null(n_to)   || !length(n_to))   NA_real_ else as.numeric(n_to)
    # AMENDMENT-01 §B: a total that does not grow does not license a guess at
    # how much new evidence arrived.
    barrowman_ok <- is.finite(n_prev) && is.finite(n_to) && n_to > n_prev
    n_new <- if (barrowman_ok) n_to - n_prev else NA_real_

    have_both <- !is.null(prev) && !is.null(new)

    rows[[i]] <- list(
      review_id = p$review_id,
      measure   = measure,
      event     = p$screen$score >= threshold,
      score     = p$screen$score,
      rcma      = ask(function() rcma(prev, new), have_both),
      ottawa    = ask(function() ottawa(prev, new),
                      have_both && measure %in% RATIO_MEASURES),
      barrowman = ask(function() barrowman(prev, n_prev = n_prev, n_new = n_new),
                      have_both && barrowman_ok)
    )
  }
  rows
}

# ------------------------------------------------------------------ metrics
#
# PROTOCOL §1 and §6. Eligibility and conditional discrimination, reported
# separately and never combined into one score.

rates_for <- function(rows, detector) {
  st    <- vapply(rows, function(r) r[[detector]]$state, character(1))
  fired <- vapply(rows, function(r) {
    f <- r[[detector]]$fired; if (is.null(f) || is.na(f)) NA else f
  }, logical(1))
  event <- vapply(rows, function(r) isTRUE(r$event), logical(1))

  n_total   <- length(rows)
  evaluable <- st != "not_evaluable" & st != "error"
  answered  <- st == "answered"

  ans_event    <- answered & event
  ans_no_event <- answered & !event

  list(
    n_pairs        = n_total,
    n_evaluable    = sum(evaluable),
    n_answered     = sum(answered),
    n_not_appl     = sum(st == "not_applicable"),
    n_not_eval     = sum(st == "not_evaluable"),
    n_error        = sum(st == "error"),
    eligibility    = if (n_total)          sum(evaluable) / n_total          else NA,
    na_rate        = if (sum(evaluable))   sum(st == "not_applicable") / sum(evaluable) else NA,
    error_rate     = if (n_total)          sum(st == "error") / n_total      else NA,
    sensitivity    = if (sum(ans_event))    mean(fired[ans_event])           else NA,
    specificity    = if (sum(ans_no_event)) mean(!fired[ans_no_event])       else NA,
    false_alarm    = if (sum(ans_no_event)) mean(fired[ans_no_event])        else NA,
    n_events       = sum(ans_event),
    n_non_events   = sum(ans_no_event)
  )
}

#' PROTOCOL §1: nonparametric bootstrap resampling WHOLE REVIEWS, 2,000
#' replicates, percentile intervals. Pairs are not resampled -- reviews
#' contribute more than one pair, and treating them as independent would
#' narrow every interval here by pretending 560 pairs are 560 studies.
boot_ci <- function(rows, detector, statistic, R = BOOT_R, seed = 20260821) {
  set.seed(seed)
  by_review <- split(seq_along(rows), vapply(rows, function(r) r$review_id, character(1)))
  reviews <- names(by_review)
  reps <- numeric(R)
  for (b in seq_len(R)) {
    pick <- sample(reviews, length(reviews), replace = TRUE)
    idx  <- unlist(by_review[pick], use.names = FALSE)
    reps[b] <- rates_for(rows[idx], detector)[[statistic]]
  }
  reps <- reps[is.finite(reps)]
  if (!length(reps)) return(c(lo = NA, hi = NA))
  a <- (1 - BOOT_LEVEL) / 2
  stats::setNames(unname(stats::quantile(reps, c(a, 1 - a), na.rm = TRUE)), c("lo", "hi"))
}

#' PROTOCOL §7, applied to the numbers rather than to an impression of them.
classify <- function(m, ci) {
  if (!is.finite(m$sensitivity) || !is.finite(m$specificity)) return("inconclusive")
  clears <- m$sensitivity >= POSITIVE$sens &&
            m$specificity >= POSITIVE$spec &&
            isTRUE(ci$sensitivity["lo"] > POSITIVE$lower_bound) &&
            isTRUE(ci$specificity["lo"] > POSITIVE$lower_bound) &&
            m$eligibility >= POSITIVE$eligibility
  if (clears) return("positive")
  wide <- isTRUE((ci$sensitivity["hi"] - ci$sensitivity["lo"]) > 0.40) ||
          isTRUE((ci$specificity["hi"] - ci$specificity["lo"]) > 0.40)
  if (wide) "inconclusive" else "negative"
}

summarise_detector <- function(rows, detector, contaminated) {
  m  <- rates_for(rows, detector)
  ci <- list(
    sensitivity = boot_ci(rows, detector, "sensitivity"),
    specificity = boot_ci(rows, detector, "specificity"),
    eligibility = boot_ci(rows, detector, "eligibility")
  )
  c(m, list(
    ci = ci,
    result = classify(m, ci),
    # PROTOCOL §3: above 2% the run is reported as FAILED for that detector
    # rather than patched.
    failed_on_errors = isTRUE(m$error_rate > ERROR_LIMIT),
    # PROTOCOL §6, and CONTAMINATED_PAIRS in the package: ottawa shares logic
    # with the conclusion-change outcome. Flagged, never dropped.
    contaminated = detector %in% contaminated
  ))
}

# ----------------------------------------------------------------- the runs

run_all <- function(pairs, threshold, contaminated) {
  rows <- evaluate(pairs, threshold)
  detectors <- c("rcma", "ottawa", "barrowman")
  list(
    n_pairs   = length(rows),
    n_reviews = length(unique(vapply(rows, function(r) r$review_id, character(1)))),
    n_events  = sum(vapply(rows, function(r) isTRUE(r$event), logical(1))),
    threshold = threshold,
    detectors = stats::setNames(
      lapply(detectors, function(d) summarise_detector(rows, d, contaminated)),
      detectors
    ),
    rows = rows
  )
}

#' S2, PROTOCOL §4: the two most extreme reweightings consistent with the
#' stratum intervals. Not a re-run of the analysis with a nicer outcome --
#' a bound on how far the prevalence could move and drag the rates with it.
s2_prevalences <- function(n) {
  pop <- c(n %/% 3, n %/% 3, n - 2 * (n %/% 3))
  c(low  = sum(STRATA$lo * pop) / n,
    high = sum(STRATA$hi * pop) / n)
}

# ------------------------------------------------------------------ selftest
#
# Invented pairs with a known answer, so the mechanics can be wrong here
# instead of on the evidence. Touches nothing under corpus/data/.

selftest <- function() {
  cat("selftest: invented data, corpus/data/ untouched\n\n")
  mk <- function(id, est_from, est_to, event, measure = "RR",
                 n_from = 1000, n_to = 1400) {
    list(review_id = id, same_abstract = FALSE, comparable_effects = TRUE,
         from_date = "2015-01-01", to_date = "2018-01-01",
         conclusions_from = "x", conclusions_to = "y",
         effect_from = list(measure = measure, est = est_from,
                            lo = est_from * 0.8, hi = est_from * 1.25),
         effect_to   = list(measure = measure, est = est_to,
                            lo = est_to * 0.8, hi = est_to * 1.25),
         counts_from = list(n_participants = n_from),
         counts_to   = list(n_participants = n_to),
         screen = list(score = if (event) 5 else 0.1))
  }
  pairs <- list(
    mk("A", 0.90, 0.40, TRUE),          # big move, event
    mk("A", 0.90, 0.42, TRUE),          # same review: clustering must see this
    mk("B", 0.90, 0.89, FALSE),         # no move, no event
    mk("C", 0.90, 0.91, FALSE),
    mk("D", 0.50, 0.20, TRUE, measure = "MD", n_from = 900, n_to = 800),
    mk("E", 0.90, 0.88, FALSE, measure = "SMD")
  )
  stopifnot(length(eligible_set(pairs)$pairs) == length(pairs))

  thr <- outcome_threshold(vapply(pairs, function(p) p$screen$score, numeric(1)), 0.5)
  cat("threshold at 50% prevalence:", thr, "\n")

  res <- run_all(pairs, thr, contaminated = CONTAMINATED_PAIRS$method)
  cat("pairs:", res$n_pairs, " reviews:", res$n_reviews, " events:", res$n_events, "\n\n")

  for (d in names(res$detectors)) {
    m <- res$detectors[[d]]
    cat(sprintf("%-10s eligibility %.2f  answered %d  n/a %d  not-eval %d  err %d%s\n",
                d, m$eligibility, m$n_answered, m$n_not_appl, m$n_not_eval, m$n_error,
                if (m$contaminated) "  [CONTAMINATED]" else ""))
  }

  # The checks that would have caught a wrong wiring.
  ott <- res$detectors$ottawa
  stopifnot(ott$n_not_eval == 2)                       # MD and SMD are outside ottawa's domain
  stopifnot(ott$contaminated)                          # ottawa x conclusion is circular
  bar <- res$detectors$barrowman
  stopifnot(bar$n_not_eval == 1)                       # D's total shrinks: no new evidence
  stopifnot(res$detectors$rcma$n_not_eval == 0)
  stopifnot(res$n_reviews == 5 && res$n_pairs == 6)    # A contributes two pairs, one review

  ci <- boot_ci(res$rows, "rcma", "specificity", R = 50)
  stopifnot(is.finite(ci["lo"]), ci["lo"] <= ci["hi"])
  cat("\nbootstrap clusters by review: OK\n")
  cat("selftest passed\n")
  invisible(TRUE)
}

# ---------------------------------------------------------------------- main

main <- function() {
  args <- commandArgs(TRUE)
  if ("--selftest" %in% args) { selftest(); return(invisible(TRUE)) }

  # §8. One run.
  if (file.exists(OUT)) {
    stop(sprintf(paste("%s already exists. PROTOCOL.md §8 allows one run and this",
                       "would be the second. If the first run genuinely failed,",
                       "move the file aside by hand and say so in the manuscript."),
                 OUT), call. = FALSE)
  }

  pairs_all <- read_pairs(file.path(DATA, "screened.jsonl"))
  es <- eligible_set(pairs_all)
  pairs <- es$pairs
  scores <- vapply(pairs, function(p) p$screen$score, numeric(1))

  threshold <- outcome_threshold(scores, TARGET_PREVALENCE)
  contaminated <- CONTAMINATED_PAIRS$method

  primary <- run_all(pairs, threshold, contaminated)

  # S1: the 120 pairs carrying a blind code. PROTOCOL §4.
  key   <- fromJSON(file.path(DATA, "coding", "key.json"), simplifyDataFrame = TRUE)
  codes <- fromJSON(file.path(DATA, "coding", "codes-claude.json"), simplifyVector = TRUE)
  key$code <- unname(codes[as.character(key$n)])
  coded_id <- paste(key$review_id, key$from_version, key$to_version)
  pair_id  <- vapply(pairs, function(p) paste(p$review_id, p$from_version, p$to_version),
                     character(1))
  in_s1 <- pair_id %in% coded_id
  s1_pairs <- pairs[in_s1]
  s1_event <- key$code[match(pair_id[in_s1], coded_id)] == "major"
  s1_rows  <- evaluate(s1_pairs, threshold)
  for (i in seq_along(s1_rows)) s1_rows[[i]]$event <- s1_event[i]
  s1 <- list(
    n_pairs = length(s1_rows),
    n_events = sum(s1_event),
    detectors = stats::setNames(
      lapply(c("rcma", "ottawa", "barrowman"),
             function(d) summarise_detector(s1_rows, d, contaminated)),
      c("rcma", "ottawa", "barrowman"))
  )

  # S2: the two extreme reweightings. PROTOCOL §4.
  s2p <- s2_prevalences(length(pairs))
  s2 <- lapply(s2p, function(pv) {
    th <- outcome_threshold(scores, pv)
    r  <- run_all(pairs, th, contaminated)
    r$rows <- NULL
    c(list(prevalence = pv), r)
  })

  primary$rows <- NULL
  out <- list(
    protocol = list(
      file = "corpus/PROTOCOL.md",
      sha256 = "c55593ec5ca8eb2777d02a92447f2047fd085c5225d0578247741abadb7abed7",
      amendment = "corpus/AMENDMENT-01.md"
    ),
    registration = list(
      commit = tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) NA),
      date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      human_coding_available = FALSE,
      outcome_rule = "screen.score >= threshold reproducing reweighted prevalence 0.52"
    ),
    eligible = list(n = length(pairs), dropped = as.list(es$dropped),
                    threshold = threshold,
                    prevalence_realised = primary$n_events / length(pairs)),
    primary = primary,
    S1 = s1,
    S2 = s2
  )
  write(toJSON(out, auto_unbox = TRUE, digits = 6, pretty = TRUE), OUT)

  cat("\n=== CONFIRMATORY RUN ===\n")
  cat(sprintf("pairs %d in %d reviews | threshold %.3f | events %d (%.1f%%)\n",
              primary$n_pairs, primary$n_reviews, threshold,
              primary$n_events, 100 * primary$n_events / primary$n_pairs))
  for (d in names(primary$detectors)) {
    m <- primary$detectors[[d]]
    cat(sprintf("\n%-10s %s%s\n", d, toupper(m$result),
                if (m$contaminated) "  [CONTAMINATED: shares logic with the outcome]" else ""))
    cat(sprintf("  eligibility %.2f [%.2f-%.2f]  n/a %.2f  errors %.3f%s\n",
                m$eligibility, m$ci$eligibility["lo"], m$ci$eligibility["hi"],
                m$na_rate, m$error_rate,
                if (m$failed_on_errors) "  <- FAILED, over 2%" else ""))
    cat(sprintf("  sensitivity %.2f [%.2f-%.2f] on %d events\n",
                m$sensitivity, m$ci$sensitivity["lo"], m$ci$sensitivity["hi"], m$n_events))
    cat(sprintf("  specificity %.2f [%.2f-%.2f] on %d non-events\n",
                m$specificity, m$ci$specificity["lo"], m$ci$specificity["hi"], m$n_non_events))
  }
  cat("\nwritten:", OUT, "\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) main()
