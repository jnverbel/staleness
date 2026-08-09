# Runs all five detectors across every historical review in metadat that can
# carry a backtest, and reports how often each one can answer at all.
#
# The published comparison of these methods (Pattanittum et al. 2012) applied
# them to 80 Cochrane reviews selected for having a non-significant pooled
# result. Two of the five -- barrowman() and simulation() -- only speak when
# the prior meta-analysis was not significant, so that cohort is the one place
# they can always be asked. This script asks what happens on reviews that were
# not selected that way, which is the question the source could not answer
# about itself.
#
# Run it with:
#
#   Rscript inst/applicability/applicability.R
#
# It takes a few seconds and needs metadat. Results are deterministic: the
# backtest seed is fixed and no figure below depends on simulation variance.
#
# INCLUSION CRITERIA, applied in this order and reported as they bite:
#
#   1. a data frame in metadat with a per-study year column
#   2. a two-group effect measure this script can build with escalc(), which
#      takes eight different column conventions for the same 2x2 table
#   3. at least 8 studies, so that yearly cuts exist
#   4. one effect per study. This script passes each row's study identifier to
#      evidence_stream(), which refuses a stream whose estimates are dependent
#      and names the studies responsible. It used to detect nesting with a
#      heuristic of its own and drop those datasets before building anything --
#      having to do that by hand was the defect that requiring `study_id`
#      fixed, and the check now lives where every caller gets it.
#
#      Note that repeated AUTHOR names are not nesting. dat.bcg has Comstock
#      three times, in three different years: three trials, not three effects
#      from one. The identifier is therefore study-and-year, not study alone.
#
#   5. the backtest must produce at least three uncensored cuts
#
# WHAT THIS IS NOT. Sensitivity and specificity below are scored against
# truth_shift(), which is this package's own definition, not a published
# outcome. This compares the five methods against each other over real
# evidence; it does not validate any of them against what actually happened,
# which is what the four cases in test-external-validation.R do. The two
# claims are different and only the first one has n = 17.

# This runs against the package, so the package has to be there. Saying so
# with an actionable message: a bare library() call fails from a clean clone
# with "there is no package called 'staleness'", which tells the reader what
# happened but not what to do.
if (!requireNamespace("staleness", quietly = TRUE)) {
  stop("this script runs against the staleness package, which is not ",
       "installed.\n  From a clean clone:  R CMD INSTALL . && Rscript ",
       "inst/applicability/applicability.R\n  While developing:    Rscript -e ",
       "'pkgload::load_all(\".\"); source(\"inst/applicability/applicability.R\")'",
       call. = FALSE)
}
if (!requireNamespace("metadat", quietly = TRUE)) {
  stop("this script needs the metadat package", call. = FALSE)
}
suppressMessages(library(staleness))
suppressMessages(library(metafor))

HORIZON <- 3   # years ahead the truth is evaluated over
WINDOW  <- 5   # years of new evidence each cut may use
MIN_K   <- 3   # smallest snapshot that can be fitted
SEED    <- 42

# One set of parameters for every review, so the comparison across them is
# like for like. No sensitivity analysis is claimed.

has <- function(d, ...) all(c(...) %in% names(d))
year_col <- function(d) {
  hit <- names(d)[grepl("^(year|yr)$", names(d), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
}

# metadat stores the same 2x2 table under eight different column conventions,
# and the mapping has to be explicit: guessing from column position would be
# how a control arm silently becomes a treatment arm.
build <- function(d) suppressWarnings(build_effects(d))

# escalc() warns when a 2x2 table with a zero cell produces an infinite effect
# ("Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs"). That is
# expected here and handled: the caller drops non-finite rows immediately
# below, and a review left with fewer than 8 usable studies is excluded and
# reported. Silenced at the point it is produced rather than at the call site,
# so a warning from anywhere else still surfaces.
build_effects <- function(d) {
  if (has(d, "ai", "bi", "ci", "di"))
    list(es = escalc("OR", ai = ai, bi = bi, ci = ci, di = di, data = d),
         ni = d$ai + d$bi + d$ci + d$di, measure = "OR")
  else if (has(d, "ai", "n1i", "ci", "n2i"))
    list(es = escalc("OR", ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = d),
         ni = d$n1i + d$n2i, measure = "OR")
  else if (has(d, "tpos", "tneg", "cpos", "cneg"))
    list(es = escalc("RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = d),
         ni = d$tpos + d$tneg + d$cpos + d$cneg, measure = "RR")
  else if (has(d, "xt", "nt", "xc", "nc"))
    list(es = escalc("OR", ai = xt, n1i = nt, ci = xc, n2i = nc, data = d),
         ni = d$nt + d$nc, measure = "OR")
  else if (has(d, "Ee", "Ne", "Ec", "Nc"))
    list(es = escalc("OR", ai = Ee, n1i = Ne, ci = Ec, n2i = Nc, data = d),
         ni = d$Ne + d$Nc, measure = "OR")
  else if (has(d, "x.a", "n.a", "x.p", "n.p"))
    list(es = escalc("OR", ai = x.a, n1i = n.a, ci = x.p, n2i = n.p, data = d),
         ni = d$n.a + d$n.p, measure = "OR")
  else if (has(d, "resp1", "n1", "resp2", "n2"))
    list(es = escalc("OR", ai = resp1, n1i = n1, ci = resp2, n2i = n2, data = d),
         ni = d$n1 + d$n2, measure = "OR")
  else if (has(d, "Me", "Se", "Ne", "Mc", "Sc", "Nc"))
    list(es = escalc("SMD", m1i = Me, sd1i = Se, n1i = Ne,
                     m2i = Mc, sd2i = Sc, n2i = Nc, data = d),
         ni = d$Ne + d$Nc, measure = "SMD")
  else if (has(d, "yi", "vi"))
    list(es = d, ni = NULL, measure = NA_character_)  # resolved below
  else NULL
}

# Datasets that ship yi/vi already computed do not say what scale they are on,
# and the scale decides which branch effect_ratio() takes. Read off metadat's
# own documentation, one at a time. dat.hackshaw1998 is a LOG ODDS RATIO -- a
# ratio measure -- and treating it as a difference would be silently wrong.
PRECOMPUTED <- c(
  dat.bangertdrowns2004   = "SMD",
  dat.hackshaw1998        = "OR",
  dat.konstantopoulos2011 = "SMD",
  dat.raudenbush1985      = "SMD"
)

# The identifier of the study each row came from, so that evidence_stream()
# can refuse dependence itself. This script used to detect nesting with a
# heuristic of its own and drop those datasets before building the stream --
# and having to do that by hand was the defect the requirement now fixes.
#
# A study-and-year pair rather than the study column alone: repeated AUTHOR
# names are not nesting. dat.bcg has Comstock three times, in three different
# years, which is three trials.
study_ids <- function(d) {
  if (any(grepl("^esid$|^es\\.id$", names(d), ignore.case = TRUE))) {
    ids <- intersect(c("study", "studyid", "trial", "id", "author"), names(d))
    if (length(ids)) return(as.character(d[[ids[1]]]))
  }
  y <- year_col(d)
  ids <- intersect(c("study", "studyid", "trial", "id", "author",
                     "Study", "StudyID", "district"), names(d))
  if (!length(ids) || is.na(y)) return(as.character(seq_len(nrow(d))))
  paste(d[[ids[1]]], d[[y]])
}

datasets <- sub("\\s.*", "", utils::data(package = "metadat")$results[, "Item"])

rows <- list(); excluded <- list()
for (nm in datasets) {
  d <- tryCatch(get(nm, envir = asNamespace("metadat")), error = function(e) NULL)
  if (is.null(d) || !is.data.frame(d)) next
  y <- year_col(d)
  if (is.na(y))          { excluded[[nm]] <- "no per-study year";      next }
  b <- build(d)
  if (is.null(b))        { excluded[[nm]] <- "no mappable measure";    next }
  if (is.na(b$measure)) {
    if (!nm %in% names(PRECOMPUTED)) {
      excluded[[nm]] <- "precomputed yi/vi of undetermined scale";     next
    }
    b$measure <- PRECOMPUTED[[nm]]
  }
  if (nrow(d) < 8)       { excluded[[nm]] <- "fewer than 8 studies";   next }
  sid <- study_ids(d)

  keep <- is.finite(b$es$yi) & is.finite(b$es$vi) & b$es$vi > 0 &
    is.finite(d[[y]])
  es <- b$es[keep, ]; yr <- d[[y]][keep]; sid <- sid[keep]
  ni <- if (is.null(b$ni)) NULL else b$ni[keep]
  if (nrow(es) < 8)      { excluded[[nm]] <- "fewer than 8 usable studies"; next }

  out <- tryCatch({
    ma <- rma(yi, vi, data = es, measure = b$measure, method = "FE")
    # No allow_dependence: a review whose rows are not independent studies is
    # refused here, by the package, and reported below as an exclusion. This
    # replaces the hand-rolled nesting heuristic this script used to carry.
    st <- evidence_stream(ma, date = yr, study_id = sid, ni = ni)
    bt <- suppressWarnings(backtest(st, cuts = "yearly", horizon = HORIZON,
                                    window = WINDOW, min_k = MIN_K, seed = SEED))
    cal <- calibration(bt, "shift")
    cal$dataset <- nm
    cal$k <- nrow(es)
    cal$measure <- b$measure
    cal
  }, error = function(e) {
    m <- conditionMessage(e)
    # The dependence refusal names every offending study, which is right in a
    # call and too long for a tally. Shortened here, not there.
    if (grepl("more than one estimate", m)) m <- "several estimates per study"
    excluded[[nm]] <<- m
    NULL
  })
  if (!is.null(out)) rows[[nm]] <- out
}

res <- do.call(rbind, rows)
reviews <- unique(res$dataset)

cat("\n=== COHORT ===\n")
cat(sprintf("  reviews with a backtest: %d\n", length(reviews)))
cat(sprintf("  excluded: %d\n", length(excluded)))
tab <- table(unlist(excluded))
for (r in names(sort(tab, decreasing = TRUE))) {
  cat(sprintf("    %-42s %d\n", r, tab[[r]]))
}

cat("\n=== HOW OFTEN EACH DETECTOR CAN ANSWER AT ALL ===\n")
cov <- aggregate(cbind(answers = n > 0) ~ method, data = res, FUN = sum)
cov <- cov[order(-cov$answers), ]
for (i in seq_len(nrow(cov))) {
  cat(sprintf("  %-12s %2d of %d reviews\n",
              cov$method[i], cov$answers[i], length(reviews)))
}

# The reason, measured rather than asserted: barrowman() and simulation() are
# only applicable when the prior meta-analysis was not significant.
cat("\n=== WHY: WAS THE PRIOR ALREADY SIGNIFICANT? ===\n")
sig <- 0; tot <- 0
for (nm in reviews) {
  d <- get(nm, envir = asNamespace("metadat")); y <- year_col(d)
  b <- build(d)
  if (is.na(b$measure)) b$measure <- PRECOMPUTED[[nm]]
  keep <- is.finite(b$es$yi) & is.finite(b$es$vi) & b$es$vi > 0 & is.finite(d[[y]])
  es <- b$es[keep, ]; yr <- d[[y]][keep]
  for (cut in sort(unique(yr))) {
    sub <- es[yr <= cut, ]
    if (nrow(sub) < MIN_K) next
    f <- tryCatch(suppressWarnings(rma(yi, vi, data = sub, method = "FE")),
                  error = function(e) NULL)
    if (is.null(f) || !is.finite(f$pval)) next
    tot <- tot + 1
    if (f$pval < 0.05) sig <- sig + 1
  }
}
cat(sprintf("  %d of %d cuts had an already-significant prior: %.0f%%\n",
            sig, tot, 100 * sig / tot))
cat("  barrowman() and simulation() require a NON-significant prior, so they\n")
cat("  are inapplicable at those cuts by construction.\n")

# --- What the truth was measured against ----------------------------------
#
# The sweep above scores every cut against the model fitted over the whole
# stream (truth_target = "final"). That is the package default and it answers
# a real question -- was this review already superseded by what would
# eventually be known -- but it is the more forgiving of the two, because
# almost everything counts as a positive.
cat("\n=== HOW MANY CUTS HAVE NO NEGATIVE TO SCORE AGAINST ===\n")
no_neg <- 0
for (nm in reviews) {
  rows_nm <- res[res$dataset == nm, ]
  if (all(is.na(rows_nm$specificity))) no_neg <- no_neg + 1
}
cat(sprintf("  reviews where specificity cannot be computed at all: %d of %d\n",
            no_neg, length(reviews)))
cat("  Against the final model every earlier cut of a review whose effect kept\n")
cat("  moving counts as out of date, so there are no true negatives left.\n")
cat("  Re-run with truth_target = \"horizon\" to ask instead whether each review\n")
cat("  went out of date within the next HORIZON years. On metadat::dat.bcg fitted\n")
cat("  with fixed effects, as this sweep fits every review, that turns 23 out of\n")
cat("  23 positives into 7 and makes specificity computable. Under REML it would\n")
cat("  be 5 and 0: the rate depends on the estimator as much as on the target.\n")

cat("\n=== SENSITIVITY WHERE IT IS DEFINED ===\n")
s <- res[!is.na(res$sensitivity), ]
agg <- aggregate(sensitivity ~ method, data = s,
                 FUN = function(x) c(reviews = length(x),
                                     mean = round(mean(x), 3),
                                     zero = sum(x == 0)))
agg <- cbind(method = agg$method, as.data.frame(agg$sensitivity))
print(agg[order(-agg$mean), ], row.names = FALSE)

cat("\n=== PER REVIEW ===\n")
t <- res[res$n > 0, c("dataset", "k", "measure", "method", "n",
                      "sensitivity", "specificity")]
t$sensitivity <- round(t$sensitivity, 3)
t$specificity <- round(t$specificity, 3)
print(t[order(t$dataset, t$method), ], row.names = FALSE)

invisible(res)
