# Searches CRAN for existing implementations of the five updating detectors,
# and writes a dated snapshot of what it found.
#
# The README and paper.md both claim that none of the five methods had a
# reusable implementation. That is a claim of ABSENCE across a whole
# repository, which is the kind that needs its evidence attached: without the
# terms, the fields and the date, a reader cannot tell a thorough search from
# a cursory one, and neither can the author a year later.
#
# Run it with:
#
#   R CMD INSTALL .
#   Rscript inst/cran-search/cran-search.R
#
# or, while developing, without installing:
#
#   Rscript -e 'pkgload::load_all("."); source("inst/cran-search/cran-search.R")'
#
# It needs network access: tools::CRAN_package_db() downloads the full package
# database (a few MB). Nothing else in the package needs the network, and no
# test runs this -- the tests read the snapshot this script writes.
#
# WHAT THIS CAN AND CANNOT ESTABLISH. It searches the Package name, Title and
# Description of every package on CRAN. It does NOT search source code, help
# pages, or vignettes: doing so would mean downloading every package. So a
# package implementing one of these criteria without naming it in its metadata
# would not appear here.
#
# That is why the claim this supports is "we did not find one", not "none
# exists". The distinction is not pedantic -- it is the difference between a
# statement the evidence carries and one it does not.

SNAPSHOT <- "cran-search-snapshot.csv"

# Terms in three groups, because they answer different questions.
#
#   method     -- the detectors themselves, by the names their sources use
#   problem    -- the task, however a package might phrase it
#   component  -- the building blocks, which we expect to FIND and whose
#                 presence is part of the claim: the pieces exist, the
#                 assembled detectors do not
TERMS <- list(
  list(group = "method",    label = "Ottawa method",             pattern = "ottawa"),
  list(group = "method",    label = "Barrowman",                 pattern = "barrowman"),
  list(group = "method",    label = "recursive cumulative MA",   pattern = "recursive cumulative"),
  list(group = "method",    label = "Shojania",                  pattern = "shojania"),
  list(group = "method",    label = "Pattanittum",               pattern = "pattanittum"),
  list(group = "problem",   label = "out-of-date review",        pattern = "out[- ]of[- ]date"),
  list(group = "problem",   label = "updating a review",         pattern = "updat\\w* (a |the )?(systematic )?(review|meta-?analys)"),
  list(group = "problem",   label = "currency / staleness",      pattern = "\\bstaleness\\b|currency of (the )?evidence"),
  list(group = "problem",   label = "when to update",            pattern = "when to update"),
  list(group = "component", label = "fail-safe N",               pattern = "fail[- ]safe"),
  list(group = "component", label = "cumulative meta-analysis",  pattern = "cumulative meta-?analys"),
  list(group = "component", label = "trial sequential analysis", pattern = "trial sequential")
)

# Every hit is adjudicated by hand and the verdict recorded, because a count
# alone cannot tell a real implementation from a package that happens to be
# written at the University of Ottawa. Anything not listed here is reported as
# UNADJUDICATED, so a new hit on a later run is visible rather than absorbed
# into a number.
VERDICTS <- c(
  SAiVE        = "false positive: the University of Ottawa, not the Ottawa method",
  themis        = "false positive: class imbalance, unrelated",
  updateme      = "false positive: warns about out-of-date R PACKAGES, not reviews",
  VegSpecIndex  = "false positive: spectral indices for vegetation",
  CRTSize       = "false positive: 'updated techniques' for cluster-trial sample size",
  metagear      = "adjacent: screens abstracts for systematic reviews; does not assess whether a completed review is out of date",
  fsn           = "component: Rosenthal's fail-safe N, one half of sufficiency_changepoint()",
  meta          = "component: general meta-analysis, computes cumulative meta-analyses",
  metafor       = "component: cumulative meta-analyses and fail-safe N; this package depends on it",
  RTSA          = "component: trial sequential analysis, a different sequential question"
)

if (!requireNamespace("tools", quietly = TRUE)) {
  stop("this script needs the tools package", call. = FALSE)
}

message("downloading the CRAN package database ...")
db <- tools::CRAN_package_db()
haystack <- paste(db$Package, db$Title, db$Description)
today <- as.character(Sys.Date())

rows <- list()
for (t in TERMS) {
  hit <- grep(t$pattern, haystack, ignore.case = TRUE)
  if (!length(hit)) {
    rows[[length(rows) + 1L]] <- data.frame(
      date = today, n_cran = nrow(db), group = t$group, term = t$label,
      pattern = t$pattern, package = NA_character_,
      verdict = "no hit", stringsAsFactors = FALSE)
    next
  }
  for (i in hit) {
    p <- db$Package[i]
    rows[[length(rows) + 1L]] <- data.frame(
      date = today, n_cran = nrow(db), group = t$group, term = t$label,
      pattern = t$pattern, package = p,
      verdict = if (p %in% names(VERDICTS)) VERDICTS[[p]] else "UNADJUDICATED",
      stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)

cat("\n=== CRAN search,", today, "===\n")
cat("packages searched:", nrow(db), "(fields: Package, Title, Description)\n\n")

for (g in c("method", "problem", "component")) {
  cat(sprintf("-- %s --\n", g))
  sub <- res[res$group == g, ]
  for (term in unique(sub$term)) {
    s <- sub[sub$term == term, ]
    if (all(s$verdict == "no hit")) {
      cat(sprintf("  %-28s no hit\n", term))
    } else {
      cat(sprintf("  %-28s %d hit(s)\n", term, nrow(s)))
      for (i in seq_len(nrow(s))) {
        cat(sprintf("      %-14s %s\n", s$package[i], s$verdict[i]))
      }
    }
  }
  cat("\n")
}

unadj <- res$package[res$verdict == "UNADJUDICATED"]
if (length(unadj)) {
  cat("!! UNADJUDICATED HITS:", paste(unique(unadj), collapse = ", "), "\n")
  cat("   A new package matched. Read it and add a verdict to VERDICTS above,\n")
  cat("   or the claim in README.md and paper.md is no longer supported.\n\n")
} else {
  cat("Every hit is adjudicated.\n\n")
}

# The claim the search supports, stated as the evidence allows.
implementations <- res[res$group == "method" &
                         res$verdict != "no hit" &
                         !grepl("^false positive", res$verdict), ]
cat("Implementations of any of the five detectors found:", nrow(implementations), "\n")
cat("Components found:",
    paste(sort(unique(res$package[res$group == "component" &
                                    !is.na(res$package)])), collapse = ", "), "\n")

# Written next to this script when run from the source tree, which is where it
# is meant to be re-run and committed; otherwise into the working directory.
# `%||%` is not used: the package declares R (>= 4.2) and it arrived in 4.4.
dest <- if (dir.exists("inst/cran-search")) {
  file.path("inst/cran-search", SNAPSHOT)
} else {
  SNAPSHOT
}
utils::write.csv(res, dest, row.names = FALSE)
cat("\nsnapshot written to", dest, "\n")

invisible(res)
