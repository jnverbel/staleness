# README.md and paper.md rest on a claim of absence: that no CRAN package
# implements any of the five detectors. inst/cran-search/ holds the evidence
# as a dated snapshot. This holds the prose to that snapshot.
#
# It does NOT touch the network. Re-running the search is a deliberate act
# (Rscript inst/cran-search/cran-search.R) that rewrites the snapshot; the
# tests only check that what the documents say matches what was found.

test_that("the CRAN-search claim matches the snapshot it rests on", {
  # Located with system.file(), which resolves under load_all(), under R CMD
  # check and from an installed library alike -- unlike a relative path, which
  # would silently skip this in the very distribution it exists to check.
  snap <- system.file("cran-search/cran-search-snapshot.csv",
                      package = "staleness")
  skip_if_not(nzchar(snap), "cran-search snapshot not found in the package")

  s <- utils::read.csv(snap, stringsAsFactors = FALSE)

  # 1. Every hit is adjudicated. An unadjudicated one means a package matched
  #    that nobody has read, and the claim is unsupported until someone does.
  expect_equal(sum(s$verdict == "UNADJUDICATED"), 0)

  # 2. The claim itself: no hit on a method name survives adjudication.
  methods <- s[s$group == "method" & s$verdict != "no hit", ]
  expect_true(all(grepl("^false positive", methods$verdict)))

  # 3. The snapshot carries a single date and package count, so the documents
  #    can quote them.
  expect_length(unique(s$date), 1)
  expect_length(unique(s$n_cran), 1)
  n_cran <- unique(s$n_cran)
  searched_on <- unique(s$date)

  # 4. And README.md and paper.md must quote those, not older ones. Both are
  #    hand-written and regenerate nothing, which is how the main README came
  #    to promise an ottawa signal of 1.21 while the package returned 0.86.
  pretty <- formatC(n_cran, format = "d", big.mark = ",")
  for (f in c("../../README.md", "../../paper.md")) {
    if (!file.exists(f)) next   # absent from an installed package; see below
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_true(grepl(pretty, txt, fixed = TRUE),
                info = paste(f, "should quote", pretty, "packages"))
    expect_true(grepl(searched_on, txt, fixed = TRUE),
                info = paste(f, "should quote the search date", searched_on))
    # The claim must be the one the evidence carries.
    expect_true(grepl("could not find", txt, fixed = TRUE),
                info = paste(f, "should say 'we could not find one', not",
                             "'none exists': the search reads metadata, not code"))
  }

  # 5. The demonstrated limit: metafor exports cumul() and fsn() and matches
  #    neither phrase in its metadata. If that ever stops being true the
  #    caveat in both documents needs rewording, so it is pinned here rather
  #    than left as a sentence nobody rechecks.
  skip_if_not_installed("metafor")
  expect_true(exists("cumul", where = asNamespace("metafor")))
  expect_true(exists("fsn", where = asNamespace("metafor")))
})
