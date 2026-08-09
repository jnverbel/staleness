# inst/applicability/README.md quotes figures produced by applicability.R.
# Like the top-level README, it is written by hand and regenerates nothing, so
# it is the kind of document that drifts without anything noticing -- which is
# exactly how the main README came to promise an ottawa signal of 1.21 while
# the package returned 0.86.
#
# This runs the script and holds its headline figures to what the README says.

test_that("the applicability sweep still produces the figures its README quotes", {
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  # Skipped on the file's absence rather than skip_on_cran(): the installed
  # package has no inst/ above the test directory under R CMD check, but
  # skip_on_cran() would also silence this in local development, where the
  # script is actually edited.
  script <- "../../inst/applicability/applicability.R"
  readme <- "../../inst/applicability/README.md"
  skip_if_not(file.exists(script) && file.exists(readme),
              "applicability files not reachable from the test dir")

  out <- utils::capture.output(source(script, local = new.env()))
  txt <- paste(out, collapse = "\n")

  # 1. The cohort size.
  n_reviews <- as.integer(sub(".*reviews with a backtest:\\s*(\\d+).*", "\\1",
                              grep("reviews with a backtest", out, value = TRUE)))
  expect_equal(n_reviews, 17)

  # 2. The headline: the share of cuts whose prior was already significant.
  #    This is the figure the whole finding rests on, so it is pinned exactly
  #    rather than as "most".
  line <- grep("cuts had an already-significant prior", out, value = TRUE)
  expect_length(line, 1)
  counts <- as.integer(regmatches(line, gregexpr("\\d+", line))[[1]])
  expect_equal(counts[1], 168)   # significant priors
  expect_equal(counts[2], 185)   # cuts in total
  expect_equal(counts[3], 91)    # per cent

  # 3. The two detectors that cannot answer, and the three that always can.
  answers <- function(m) {
    l <- grep(paste0("^\\s+", m, "\\s+\\d+ of \\d+ reviews"), out, value = TRUE)
    as.integer(sub(".*?(\\d+) of \\d+ reviews.*", "\\1", l))
  }
  expect_equal(answers("barrowman"), 4)
  expect_equal(answers("simulation"), 5)
  expect_equal(answers("ottawa"), 17)
  expect_equal(answers("rcma"), 17)
  expect_equal(answers("sufficiency"), 17)

  # 4. And the README must say the same numbers as the run just did.
  md <- readLines(readme, warn = FALSE)
  expect_true(any(grepl("168 of the 185 cuts \\(91%\\)", md)))
  expect_true(any(grepl("`simulation` \\| \\*\\*5 of 17\\*\\*", md)))
  expect_true(any(grepl("`barrowman` \\| \\*\\*4 of 17\\*\\*", md)))
})
