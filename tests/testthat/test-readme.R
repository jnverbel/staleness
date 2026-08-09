# The README is the first thing anyone reads, and its output block is written
# by hand -- unlike the vignettes, which are regenerated on every build. That
# makes it the one piece of documentation that can drift without anything
# noticing. It did twice, in two different ways, and the second is why this
# file now runs the README's own code:
#
#   1. Its FIGURES drifted. When ottawa's effect criterion moved to risk
#      reductions its signal went from 1.21 to 0.86 and the README kept
#      promising 1.21. The check below catches that.
#
#   2. Its CODE broke. Requiring `study_id` invalidated the very first call in
#      the worked example, and the test that was supposed to guard the README
#      did not notice -- because it re-implemented the example in its own
#      words and compared only the numbers. A copy of the code cannot go
#      stale the way the original does.
#
# So the block is now extracted from the file and evaluated. Comparing output
# is not enough if the input is never run.
#
# On the eval(): the text being parsed is this package's own README, read from
# its source tree during its own test suite. There is no external or
# user-supplied input anywhere in the path, and the alternative -- a
# paraphrase of the code -- is precisely what failed to catch the breakage.
# Each evaluation gets a fresh environment so nothing leaks between tests.
# This is what knitr does to the same file when the vignettes are built.

readme_path <- function() "../../README.md"

# The worked example: the first ```r block that actually calls the package.
readme_example <- function(lines) {
  fence <- grep("^```", lines)
  for (i in seq(1, length(fence) - 1, by = 2)) {
    body <- lines[(fence[i] + 1):(fence[i + 1] - 1)]
    if (any(grepl("evidence_stream\\(", body))) {
      # `#>` lines are the promised output, not input.
      return(body[!grepl("^#>", body)])
    }
  }
  character(0)
}

test_that("the README's worked example actually runs", {
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  # Skipped on the file's absence, not on skip_on_cran(). The installed
  # package has no README.md above the test directory, so the check cannot
  # run there -- but skip_on_cran() also silences it in local development,
  # where NOT_CRAN is unset, and a test that never runs protects nothing.
  skip_if_not(file.exists(readme_path()),
              "README.md not reachable from the test dir")

  code <- readme_example(readLines(readme_path(), warn = FALSE))
  expect_gt(length(code), 0)

  # Evaluated in a fresh environment, exactly as written. If a reader copying
  # the README would hit an error, so does this.
  env <- new.env(parent = globalenv())
  expect_no_error(eval(parse(text = paste(code, collapse = "\n")), envir = env))

  # And the objects it claims to build must exist afterwards.
  expect_true(all(c("stream", "prev", "new") %in% ls(env)))
  expect_s3_class(get("stream", envir = env), "staleness_stream")
})

test_that("the README's worked example still produces what it claims", {
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  skip_if_not(file.exists(readme_path()),
              "README.md not reachable from the test dir")

  readme <- readLines(readme_path(), warn = FALSE)
  claimed <- as.numeric(sub("^#>\\s+signal:\\s*", "",
                            grep("^#>\\s+signal:", readme, value = TRUE)))
  expect_length(claimed, 3)

  # Re-run the README's own code rather than a paraphrase of it, then ask the
  # resulting objects for the figures the README prints.
  env <- new.env(parent = globalenv())
  eval(parse(text = paste(readme_example(readme), collapse = "\n")), envir = env)
  chk <- check_currency(get("prev", envir = env), get("new", envir = env),
                        methods = c("rcma", "ottawa", "sufficiency_changepoint"))
  actual <- vapply(chk$verdicts, function(v) round(as.numeric(v$signal), 2),
                   numeric(1))
  expect_equal(unname(actual), claimed, tolerance = 1e-8)

  # The prose figures in the same block.
  expect_true(any(grepl("studies: 7 prior \\+ 6 new", readme)))
  expect_equal(chk$prev$k, 7)
  expect_equal(chk$updated$k - chk$prev$k, 6)
  expect_equal(round(chk$i2, 1), 92.2)
})
