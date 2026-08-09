# The README is the first thing anyone reads, and its output block is written
# by hand -- unlike the vignettes, which are regenerated on every build. That
# makes it the one piece of documentation that can drift without anything
# noticing. It did: when ottawa's effect criterion moved to risk reductions,
# its signal went from 1.21 to 0.86 and the README kept promising 1.21.

test_that("the README's worked example still produces what it claims", {
  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  # Skipped on the file's absence, not on skip_on_cran(). The installed
  # package has no README.md above the test directory, so the check cannot
  # run there -- but skip_on_cran() also silences it in local development,
  # where NOT_CRAN is unset, and a test that never runs protects nothing.
  path <- "../../README.md"
  skip_if_not(file.exists(path), "README.md not reachable from the test dir")

  readme <- readLines(path, warn = FALSE)
  block <- grep("^#>\\s+signal:", readme, value = TRUE)
  claimed <- as.numeric(sub("^#>\\s+signal:\\s*", "", block))
  expect_length(claimed, 3)

  es <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos,
                        di = cneg, data = metadat::dat.bcg)
  st <- evidence_stream(metafor::rma(yi, vi, data = es), date = es$year)
  chk <- check_currency(snapshot_at(st, 1970), window_between(st, 1970, 1980),
                        methods = c("rcma", "ottawa", "sufficiency"))
  actual <- vapply(chk$verdicts, function(v) round(as.numeric(v$signal), 2),
                   numeric(1))

  expect_equal(unname(actual), claimed, tolerance = 1e-8)

  # The prose figures in the same block.
  expect_true(any(grepl("studies: 7 prior \\+ 6 new", readme)))
  expect_equal(chk$prev$k, 7)
  expect_equal(chk$updated$k - chk$prev$k, 6)
  expect_equal(round(chk$i2, 1), 92.2)
})
