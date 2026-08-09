# DESCRIPTION, CITATION.cff and .zenodo.json each state the version, the title
# and the licence by hand. Three hand-written copies of the same fact drift --
# that is how the README came to promise an ottawa signal of 1.21 while the
# package returned 0.86 -- and these three drift at the worst moment, because
# `version` and `date-released` are exactly what a release changes.
#
# Parsed with regular expressions rather than yaml/jsonlite: both files are
# small and flat, and adding two dependencies to Suggests so that a metadata
# test can run is a worse trade than a strict pattern.

field <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(gsub(pattern, "\\1", hit[1]))
}

test_that("DESCRIPTION, CITATION.cff and .zenodo.json agree", {
  # These three live at the package root and are excluded from the tarball via
  # .Rbuildignore, so they are genuinely absent under R CMD check -- unlike
  # inst/, which is installed. Skipping on the file's absence is right here
  # for the same reason it is right in test-readme.R and was wrong in
  # test-applicability.R.
  desc_path <- "../../DESCRIPTION"
  cff_path  <- "../../CITATION.cff"
  zen_path  <- "../../.zenodo.json"
  skip_if_not(all(file.exists(desc_path, cff_path, zen_path)),
              "metadata files not reachable from the test dir")

  desc <- read.dcf(desc_path)
  cff  <- readLines(cff_path, warn = FALSE)
  zen  <- readLines(zen_path, warn = FALSE)

  # 1. Version.
  cff_version <- field(cff, '^version:\\s*"?([0-9.]+)"?\\s*$')
  zen_version <- field(zen, '^\\s*"version":\\s*"([^"]+)".*$')
  expect_equal(cff_version, unname(desc[1, "Version"]))
  expect_equal(zen_version, unname(desc[1, "Version"]))

  # 2. Title. DESCRIPTION may wrap it across lines; the other two carry it on
  #    one, prefixed by the package name.
  desc_title <- gsub("\\s+", " ", unname(desc[1, "Title"]))
  cff_title  <- field(cff, '^title:\\s*"staleness: (.+)"\\s*$')
  zen_title  <- field(zen, '^\\s*"title":\\s*"staleness: (.+)",?\\s*$')
  expect_equal(cff_title, desc_title)
  expect_equal(zen_title, desc_title)

  # 3. Licence. DESCRIPTION says "MIT + file LICENSE"; the other two say MIT.
  expect_match(unname(desc[1, "License"]), "^MIT")
  expect_equal(field(cff, "^license:\\s*(\\S+)\\s*$"), "MIT")
  expect_equal(field(zen, '^\\s*"license":\\s*"([^"]+)".*$'), "MIT")

  # 4. The release date must be a real date, not a placeholder. It is the
  #    field most likely to be left behind at a tag.
  released <- field(cff, '^date-released:\\s*"?([0-9-]+)"?\\s*$')
  expect_match(released, "^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
  expect_false(is.na(as.Date(released, format = "%Y-%m-%d")))

  # 5. Both must name the maintainer as DESCRIPTION does, so a citation
  #    generated from either credits the same person.
  expect_true(any(grepl("Núñez", cff, fixed = TRUE)))
  expect_true(any(grepl("Núñez", zen, fixed = TRUE)))
})
