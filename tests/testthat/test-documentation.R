# The two most-repeated requests from CRAN reviewers are that every exported
# function document its \value and every exported dataset its \format. Neither
# is checked by `R CMD check`: man/truth.Rd documented four exported functions
# with no \value at all, and CONTAMINATED_PAIRS had no \format, through eight
# green platforms and a win-builder run.
#
# So the only thing that can catch the next one is a test. These read the
# package's own Rd database, which means they work against the installed
# package -- the same form a reviewer receives -- rather than against files in
# the source tree.

rd_tags <- function(rd) vapply(rd, function(x) attr(x, "Rd_tag"), character(1))

rd_database <- function() {
  db <- tryCatch(tools::Rd_db("staleness"), error = function(e) NULL)
  # Under pkgload::load_all() there is no installed Rd database, so fall back
  # to the source tree. Skipping instead would silence this in exactly the
  # place the documentation gets written -- the mistake this very file exists
  # to stop being repeatable. Rd_db(dir=) wants the package root, not man/.
  if (is.null(db) || !length(db)) {
    db <- tryCatch(tools::Rd_db(dir = "../.."), error = function(e) NULL)
  }
  db
}

test_that("every documented function states what it returns", {
  db <- rd_database()
  skip_if_not(!is.null(db) && length(db) > 0, "no Rd database available")

  missing <- names(db)[vapply(db, function(rd) {
    t <- rd_tags(rd)
    "\\usage" %in% t && !("\\value" %in% t)
  }, logical(1))]

  expect_equal(missing, character(0),
               info = paste("Rd pages with \\usage and no \\value:",
                            paste(missing, collapse = ", ")))
})

test_that("every exported dataset states its format", {
  db <- rd_database()
  skip_if_not(!is.null(db) && length(db) > 0, "no Rd database available")

  # A data object's \usage has no call syntax: it is just the object's name.
  # That is how an Rd page for data is told apart from one for a function.
  is_data_page <- function(rd) {
    t <- rd_tags(rd)
    if (!("\\usage" %in% t)) return(FALSE)
    usage <- paste(as.character(rd[[which(t == "\\usage")[1]]]), collapse = "")
    !grepl("(", usage, fixed = TRUE)
  }

  data_pages <- names(db)[vapply(db, is_data_page, logical(1))]
  # The suite would quietly pass if this ever found nothing to check.
  expect_true(length(data_pages) >= 1)

  missing <- data_pages[vapply(db[data_pages], function(rd) {
    !("\\format" %in% rd_tags(rd))
  }, logical(1))]

  expect_equal(missing, character(0),
               info = paste("data pages with no \\format:",
                            paste(missing, collapse = ", ")))
})
