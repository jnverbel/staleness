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

# A third class, found by a reviewer rather than by anything here. lead_time()
# changed its `within` default from Inf to the backtest's own horizon; the
# \arguments entry was updated and a \section three paragraphs above still
# said "The default, `Inf`, keeps the definition above". Both sentences were
# in the same help page, and neither R CMD check nor any test could see that
# they contradicted each other.
#
# This reads the actual formals and flags a page that names a default the
# function does not have. Two exemptions, both principled rather than
# convenient:
#
#   - A claim qualified as old/previous/earlier is a historical note, not a
#     claim about now. Several of these are deliberate, so the exemption is
#     matched explicitly.
#   - Only CONSTANTS are checked. "Defaults to `ma$ni`" and "defaults to the
#     backtest's own `horizon`" name where the effective value comes from when
#     the formal is NULL, which is a correct description of a real pattern
#     here. A bare literal -- `Inf`, a number, a string -- asserts what the
#     value IS, and that is the claim that can be flatly false. `Inf` is what
#     lead_time() wrongly claimed.
test_that("no help page claims a default the function does not have", {
  db <- rd_database()
  skip_if_not(!is.null(db) && length(db) > 0, "no Rd database available")

  literal <- function(x) {
    if (missing(x) || identical(x, quote(expr = ))) return(NULL)
    if (is.call(x) || is.name(x)) return(NULL)   # match.arg(), c(...), a symbol
    d <- tryCatch(deparse(x), error = function(e) NULL)
    if (length(d) == 1L) d else NULL
  }

  offenders <- character()
  for (nm in names(db)) {
    txt <- paste(as.character(db[[nm]]), collapse = "")
    fname <- sub("\\.Rd$", "", nm)
    f <- tryCatch(get(fname, envir = asNamespace("staleness")),
                  error = function(e) NULL)
    if (!is.function(f)) next

    have <- unlist(lapply(formals(f), literal))
    have <- c(unname(have), "NULL")   # an argument may document NULL meaning

    # "The default, `X`" / "the default is `X`" / "defaults to `X`", but not
    # "The old default was `X`".
    hits <- gregexpr(
      "(?i)(old |previous |earlier )?defaults?(?: is| was| to|,)?\\s*\\\\code\\{([^}]+)\\}",
      txt, perl = TRUE)
    m <- regmatches(txt, hits)[[1]]
    for (one in m) {
      if (grepl("(?i)^(old|previous|earlier) ", one)) next
      claimed <- sub(".*\\\\code\\{([^}]+)\\}.*", "\\1", one)
      # A symbol or an expression describes where the value comes from; only a
      # constant asserts what it is. See the note above.
      parsed <- tryCatch(str2lang(claimed), error = function(e) NULL)
      if (is.null(parsed) || is.call(parsed) || is.name(parsed)) next
      if (!claimed %in% have) {
        offenders <- c(offenders, paste0(nm, ": claims default ", claimed,
                                         "; formals offer ",
                                         paste(have, collapse = ", ")))
      }
    }
  }
  expect_equal(offenders, character())
})
