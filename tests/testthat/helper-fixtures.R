# Fixtures shared across test files. testthat sources helper-*.R before any
# test, so these are available everywhere without each file rebuilding them.
#
# The BCG escalc() call appeared verbatim in ten places, which the decision
# ledger flagged as deferred duplication. Repeated setup is not just noise: it
# is ten copies to keep in step when the fixture has to change.
#
# Note the asymmetry with the roxygen `\examples{}` blocks, which repeat the
# same data on purpose. R CMD check runs each example independently, so an
# example that leans on a shared object is an example that does not run.

bcg_es <- function() {
  metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                  data = metadat::dat.bcg)
}

bcg_ni <- function(es = bcg_es()) {
  es$tpos + es$tneg + es$cpos + es$cneg
}
