# Contributing to staleness

Bug reports, questions and pull requests are all welcome, at
<https://github.com/jnverbel/staleness/issues>.

## Reporting a bug

The most useful report is a small, self-contained example that reproduces the
problem — ideally one that builds its own `yi` and `vi` rather than depending
on an external dataset. Please include the output of `sessionInfo()`.

Statistical behaviour counts as a bug. If a detector fires where you can show
it should not, or stays silent where it should speak, that is worth reporting
even if nothing errors: several of the defects fixed so far produced a
confident, wrong number rather than a failure.

## Pull requests

* Every change needs a test, and the test should fail without the change.
  Please check that: revert your fix, confirm the test goes red, restore it.
  A test that passes both ways documents an intention rather than protecting
  one.
* Code and documentation are in English; the prose is British spelling
  (`Language: en-GB`).
* Run `devtools::document()` if you touched roxygen comments, and
  `R CMD check --as-cran` before opening the request. CI runs the check on
  macOS, Windows and Linux, from R 4.2 to R devel.
* If a change alters a number that appears in a vignette, update the vignette
  in the same commit. Prose figures are not regenerated automatically and will
  otherwise go quietly stale.

## Scope

This package applies published methods and calibrates them. It deliberately
does not search literature, screen studies, or fit models of its own — all
estimation is delegated to `metafor`. Proposals that add a *new* detector are
welcome if the method is published and citable; the point of the package is to
make existing methods comparable, not to invent another one.

Where an implementation departs from its published description, the deviation
is documented in `?the-function` and in `vignette("methods")`, along with the
evidence that motivated it. Please hold new contributions to the same standard.
