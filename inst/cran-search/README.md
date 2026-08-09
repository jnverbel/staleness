# Is any of this already on CRAN?

`README.md` and `paper.md` both rest on a claim of absence: that none of the
five updating detectors had a reusable implementation. A claim about an entire
repository needs its evidence attached — without the terms, the fields and the
date, a reader cannot tell a thorough search from a cursory one, and neither
can the author a year later.

`cran-search.R` is that evidence. It needs network access; the tests do not
run it, they read the snapshot it writes.

    R CMD INSTALL .
    Rscript inst/cran-search/cran-search.R

or, while developing, without installing:

    Rscript -e 'pkgload::load_all("."); source("inst/cran-search/cran-search.R")'

## What it found (2026-08-09, 24,734 packages)

Searching the `Package`, `Title` and `Description` fields:

**No implementation of any of the five detectors.** Every hit on a method name
is a false positive, and they are listed rather than counted away:

| Package | Why it matched | Verdict |
|---|---|---|
| `SAiVE` | the University of **Ottawa**'s research group | false positive |
| `themis` | unrelated (class imbalance) | false positive |
| `updateme` | warns about **out-of-date R packages**, not reviews | false positive |
| `VegSpecIndex` | spectral indices for vegetation | false positive |
| `CRTSize` | "updated techniques" for cluster-trial sample size | false positive |

`barrowman`, `recursive cumulative`, `shojania`, `pattanittum`, `staleness`,
`currency of evidence` and `when to update` return nothing at all.

**The components do exist**, which is the other half of the claim and the part
that makes it defensible rather than merely bold: `metafor` (cumulative
meta-analyses, fail-safe N), `fsn` (Rosenthal's fail-safe number with
confidence intervals), `meta` (general meta-analysis including cumulative),
and `RTSA` (trial sequential analysis — a different sequential question).

**One neighbour is worth naming**: `metagear` provides tools for systematic
reviews including abstract screening. It does the half of the problem this
package explicitly does not do — finding and screening literature — and does
not assess whether a completed review has gone out of date.

## What this can and cannot establish

It searches package metadata. It does **not** search source code, help pages
or vignettes; doing so would mean downloading every package on CRAN. A package
implementing one of these criteria without naming it in its `Description`
would not appear.

That is not a hypothetical caveat, and the search demonstrates its own limit:

> **`metafor` matches neither `cumulative meta-analys` nor `fail-safe` in its
> metadata**, and yet it exports `cumul()` and `fsn()` — the two components
> this package builds on. Its `Description` simply does not use those words.

So the claim the evidence supports is **"we did not find one"**, not "none
exists", and `README.md` and `paper.md` say it that way.

## Keeping it honest over time

`cran-search-snapshot.csv` holds one row per term-and-hit, with the date and
the package count. Re-running the script overwrites it, and any package that
matches without a recorded verdict is reported as `UNADJUDICATED` rather than
absorbed into a total — so a new arrival is visible, and the claim can be
withdrawn if one turns up.

`tests/testthat/test-cran-search.R` holds the prose in `README.md` and
`paper.md` to the snapshot. It does not touch the network.
