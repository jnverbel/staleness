---
title: 'staleness: applying and calibrating published signals for updating meta-analyses in R'
tags:
  - R
  - meta-analysis
  - systematic reviews
  - evidence synthesis
  - research methodology
authors:
  - name: Javier Núñez
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 10 August 2026
bibliography: paper.bib
---

# Summary

A meta-analysis is a photograph, not a standing fact. About a quarter of
systematic reviews are out of date within two years of publication and half
within five and a half [@shojania2007]. Five statistical methods for detecting
this were published between 1999 and 2007 [@ioannidis1999; @mullen2001;
@barrowman2003; @shojania2007; @sutton2007], and the only study to run all five
side by side found that two flagged nothing at all and that the three which did
discriminate agreed at Kappa = 0.14, essentially chance [@pattanittum2012].

`staleness` applies those five detectors to an existing meta-analysis and the
evidence published since, and backtests them against the history of any body of
evidence, so their sensitivity, specificity and lead time can be measured
rather than assumed. It performs no literature searching and no screening, and
fits no models of its own: all estimation is delegated to `metafor`
[@viechtbauer2010].

# Statement of need

Each method was published as a description in a paper, and we could not find a
reusable implementation of any of them. Searching the metadata of all 24,734
CRAN packages on 2026-08-09 returns no implementation of the Ottawa method, of
Barrowman, or of recursive cumulative meta-analysis as an updating diagnostic;
every hit on those names is a false positive, and each is adjudicated in
`inst/cran-search/`, which holds the search as a runnable script and a dated
snapshot. The components exist — `metafor` computes cumulative meta-analyses
and Rosenthal's fail-safe N, `RTSA` covers trial sequential analysis, and
`metagear` screens literature, the half of the problem this package does not
attempt — but nothing assembles them into a decision about whether a review has
gone out of date.

The claim is *we did not find one*, not *none exists*: the search reads package
metadata, not source code, and `metafor` matches neither "cumulative
meta-analysis" nor "fail-safe" in its own metadata while exporting `cumul()`
and `fsn()`. The consequence is that the comparison question has stayed open. A
team wanting to know which updating signal to trust must first reimplement five
methods from prose, and each reimplementation is a fresh opportunity to differ
from every other.

The package is built around the failure mode such a comparison invites.
Defining "was this review really out of date?" with the same rule a detector
uses makes that detector correct by construction, so three evaluation targets
are implemented, and the one detector–target pair that remains circular is
named in an exported object, `CONTAMINATED_PAIRS`, whose flag travels on every
row of the results rather than as a footnote. Lead time — the interval between
a detector firing and the evidence actually moving — is reported, which no
methods paper does and which decides whether a signal is useful or merely
eventually correct. A detector that cannot answer returns `not_applicable`
rather than `current`, so declining is never scored as a correct call. Every
snapshot is refit from the studies available at that date, so no detector sees
evidence that did not yet exist, and cuts too close to the end of the series
are marked censored and excluded rather than scored against a target that
cannot be known.

# Declared substitutions

Two detectors are not literal reproductions, and both are named for it.
`sufficiency_changepoint()` replaces the stability half of the sufficiency
method [@mullen2001] — described as the slope of a least-squares fit to the
cumulative effect series — with a change-point statistic assessed against an
order-permutation null. That slope has no valid null distribution, since a
cumulative mean is autocorrelated by construction and convergent by the law of
large numbers, and on simulated evidence containing no change at all it fired
on 209 of 300 samples. It is still computed and returned as a diagnostic that
decides nothing. `simulation()` simulates effects rather than participants,
because the package never sees participant-level data. Both departures, their
calibration across variance regimes, and the regimes where they lose power are
measured in the `methods` vignette; `inst/calibration/` regenerates every
figure from the original seeds and reconstructs the replaced statistics, so the
before-and-after can be re-derived rather than taken on trust.

# Scope, and what the numbers mean

Deciding when to update a review is not a statistical question alone. The
consensus checklist [@garner2016] treats the pooled estimate as one input among
several, and living systematic reviews [@elliott2014] replace the discrete
update with continuous surveillance. `staleness` does none of that: it takes
evidence already found and screened and asks what five published signals say
about it, so it is not an alternative to those strategies and cannot be
benchmarked against them.

The three evaluation targets are operational choices, not ground truth. All
observe the pooled estimate moving; none observes what a review team did,
whether a recommendation changed, or whether anyone was harmed by acting on the
old estimate. Two share an identical numerator and differ only in which
standard error divides it, so they are one distance on two scales rather than
independent checks. A sensitivity reported here means "agreed with a stated
criterion", not "was right". That bound exists because no published corpus of
reviews with recorded update decisions exists to score against; building one is
the obvious next step, and it is not this package.

# Acknowledgements

`staleness` builds entirely on `metafor` [@viechtbauer2010] for model fitting.

# References
