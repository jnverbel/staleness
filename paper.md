---
title: 'staleness: Detecting and calibrating out-of-date meta-analyses in R'
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
date: 8 August 2026
bibliography: paper.bib
---

# Summary

A meta-analysis is a photograph, not a standing fact. About a quarter of
systematic reviews are out of date within two years of publication and half
within five and a half [@shojania2007], yet deciding *when* a particular review
needs updating is still largely a matter of judgement. Five statistical methods
for making that decision have been published since 2001
[@ioannidis2001; @mullen2001; @barrowman2003], and the only study to run them
side by side found they agreed at Kappa = 0.14 — essentially chance
[@pattanittum2012].

`staleness` is an R package that does two things. It applies all five published
detectors to an existing meta-analysis and the evidence published since,
returning a verdict and a signal for each. And it backtests those detectors
against the history of any body of evidence, so that their sensitivity,
specificity and lead time can be measured rather than assumed. The package
performs no literature searching and no screening, and fits no models of its
own: all estimation is delegated to `metafor` [@viechtbauer2010].

# Statement of need

Every one of the five methods was published as a description in a paper, and
none had a reusable software implementation. That absence has a consequence
beyond inconvenience: it is why the comparison question has stayed open. A
research team that wants to know which updating signal to trust must first
reimplement five methods from prose, and each reimplementation is a fresh
opportunity to differ from every other.

`staleness` closes that gap, and is built around the failure mode such a
comparison invites. Defining "was this review really out of date?" with the
same rule a detector uses makes that detector correct by construction, so the
package implements three independent truth definitions — a shift in the pooled
estimate, a shift that would have surprised an analyst standing at the time,
and a change in the practical conclusion — and names the one detector–truth
pair that remains circular in an exported object, `CONTAMINATED_PAIRS`, whose
flag travels on every row of the results rather than as a footnote. It reports
lead time, the interval between a detector firing and the evidence actually
moving, which no methods paper reports and which decides whether a signal is
useful or merely eventually correct. And a detector that cannot answer returns
`not_applicable` rather than `current`, so that declining to answer is never
scored as a correct call.

Backtesting is done by refitting the meta-analysis at every cut point from the
studies available at that date, so no detector can see evidence that did not
yet exist. Cuts too close to the end of the series to have an observable future
are marked censored and excluded from the metrics instead of being scored
against a truth that cannot be known.

# Implementation note

One method is implemented with a declared deviation. The stability half of
the sufficiency method [@mullen2001] is described as the slope of an ordinary
least squares fit to the cumulative effect series. That test has no valid null
distribution — a cumulative mean is autocorrelated by construction and
convergent by the law of large numbers — and on simulated evidence containing
no change at all it fired on 209 of 300 samples. `staleness` replaces it with a
change-point statistic, the largest standardised split in the cumulative
series, assessed against an order-permutation null. The published slope is
still computed and returned alongside. The substitution, its measured
calibration across nine variance regimes, and the regimes in which it loses
power are documented in the package's `methods` vignette. Rosenthal's fail-safe
N, on which the sufficiency half rests, has been discredited [@becker2005] and
is implemented for fidelity to the published method rather than as an
endorsement.

# Acknowledgements

`staleness` builds entirely on `metafor` [@viechtbauer2010] for model fitting.

# References
