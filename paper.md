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
date: 8 August 2026
bibliography: paper.bib
---

# Summary

A meta-analysis is a photograph, not a standing fact. About a quarter of
systematic reviews are out of date within two years of publication and half
within five and a half [@shojania2007], yet deciding *when* a particular review
needs updating is still largely a matter of judgement. Five statistical methods
for making that decision were published between 1999 and 2007
[@ioannidis1999; @mullen2001; @barrowman2003; @shojania2007; @sutton2007]. The
only study to run all five
side by side, on 80 Cochrane reviews, found that two of them flagged nothing
at all, and that the three which did discriminate agreed at Kappa = 0.14 —
essentially chance [@pattanittum2012].

`staleness` is an R package that does two things. It applies all five published
detectors to an existing meta-analysis and the evidence published since,
returning a verdict and a signal for each — three of them as their sources
specify, two with declared substitutions described below. And it backtests
those detectors against the history of any body of evidence, so that their
sensitivity, specificity and lead time can be measured rather than assumed.
The package performs no literature searching and no screening, and fits no
models of its own: all estimation is delegated to `metafor`
[@viechtbauer2010].

# Statement of need

Every one of the five methods was published as a description in a paper, and
we could not find a reusable software implementation of any of them. Searching
the metadata of all 24,734 packages on CRAN on 2026-08-09 returns no
implementation of the Ottawa method, of Barrowman, or of recursive cumulative
meta-analysis as an updating diagnostic: every hit on those names is a false
positive, and each is listed and adjudicated in `inst/cran-search/`, which
holds the search as a runnable script and a dated snapshot. The components
exist — `metafor` [@viechtbauer2010] computes cumulative meta-analyses and
Rosenthal's fail-safe N, `fsn` and `meta` cover parts of the same ground, and
`RTSA` covers trial sequential analysis — and `metagear` screens literature,
the half of the problem this package does not attempt. But no package
assembles any of them into a decision about whether a review has gone out of
date.

The claim is stated as *we did not find one*, not *none exists*, because the
search reads package metadata rather than source code. The limit is
demonstrable rather than hypothetical: `metafor` matches neither "cumulative
meta-analysis" nor "fail-safe" in its own metadata, yet exports `cumul()` and
`fsn()`. That absence has a consequence
beyond inconvenience: it is why the comparison question has stayed open. A
research team that wants to know which updating signal to trust must first
reimplement five methods from prose, and each reimplementation is a fresh
opportunity to differ from every other.

`staleness` closes that gap, and is built around the failure mode such a
comparison invites. Defining "was this review really out of date?" with the
same rule a detector uses makes that detector correct by construction, so the
package implements three evaluation targets — a shift in the pooled estimate,
a shift that would have surprised an analyst standing at the time, and a
change in the practical conclusion — and names the one detector–target pair
that remains circular in an exported object, `CONTAMINATED_PAIRS`, whose flag
travels on every row of the results rather than as a footnote. It reports lead
time, the interval between a detector firing and the evidence actually moving,
which no methods paper reports and which decides whether a signal is useful or
merely eventually correct. And a detector that cannot answer returns
`not_applicable` rather than `current`, so that declining to answer is never
scored as a correct call.

Those three are **operational choices, not ground truth**, and the package
says so in the function that defines them. All three observe the pooled
estimate moving; none observes what a review team did, whether a
recommendation changed, or whether anyone was harmed by acting on the old
estimate. Two of them, `truth_shift` and `truth_surprise`, share the identical
numerator and differ only in which standard error divides it, so they are one
distance on two scales rather than two independent checks — agreement between
them is expected and disagreement is the informative case. A sensitivity
reported here therefore means "agreed with a stated criterion about the
estimate", not "was right". That bound is the price of there being no
published corpus of reviews with recorded update decisions to score against;
building one is the obvious next step, and it is not this package.

Backtesting is done by refitting the meta-analysis at every cut point from the
studies available at that date, so no detector can see evidence that did not
yet exist. Cuts too close to the end of the series to have an observable future
are marked censored and excluded from the metrics instead of being scored
against a truth that cannot be known.

# Fidelity, checked against the published application

Implementing a method from prose is where a comparison study silently goes
wrong, so each detector here is checked against the published applications of
these methods. That check found an error in this package, and because the
correction below carries a criticism of the Ottawa method, the evidence for
the reading comes first.

The original description [@shojania2007] states the criterion as "relative
changes in effect magnitude of at least 50%" without saying of *what*, and the
full text is behind a subscription. Two independent applications resolve it
the same way, on data that can be checked. Pattanittum et al. [-@pattanittum2012],
Table 1, computes the change on relative risk *reductions*,
`(1 - RR_new) / (1 - RR_prev)`; of the ten reviews with the largest Ottawa
indicator in that study's appendix, all ten fire under that reading and none
fires under the ratio of the risk ratios. Mickenautsch and Yengopal
[-@mickenautsch2013], applying the same method and tracing it to Shojania's
priority approach, work two examples — RR 2.10 to 1.51 and RR 2.61 to 1.66,
both declared in their own text to represent "a change in relative effect size
of over 50%". Neither meets that on the ratio of the effects (0.719 and
0.636), nor as a percentage change in the estimate (28% and 36%); both meet it
on the ratio of risk reductions (0.464 and 0.410, i.e. 54% and 59%). A fourth
candidate reading, the change expressed over the *new* estimate, gives 39% and
57% and so fires on one example but not the other. The ratio of risk
reductions is the only reading under which all four worked examples across the
two papers behave as their authors say. The tests reproduce all of them. A reader who takes Shojania's sentence to mean
the ratio of the effects gets a materially different detector, and the
criticism below applies to the reading above, not to that one.

The correction carries a finding. The denominator goes to zero as the prior
effect approaches no effect, so the criterion is unstable precisely on the
null meta-analyses the method targets: on evidence containing no change at
all, its effect signal fires on 64% of samples under a null effect and on 0%
once the effect is real. That is a property of the criterion, and it accounts
for the published comparison flagging 34 of 80 reviews by Ottawa against 7
each by the other two discriminating methods, on a cohort of null
meta-analyses.

# Implementation note

One method is implemented with a declared deviation. The stability half of
the sufficiency method [@mullen2001] is described as the slope of an ordinary
least squares fit to the cumulative effect series. That test has no valid null
distribution — a cumulative mean is autocorrelated by construction and
convergent by the law of large numbers — and on simulated evidence containing
no change at all it fired on 209 of 300 samples. `staleness` replaces it with a
change-point statistic, the largest standardised split in the cumulative
series, assessed against an order-permutation null. The published slope is
still computed and returned alongside, as a diagnostic that decides nothing.
The detector is named `sufficiency_changepoint()` rather than `sufficiency()`
so that the substitution is visible at the call site and in every results
table, not only in the documentation. The substitution, its measured
calibration across nine variance regimes, and the regimes in which it loses
power are documented in the package's `methods` vignette. Every figure behind
that decision is reproducible from the package itself: `inst/calibration/`
regenerates all of them from the original seeds, and reconstructs both
statistics that were replaced, so the before-and-after can be re-derived
rather than taken on trust. Rosenthal's fail-safe
N, on which the sufficiency half rests, has been discredited [@becker2005] and
is implemented for fidelity to the published method rather than as an
endorsement.

# Scope: a measuring instrument, not an updating strategy

Deciding when to update a review is not, in practice, a statistical question
alone. The consensus checklist for updating [@garner2016] treats the pooled
estimate as one input among several — whether the question is still relevant,
whether new methods or outcomes have appeared, whether the evidence would
change a recommendation — and living systematic reviews [@elliott2014;
@elliott2017] replace the discrete update with continuous surveillance,
supported increasingly by machine-assisted screening [@marshall2019]. Those
are multi-component workflows whose first step is finding and screening new
studies.

`staleness` does none of that, deliberately. It takes evidence that has
already been found and screened, and asks what five published statistical
signals say about it. So it is not an alternative to those strategies and
cannot be benchmarked against them: comparing a single statistical signal with
a surveillance workflow compares two different things, and a favourable result
either way would mean nothing.

What it can do is measure the component those strategies contain. Every
updating strategy that consults a pooled estimate is using one of these
signals or something like it, and until now no one could say how any of them
behaves on real history, because no implementation existed to run. The
applicability sweep in `inst/applicability/` is an example of what becomes
answerable: across the 17 reviews in `metadat` that carry a backtest, 168 of
185 cuts had an already-significant prior meta-analysis, so `barrowman()` and
`simulation()` — both of which require a non-significant prior — can be asked
in only 4 and 5 of the 17. That is not a finding about those methods being
wrong. It is a finding about when they can be used at all, and it was
invisible to the one published comparison because its 80 reviews were selected
for having a non-significant pooled result.

Three bounds on that figure, stated because they are not small. The 17 are
what survives from `metadat`'s 110 data frames, and 54 of the 93 exclusions
are for one reason: the dataset records no per-study publication year. Sets
that do record one skew towards the well-curated classics, so the 17 are a
convenience set whose selection may correlate with what is being measured, in
a direction nobody knows. Nothing is held out — the 91% and every other figure
here are generated and stated on the same 17. And two of the five detectors
are the substitutions declared above. These results are therefore
**exploratory**: they are strong enough to say the published comparison's
cohort was unrepresentative in a way that mattered, and not strong enough to
be read as calibration of the five methods in general.

The honest positioning, then, is that this is a reproducible platform for
studying updating signals, not an engine that decides when to update. The
package title says `Apply and Calibrate` for that reason.

# Acknowledgements

`staleness` builds entirely on `metafor` [@viechtbauer2010] for model fitting.

# References
