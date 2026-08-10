---
title: "Five published signals for updating meta-analyses, applied and evaluated: an exploratory study of when they can be used at all"
author: |
  Javier Núñez\
  \small Independent researcher · ORCID [0009-0003-9770-4986](https://orcid.org/0009-0003-9770-4986) · jnverbel@gmail.com
date: 2026-08-10
keywords: [systematic review updating, meta-analysis, evidence synthesis, research methodology, reproducibility]
---

# Abstract

**Background.** Five statistical methods for deciding whether a systematic
review has gone out of date were published between 1999 and 2007. They have
been run side by side exactly once, on 80 Cochrane reviews, where two flagged
nothing and the three that discriminated agreed at Kappa = 0.14. No reusable
software implementation of any of them appears to exist, so the comparison has
not been repeated.

**Methods.** We implemented all five in an open-source R package and applied
them to historical evidence in two arms. Arm A sweeps the 17 meta-analyses in
the `metadat` collection that carry per-study data and publication years,
running all five detectors against three operational definitions of the pooled
estimate having moved. Arm B assembles 6,686 consecutive version pairs
from the 4,132 Cochrane reviews with more than one version, harvested from
Europe PMC, and attaches the editorial outcome recorded in each version's
authors' conclusions, supporting the three detectors that need
only pooled estimates.

**Results.** In Arm A, 168 of 185 cuts (91%) had an already-significant prior
meta-analysis, so `barrowman` and `simulation` — both of which require a
non-significant prior — could be asked in only 4 and 5 of the 17 reviews. This
was invisible to the published comparison because its cohort was selected for
non-significance. The Ottawa method's effect criterion is a ratio of relative
risk reductions whose denominator approaches zero as the prior effect
approaches the null, so it fires on 64% of samples containing no change at all;
its specificity falls to 0.14 on a null review where the other four hold at
1.00. The stability half of the sufficiency method is defined as the slope of a
least-squares fit to the cumulative effect series; that slope has no valid null
distribution and fired on 209 of 300 samples of unchanging evidence. In Arm B,
560 pairs carry both a comparable pooled estimate and authors' conclusions at
each end; an automated screen separates likely conclusion changes 11.3-fold
across strata.

**Conclusions.** The five divide into two kinds of problem, and the distinction
matters. Two are simply not defined once the prior meta-analysis is
significant, which is the case at 91% of cuts here: that is a domain
restriction rather than a defect, but it means they usually cannot be asked at
all, and the one published comparison could not see this because its cohort was
selected for the single condition under which they can. The other two failures
are defects: one criterion is unstable by construction on exactly the reviews
its method targets, and one statistic has no valid null distribution. These are
exploratory findings on 17 reviews with no held-out set, scored against
operational targets that observe the pooled estimate moving rather than what
any review team did. The software, the sweeps and the corpus
are reproducible from public sources.

---

# 1. Introduction

A meta-analysis is a photograph, not a standing fact. About a quarter of
systematic reviews are out of date within two years of publication and half
within five and a half [@shojania2007]. Deciding *when* a particular review
needs updating remains largely a matter of judgement, and the consensus
checklist for that judgement [@garner2016] treats the pooled estimate as one
input among several.

Five statistical methods have been proposed to inform the quantitative part of
that decision, published between 1999 and 2007: recursive cumulative
meta-analysis [@ioannidis1999], the sufficiency-and-stability indicators
[@mullen2001], Barrowman's participant-ratio criterion [@barrowman2003], the
Ottawa method [@shojania2007], and prospective-power simulation [@sutton2007].

They have been run side by side exactly once. Pattanittum and colleagues
[@pattanittum2012] applied all five to 80 Cochrane reviews and reported that
two flagged nothing at all, while the three that did discriminate — Ottawa
flagging 34 reviews, recursive CMA 7, Barrowman 7 — agreed at Kappa = 0.14,
indistinguishable from chance. Only one review was flagged by all three.

That result has stood unexamined for over a decade, and the reason is
mundane. Each method was published as a description in a paper, and we could
not find a reusable software implementation of any of them. A search of the
metadata of all 24,734 CRAN packages on 2026-08-09 returns no implementation of
the Ottawa method, of Barrowman, or of recursive cumulative meta-analysis as an
updating diagnostic; every hit on those names is a false positive, and each is
adjudicated in the package's `inst/cran-search/`. The claim is *we did not find
one*, not *none exists*: the search reads package metadata rather than source
code, and `metafor` itself matches neither "cumulative meta-analysis" nor
"fail-safe" in its own metadata while exporting `cumul()` and `fsn()`.

The consequence is that a team wanting to know which updating signal to trust
must first reimplement five methods from prose, and each reimplementation is a
fresh opportunity to differ from every other.

This paper reports what happens when the five are implemented once, in the
open, and run against real history. The findings are mostly negative, and the
most useful of them is about **when the methods can be used at all** rather
than about whether they are right.

# 2. What this study can and cannot establish

We state the bounds before the results, because they qualify every number that
follows.

**The evaluation targets are not outcomes.** Arm A scores each detector against
three operational definitions of the pooled estimate having moved between a cut
point and a later target. None of the three observes what a review team did,
whether a recommendation changed, or whether anyone was harmed by acting on the
old estimate. A sensitivity of 0.32 therefore means "agreed with a stated
criterion about the estimate 32% of the time", not "was right 32% of the time".
Two of the three share an identical numerator and differ only in which standard
error divides it, so they are one distance on two scales rather than
independent checks; agreement between them is expected and disagreement is the
informative case.

**Arm A is 17 reviews and nothing is held out.** The 17 are what survives from
the 110 data frames in `metadat`, and 54 of the 93 exclusions are for a single
reason: the dataset records no per-study publication year. Sets that do record
one skew towards the well-curated classics, so the 17 are a convenience sample
whose selection may correlate with what is being measured, in a direction we do
not know. Every figure reported from Arm A is generated and stated on the same
17.

**Two of the five detectors are not literal reproductions.** The stability half
of the sufficiency method is replaced by a change-point statistic, for reasons
given in §5.3, and the substitution is carried in the function's name. The
simulation detector simulates effects rather than participants, because no
participant-level data is available. Both departures are measured; neither is
hidden.

**Arm B does not reach the per-study data.** Cochrane full text is behind a
paywall, and the analysis tables are absent from the openly available full text
in Europe PMC (checked on five reviews: no data tables in any). The finest
grain available at scale is therefore the pooled estimate, which supports
`rcma`, `ottawa` and `barrowman`, and not `sufficiency_changepoint` or
`simulation`. That is a property of data access, not of the design.

# 3. Methods

## 3.1 Implementation

All five detectors are implemented in `staleness`, an MIT-licensed R package
built entirely on `metafor` [@viechtbauer2010] for model fitting. Each detector
is a pure function with a common signature and returns a verdict, a signal and
a detail record. A detector that cannot answer returns `not_applicable` rather
than `current`, so declining to answer is never scored as a correct call.

Backtesting refits the meta-analysis at every cut point from the studies
available at that date, so no detector sees evidence that did not yet exist.
Cuts too close to the end of the series to have an observable future are marked
censored and excluded from the metrics rather than scored against a target that
cannot be known.

Where a published procedure is ambiguous, the reading adopted is stated and
justified from worked examples in the applying literature rather than chosen
for convenience. §5.2 gives the case where this mattered most.

## 3.2 Arm A: historical sweep with per-study data

Every data frame in `metadat` was screened for a per-study publication year, a
two-group effect measure constructible with `escalc()`, at least eight studies,
one effect per study-year, and at least three uncensored cuts. Seventeen
survived. Each was backtested with one set of parameters (`horizon = 3`,
`window = 5`, `min_k = 3`) so the comparison is like for like, and all five
detectors were run at every yearly cut.

## 3.3 Arm B: version chains with an editorial outcome

Cochrane assigns each review a stable identifier and each version a suffix
within the DOI: `10.1002/14651858.CD005343.pub7` is review `CD005343`, version
7. Grouping by identifier and sorting by suffix therefore reconstructs a
review's full version history from metadata alone, with no additional lookup;
Crossref independently exposes `update-to` and `updated-by` relations with
dates.

We harvested 17,247 review versions from Europe PMC's open REST API, covering
9,862 distinct reviews. Of these, 4,132 have more than one version and so
contribute the 6,686 consecutive pairs that carry an update; the rest have a
single version and cannot. The authors'
conclusions section is present in the free abstract of each version, giving the
editorially recorded position at both ends of a pair.

Comparing pooled estimates across a pair requires that both quote the same
outcome. We require the same effect measure and a text-similarity of at least
0.80 between the words preceding each estimate; 746 pairs meet this, of which
560 also carry conclusions at both ends.

Nothing in Arm B touches cochranelibrary.com. Every input is served by Europe
PMC or Crossref under their normal terms, which is what allows the corpus to be
deposited alongside this paper rather than described and withheld.

# 4. Reproducibility

Three sweeps regenerate every quantitative claim below from public sources:
`inst/applicability/` produces Arm A, `inst/calibration/` regenerates each
calibration figure from its original seeds **and reconstructs the statistics
that were replaced**, so the before-and-after can be re-derived rather than
taken on trust, and `corpus/` rebuilds Arm B in about ten minutes.

The package carries 1,176 tests and the corpus 35, with mutation testing as the
standard: for each guard, the code is changed and the test is required to turn
red. This is reported because it caught errors that reading did not, including
four in the corpus pipeline described in §5.5.

# 5. Results

## 5.1 Two of five detectors are structurally unable to answer

Across the 17 reviews of Arm A and their 185 yearly cuts, **168 cuts (91%) had
an already-significant prior meta-analysis**. In 11 of the 17 reviews, every
cut did.

`barrowman` and `simulation` both require a non-significant prior: the former
asks how many participants would be needed to reach significance, the latter
simulates the power of the next batch to achieve it. Neither question is
defined once the prior is already significant. They can therefore be asked in
only **4 and 5 of the 17 reviews** respectively.

| Detector | Published? | Can answer | Mean sensitivity | Reviews scoring zero |
|---|---|---|---|---|
| `ottawa` | as published | 17 of 17 | 0.320 | 8 |
| `rcma` | as published | 17 of 17 | 0.309 | 9 |
| `barrowman` | as published | **4 of 17** | 0.167 | 2 of 3 |
| `simulation` | effect-level | **5 of 17** | 0.167 | 2 of 3 |
| `sufficiency_changepoint` | **substitute** | 17 of 17 | 0.041 | 11 |

The second column is not decoration. `sufficiency_changepoint` runs a statistic
its source never described, so **0.041 is a property of our substitute and
carries no information about the published sufficiency method** — whose own
stability statistic is shown in §5.3 to have no valid null distribution and
therefore no meaningful sensitivity to report. `simulation` simulates effects
rather than participants. Neither row may be read as evidence about the
procedure its name refers to, and no conclusion below does so.

This is not a finding about those methods being wrong. It is a finding about
when they can be used, and it was **invisible to the published comparison
because its 80 reviews were selected for having a non-significant pooled
result** — the one cohort in which these two detectors can always be asked.

The figure does not depend on the backtest configuration: run at `horizon = 6`
the sweep returns the same 17 reviews, the same 168 of 185 cuts, and the same
coverage. That is expected rather than reassuring — whether a prior was already
significant is a fact about the evidence at a cut — but it means the finding is
a property of the reviews and not of this parameterisation.

## 5.2 The Ottawa effect criterion is unstable where the method is meant to be used

The Ottawa method states its quantitative signals as "changes in statistical
significance or relative changes in effect magnitude of at least 50%"
[@shojania2007]. The phrasing does not say *of what*, and the choice matters:
the two readings give materially different detectors.

Two independent applications resolve it the same way, on data that can be
checked. Pattanittum et al. [@pattanittum2012], Table 1, computes the change on
relative risk *reductions*, `(1 - RR_new) / (1 - RR_prev)`; of the ten reviews
with the largest Ottawa indicator in that study's appendix, all ten fire under
that reading and none fires under the ratio of the risk ratios. Mickenautsch
and Yengopal [@mickenautsch2013] work two further examples — RR 2.10 to 1.51
and RR 2.61 to 1.66 — describing both in their own words as "a change in
relative effect size of over 50%". Neither meets that on the ratio of the
effects (0.719, 0.636) nor as a percentage change in the estimate (28%, 36%);
a fourth candidate, the change expressed over the *new* estimate, gives 39% and
57% and so clears the bar on one example and not the other. Only the ratio of
risk reductions fits all four worked examples across the two papers.

The correction carries a finding. The denominator `1 - RR_prev` approaches zero
as the prior effect approaches no effect, so the criterion is unstable
precisely on the null meta-analyses the method targets. On simulated evidence
containing **no change at all**, the effect signal fires on **64% of samples
under a null effect** and on **0%** once the effect is real and precise.
Confirmed outside simulation on `metadat::dat.laopaiboon2015`, a null review
where this detector's specificity falls to **0.14** while the other four hold
at 1.00.

This is not an implementation artefact; it follows from the criterion as
written. It also accounts for the published comparison: Ottawa flagged 34 of 80
reviews where recursive CMA and Barrowman each flagged 7, on a cohort of null
meta-analyses by inclusion criterion.

We did not correct it. Correcting it would mean implementing a different
method, and the point of the exercise is to find out how the published ones
behave. What we corrected is the silence: the verdict now reports the
denominator it divided by, and flags when that denominator is near zero.

## 5.3 The sufficiency stability statistic has no valid null distribution

The stability half of the sufficiency method is described as the "absolute
slope of the linear regression fitted across the cumulative treatment effects
versus information increment" [@pattanittum2012]. Taken literally the rule is
degenerate — on continuous data that slope is never exactly zero — so a
significance rule has to stand in for it.

Three implementations were tried and measured:

1. **The ordinary least-squares t-test.** A cumulative mean is
   near-perfectly autocorrelated by construction and converges on the pooled
   effect by the law of large numbers, so the test detects *convergence* and
   reports *instability*. Over 300 samples of genuinely unchanging evidence it
   returned `out_of_date` **209 times**, where `rcma` and `ottawa` returned it
   none.

2. **Permuting study order, keeping the slope.** This brings the false-alarm
   rate to 16 of 300, nominal on average, but the slope of a cumulative series
   is dominated by its first few points, so a *late* change cannot clear the
   permutation null: power against ten new studies at RR 0.30 was **1 in 200**,
   falling to zero as the shift grew larger. Worse, the permutation null is
   itself invalid when study variances change over calendar time: with early
   small trials, later large ones and no drift at all, it fired on **28%** of
   samples, and on a schedule of 20 small then 10 large trials, on **42%**.

3. **A change-point statistic** — the largest standardised difference between
   the studies before and after any split — assessed against the same
   order-permutation null. False alarms 15 of 300 (5.0%); power 200 of 200
   against ten new studies at each of four effect sizes; and the
   heteroscedastic no-drift false-alarm rate falls to 5.3%, 6.7% and 6.7% on
   the schedules that produced 28% and 42% above.

The third is what the package runs, and the detector is named
`sufficiency_changepoint()` rather than `sufficiency()` so the substitution is
visible at the call site and in every results table. The published slope is
still computed and returned as a diagnostic; it decides nothing.

The substitute is not uniformly calibrated either, and we report where it
fails as carefully as where the original does. Across eleven variance regimes
the false-alarm rate runs from 2.4% (conservative) to 11.1%, roughly twice
nominal, with the worst case a smooth 650:1 monotone precision ramp. The
pattern is not the direction of the trend — growth and decay are equally bad —
but its smoothness, monotonicity and range: the same `dat.bcg` variances
rearranged into sorted order move the rate from 5.0% to 9.4%, which is the
cleanest demonstration that it is the arrangement and not the numbers.

A second limit matters more for interpretation, and it is the one place where
this paper corrects a figure the package itself had published. The statistic is
a maximum over split points, and the split isolating a single study is among
them, so a lone discordant study arriving last can carry it. Measured by
`inst/persistence/persistence.R` over 400 replicates per cell, against a
baseline with no outlier at all:

| k | no outlier | one 5-SE study, last | one 5-SE study, mid-series |
|---|---|---|---|
| 20 | 0.050 | 0.035 | 0.005 |
| 30 | 0.052 | 0.207 | 0.003 |
| 40 | 0.043 | 0.395 | 0.000 |

So the effect is real and grows with `k` — at 40 studies a single discordant
final study forces `unstable` about eight times more often than chance — but it
is confined to studies arriving **last**. One in mid-series does essentially
nothing, because the split isolating it leaves large blocks on both sides.

An earlier version of this claim, in the package documentation, reported far
higher rates (0.995 at `k` = 40) and stated that such a study forces `unstable`
"almost always from about `k` = 25". We could not reproduce those figures. An
independent reimplementation that agrees with the shipped statistic exactly on
200 random inputs, and that reproduces the no-outlier control row, returns
0.44–0.48 at `k` = 40 under every reading of "five standard errors" we tried,
and saturates near 0.55 even at twenty. The documented table had no generating
script in the repository, which is why it survived; the figures above do, and
the qualitative claim — that this is substantially a single-outlier detector at
the tail — survives with it while the magnitude does not.

### A persistence requirement fixes most of it

If a shift has to persist to count, the split isolating one study stops being
available. Requiring both sides of a split to hold at least `r` studies gives,
averaged over `k` = 20, 30 and 40:

| r | drift detected | one outlier, last | one outlier, mid | no change |
|---|---|---|---|---|
| 1 (shipped) | 0.802 | 0.212 | 0.003 | 0.048 |
| 2 | 0.815 | 0.226 | 0.003 | 0.051 |
| 3 | 0.825 | 0.147 | 0.006 | 0.054 |
| 5 | 0.845 | **0.083** | 0.013 | 0.052 |

Requiring five is strictly better on this evidence: power against a real late
shift rises slightly (0.802 to 0.845), the single-outlier false positive falls
by more than half (0.212 to 0.083), and calibration under no change is
untouched (0.048 to 0.052). At `k` = 40 the outlier rate falls from 0.395 to
0.090.

We have not adopted it. The measurement is four regimes on simulated evidence
with one variance schedule, and the package ships the statistic these results
question rather than the one they favour, because changing it on the strength
of a single experiment would repeat the error this section documents. It is
stated as the next test to run, with the script that runs it.

## 5.4 An editorial outcome is reachable, for three of the five detectors

Arm B attaches to each version pair the authors' conclusions as recorded at
both ends. Of 6,686 consecutive pairs, 4,530 carry conclusions at both ends once
duplicated records are removed, and 560 additionally carry pooled estimates
that can be compared.

A screen combining text similarity with Cochrane's own controlled vocabulary —
GRADE certainty tiers, the standardised hedging statements, and the appearance
or disappearance of "insufficient evidence" — was applied to all pairs, and a
sample of 120 was drawn stratified across the score range and coded blind
against the definition of French et al. [@french2005], for whom a changed
conclusion is one that "alter[s] the substance or meaning of a section or
alter[s] the interpretation", with style rewrites explicitly excluded.

| Stratum | n | Major | Minor | None | % major (95% CI) |
|---|---|---|---|---|---|
| High | 40 | 34 | 6 | 0 | 85% (71–93) |
| Medium | 40 | 26 | 14 | 0 | 65% (50–78) |
| Low | 40 | 3 | 20 | 17 | 8% (3–20) |

The screen separates 11.3-fold, which is what makes it usable to order the
remaining pairs. The sample was stratified rather than taken from the top of
the ranking precisely so that this could be measured: a top-N sample would have
shown a high rate and said nothing about what the screen misses.

**The resulting rate of 52% is six times the 9% French et al. measured over 254
updated reviews, and we report it as a working figure rather than a finding.**
Four explanations are available and none has been tested. The largest is
probably selection: French's denominator is updated reviews, while ours is
pairs carrying a comparable quantitative estimate at both ends — a review that
reports a poolable effect twice has something that can move, and one that says
"insufficient evidence" throughout does not. Cochrane's restructuring of how
conclusions are written around 2011 is a second candidate, the greater volume
of modern per-outcome reporting a third, and a single coder rather than two
independent human investigators a fourth. Distinguishing "the corpora differ"
from "the coder is liberal" requires either human coding of a subsample or
access to the original French coding as an external standard; neither is yet
available.

## 5.5 What the tests caught that reading did not

Four defects in the corpus pipeline were silent, and three of them produced
well-formed output that would have entered the analysis unnoticed.

An unsuffixed DOI was read as version 1; because Europe PMC carries several
index records under the same bare DOI, 551 reviews had multiple "version 1"
records and 681 pairs (9%) had the same version at both ends. What surfaced it
was a **range check**: the gap between versions came back as −13 years, and a
negative gap cannot exist.

142 pairs had a byte-identical abstract at both ends. This is not the same as
an update whose conclusions stood — that case is real and common, 662 pairs
have identical conclusions with the Main results section visibly rewritten —
but the same record reaching the index twice.

The count of pairs usable for the detector analysis was initially reported as
1,907 and is 746, because nothing checked that two versions quoted the same
outcome; 42% name clearly different ones and 12% do not share an effect
measure. It surfaced from one case read by eye.

And 14 pairs carried an interval that does not contain its own point estimate.
Here the parser was correct and **the published abstracts are wrong**:
`MD -2.76, 95% CI 3.57 to 1.96` with the minus signs dropped, `RD 0.03, 95% CI
-0.01 to -0.07` running backwards. These are refused rather than repaired,
since guessing which digit is wrong invents data, and not reordered either,
since a sign lost on both bounds is not fixed by swapping them.

We report these because a corpus assembled from bibliographic metadata invites
exactly this class of error, and because the checks that caught them were
mostly range checks rather than sophisticated tests.

# 6. Discussion

The results are of two kinds and should not be read as one. `barrowman` and
`simulation` are **not defined** once the prior meta-analysis is significant,
which is the case at 91% of cuts in this sample. That is a property of the
questions they ask — how many participants would be needed to reach
significance, and what power the next batch has to achieve it — and neither is
a defect. It does mean that on evidence taken as it comes, they usually cannot
be asked.

The Ottawa effect criterion and the sufficiency stability slope are a different
matter, because both fail where they do apply: the first has a denominator that
approaches zero precisely on the null meta-analyses its method targets, and the
second has no valid null distribution at all. `rcma` is implementable as
described and behaves unremarkably.

These are not new methods and this is not a proposal to replace them. The
contribution is that the question can now be asked at all, repeatedly, by
anyone, on any body of evidence — and that when it is asked, the answers are
mostly unflattering.

**Coverage belongs beside sensitivity and specificity, not inside a
footnote.** The conventional way to compare detectors is a contingency table,
and a contingency table has no cell for "this method had no right to answer".
Our results suggest the first question is not which detector is most accurate
but *in what fraction of situations a detector is defined at all* — and that
those two questions can have opposite answers, since a method that answers
rarely can look excellent on the few occasions it does. The package therefore
returns a row with `n = 0` and `NA` metrics for a detector that never applied,
rather than omitting it, so that an absence is a value in the table instead of
something a reader has to notice.

The applicability result deserves emphasis because it changes how the one
existing comparison should be read. Pattanittum et al. did not find that
`barrowman` and `simulation` flag nothing because those methods are
insensitive; they found it on a cohort selected for the one condition under
which those methods are defined. On evidence taken as it comes, the problem is
more basic: they usually cannot be asked.

The route to a stronger claim is visible and partly built. Arm B attaches an
outcome that a human recorded, which is what the operational targets of Arm A
are not, and it reaches thousands of pairs rather than seventeen. What it
cannot reach, because the analysis tables are paywalled, is the per-study data
that two of the five detectors require. An outcome-anchored evaluation of all
five would need either those tables or an equivalent corpus with per-study
effects and dates.

# 7. Limitations

Restated compactly, because each has already qualified a result above.

Arm A is 17 reviews, selected by data availability rather than sampled, with
nothing held out; its figures are exploratory. Its evaluation targets observe
the pooled estimate moving and not any decision, so its rates measure agreement
with a stated criterion. Two of the five detectors depart from their published
procedure, and both are named for it. Arm B reaches only the pooled estimate,
supporting three detectors and not five; its outcome coding was performed once,
by an automated coder, and its 52% rate against a published 9% remains
unexplained and is reported as provisional. The comparability criterion for
pooled estimates is a text-similarity proxy for outcome identity, not a
verified match.

# 8. Availability

The package is public and MIT-licensed at
<https://github.com/jnverbel/staleness>. Arm A is regenerated by
`inst/applicability/applicability.R`, the calibration figures by
`inst/calibration/calibration.R`, and Arm B by the three scripts in `corpus/`.
The dated CRAN search behind the claim in §1 is in `inst/cran-search/`.

# References
