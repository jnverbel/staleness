# Outcome corpus — tier 1

Version chains of Cochrane reviews with the editorial outcome attached: did
the authors' conclusions change between one version and the next?

This is what turns the package's evaluation from "agreed with a stated
criterion" into "predicted what the review team did". Everything here comes
from Europe PMC's open REST API. Nothing touches cochranelibrary.com, so the
corpus does not depend on anyone's subscription and can be deposited with the
paper.

## Why the version chain needs no lookup

A Cochrane DOI carries it: `10.1002/14651858.CD005343.pub7` is review
`CD005343`, version 7. Grouping by review and sorting by version reconstructs
the whole history from metadata already in hand.

## Run

    python3 01-harvest.py     # -> data/versions.jsonl   (~10 min, 21 requests)
    python3 02-chains.py      # -> data/pairs.jsonl

## What it produced, 2026-08-10

| | |
|---|---|
| versions harvested | 17,247 |
| distinct reviews harvested | 9,862 |
| of those, with more than one version | 4,132 |
| consecutive version pairs | 6,686 |
| pairs with authors' conclusions at **both** ends | 4,661 (70%) |
| minus pairs whose entire abstract is identical (see below) | 4,530 |
| pairs also with a parseable pooled effect at both ends | 2,236 |
| **of those, with COMPARABLE effects** | **746 (11%)** |
| median gap between versions | 4 years |

**There are two denominators and they are far apart.** The outcome — did the
authors' conclusions change — needs only the conclusions text, so it rests on
**4,530 pairs**. The detector analysis needs two comparable pooled estimates,
and that is a much smaller set.

At the 9% base rate French et al. (2005) measured over 254 updated reviews,
4,530 pairs carry roughly 400 conclusion-change events and 746 carry about
67. Size the study from those, not from the harvest total.

## The correction that halved the detector denominator

An earlier version of this file said 1,907 pairs were usable for the detector
analysis. That was wrong, and wrong in the direction that flatters.

The effect parser takes the first `measure = value, 95% CI lo to hi` it finds
in an abstract. Nothing checked that version N and version N+1 were quoting
**the same outcome**. They frequently are not: measured on the context
preceding each match, only 31% of pairs look like the same outcome at both
ends, 42% look clearly different, and 12% do not even use the same effect
measure — an RR against an OR is a definite mismatch, not a borderline one.

It surfaced from a single case read by eye. `CD007145` went from
RR 0.72 [0.54, 0.95] to RR 1.14 [1.02, 1.27], which reads as a dramatic
reversal and is almost certainly two different outcomes. Comparing them would
have fed the detectors a signal computed from unrelated quantities, and every
detector would have returned a confident number.

French et al. avoided this with a pre-specified rule: the primary outcome is
the one the authors state, else the first listed under Objectives, else
mortality. Applying that rule needs the outcome to be identifiable, which at
abstract level it often is not. Comparability is now decided in the pipeline rather than eyeballed
afterwards: `comparable_effects` requires the same measure and a context
similarity of at least 0.80, and 746 pairs clear it. That is a proxy for
outcome identity, not a verified count, and the threshold is written where it
can be argued with.

The effect parser is deliberately the same naive one the feasibility spike
measured, so the figure quoted in the write-up is the figure this produces. A
tuned parser raises 29%; it has not been tuned yet, on purpose.

## Defects found here, recorded because each was silent

An unsuffixed DOI was read as version 1. That is wrong: Europe PMC carries
several index records under the same bare DOI for one review, so 551 reviews
had multiple "version 1" records, and zipping the sorted list paired those
with **each other**. 681 pairs (9%) had the same version at both ends.

Nothing about such a pair looks wrong — it has two abstracts, two conclusions,
two effects — and it would have entered the analysis as an update that never
happened. What surfaced it was a range check: the gap between versions came
back as −13 to 18 years, and a negative gap is impossible. 542 of the 681 had
inverted dates.

Fixed with two guards rather than one, since either alone leaves a hole:
records sharing (review, version) are collapsed keeping the earliest date, and
only ends whose versions actually differ are paired. Five pairs with genuinely
distinct versions still carry inverted dates (0.07%); those are real
indexing oddities and are left visible rather than dropped.

A second one, smaller: 142 pairs have a byte-identical abstract at both ends.
Those are not an update whose conclusions stood — that case is real and common
(662 pairs have identical conclusions while Main results is visibly rewritten)
— they are the same record reaching the index twice. Flagged in `pairs.jsonl`
as `same_abstract` and excluded by the screen, rather than dropped silently.

## Screening, and what it can and cannot do

`03-screen.py` produces interpretable flags and a score, then draws a sample
stratified across the score range. It orders the work; it does not decide it.
French's definition turns on whether a change "alter[s] the substance or
meaning", with style rewrites explicitly excluded, and no text statistic makes
that call.

The sample is stratified rather than top-N on purpose: adjudicating only the
highest scores would measure nothing about what the screen misses, and the
miss rate is exactly what decides whether the screen can be trusted to order
the rest.

One flag family is weaker than its overall rate suggests. The GRADE certainty
and hedging terms only bite on recent pairs, because the vocabulary is recent:
0.0% of pairs ending by 2010 move a certainty tier, 0.5% for 2011-2016, and
3.5% from 2017. Their low overall rate is a fact about when Cochrane
standardised its phrasing, not about how often certainty changes.

A third, which the tests found on their first run rather than a reader: 14
pairs carried an interval that does not contain its own estimate. The parser
was right and **the published abstracts are wrong** — `MD -2.76, 95% CI 3.57
to 1.96` with the minus signs dropped, `RD 0.03, 95% CI -0.01 to -0.07`
running backwards, `RR 1.07, 95% CI 1.98 to 1.18`. Not repaired, because
guessing which digit is wrong invents data, and not reordered either, since a
sign lost on both bounds is not fixed by swapping them. Refused at parse time,
so the invariant holds by construction.

## Tests

    python3 -m unittest discover -s corpus -p "test_*.py"

35 of them. Unit tests on the parsing run anywhere; invariants over the
produced files skip when `data/` is absent.

Every defect above has a test that fails without its fix, verified by mutation
rather than assumed:

| mutation | tests that go red |
|---|---|
| `comparable` stops requiring the same measure | 1 |
| `comparable` stops checking the outcome context | 2 |
| intervals that exclude their estimate are accepted | 1 |
| one null value for ratios and differences alike | 2 |
| `tier()` takes the first match instead of the lowest | 2 |
| the sample takes the top 120 instead of stratifying | 1 |

The invariants are mostly range checks, because a range check is what caught
the first defect: a gap between versions came back as −13 years, and a
negative gap cannot exist. The nesting invariant — comparable ⊆ has-effect ⊆
all pairs — exists so no future edit can quietly report the wider set as the
usable one, which is exactly the mistake that was made once already.

## The coding, and a rate that does not match the one it was planned around

All 120 sampled pairs were coded blind — no score, no flags, no stratum
visible — against French's three categories, with major defined as a change
that alters the substance, meaning or interpretation.

| stratum | n | major | minor | none | % major (95% CI) |
|---|---:|---:|---:|---:|---|
| high | 40 | 34 | 6 | 0 | 85% (71–93) |
| medium | 40 | 26 | 14 | 0 | 65% (50–78) |
| low | 40 | 3 | 20 | 17 | 8% (3–20) |

**The screen works.** 85% against 8% is an 11.3-fold separation, which is what
makes it worth using to order the remaining 440 pairs. The reweighted
population estimate equals the raw sample proportion here only because the
strata are equal thirds and 40 were drawn from each; the reweighting is in
`04-analyse-coding.py` so it stays correct if either changes.

**And the rate is 52%, where French measured 9%.** That is a sixfold
discrepancy and it is not a detail. Four explanations, none yet tested:

1. **Selection.** French's 9% is over *updated reviews*. This 52% is over
   pairs that carry a comparable quantitative estimate at both ends. A review
   that reports a poolable effect twice is a review with something that could
   move; one that says "insufficient evidence" throughout cannot change its
   conclusion much. This is probably the largest of the four.
2. **The GRADE transition.** Cochrane restructured how conclusions are written
   around 2011. Pairs spanning that boundary had their conclusions rewritten
   wholesale, and telling a template change from a substantive one is exactly
   the judgement French's "minor change" category exists to make.
3. **Era.** French compared 1998 against 2002. Modern reviews state findings
   per outcome with certainty grades, so there is simply more text in which a
   substantive change can occur.
4. **One coder, not two.** French used two investigators independently. This
   was coded once, and the borderline cases went to major more often than not.

**This makes human validation more important, not less.** A 52% rate produced
by a single automated coder, against a published 9% from two humans, is the
first thing a reviewer will attack — and correctly, because nothing here
distinguishes "the corpora genuinely differ" from "the coder is liberal". A
sub-sample coded by a person, or French's own coding as an external standard,
is what settles it. Until one of those exists, **52% is a working figure and
not a finding.**

The practical consequence is the opposite of bad news: ~290 expected events in
560 pairs rather than ~50, which is a well-powered study. It just cannot be
claimed yet.

## Not in here

Per-study forest plot data. Cochrane full text is behind a paywall and
cochranelibrary.com returns 403 to any automated request, so the pooled
estimate is the finest grain available at scale. That supports `rcma`,
`ottawa` and `barrowman`, and not `sufficiency_changepoint` or `simulation`,
which need the per-study series. The validation substudy addresses that
separately, by hand, on a few dozen reviews.
