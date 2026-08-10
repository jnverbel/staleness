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
| distinct reviews | 9,862 |
| consecutive version pairs | 6,686 |
| pairs with authors' conclusions at **both** ends | 4,661 (70%) |
| minus pairs whose entire abstract is identical (see below) | 4,519 |
| pairs also with a parseable pooled effect at both ends | 1,837 |
| **of those, plausibly reporting the SAME outcome** | **~570 (31%)** |
| median gap between versions | 4 years |

**There are two denominators and they are far apart.** The outcome — did the
authors' conclusions change — needs only the conclusions text, so it rests on
**4,519 pairs**. The detector analysis needs two comparable pooled estimates,
and that is a much smaller set.

At the 9% base rate French et al. (2005) measured over 254 updated reviews,
4,519 pairs carry roughly 400 conclusion-change events and ~570 carry about
50. Size the study from those, not from the harvest total.

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
abstract level it often is not. Until it is applied, ~570 is the honest
ceiling and it is a proxy, not a verified count.

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

## Not in here

Per-study forest plot data. Cochrane full text is behind a paywall and
cochranelibrary.com returns 403 to any automated request, so the pooled
estimate is the finest grain available at scale. That supports `rcma`,
`ottawa` and `barrowman`, and not `sufficiency_changepoint` or `simulation`,
which need the per-study series. The validation substudy addresses that
separately, by hand, on a few dozen reviews.
