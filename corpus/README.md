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
| pairs also with a parseable pooled effect at both ends | **1,907 (29%)** |
| median gap between versions | 4 years |

1,907 is the denominator for the detector analysis. At the 9% base rate
French et al. (2005) measured over 254 updated reviews, that is roughly **170
conclusion-change events** — a usable positive class, and the reason to size
the study from this number rather than discover it later.

The effect parser is deliberately the same naive one the feasibility spike
measured, so the figure quoted in the write-up is the figure this produces. A
tuned parser raises 29%; it has not been tuned yet, on purpose.

## One defect found and fixed here, recorded because it was silent

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

## Not in here

Per-study forest plot data. Cochrane full text is behind a paywall and
cochranelibrary.com returns 403 to any automated request, so the pooled
estimate is the finest grain available at scale. That supports `rcma`,
`ottawa` and `barrowman`, and not `sufficiency_changepoint` or `simulation`,
which need the per-study series. The validation substudy addresses that
separately, by hand, on a few dozen reviews.
