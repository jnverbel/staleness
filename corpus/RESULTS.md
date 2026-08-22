# Confirmatory evaluation on Arm B — result, and what is wrong with it

Executed **2026-08-22 07:09 −0500**, commit `347a281`, against
[`PROTOCOL.md`](PROTOCOL.md) (SHA-256 `c55593ec5ca8eb27…`, frozen 2026-08-10 and
quoted in the preprint) and [`AMENDMENT-01.md`](AMENDMENT-01.md) (2026-08-21).
Output: `data/confirmatory.json`. One run, per §8, and the script refuses to
overwrite it.

## The result

**Negative for all three detectors**, which §7 fixed in advance as a publishable
outcome: after two decades, no published updating signal shows confirmatory
utility against a recorded editorial outcome on held-out evidence.

560 pairs in 483 reviews. Threshold 0.104, 291 events (52.0%).

| detector | eligibility | sensitivity | specificity | n answered |
|---|---|---|---|---|
| `rcma` | 1.00 [1.00–1.00] | 0.03 [0.01–0.05] | 1.00 [0.99–1.00] | 541 |
| `ottawa` | 0.89 [0.86–0.92] | 0.11 [0.07–0.15] | 0.98 [0.97–1.00] | 499 |
| `barrowman` | 0.33 [0.29–0.37] | 0.04 [0.00–0.11] | 0.94 [0.80–1.00] | 62 |

No detector raised an error on any pair, so the 2% rule in §3 never bit.

`ottawa` is flagged **contaminated** throughout, as `CONTAMINATED_PAIRS` in the
package requires: it shares logic with the conclusion-change outcome, and that
is the outcome here. Its 0.11 is the most favourable number in the table and
the least trustworthy.

## What is wrong with it

**The outcome threshold was fixed on a prevalence belonging to a different
population.** `AMENDMENT-01 §A` set the threshold at the score reproducing a
reweighted prevalence of 52%, taken from `04-analyse-coding.py`, which reported
that figure as the prevalence of "pairs with comparable effects", i.e. the 560
the detectors run on.

That was wrong, and the error was in `04-analyse-coding.py`, not in the coding
or the sample. The stratified sample of 120 was drawn in `03-screen.py` from
pairs parsing an effect **at both ends** — 1,825 of them — and `04` reweighted
the strata by the size of a different and much smaller set. The 52% estimates
the prevalence of the 1,825.

The two sets are not interchangeable, because requiring the two effects to be
*comparable* selects pairs whose abstracts repeat the same outcome, and those
score low:

| stratum | share of the 1,825 sampled | share of the 560 analysed | measured event rate |
|---|---|---|---|
| alto | 33% | 6% | 85% [71–93] |
| medio | 33% | 19% | 65% [50–78] |
| bajo | 33% | 74% | 8% [3–20] |

Rebuilt from its own composition, the prevalence of the analysed set is **23%**,
about 131 events, not 52% and 291. The confirmatory run therefore labelled
roughly 291 pairs as events, and **147 of those sat in the bottom stratum**,
where 8% of pairs are real events. The primary outcome carries a large and
one-directional false-positive load, which depresses every sensitivity in the
table above.

Two smaller defects, both stated rather than discovered later:

- **S1 collapsed to 42 pairs of the 120**, for the same reason: only 42 of the
  coded pairs survive the comparability requirement. §4 expected S1 to carry
  the comparison against a code not produced by a threshold, and 42 pairs with
  8 events cannot carry it.
- **`classify()` treats an interval wider than 0.40 as inconclusive.** That
  0.40 appears in neither the protocol nor the amendment. It was written before
  any result was seen, but it was not pre-specified, and it decides the label
  printed next to each detector.

`04-analyse-coding.py` is fixed as of this commit: it now takes its population
from the same definition `03-screen.py` uses, and reports the analysed subset's
prevalence separately instead of conflating the two.

## What survives, and what does not

**Does not survive:** the individual sensitivities as point estimates. They are
computed against a labelling that is wrong often enough to matter, and they
should not be quoted as the measured sensitivity of these detectors.

**Survives:** the firing rate, which does not depend on the labelling at all.

| detector | fires on | of pairs answered | rate |
|---|---|---|---|
| `rcma` | 8 | 541 | **1.5%** |
| `ottawa` | 31 | 499 | **6.2%** |
| `barrowman` | 3 | 62 | **4.8%** |

A detector that says `out_of_date` on 1.5% of the pairs it answers cannot
detect a quarter of them, nor a half, whatever the labelling. That is the
finding, and it is arithmetic rather than inference.

The pre-specified **S2** already bounded this across the extremes the stratum
intervals allow — prevalence 0.41 and 0.63 — and the sensitivities moved from
0.02 to 0.12 across that whole range. The corrected 23% falls outside it, which
is why the post-hoc re-analysis below was run rather than argued.

## Post-hoc re-analysis at the corrected prevalence

`06-posthoc-prevalence.R`, run after the results were seen and **labelled post
hoc**, per the protocol's own instruction. It does not replace
`data/confirmatory.json`. Threshold 1.346, 128 events (22.9%):

| detector | sensitivity | specificity | firing rate |
|---|---|---|---|
| `rcma` | 0.032 [0.007–0.067] | 0.990 | 1.5% |
| `ottawa` | 0.158 [0.095–0.227] | 0.966 | 6.2% |
| `barrowman` | 0.071 [0.000–0.179] | 0.971 | 4.8% |

**Negative for all three, again.** Correcting the prevalence roughly halves the
event count and moves `ottawa` from 0.11 to 0.16; §7 asked for 0.60 with a
lower bound above 0.50. Nothing here comes near it, and the firing rates are
identical to three significant figures, because they never depended on the
labelling.

## How to report this

The defensible sentence is not "the sensitivity of `rcma` is 0.03". It is:

> On 560 consecutive Cochrane review pairs, the three detectors that can be
> evaluated from pooled estimates alone fire on 1.5%, 6.2% and 4.8% of the
> pairs they answer, against an event prevalence of between 23% and 52%
> depending on the labelling. None approaches the pre-specified threshold for
> confirmatory utility under either.

That claim is robust to the defect. The per-detector sensitivities should be
reported with the labelling error attached, or not reported as measurements at
all.
