# Pre-specified protocol — confirmatory evaluation of the detectors on Arm B

**Status: frozen and not yet executed.**

This document fixes every analytic decision before any detector is run against
the Arm B corpus. It exists because Arm B is the only held-out evidence this
project has, and it can only be spent once.

Nothing below may be revised after results are seen. Anything run afterwards is
post hoc and will be labelled so in the manuscript, including any decision to
report a subgroup, drop an exclusion, or change a metric.

## 0. Provenance of the held-out status

Detector behaviour frozen at `f59bb41` (2026-08-09 19:30 −0500). Arm B corpus
first introduced at `c7e9505` (2026-08-10 07:56 −0500). Since `c7e9505`,
`git diff` reports zero non-comment lines changed in any detector source.

At the time of writing this protocol, **no detector has ever been run against
Arm B**: nothing under `corpus/` invokes one, and no harvested record carries a
verdict field. Both are checkable from the repository.

Execution commit and this file's hash are recorded in §9 before the run.

## 1. Unit and estimand

The unit is a **consecutive version pair** of one Cochrane review.

For each detector *d* we estimate two quantities and report them separately,
because a detector that answers rarely can look excellent on the occasions it
does:

- **Eligibility.** `P(d can be asked | pair)` — the fraction of pairs at which
  the detector is defined at all.
- **Conditional discrimination.** `P(d fires | conclusion changed, d eligible)`
  and `P(d does not fire | conclusion unchanged, d eligible)` — sensitivity and
  specificity *within* the eligible subset.

No single composite score is computed. Combining coverage with accuracy into
one number is the error §6 of the manuscript argues against.

Reviews contribute more than one pair. All uncertainty is clustered by
`review_id`: a nonparametric bootstrap resampling **whole reviews**, 2,000
replicates, percentile intervals. Pairs are not resampled.

## 2. Eligible population

The analysis set is pairs with `comparable_effects == TRUE` and authors'
conclusions at both ends: **n = 560** at the time of writing.

Per detector, within that set:

| Detector | Evaluable | Why |
|---|---|---|
| `rcma` | yes | needs prior and updated pooled effect only |
| `ottawa` | yes, on ratio measures | its effect criterion is a ratio of risk reductions and is undefined elsewhere; the significance signal alone is not the published detector |
| `barrowman` | only where participant counts parse at both ends | it sums participants across a snapshot |
| `sufficiency_changepoint` | **no** | needs the per-study cumulative series |
| `simulation` | **no** | needs per-study variances of the new evidence |

The two excluded detectors are excluded for lack of data, not for anticipated
performance, and this is fixed before execution.

## 3. Handling of `not_applicable`

Three states are kept apart and never merged:

1. **Eligible and answered** — enters sensitivity and specificity.
2. **Eligible but `not_applicable`** — the detector was asked and declined.
   Counted in the eligibility denominator, excluded from discrimination, and
   reported as its own column.
3. **Not evaluable** — the pair lacks the input the detector requires. Excluded
   from both, reported as a count.

A `not_applicable` verdict is never scored as `current`. An error is not a
verdict: any pair where a detector raises is logged, excluded, and its count
reported. If errors exceed 2% for any detector the run is reported as failed
for that detector rather than patched.

## 4. Outcome

The outcome is **a major change in the authors' conclusions** between the two
versions, per French et al. (2005): a change that alters the substance or
meaning, with style rewrites excluded.

Its validation status is partial and that is declared now, not discovered
later. The coding available at execution time is one automated pass over the
whole corpus, calibrated against 120 blind-coded pairs stratified by screen
score. The primary analysis uses that coding. **Two sensitivity analyses are
pre-specified:**

- **S1.** Restricted to the 120 pairs with a blind human-comparable code.
- **S2.** Under the two most extreme reweightings consistent with the stratum
  confidence intervals in §5.4 of the manuscript, to bound how far the outcome
  prevalence could move.

If independent human coding arrives before execution, it becomes the primary
outcome and the automated coding becomes S1. If it arrives after, it is a
labelled post-hoc re-analysis.

## 5. Exclusions, fixed in advance

- Pairs where either version is a protocol.
- Pairs flagged `same_abstract` (duplicate index records).
- Pairs whose two versions do not both parse a comparable effect.
- Pairs with a non-positive interval between versions.

No exclusion may be added after results are seen.

## 6. Metrics reported

Per detector: eligibility, sensitivity, specificity, false-alarm rate, the
`not_applicable` rate within eligible pairs, and the number of pairs behind
each. Every rate with a review-clustered bootstrap interval. Contaminated
detector–outcome pairs are flagged as they are elsewhere in the package.

No p-values. Nothing here is a hypothesis test.

## 7. What counts as a result

Fixed before execution so that neither direction can be presented as the
expected one.

- **Positive.** A detector reaches sensitivity ≥ 0.60 and specificity ≥ 0.60,
  with the lower bound of both intervals above 0.50, on eligibility ≥ 0.50.
- **Negative.** No detector clears the above. This is a publishable result and
  the manuscript is written to accommodate it: after two decades, no published
  updating signal shows confirmatory utility against a recorded editorial
  outcome on held-out evidence.
- **Inconclusive.** Intervals too wide to separate the two, which given ~560
  pairs and an expected event rate above 0.30 would itself be informative about
  the ceiling of this design.

The manuscript's argument does not depend on which of the three occurs. That is
the condition under which a held-out set is worth opening.

## 8. Execution

One run. A single script, `corpus/05-confirmatory.R`, written after this
protocol and reviewed against it before execution. Output written once to
`corpus/data/confirmatory.json` and not regenerated with altered options.

## 9. Registration

To be completed immediately before the run and not after:

- Protocol file SHA-256: _______
- Repository commit at execution: _______
- Date of execution: _______
- Human coding available at execution: yes / no

## 10. What this protocol does not cover

The persistence variant of the change-point statistic (`r` = 5) is **not**
evaluated here. It is an exploratory result from `inst/persistence/`, the
package ships `r` = 1, and adopting a modification on the strength of the
experiment designed to study it is the error this project has documented
elsewhere. Any evaluation of it needs its own held-out evidence.
