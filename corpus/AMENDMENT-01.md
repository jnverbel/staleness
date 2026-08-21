# Amendment 1 to the confirmatory protocol — 2026-08-21

**Written before execution. Nothing here was decided with a result in view:
at the time of writing no detector has been run against Arm B, which §0 of the
protocol makes checkable and which still holds.**

`PROTOCOL.md` is deliberately not edited. Its SHA-256 `c55593ec5ca8eb27…` is
quoted in the preprint published today under
[10.5281/zenodo.22050352](https://doi.org/10.5281/zenodo.22050352), and that
hash is the evidence that the analysis was specified before the evidence was
opened. Changing one byte of that file would destroy the thing it exists to
prove. So the two decisions below live here, dated, as an addition rather than
a revision.

Both fill gaps that the protocol assumed were already filled. Neither relaxes
anything it fixed.

## A. The outcome had no per-pair definition

§4 names the outcome — a major change in the authors' conclusions, per French
et al. (2005) — and says the coding available at execution time is *"one
automated pass over the whole corpus, calibrated against 120 blind-coded
pairs"*.

That pass does not exist. What exists is:

- `screen.score`, **continuous**, on all 560 eligible pairs;
- **120 pairs** carrying a blind code: 63 `major`, 40 `minor`, 17 `none`;
- a stratum-reweighted prevalence of **52%** — about 293 events in 560 — with
  the screen separating 85% in the top stratum from 8% in the bottom, a
  ratio of 11.3, which is the figure the preprint reports.

A continuous score is not a label, and without a label per pair there is no
sensitivity and no specificity. §1 and §6 cannot be computed as written.

**Rule adopted.** A pair counts as an event when its `screen.score` falls at or
above the threshold that reproduces the reweighted prevalence: the threshold
`t` such that `#{score >= t} / 560 = 0.52`, ties resolved in the direction that
keeps the count nearest 293 from below.

Why this rule and not another. The threshold is pinned by a quantity that was
estimated and published before this decision, from a sample drawn and coded
before it. It cannot be tuned toward a detector because it never looks at one.
The obvious alternative — the cut that maximises Youden's J against the 120
coded pairs — optimises a criterion computed on the outcome, which is a
defensible choice in general and the wrong one here, because it would let the
labelling absorb some of the discrimination the run is meant to measure.

What it costs, stated now. The labels are automatic and carry error: the
calibration says roughly 85% of top-stratum pairs and 8% of bottom-stratum
pairs are real events, so the label is wrong often enough to matter, and the
error is not symmetric. Every rate this run reports is a rate against an
automated label, not against adjudicated truth, and must be written that way.
The pre-specified S1 — the 120 pairs with a blind code — is what carries the
comparison against a code that was not produced by a threshold, and it stays
exactly as §4 specifies. Independent human coding, if it ever arrives, is a
labelled post-hoc re-analysis, per §4, and this amendment does not change that.

## B. barrowman's inputs were never parsed

§2 lists `barrowman` as evaluable *"only where participant counts parse at both
ends"*. No participant count was parsed anywhere: `N_PART` and `K_STUDIES` were
written in `02-chains.py` and never called, and no output record carried the
field. The detector would have fallen to *not evaluable* on all 560 pairs — not
because the corpus lacks the information, but because nobody had read it out.

Repaired in commit `a79bdbe`, before execution. `counts()` now records
`n_participants`, every participant match found (`n_participants_all`) and
`k_studies` for both ends of every pair. The eligible set is unchanged: the
same 560 pairs, the same effect values, and `sample.jsonl` still matches
`key.json` pair for pair.

**Rule adopted.** `n_prev` is the participant count parsed from the earlier
version. `n_new` is `n_to - n_prev`, and the pair is evaluable for `barrowman`
only where both counts parse and `n_to > n_prev`. A total that does not grow
does not license a guess at how much new evidence arrived; those pairs are
*not evaluable* under §3.3 and are reported as a count, never as `current`.

What it costs, stated now. 441 of 560 pairs parse a count at both ends, but
only **185** have a total that grows, so `barrowman` is asked on at most a
third of the corpus. And 282 of those 441 cite more than one participant
number in the same abstract, so "the first one" is a judgement rather than a
reading. Every match is kept in `n_participants_all` so the choice can be
audited instead of trusted. This is weaker evidence than the other two
detectors get, and the run reports it as its own row rather than blending it in.

## C. What this amendment does not touch

The estimand (§1), the eligible population (§2), the three-way handling of
`not_applicable` (§3), the exclusions (§5), the metrics (§6), the definition of
a positive, negative or inconclusive result (§7), the one-run rule (§8), and
the exclusion of the persistence variant (§10) are all unchanged.
