# staleness 0.1.0

First public release.

## What the package does

* `check_currency()` applies five published detectors — `rcma()`, `ottawa()`,
  `barrowman()`, `sufficiency()` and `simulation()` — to decide whether an
  existing meta-analysis is still current given the evidence published since.

* `backtest()` replays those detectors over the history of a body of evidence,
  refitting at every cut so no detector ever sees the future, and
  `calibration()` and `lead_time()` score them against three independent
  definitions of ground truth (`truth_shift()`, `truth_surprise()`,
  `truth_conclusion()`).

* Detector–truth pairs that are circular by construction are named in
  `CONTAMINATED_PAIRS` and flagged on every row of `calibration()`, rather
  than being dropped where a reader might not notice the absence.

* A detector that cannot answer returns `not_applicable`, never `current`.
  Scoring inapplicability as a correct call would flatter the detectors that
  decline most often.

## Metric semantics

* `lead_time()` no longer discards a detector's firing because the truth of
  that firing's own cut could not be determined. It is the only metric that
  relates different rows — the event in one, the firing that preceded it in
  another — and the two sides do not carry the same requirement. An event
  needs a known truth, because unknown is not true; a firing is an observed
  fact whose role is to precede a later event. Reusing `calibration()`'s
  row-wise eligibility filter on both sides deleted real early warnings, and
  only ever against the detector: dropping a firing can lengthen or erase a
  lead, never shorten it.

## Fidelity

* `ottawa()`'s effect criterion now compares relative risk **reductions**,
  `(1 - RR_new) / (1 - RR_prev)`, which is what the method's published
  application specifies for ratio measures; mean differences still defer to
  the `rcma()` rule. The package previously compared the effects themselves.
  Of the ten reviews with the largest Ottawa indicator in Pattanittum et al.
  (2012), Appendix S3, all ten fire under the corrected definition and none
  fired under the old one. Those ten are now a test fixture.

* Consequently `ottawa()` no longer contains `rcma()`. The two use different
  quantities on ratio measures, so an `rcma` firing is not a subset of the
  `ottawa` firings and the containment declared in earlier documentation was
  wrong. On difference measures they still coincide.

* The corrected criterion is unstable where the method is meant to be used:
  `1 - RR_prev` approaches zero on a null meta-analysis, and on evidence with
  no change at all the effect signal fires on 64% of samples under a null
  effect versus 0% under a real one. Measured, pinned by a test, and
  documented in `?ottawa` and `vignette("methods")` rather than smoothed over.

* `simulation()` now follows the published procedure step for step
  (Pattanittum et al. 2012, Appendix S1): the effect of the new study is drawn
  from a **t** distribution rather than a normal one; **one** study is
  simulated carrying the combined precision of the recent studies, rather than
  one per recent study; and the threshold is strict, matching the source's
  "Power >80%". The one deviation that cannot be removed — simulating at the
  level of effects rather than of participants, because the package never sees
  2x2 tables — is now declared in `?simulation` instead of going unmentioned.

* `ottawa()` now reports `detail$effect_unstable` and `detail$rrr_prev`. The
  effect criterion divides by `1 - RR_prev`, which approaches zero on exactly
  the null reviews the method targets, so the ratio can be arbitrarily large.
  That behaviour is the published method's and is deliberately left alone —
  repairing it would mean implementing something else — but a caller used to
  receive a ratio of -19 with an empty `reason` and no indication that it came
  from dividing by -0.005. Now the flag and the denominator travel with the
  verdict, which is unchanged.

## Reproducibility

* `inst/calibration/calibration.R` regenerates every calibration figure quoted
  in `?sufficiency`, `vignette("methods")` and the paper, from the seeds the
  original measurements used, and reconstructs both of the statistics that
  were replaced so the before-and-after comparison can be re-derived. Two
  published figures did not survive that check and were corrected: the
  20-small-then-10-large schedule reads 20/300 (6.7%), not 19/300, and the
  schedule itself is 50:1.

## Notes on two published methods

* `sufficiency()` tests stability with a change-point statistic
  (`max_m |Z_m|`) under an order-permutation null, not with the ordinary least
  squares slope of the cumulative series that the source describes. A *t*-test
  on a cumulative mean has no valid null distribution — it is autocorrelated
  by construction and convergent by the law of large numbers — and fired on
  209 of 300 samples containing no change at all. The published slope is still
  computed and reported in `detail$slope`. `vignette("methods")` gives the
  measured calibration of the replacement across nine variance regimes, and
  where it degrades.

* `rcma()`'s rule is contained verbatim in `ottawa()`: same effect ratio, same
  thresholds. Whenever `rcma()` fires, `ottawa()` fires, by arithmetic rather
  than by agreement. Declared in both help pages so that nobody recomputes
  inter-method agreement treating one of the ten detector pairs as data.
