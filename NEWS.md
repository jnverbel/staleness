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
