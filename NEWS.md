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

* `evidence_stream()` now carries `test`, `weighted` and a user-fixed `tau2`
  into every snapshot, alongside `method`. Previously only `method` survived,
  so a model fitted with `test = "knha"` was scored with the default z test --
  a p-value differing by a factor of 178 on `metadat::dat.bcg`, in a detector
  that decides on p-values. Snapshots now match a direct refit exactly for
  each of those options.

* `evidence_stream()` refuses, with an explanation, models it cannot honour:
  meta-regressions, whose `beta` is a vector of coefficients rather than a
  pooled effect, and models with custom per-study `weights`, which cannot
  follow the subsetting each snapshot does. Both were previously accepted and
  quietly flattened into a plain pooled analysis.

* `backtest()` validates its arguments. A negative `horizon`, a zero or
  negative `window`, a `min_k` below 2, `cuts` containing `NA` and an empty
  `methods` were all accepted or failed with internal R messages such as
  "missing value where TRUE/FALSE needed". Each now says what is wrong.

* `cuts` are sorted and de-duplicated. A repeated cut used to enter the
  results once per repetition, so passing the same cut three times raised `n`
  from 3 to 5 in `calibration()` — the denominator of every rate, inflated
  without a warning.

* The five detectors validate their scalar arguments too — thresholds,
  significance levels, replicate and permutation counts, sample sizes and
  seeds. The engine had been swept and they had not, and the same defect was
  in every one of them. Two were silent, which is the worse failure:
  `ottawa(alpha = NA)` returned a verdict of `out_of_date`, because both
  significance comparisons became `NA` and the code read that as a change, and
  `sufficiency(alpha_stability = NA)` returned `current`. Neither warned.

* `seed` is validated in `simulation()`, `sufficiency()` and `backtest()`, and
  is the reason the sentence above says *scalar arguments* rather than *all
  arguments*. `set.seed()` truncates towards zero, so `seed = 1.5` was accepted
  in silence and produced the identical stream to `seed = 1` — two values a
  reader would record as different runs, giving the same numbers. A vector was
  accepted too, silently using its first element. `NULL` remains valid and
  still means an unseeded run. The check sits at each entry point rather than
  only at `with_preserved_seed()`, which the detectors reach after their early
  returns: a `not_applicable` verdict must not swallow a malformed seed.

  The accepted range is `[-.Machine$integer.max, .Machine$integer.max]`, which
  is what `set.seed()` takes — note the lower end, since R reserves
  `-2147483648` for `NA_integer_` and rejects it as a seed. A whole number
  past that range used to satisfy the check and then fail inside `set.seed()`
  with a coercion warning, so the promise was wider than the code and the
  failure arrived without an explanation.

* `check_currency()` and `simulation()` require the new-evidence object to be
  internally consistent: `yi` and `vi` of the same length, and `k` equal to
  that length. `k` decided whether there was anything to assess and `yi`/`vi`
  were what got fitted, with nothing tying them together, so the object could
  say one thing and carry another — and it failed silently in both directions.
  `k = 1` with no studies returned a verdict of `current` from an updated model
  refitted on the prior evidence alone, walking straight past the guard written
  to prevent exactly that ("absence of new evidence is not evidence of
  currency"); `k = 0` with real studies discarded them; `k = 99` with one study
  left the verdict unchanged but reported `detail$k_new` as 99. Mismatched
  lengths used to surface as metafor's "Length of 'yi' and 'vi' (or 'sei') are
  not the same", a message about another package's arguments.
  [window_between()] has always produced `k` as the size of the same subset it
  takes `yi` and `vi` from, so the contract is what the canonical source
  already builds.

* The same objects must also be statistically usable, not merely well shaped:
  `yi` finite, `vi` finite and strictly positive. Shape alone was not enough,
  and the gap reached the same wrong answer by another route — `yi = NA`,
  `vi = NA`, `vi = Inf`, `vi = 0` and `vi = -1` each returned `current` from
  both `rcma()` and `ottawa()`, because metafor drops or ignores such a study,
  the updated model comes back identical to the prior one, and every ratio is
  1. Zero-length vectors stay valid, being vacuously finite and positive, so
  the empty case still returns its own class.

* `check_currency()` refuses an empty `methods`, as `backtest()` already did.
  It used to return a `staleness_check` holding zero verdicts: an object that
  answers no question.

* `barrowman()` refuses impossible sample sizes. `n_prev = 0` drove the
  required sample size to zero, so the participant ratio was `Inf` and the
  verdict was `out_of_date`; a negative size made the ratio negative, which
  can never exceed 1, so it always read `current` — a failure biased in one
  direction. `n_new = 0` remains valid and means what it says: no new
  participants, so the ratio is 0 and the review reads as current. An `NA` or
  infinite size is still a fact about the evidence, not the call, and still
  yields `"not_applicable"` with its reason.

## Units, scales and unknowns

* `evidence_stream()` requires publication **years** and refuses a `Date`.
  `date` had been checked for length and for missing values but never for
  type, so a `Date` went through `as.numeric()` into days since 1970 — and
  every window here is denominated in years. Ten studies three days apart
  produced **21 cuts numbered 11329, 11330, …** under `cuts = "yearly"`, and
  `horizon` and `window` silently became days as well. The documentation had
  promised `Date` support outright. Convert with
  `as.numeric(format(date, "%Y"))`; a factor or character vector is refused
  too, rather than failing later inside a comparison.

* `rcma()` and `ottawa()` require both models to be on the same `measure`.
  Each read `$measure` off `prev` to choose its branch and neither checked
  `new_ma`, so a risk ratio compared against a mean difference returned a
  signal of 2.98 and a verdict of `out_of_date` — a number with no meaning,
  presented like any other. `check_currency()` never produced this because it
  refits both models itself, but both detectors are exported.

* `ottawa()` refuses an `NA` in `qualitative`. `nzchar(NA_character_)` is
  `TRUE`, so an analyst recording "unknown" for one of the four qualitative
  signals received `"out_of_date"`. Unknown is not present. Pass `""` for a
  signal that was checked and found absent.

* `truth_shift()`, `truth_surprise()` and `truth_conclusion()` require scalar
  inputs, as `?truth` has always said they take. `truth_shift()` returned one
  logical per element and `truth_conclusion()` died inside R's own coercion —
  the documentation was stronger than the code, and the `\value` section
  added a commit earlier had just restated the claim without enforcing it.

## Validation against real reviews, with the criterion matched to the detector

* Two further historical cases, both scored against what Cochrane went on to
  publish, read from the raw abstracts of the successive versions rather than
  from a summarising fetch.

  `metadat::dat.damico2009` is the evidence behind CD000022.pub3 (2009). The
  review was updated as .pub4 in 2021 and its conclusion did not change. No
  detector fires on any of its 12 cuts, and the pooled effect moves 7% across
  the whole series. Twelve years and a full update later, silence was right.

  `metadat::dat.lee2004` is CD003281.pub2 (2004). Its conclusion never changed
  either — P6 still works in .pub3 and .pub5 — but the pooled odds ratio moves
  from 0.38 to 0.60 inside the window the dataset covers, a 55% shift past
  `rcma()`'s threshold, and `rcma()` fires. `ottawa()` stays quiet on the same
  evidence, because significance never changed. **Both are correct.**

* That pair carries a lesson worth stating, because it invalidates the obvious
  way to score these detectors: **a criterion only validates the detector whose
  object it measures.** Cochrane declares whether CONCLUSIONS changed, which
  judges `ottawa()` and `truth_conclusion()`. It cannot judge `rcma()`, which
  measures magnitude — scoring `dat.lee2004` against the conclusion would have
  recorded a correct firing as a false positive. The earlier `dat.lau1992` and
  `dat.li2007` cases did not expose this because magnitude and conclusion moved
  together in both.

## The claim that nothing like this existed

* `inst/cran-search/` turns that claim into evidence: a runnable search, the
  terms and fields it uses, a dated snapshot, and a hand-written verdict for
  every hit. Searching the metadata of all 24,734 CRAN packages on 2026-08-09
  finds no implementation of any of the five detectors — every hit on a method
  name is a false positive, and each is named rather than counted away, from
  the University of **Ottawa**'s research group to a package that warns about
  out-of-date **R packages**.

* The claim is now stated as *we did not find one* rather than *none exists*,
  because the search reads package metadata and not source code. The limit is
  demonstrated rather than asserted: `metafor` matches neither "cumulative
  meta-analysis" nor "fail-safe" in its own metadata, and exports `cumul()`
  and `fsn()` regardless.

* The search also named components the documentation had not: `fsn` and `meta`
  alongside `metafor` and `RTSA`, and `metagear`, which screens literature —
  the half of the problem this package explicitly does not attempt.

* A hit with no recorded verdict is reported as `UNADJUDICATED` rather than
  absorbed into a total, so a package arriving later is visible and the claim
  can be withdrawn if one turns up.

## Applicability across historical reviews

* `inst/applicability/` sweeps every dataset in `metadat`, keeps the 17 that
  can carry a backtest, and runs all five detectors over each. **168 of those
  reviews' 185 yearly cuts (91%) had an already-significant prior
  meta-analysis**, and in 11 of the 17 every cut did. `barrowman()` and
  `simulation()` require a non-significant prior, so they can answer in only 4
  and 5 of the 17 — not wrongly, but unable to be asked at all. The published
  comparison of these methods could not see this: its 80 reviews were selected
  for having a non-significant pooled result, which is precisely the case
  where those two are applicable. `sufficiency()` fails the other way,
  answering in all 17 with a mean sensitivity of 0.041 and zero in 11 of them.

  Sensitivity and specificity there are scored against `truth_shift()`, this
  package's own definition, so the sweep compares the five methods against
  each other over real evidence. It does not validate them against published
  outcomes; that claim belongs to the four cases in the test suite and still
  numbers four.

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

* `rcma()` and `ottawa()` are two criteria, not one counted twice. Both compare
  an updated quantity against a prior one at thresholds of 0.5 and 1.5, but
  `rcma()` takes the ratio of the pooled effects and `ottawa()` the ratio of
  the relative risk reductions, so on ratio measures they disagree in both
  directions; on difference measures they coincide, because the source defers
  to the rCMA rule there. Earlier documentation declared the first to be
  contained in the second — see the correction above — and both help pages now
  say otherwise, so that inter-method agreement is computed on what the
  detectors actually do.
