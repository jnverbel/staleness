#' Rosenthal's fail-safe N
#'
#' Included because the sufficiency and stability method specifies it. Note that
#' Rosenthal's fail-safe N has been discredited as a measure of publication bias
#' since Becker (2005). It is implemented faithfully so that the backtesting
#' engine can settle the question with data rather than opinion.
#'
#' @param yi,vi Effect sizes and their variances.
#' @param z_crit One-sided critical value, 1.645 for alpha = 0.05.
#' @return The fail-safe N.
#' @keywords internal
failsafe_n <- function(yi, vi, z_crit = 1.645) {
  # Empty input would otherwise evaluate to exactly 0, and 0 is a meaningful
  # fail-safe N: it flows on into index = 0 / (5 * 0 + 10) and reports
  # `sufficient = FALSE`. "No studies" and "not enough unpublished nulls to
  # matter" are different facts, and only one of them is an answer.
  if (!length(yi)) return(NA_real_)
  z <- yi / sqrt(vi)
  (sum(z)^2) / (z_crit^2) - length(z)
}

#' Cumulative fixed-effect estimate after each study
#'
#' @param yi,vi Effect sizes and their variances, in the order studies are to
#'   be accumulated.
#' @return Numeric vector of length `length(yi)`.
#' @keywords internal
cumulative_effect <- function(yi, vi) {
  w <- 1 / vi
  cumsum(w * yi) / cumsum(w)
}

#' Slope of a cumulative series against accumulated information
#'
#' The x-axis is accumulated Fisher information, `cumsum(1 / vi)`, which is
#' what the primary source means by *"versus information increment"*
#' (Pattanittum et al. 2012, Table 1) and what trial-sequential-analysis
#' methods use. An earlier implementation regressed on the study index
#' `1, 2, 3, ...` instead; the two coincide only when every study carries the
#' same variance.
#'
#' The first element of a cumulative series is just the first study, carrying
#' no information about accumulation, so it is dropped before the slope is
#' taken. Computed in closed form rather than through [stats::lm()]: the same
#' number, but cheap.
#'
#' This slope is **reported** in `detail$slope` for fidelity to the published
#' rule, so that anyone who wants to apply `|slope| > 0` literally can read it
#' straight off the verdict. It is not what decides stability — see
#' [stability_shift_z()] and the roxygen block of
#' [sufficiency_changepoint()] for why a slope over a cumulative series turned
#' out to be an unusable instrument.
#'
#' @param cum_theta Numeric vector, as returned by [cumulative_effect()].
#' @param info Numeric vector of the same length, the accumulated information
#'   `cumsum(1 / vi)` in the same order.
#' @return The ordinary-least-squares slope, or `NA_real_` if the series is
#'   too short to have one.
#' @keywords internal
cum_drift_slope <- function(cum_theta, info) {
  y <- cum_theta[-1]
  if (length(y) < 2) return(NA_real_)
  x  <- info[-1]
  xc <- x - mean(x)
  denom <- sum(xc^2)
  if (denom <= 0) return(NA_real_)
  sum(xc * (y - mean(y))) / denom
}

#' Largest standardised movement still left in a cumulative series
#'
#' The stability statistic actually tested by [sufficiency_changepoint()]. For
#' every split point `m` it compares the fixed-effect pooled estimate of the
#' studies accumulated up to `m` against the pooled estimate of everything
#' after `m`,
#' standardised by the standard error of that difference:
#'
#' \deqn{Z_m = \frac{\hat\theta_{1:m} - \hat\theta_{(m+1):k}}{\sqrt{1/I_m +
#'   1/(I_k - I_m)}}}{Z_m = (theta(1..m) - theta(m+1..k)) /
#'   sqrt(1/I_m + 1/(I_k - I_m))}
#'
#' where \eqn{I_m = \sum_{i \le m} 1/v_i}{I_m = sum of 1/v_i over the first m
#' studies} is the information accumulated by step `m`. The statistic is
#' `max |Z_m|` over `m = 1, ..., k - 1`.
#'
#' It measures the same thing the published slope is reaching for. The
#' cumulative estimate's remaining movement after step `m` is exactly
#' `(1 - I_m / I_k)` times `theta_after - theta_before`, so `Z_m` is that
#' movement expressed in units of its own standard error, and the maximum over
#' `m` is the largest step the cumulative curve has left in it.
#'
#' @section How far the pivotality goes, and where it stops:
#' Two things are exactly true. **Scale pivotality**: rescaling the effects and
#' their variances together — `yi` to `theta + a * (yi - theta)` and `vi` to
#' `a^2 * vi` — leaves the statistic bit-identical, since it is a ratio of a
#' movement to its own standard error and both halves move by `a`. Note what
#' this does *not* say: multiplying `vi` alone by `c` while holding `yi` fixed
#' divides every `Z_m` by `sqrt(c)` (measured: 1.2247 at `vi = 1`, 0.6124 at
#' `vi = 4`, for `yi = c(1, 2, 3)`). That weaker statement is what an earlier
#' draft of this section claimed, and it is false. The permutation *p*-value is
#' unaffected either way, because every permuted statistic rescales by the same
#' factor. And **marginal pivotality**: under no drift each individual `Z_m` is
#' standard normal whatever the variance schedule is (measured standard
#' deviations 0.998–1.002 across schedules spanning 650:1).
#'
#' The stronger claim — that the statistic therefore carries no imprint of the
#' variance schedule at all — does **not** follow, and an earlier draft of this
#' documentation asserted it wrongly. `max_m |Z_m|` is the discrete supremum of
#' a standardised Brownian bridge sampled at the information fractions
#' `I_m / I_k`, and those fractions are themselves a function of the schedule.
#' Measured null quantiles at `k = 30` over 60,000 draws:
#'
#' \preformatted{
#'   schedule (variance)          ratio    median   95th pct
#'   flat                           1:1     1.817      2.881
#'   V-shaped,      big-small-big      16:1     1.923      2.950
#'   Lambda-shaped, small-big-small    16:1     1.599      2.732
#'   geometric decay                  646:1     1.930      2.956
#' }
#'
#' The V and Lambda nulls differ with a Kolmogorov-Smirnov `D` of 0.25, and the
#' flat-schedule 95th percentile delivers 6.1% under V and 3.4% under Lambda.
#'
#' So the justification for the order-permutation null in
#' [sufficiency_changepoint()] is that the dependence is **weak**, not that it
#' is absent — quantiles move by a few percent where the old slope statistic
#' moved by orders of magnitude
#' (0% to 42% false alarms across the same schedules). Since the null is
#' rebuilt from each dataset's own studies rather than read off a reference
#' distribution, most of that residual shows up only where the schedule is
#' smooth and strongly monotone. See the calibration section of
#' [sufficiency_changepoint()] for what it costs in practice, and where.
#'
#' @param yi,vi Effect sizes and their variances, in the order studies are to
#'   be accumulated.
#' @return A non-negative number, or `NA_real_` for a series too short to
#'   split.
#' @keywords internal
stability_shift_z <- function(yi, vi) {
  z <- split_z(yi, vi)
  # A split can be undefined without the series being: a study with infinite
  # variance at the end contributes no weight, so the "after" block of the last
  # split holds no information and its Z is 0/0. Drop those splits rather than
  # letting one NaN swallow the maximum; if none survive there is nothing to
  # test and the caller resolves it as stable.
  z <- z[is.finite(z)]
  if (!length(z)) return(NA_real_)
  max(abs(z))
}

#' Standardised split differences, one per possible changepoint
#'
#' @inheritParams stability_shift_z
#' @return Numeric vector of length `length(yi) - 1`, possibly containing
#'   non-finite entries for splits that carry no information.
#' @keywords internal
split_z <- function(yi, vi) {
  k <- length(yi)
  if (k < 2) return(numeric(0))
  w   <- 1 / vi
  sw  <- cumsum(w)
  swy <- cumsum(w * yi)
  m   <- seq_len(k - 1L)
  before <- sw[m]
  after  <- sw[k] - before
  (swy[m] / before - (swy[k] - swy[m]) / after) /
    sqrt(1 / before + 1 / after)
}

#' Split point at which [stability_shift_z()] attains its maximum
#'
#' Reported in `detail$split` purely as a diagnostic: it says *where* the
#' cumulative estimate was still moving most, which is the difference between
#' a review whose conclusion shifted years ago and one that is moving now.
#'
#' @inheritParams stability_shift_z
#' @return An integer in `1:(k - 1)`, or `NA_integer_`.
#' @keywords internal
stability_shift_at <- function(yi, vi) {
  z <- split_z(yi, vi)
  if (!any(is.finite(z))) return(NA_integer_)
  z[!is.finite(z)] <- NA_real_       # which.max() skips NA; -Inf would win
  which.max(abs(z))
}

#' Sufficiency and stability, with a change-point test for stability
#'
#' The name carries the substitution rather than hiding it. This is the
#' sufficiency-and-stability method, and its sufficiency half is the published
#' one; its stability half is **not** the published statistic. The source tests
#' stability with the slope of a least-squares line through the cumulative
#' effects, and that slope proved to be an unusable instrument — the section
#' *How stability is tested, and why not as published* below gives the three
#' measurements that established it. What runs instead is a change-point
#' statistic: the largest standardised difference between the studies before
#' and after any split of the series. It asks the method's own question and it
#' is calibrated; it is still a substitution, so the function is not called
#' `sufficiency()`.
#'
#' The published slope is still computed and returned in `detail$slope`, as a
#' diagnostic and for anyone who wants to apply the source's rule literally. It
#' decides nothing.
#'
#' Sufficiency is the fail-safe N scaled by `5k + 10`; a review is sufficient
#' when this index exceeds 1, Rosenthal's own rule of thumb for a pooled effect
#' being robust to unpublished null studies. Stability asks whether the
#' cumulative pooled effect has stopped moving; a review is stable when the
#' largest standardised movement left in its cumulative series is no larger
#' than the same studies would produce in a random order.
#'
#' The two indicators deliberately read from different bodies of evidence.
#' Per the primary source that operationalises this method as a two-snapshot
#' comparison (Pattanittum et al. 2012, Table 1), sufficiency — both the
#' fail-safe N and the `k` in `5k + 10` — is computed on the meta-analysis
#' **as previously published** (`prev`), while stability is computed on the
#' cumulative series of the **updated** meta-analysis (`new_ma`), which is
#' where new studies show up as drift or lack of it. Confirmed against a
#' second, independent secondary source; the original method paper (Mullen,
#' Muellerleile & Bryant 2001) was not reachable in full text, and does not
#' itself define a two-snapshot variant to compare against — see the design
#' doc (section 4.4) for the full trail of evidence.
#'
#' Per the same source, an out-of-date review is one that is BOTH sufficient
#' and unstable: enough evidence had already accumulated, as of the prior
#' review, to be confident the effect is real, but the pooled estimate
#' (including what came after) is still drifting, so its magnitude is not yet
#' settled. Insufficient evidence alone is never grounds for "out of date" —
#' the opposite combination from what a first reading of the secondary source
#' suggests, and also confirmed in the design doc.
#'
#' @section How stability is tested, and why not as published:
#' This is the largest declared departure in the package, and it is a
#' departure from the *statistic*, not only from the decision rule. It was
#' arrived at by measurement, over three implementations; the numbers below
#' are all reproducible from `vignette("methods")`.
#'
#' **What the source says.** Pattanittum et al. (2012, Table 1) states the
#' instability criterion as an *"absolute slope of the linear regression
#' \[fitted across the cumulative treatment effects versus information
#' increment\] >0"*. Taken literally the rule is degenerate: on continuous
#' data that slope is never exactly zero, so every review with `index > 1`
#' would be flagged and the indicator would carry no information. Some
#' significance rule has to stand in for it.
#'
#' **Attempt 1, the OLS t-test: invalid.** A cumulative mean is
#' near-perfectly autocorrelated by construction and converges on the pooled
#' effect by the law of large numbers, so the `lm()` t-test detects
#' *convergence* and reports *instability*. Over 300 samples of genuinely
#' unchanging evidence (20 prior plus 10 new studies from one distribution)
#' it returned `out_of_date` 209 times, where `rcma()` and `ottawa()`
#' returned it none.
#'
#' **Attempt 2, permuting study order but keeping the slope: calibrated on
#' average, blind where it matters.** Comparing the observed slope against
#' the slopes of the same studies reshuffled brought the false-alarm rate on
#' that experiment to 16/300, nominal. But the slope of a cumulative series
#' is dominated by its first few points, where little information has
#' accumulated and the series swings hardest, so the permutation null is
#' wide and a *late* change cannot clear it. Measured: with 20 prior studies
#' at RR 0.5 and 10 new ones at RR 0.30, power was **1 in 200** — and it fell
#' to 0 as the new studies got further away, because a larger shift makes the
#' permuted series swing harder too. On a 20-study series that sits at RR 0.5
#' and then drops to RR 0.05, it fired when the drop came after study
#' 2, 5, 8 or 10 and was silent when it came after study 12, 15 or 18. Worse,
#' the permutation null itself is not valid when study variances change
#' systematically over calendar time: with early small trials and later large
#' ones and *no drift at all*, it fired on 83/300 samples (28%), and on a
#' schedule of 20 small trials followed by 10 large ones, on 127/300 (42%).
#' Study order is exchangeable under no change only if the variances are
#' exchangeable too, and in a real evidence stream they are not.
#'
#' Fitting the regression against accumulated information rather than the
#' study index — which is what the source's *"versus information increment"*
#' actually specifies, and what [cum_drift_slope()] now reports — is a real
#' fidelity fix but **not** a fix for any of this. When variances are equal,
#' `cumsum(1 / vi)` is an affine function of the study index, so every slope
#' is rescaled by the same constant and the permutation p-value is *identical
#' to the last bit*. It was measured: E1, E2 and E3 above were unchanged, and
#' the 28% false-alarm rate moved only to 25%.
#'
#' **What runs instead.** The statistic is replaced with the largest
#' standardised movement the cumulative series has left in it
#' ([stability_shift_z()]): for each split point `m`, the fixed-effect pooled
#' estimate of the first `m` studies is compared with the pooled estimate of
#' the rest, divided by the standard error of that difference, and the
#' statistic is the maximum absolute value over `m`. This still asks the
#' method's own question — *has the pooled estimate stopped moving?* — and it
#' is monotone in the same quantity the published slope is a proxy for, since
#' the movement remaining in the cumulative curve after step `m` is exactly
#' `(1 - I_m / I_k)` times the difference it standardises. What it adds is
#' that its dependence on the variance schedule is **weak instead of
#' structural**: it is exactly scale-pivotal and each `Z_m` is marginally
#' standard normal under any schedule, so the schedule no longer drives the
#' statistic the way it drove the slope. It does not vanish from it entirely —
#' see the pivotality section of [stability_shift_z()] for the measured
#' residual and the calibration table below for what it costs. The permutation
#' machinery is otherwise unchanged — `n_perm` draws, the two-sided
#' `(1 + count) / (n_perm + 1)` estimator so the p-value is never zero, a fixed
#' seed, the caller's random stream preserved, and a cumulative series that is
#' constant to floating-point rounding short-circuited as maximally stable
#' before anything is computed.
#'
#' Measured against the same four experiments: false alarms 15/300 (5.0%);
#' power against 10 new studies at RR 0.40/0.30/0.15/0.02, 200/200 at every
#' level; the shift-position scan fires at all seven positions including
#' 12, 15 and 18; and the heteroscedastic no-drift false-alarm rate is 16/300
#' (5.3%) with variance falling over time, 20/300 (6.7%) with it rising, and
#' 20/300 (6.7%) on the 20-small-then-10-large schedule that produced 127/300
#' (42%) before. Every figure in this section is regenerated by
#' \code{Rscript} on the script at
#' \code{system.file("calibration", "calibration.R", package = "staleness")},
#' which also reconstructs the replaced statistic so the before-and-after is
#' checkable rather than asserted.
#'
#' @section Calibration, and the variance schedules that break it:
#' An earlier draft of this documentation reported that calibration "stays
#' between 4.0% and 6.7% across nine variance and heterogeneity regimes". That
#' was true of those nine regimes and false as a general claim: all nine were
#' flat, mildly linear or single-step. Measured over 1000 no-drift samples at
#' `k = 30` (2000 for the worst case), against the actual permutation test:
#'
#' \preformatted{
#'   variance schedule                              ratio   false alarms
#'   flat                                             1:1       5.3%
#'   linear decay (the E4 schedule)                  16:1       5.7%
#'   dat.bcg, real variances in year order           60:1       5.0%
#'   dat.egger2001, real variances                   90:1       3.0%
#'   L-shaped, small-big-small                       16:1       2.4%  <- conservative
#'   step, 20 small then 10 large trials             50:1       8.1%
#'   V-shaped, big-small-big                         16:1       8.4%
#'   dat.bcg, the same variances sorted decreasing   60:1       9.4%
#'   geometric decay, vi = 0.5 * 0.845^j            114:1       9.3%
#'   geometric decay, vi = 0.5 * 0.80^j             650:1      11.1%  <- worst measured
#'   geometric growth, the same reversed            650:1      10.9%
#' }
#'
#' The pattern is not the direction of the trend — growth and decay are
#' equally bad — but its **smoothness, monotonicity and range**. Irregular
#' real-world schedules are nominal; the same `dat.bcg` variances rearranged
#' into sorted order go from 5.0% to 9.4%, which is the cleanest demonstration
#' that it is the arrangement and not the numbers. A smooth 650:1 monotone
#' precision ramp is the measured worst case at about **11%**, roughly twice
#' nominal. That is a disclosed limit, not a claim of nominality; it is pinned
#' by a test in `test-invariants.R` so it cannot silently get worse.
#'
#' `detail$slope` still reports the literal published quantity — the
#' least-squares slope of the cumulative effects against accumulated
#' information — so a reader who wants the source's own rule can apply it.
#' It no longer decides anything.
#'
#' @section Known limits of the stability test:
#' Four, all measured, none of them hidden. Three of the four were understated
#' in the first draft of this documentation and are restated here with the
#' numbers that falsified the original wording.
#'
#' **It is conservative well past `k = 5`.** The false-alarm rate on unchanging
#' evidence is 0.9% at `k = 5`, 1.7% at `k = 6`, 3.4% at `k = 8`, 3.3% at
#' `k = 10` and 5.5% at `k = 12` — still half-nominal at twice the minimum, and
#' nominal only from about `k = 12`. (The first draft said "conservative at
#' `k = 5`" and blamed the 120 distinct orderings, which cannot explain `k = 8`
#' and its 40,320.) The real mechanism is the same one that produces the next
#' limitation, with its sign reversed: the maximum is often driven by the split
#' that isolates the single most extreme study, and a given study lands at an
#' end in `2 / k` of all orderings — 40% at `k = 5`, 25% at `k = 8`, 6.7% at
#' `k = 30`. At small `k` a large fraction of permutations therefore reproduce
#' the observed extreme, the null is diffuse, and the test under-fires. It
#' under-fires rather than over-fires, which is the safe direction, but a
#' `current` verdict from a short update is weak evidence of stability.
#'
#' **A single discordant study arriving LAST can carry the statistic.** The
#' statistic is a maximum over split points, and the split that isolates one
#' study is one of them. Measured by `inst/persistence/persistence.R` over 400
#' replicates per cell, against a baseline with no outlier:
#'
#' \preformatted{
#'    k     no outlier   one 5-SE study, last   the same study, mid-series
#'   20         0.050            0.035                    0.005
#'   30         0.052            0.207                    0.003
#'   40         0.043            0.395                    0.000
#' }
#'
#' The effect is real, grows with `k`, and is confined to the tail: a study in
#' mid-series does essentially nothing, because the split isolating it leaves
#' large blocks on both sides. `detail$split` tells the two apart, since it
#' reports where the maximum was attained.
#'
#' An earlier version of this block reported far higher rates -- 0.995 at
#' `k` = 40, and "almost always from about `k` = 25" -- and also a superseded
#' claim before that. **Neither is reproducible.** An independent
#' reimplementation that agrees with `stability_shift_z()` exactly on 200
#' random inputs, and that reproduces the no-outlier control row above, returns
#' 0.44-0.48 at `k` = 40 under every reading of "five standard errors" tried,
#' and saturates near 0.55 even at twenty. The old table had no generating
#' script anywhere in the package, which is exactly how it survived two
#' revisions: nothing could contradict it. The figures above have one, and it
#' is named.
#'
#' Requiring a shift to PERSIST removes most of the effect. Demanding that both
#' sides of a split hold at least `r` studies gives, averaged over
#' `k` = 20, 30, 40:
#'
#' \preformatted{
#'    r        drift   outlier(last)   outlier(mid)   no change
#'    1 (this) 0.802       0.212           0.003         0.048
#'    2        0.815       0.226           0.003         0.051
#'    3        0.825       0.147           0.006         0.054
#'    5        0.845       0.083           0.013         0.052
#' }
#'
#' `r` = 5 is better on every axis measured: more power against a real late
#' shift, less than half the single-outlier rate, and calibration untouched.
#' It is **not** what this function runs. Changing the shipped statistic on one
#' experiment over four regimes and a single variance schedule would repeat the
#' mistake documented two paragraphs above, so it is recorded as the next test
#' rather than adopted.
#'
#' **A smooth, strongly monotone variance schedule runs anti-conservative**, up
#' to about 11% against a nominal 5%. See the calibration section above for the
#' full table and for the regimes where it does not happen, which include every
#' real variance schedule tested.
#'
#' **Power falls away under strong between-study heterogeneity.** Against 20
#' prior studies at RR 0.5 and 10 new ones at RR 0.30, power is 200/200 at
#' `tau = 0` and `tau = 0.10`, 199/200 at `tau = 0.20`, and 160/200 at
#' `tau = 0.35`. The permutation null is built from the same heterogeneous
#' studies, so calibration holds throughout — it is power, not the false-alarm
#' rate, that degrades.
#'
#' @section Monte Carlo error, reported and never acted on:
#' The permutation p-value `p_stability` is estimated from a finite
#' simulation, so it moves with the seed.
#' Three fields report that movement and one interprets it, and none of them
#' touches the verdict:
#'
#' * `detail$mc_se` -- the Monte Carlo standard error of `p_stability`.
#' * `detail$mc_lo`, `detail$mc_hi` -- a 95% Wilson interval for it. Wilson
#'   rather than Wald, because the Wald interval has width exactly zero at 0
#'   and 1, which is where a reader most needs to be told the estimate is not
#'   certain.
#' * `detail$near_threshold` -- `TRUE` when `alpha_stability` falls inside that
#'   interval, meaning a rerun under a different seed could plausibly have
#'   returned the other verdict. `FALSE` means the simulation is resolved on
#'   this question; it does not mean the detector is right.
#'
#' This is Monte Carlo error alone -- the variability of the estimate around
#' what infinitely many draws would give. It says nothing about sampling error
#' in the underlying studies, which is much larger.
#'
#' The interval is built on the exceedance count out of `n_perm` draws and then
#' carried through the same `(1 + x) / (n_perm + 1)` map the reported p-value
#' uses, so the bounds sit on exactly the scale of the number they bracket. The
#' degenerate branches -- a cumulative series constant to rounding, or a series
#' where no split carries information -- draw no permutations at all, so their
#' `mc_*` fields are `NA` and `near_threshold` is `FALSE`: those p-values are
#' determined, not estimated.
#'
#' @param prev An `rma.uni` object, the meta-analysis as previously published.
#' @param new_ma An `rma.uni` object refitted with the new evidence included.
#' @param min_k Minimum number of studies in `new_ma`. Below this the
#'   stability test is meaningless: an order-permutation null needs enough
#'   orderings to have a distribution.
#'
#'   It gates the **updated** evidence only. The two halves of this method
#'   read different objects — sufficiency from `prev`, stability from
#'   `new_ma` — so nothing here constrains the size of `prev`, and a
#'   two-study prior review can drive the sufficiency half on its own. That
#'   is deliberate rather than overlooked: the fail-safe N is defined for any
#'   number of studies, and the source sets no floor. Read `detail$k` if the
#'   size of the prior review matters to you; it is reported for that reason.
#' @param alpha_stability Cutoff for the stability permutation p-value: the
#'   review is unstable when `p_stability < alpha_stability`.
#' @param n_perm Number of order permutations used to build the null
#'   distribution of the stability statistic.
#' @param seed Integer seed for the permutation draw. Fixed by default so
#'   verdicts are reproducible; the caller's own random stream is saved and
#'   restored around the draw, so calling this detector never changes it.
#' @return A `staleness_verdict`. Its `detail` carries `index`, `sufficient`,
#'   `stable`, `slope` (the published slope, reported only), `z_shift` and
#'   `split` (the statistic actually tested and where it peaked),
#'   `p_stability`, its Monte Carlo error as `mc_se`, `mc_lo`, `mc_hi` and
#'   `near_threshold` (see the section above), `k` (of `prev`) and `k_new` (of
#'   `new_ma`).
#' @examples
#' library(metafor)
#' # Sufficiency is judged on the PRIOR review, stability on the UPDATED one.
#' # The two arguments are not interchangeable.
#' prev <- rma(yi = rep(log(0.50), 6), vi = rep(0.02, 6), measure = "RR")
#'
#' # Six more studies telling the same story: sufficient and stable.
#' steady <- rma(yi = rep(log(0.50), 12), vi = rep(0.02, 12), measure = "RR")
#' sufficiency_changepoint(prev, steady)
#'
#' # Now the later studies break away from the earlier ones. The statistic is
#' # the largest standardised split in the cumulative series, so a late shift
#' # is exactly what it is built to catch.
#' shifted <- rma(yi = c(rep(log(0.50), 6), rep(log(1.30), 6)),
#'                vi = rep(0.02, 12), measure = "RR")
#' res <- sufficiency_changepoint(prev, shifted)
#' res
#' res$detail[c("sufficient", "stable", "p_stability", "split")]
#' @export
sufficiency_changepoint <- function(prev, new_ma, min_k = 5,
                                    alpha_stability = 0.05,
                                    n_perm = 999, seed = 20260807) {
  check_rma_uni(prev, "prev")
  check_rma_uni(new_ma, "new_ma")
  check_count(min_k, "min_k")
  check_probability(alpha_stability, "alpha_stability")
  check_count(n_perm, "n_perm")
  check_seed(seed)
  yi <- as.numeric(new_ma$yi)
  vi <- as.numeric(new_ma$vi)
  k_new <- length(yi)
  if (k_new < min_k) {
    return(verdict_na("sufficiency_changepoint",
      paste0("needs at least ", min_k, " studies; found ", k_new)))
  }

  # Sufficiency is read off the meta-analysis AS PREVIOUSLY PUBLISHED, not the
  # updated one (see the roxygen note above and design doc section 4.4).
  yi_prev <- as.numeric(prev$yi)
  vi_prev <- as.numeric(prev$vi)
  k_prev  <- length(yi_prev)
  index <- failsafe_n(yi_prev, vi_prev) / (5 * k_prev + 10)
  # Resolved here for the same reason p_stability is resolved below: an index
  # that is not finite makes `sufficient` NA, and `NA && !stable` is NA
  # whenever the evidence is unstable, which turns the verdict's `if ()` into
  # an error. Half the method cannot be evaluated, so the answer is that it
  # does not apply -- not a "current" arrived at by NA && FALSE being FALSE.
  if (!is.finite(index)) {
    return(verdict_na("sufficiency_changepoint",
      paste0("the sufficiency index is not finite (a degenerate prior ",
             "meta-analysis?); sufficiency cannot be assessed")))
  }
  sufficient <- index > 1

  # Stability, by contrast, is read off the updated evidence: the cumulative
  # fixed-effect estimate after each study, in input order.
  cum_theta <- cumulative_effect(yi, vi)
  info      <- cumsum(1 / vi)

  if (!all(is.finite(cum_theta)) || !all(is.finite(info))) {
    return(verdict_na("sufficiency_changepoint",
      paste0("the cumulative effect is not finite (a study variance of ",
             "zero?); stability cannot be assessed")))
  }

  # Degenerate short circuit, BEFORE anything is computed. Byte-identical
  # studies span only rounding noise -- 12 identical studies span 2.2e-16 --
  # and there is no split to find in them: every study sits on the same
  # effect, which is maximal stability. Never test rounding noise.
  #
  # Keyed to the spread of `yi`, because that is what the change-point
  # statistic reads. It used to be keyed to `cum_theta[-1]`, the series the
  # old OLS slope was fitted to, and that is not the same condition: a
  # cumulative series can be flat from its second entry onwards while the
  # split isolating the FIRST study is enormous. `yi = c(10, -10, 0, 0, ...)`
  # is exactly that shape -- constant cumulative tail, max|Z_m| = 10.2 -- and
  # the old condition reported it as perfectly stable, p = 1.
  scale <- max(abs(yi))
  tol <- 4 * k_new * .Machine$double.eps * max(scale, .Machine$double.eps)
  if (diff(range(yi)) <= tol) {
    slope       <- 0
    z_shift     <- 0
    split       <- NA_integer_
    p_stability <- 1
    stable      <- TRUE
    # No permutation was drawn, so there is no Monte Carlo error to report.
    # This p is exactly 1 by construction, not estimated.
    perm_exceed <- NA_integer_
  } else {
    # The published slope, on the source's own x-axis (accumulated
    # information). Reported, not tested -- see the roxygen section above.
    slope   <- cum_drift_slope(cum_theta, info)
    z_shift <- stability_shift_z(yi, vi)
    split   <- stability_shift_at(yi, vi)
    if (!is.finite(z_shift)) {
      # No split carries information -- every study after the first is
      # weightless. There is nothing to test. Resolved here rather than left to
      # the p-value below, because `abs(perm) >= NA` is all NA and `na.rm` would
      # then count zero exceedances and report p = 1/1000, i.e. "unstable", from
      # an absence of evidence. The degenerate short circuit above normally
      # catches this shape first (a weightless tail leaves the cumulative series
      # constant); this is defence in depth behind it.
      p_stability <- NA_real_
      stable      <- TRUE
      perm_exceed <- NA_integer_
    } else {
      # Permutation test over study order (see the roxygen section above). The
      # seed is applied through with_preserved_seed() so that the caller's own
      # random stream is left exactly as it was found.
      perm <- with_preserved_seed(seed = seed, {
        vapply(seq_len(n_perm), function(i) {
          o <- sample.int(k_new)
          stability_shift_z(yi[o], vi[o])
        }, numeric(1))
      })
      # stability_shift_z() is already a maximum absolute value, so abs() here
      # is a no-op kept for symmetry with the two-sided reading of the source's
      # "absolute slope".
      perm_exceed <- sum(abs(perm) >= abs(z_shift), na.rm = TRUE)
      p_stability <- (1 + perm_exceed) / (n_perm + 1)
      # Defensive: a p-value that is not finite must be resolved deliberately,
      # never left to turn `if (sufficient && !stable)` into `if (NA)`. An
      # undeterminable drift is not evidence of drift, so it reads as stable.
      stable <- if (is.finite(p_stability)) {
        p_stability >= alpha_stability
      } else {
        TRUE
      }
    }
  }

  # Out-of-date requires sufficiency AND instability together. Insufficient
  # evidence (index <= 1) never triggers "out_of_date" by itself, regardless
  # of stability.
  out <- if (sufficient && !stable) "out_of_date" else "current"

  # `p_stability` is estimated from `n_perm` draws, so it moves with the seed.
  # The interval is built on the exceedance count and then carried through the
  # same (1 + x) / (n_perm + 1) map the reported p-value uses, so the bounds
  # live on exactly the scale of the number they bracket. Reported only: the
  # verdict above is unchanged, and `near_threshold` says whether a rerun
  # could plausibly have crossed `alpha_stability`.
  #
  # The degenerate branches leave `perm_exceed` as NA because they drew no
  # permutations at all, and both helpers return NA for a non-finite count, so
  # those verdicts come back with no interval rather than with one invented
  # around a p-value that was determined instead of estimated. No guard here:
  # an explicit `if (is.na(perm_exceed))` was written first and turned out to
  # be inert -- a mutation replacing it with `if (FALSE)` changed nothing --
  # so the NA propagation is the mechanism, and it is tested as such.
  shift   <- function(x) (1 + x * n_perm) / (n_perm + 1)
  # p = (1 + count) / (n + 1) is an affine map of count/n with slope
  # n / (n + 1), so the standard error of p is the standard error of the raw
  # proportion scaled by the same factor.
  p_mc_se <- mc_se(perm_exceed, n_perm) * n_perm / (n_perm + 1)
  p_ci    <- shift(mc_interval(perm_exceed, n_perm))
  new_verdict("sufficiency_changepoint", out, signal = index,
              detail = list(index = index, sufficient = sufficient,
                            stable = stable, slope = slope,
                            z_shift = z_shift, split = split,
                            p_stability = p_stability,
                            mc_se = p_mc_se,
                            mc_lo = p_ci[1], mc_hi = p_ci[2],
                            near_threshold =
                              mc_near_threshold(p_ci, alpha_stability),
                            k = k_prev, k_new = k_new))
}
