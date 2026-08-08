# staleness 0.0.0.9000

## Bug fixes

* `sufficiency()` no longer reports perfect stability over a real change
  point. The degenerate short circuit was still keyed to `cum_theta[-1]`, the
  series the *previous* OLS statistic was fitted to, and a cumulative series
  can be flat from its second entry onwards while the split isolating the
  first study is enormous. On `yi = c(10, -10, 0, ...)` the detector reported
  `z_shift = 0`, `p_stability = 1`, `stable = TRUE` where the statistic is
  10.17 and the permutation *p* is 0.137. The condition now reads the spread
  of `yi`, which is what the change-point statistic actually measures.

* `simulation()` without a `seed` now returns different draws on each call.
  The caller's random stream was being restored unconditionally, so every
  call restarted from the same state and returned byte-identical results
  while the documentation promised the opposite. Inside `backtest(seed =
  NULL)` this made the Monte Carlo error perfectly correlated across cuts
  rather than independent. The caller's stream is still left exactly as it
  was found, with a seed or without one.

* `with_preserved_seed()` no longer creates `.Random.seed` in a session that
  had none — a global footprint in the function whose purpose is to leave
  none.

* `evidence_stream()` no longer refuses to build when the `ni` that `metafor`
  derived on its own contains `NA`. An `NA` in a sample size the caller typed
  is still an error, but one in an auto-filled `ni` used to abort
  construction and take down the four detectors that never look at `ni`.
  `barrowman()` already answers `not_applicable` naming the non-finite `n`.

* `sufficiency()` with `min_k = 1` no longer emits three "no non-missing
  arguments" warnings before quietly reporting `stable`.

* `check_currency()` now rejects a fitted `rma` passed as `new`. That
  argument wants the new evidence alone, while `rcma()` and `ottawa()` want
  the updated model, so confusing the two is the natural mistake — and an
  `rma.uni` carries `$yi`, `$vi` and `$k`, so duck typing let it through
  every guard and pooled the prior studies twice without a word.

## Documentation

* Corrected the scale-pivotality claim for `stability_shift_z()`, in both the
  roxygen and `vignette("methods")`. The statistic is invariant when the
  effects and their variances are rescaled *together*; multiplying `vi` alone
  by `c` divides every `Z_m` by `sqrt(c)`. The *p*-value is unaffected either
  way, since the permuted statistics rescale with it.

* Every exported function and dataset now carries a runnable `\examples{}`
  block (17 of 17).

## Infrastructure

* `R CMD check` now runs on GitHub Actions across macOS, Windows and Linux,
  from R 4.2 through R devel.

* `DESCRIPTION` declares `R (>= 4.2)`, which is the oldest release the checks
  actually cover. It previously claimed 4.1, which nothing verified and
  nothing in the package required.
