# Calibration measurements

`sufficiency_changepoint()` departs from its published description in one place: the
stability half is tested with a change-point statistic under a permutation
null, instead of the ordinary least squares slope of the cumulative series
that the source specifies.

A deviation from a published method needs evidence, not an assertion. This
directory holds that evidence in runnable form. Every calibration figure
quoted in `?sufficiency_changepoint`, `vignette("methods")` and the JOSS paper is produced
by `calibration.R`:

The script runs *against* the package, so the package has to be available
first. From a clean clone:

```
R CMD INSTALL .
Rscript inst/calibration/calibration.R
```

or, while developing, without installing:

```
Rscript -e 'pkgload::load_all("."); source("inst/calibration/calibration.R")'
```

It takes about a minute, needs nothing beyond the package itself, and is
deterministic — the seeds are the ones the original measurements used.

## What it measures

| | |
|---|---|
| **E1** | No drift at all. False alarms must sit near the nominal 5%. |
| **E2** | Power against a shift confined to the last 10 of 30 studies — a mature review plus a small batch of new trials, the regime this package exists for. |
| **E3** | Deterministic scan of *where* the shift sits. The replaced statistic went silent once the shift moved past study 10. |
| **E4** | No drift, but a variance schedule correlated with time. This is the failure mode that a permutation null on the raw series cannot see. |

It also reconstructs **both** statistics that were replaced — the OLS *t*-test
the source describes, and the permutation test on the slope that was tried
next — and runs them on the same series. Otherwise the "before" half of every
before-and-after claim would be a number nobody can re-derive.

## If a figure stops reproducing

That is a finding, not a nuisance. The script exits non-zero and says so.
Either the statistic changed and the documentation has not caught up, or the
published number was wrong. Both are worth knowing; neither should be
resolved by editing the expected value until it matches.
