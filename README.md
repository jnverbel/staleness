# staleness

A meta-analysis is a photograph, not a standing fact. Roughly a quarter of
systematic reviews are out of date within two years of publication, and half
within five and a half (Shojania et al., 2007). Five statistical methods for
detecting this have been published since 2003. None of them, until now, had
a reusable software implementation — which is also why nobody has ever been
able to run all five against real history and find out which of them
actually work.

`staleness` does two things:

- Applies the five published detectors to decide whether an existing
  meta-analysis is still current given evidence published since.
- Backtests those detectors against historical evidence, so their
  sensitivity, specificity and — the metric nobody has reported before —
  lead time can be measured with data rather than assumed from the
  methods papers.

The package performs no literature searching, no study screening, and no
meta-analysis fitting of its own. It is built entirely on `metafor`.

## Installation

`staleness` is not yet on CRAN. Install the development version from
source:

```r
# install.packages("remotes")
remotes::install_local("path/to/staleness")
```

## Usage

```r
library(staleness)
library(metafor)

dat <- metadat::dat.bcg
es  <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
              data = dat)
ma  <- rma(yi, vi, data = es)

stream <- evidence_stream(ma, date = es$year)
prev   <- snapshot_at(stream, 1970)
new    <- window_between(stream, 1970, 1980)

check_currency(prev, new, methods = c("rcma", "ottawa", "sufficiency"))
#> <staleness_check>
#>   studies: 7 prior + 6 new
#>   I2 (updated): 92.2 %
#>   note: heterogeneity above 75%; the pooled estimate is itself debatable
#>
#> rcma         current
#>   signal: 1.21
#> ottawa       current
#>   signal: 1.21
#> sufficiency  OUT OF DATE
#>   signal: 4.24
#>
#>   detectors disagree

bt <- backtest(stream, cuts = "yearly")
summary(bt)
plot(bt)
```

## The five detectors

| Detector | Idea | Known limitation |
|---|---|---|
| `rcma` | Signal when the pooled effect moves by 50% or more. | Undefined for difference measures with a prior effect near zero; returns `not_applicable` rather than a spurious number. |
| `ottawa` | Change in significance (p < 0.04) or effect size (≥ 50%). | Four qualitative signals in the original method are not automated; supplied by the analyst instead. |
| `barrowman` | Ratio of participants contributed by new studies to the number needed to reach significance. | Only applies when the prior meta-analysis was not significant. |
| `sufficiency` | Rosenthal's fail-safe N and the stability of the cumulative effect. | Fail-safe N has been discredited since Becker (2005); implemented for fidelity to the published method, not as an endorsement. |
| `simulation` | Simulated power of the next batch of studies. | Flagged none of 80 reviews in the one published comparison of all five methods. |

`vignette("methods", package = "staleness")` covers each of these in full,
including the formula, the original source, and the critique in more
detail.

## Backtesting

The five methods agreed with each other at Kappa = 0.14 in the one published
comparison that ran them side by side — essentially chance. `staleness`
exists to make that comparison repeatable, on any body of evidence, with a
design that avoids the trap of defining "ground truth" using the same rule a
detector is scored against: three independent truth definitions are
implemented, contaminated detector-truth pairs are marked in the data
itself, and `lead_time()` measures how far ahead of the evidence a detector
actually fires, not just whether it eventually fires at all.

See `vignette("backtesting", package = "staleness")` for the full argument
and a worked backtest.

## License

MIT. See `LICENSE`.

## Citing this package

```r
citation("staleness")
```
