# staleness

<!-- badges: start -->
[![R-CMD-check](https://github.com/jnverbel/staleness/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jnverbel/staleness/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

A meta-analysis is a photograph, not a standing fact. Roughly a quarter of
systematic reviews are out of date within two years of publication, and half
within five and a half (Shojania et al., 2007). Five statistical methods for
detecting this were published between 1999 and 2007. None of them, until now, had a
reusable software implementation: a search of all 24,708 CRAN packages returns
no hit for the Ottawa method, for Barrowman, for recursive cumulative
meta-analysis as an updating diagnostic, or for updating systematic reviews at
all. The building blocks are there — `metafor` computes cumulative
meta-analyses and Rosenthal's fail-safe N — but the detectors are not. Which is
also why nobody has ever been able to run all five against real history and
find out which of them actually work.

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
GitHub:

```r
# install.packages("remotes")
remotes::install_github("jnverbel/staleness")
```

`metafor` does all the model fitting and is installed alongside. Building the
vignettes additionally needs `knitr`, `rmarkdown` and `metadat`:

```r
remotes::install_github("jnverbel/staleness", build_vignettes = TRUE)
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
#>   signal: 0.86
#> sufficiency  current
#>   signal: 4.24

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

The two signals above differ because they measure different things:
`rcma` compares the pooled effects, `ottawa` compares the relative risk
reductions. See `?ottawa`.

### Where the implementation departs from its source

The column above lists limitations of the *methods*. Two detectors also depart
from the *procedure* their sources describe, and the results of those two
should not be read as a literal reproduction:

- **`sufficiency`** tests stability with a change-point statistic
  (`max_m |Z_m|`) under an order-permutation null, not with the ordinary least
  squares slope of the cumulative series the source specifies. That slope has
  no valid null distribution — the cumulative mean is autocorrelated by
  construction and convergent by the law of large numbers — and fired on 209
  of 300 samples containing no change at all. The published slope is still
  computed and reported in `detail$slope`.
- **`simulation`** departs in four ways: it draws from a *t* distribution
  rather than a normal, simulates one study carrying the combined precision of
  the recent ones rather than one per study, and uses a strict threshold. The
  fourth cannot be removed: the source simulates participants, and this package
  only ever sees effect sizes and their variances.

`vignette("methods")` measures each departure and says where it degrades.

`vignette("methods", package = "staleness")` covers each of these in full,
including the formula, the original source, and the critique in more
detail.

## Backtesting

In the one published comparison that ran all five side by side, on 80
Cochrane reviews, two of them flagged nothing at all: the sufficiency method
and the simulation method each identified zero out-of-date reviews. The three
that did discriminate — Ottawa (34 reviews), recursive CMA (7) and Barrowman
(7) — agreed at Kappa = 0.14, essentially chance (Pattanittum et al., 2012). `staleness`
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
