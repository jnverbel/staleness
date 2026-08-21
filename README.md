# staleness

<!-- badges: start -->
[![R-CMD-check](https://github.com/jnverbel/staleness/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jnverbel/staleness/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22050352.svg)](https://doi.org/10.5281/zenodo.22050352)
<!-- badges: end -->

A meta-analysis is a photograph, not a standing fact. Roughly a quarter of
systematic reviews are out of date within two years of publication, and half
within five and a half (Shojania et al., 2007). Five statistical methods for
detecting this were published between 1999 and 2007, and we could not find a
reusable software implementation of any of them. Searching the metadata of all
24,734 CRAN packages on 2026-08-09 returns no implementation of the Ottawa
method, of Barrowman, or of recursive cumulative meta-analysis as an updating
diagnostic; every hit on those names is a false positive, listed one by one in
`inst/cran-search/`. The building blocks are there — `metafor` computes
cumulative meta-analyses and Rosenthal's fail-safe N, `fsn` and `meta` cover
pieces of the same ground, and `metagear` screens literature, which is the half
of the problem this package does not touch — but the assembled detectors are
not. Which is also why nobody has ever been able to run all five against real
history and see how they behave.

That search is reproducible and dated, and it says *we did not find one*
rather than *none exists*: it reads package metadata, not source code, and
`metafor` itself would not turn up under "cumulative meta-analysis" because
its `Description` never uses the phrase. See `inst/cran-search/README.md`.

`staleness` does two things:

- Applies the five published detectors to decide whether an existing
  meta-analysis is still current given evidence published since.
- Backtests those detectors against historical evidence, so their
  sensitivity, specificity and — the metric nobody has reported before —
  lead time can be measured with data rather than assumed from the
  methods papers.

Applying a published method is not the same as endorsing it, and the second
job keeps finding that out. Where a source's own statistic turns out to be
invalid — not merely weak, but unable to distinguish change from no change —
the package substitutes a calibrated one, keeps computing and reporting the
original as a diagnostic, and **renames the detector** so the substitution
travels with the name instead of hiding behind it. That is why the third
detector is `sufficiency_changepoint()` and not `sufficiency()`.

The package performs no literature searching, no study screening, and no
meta-analysis fitting of its own. It is built entirely on `metafor`.

## What this establishes, and what it does not

Read this before the numbers, because it bounds every one of them.

**It is a measuring instrument, not an updating strategy.** Deciding when to
update a review is not a statistical question alone: the consensus checklist
(Garner et al., 2016) treats the pooled estimate as one input among several,
alongside whether the question is still relevant and whether a recommendation
would change. Four of the six signals in the Ottawa method itself are
qualitative and no program can infer them from `yi` and `vi`. If you want a
tool that tells a review team when to update, this is not it and cannot become
it.

**The evaluation targets are not outcomes.** `calibration()` scores a detector
against `truth_shift()`, `truth_surprise()` or `truth_conclusion()`, and all
three observe *the pooled estimate moving*. None of them observes what a
review team did, whether a recommendation changed, or whether anyone was
harmed by acting on the old estimate. A sensitivity of 0.32 means "agreed with
a stated rule about the estimate 32% of the time", not "was right 32% of the
time". And `truth_shift` and `truth_surprise` share the identical numerator —
they are one distance on two scales, not two independent checks.

**The historical evidence is 17 reviews, and they are not a sample.** The
sweep in `inst/applicability/` covers the 17 of `metadat`'s 110 data frames
that can carry a backtest. 54 of the 93 exclusions are for the single reason
that the dataset records no per-study publication year. Datasets that do
record one skew towards the well-curated classics, so the surviving 17 are a
convenience set whose selection may correlate with what is being measured.
Which direction, if any, is unknown.

**Nothing is held out.** The 91% figure, `ottawa`'s specificity of 0.14 on a
null review, the six reviews with no true negatives — all come from those same
17. The findings are **exploratory**: generated and stated on one body of
evidence, with no confirmatory set behind them.

**Two of the five are not the published procedure.** `sufficiency_changepoint()`
substitutes a calibrated statistic for one that could not distinguish change
from no change, and `simulation()` simulates effects rather than participants
because the package never sees participant-level data. Both are named and
measured; neither is a literal reproduction. See below.

What survives all of that is worth having, and it is what the package is for:
a reproducible platform for studying these signals, plus three findings about
them that hold up — two of the five are structurally unable to answer on
ordinary evidence, one is unstable by construction, and one rests on a
statistic with no valid null distribution.

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

# `study_id` is required, and asked for rather than assumed: without it,
# several outcomes or time points from one trial would enter as separate
# studies and be counted twice. Here each row is a distinct trial.
stream <- evidence_stream(ma, date = es$year, study_id = es$trial)
prev   <- snapshot_at(stream, 1970)
new    <- window_between(stream, 1970, 1980)

check_currency(prev, new, methods = c("rcma", "ottawa", "sufficiency_changepoint"))
#> <staleness_check>
#>   studies: 7 prior + 6 new
#>   I2 (updated): 92.2 %
#>   note: heterogeneity above 75%; the pooled estimate is itself debatable
#>
#> rcma                    current
#>   signal: 1.21
#> ottawa                  current
#>   signal: 0.86
#> sufficiency_changepoint current
#>   signal: 4.24

bt <- backtest(stream, cuts = "yearly")
summary(bt)
plot(bt)

# What series those rates are rates over: how many estimates from how many
# studies, over what span, on what measure, under what model, and whether
# sample sizes are there for barrowman() to use.
eligibility(bt)
```

## The five detectors

| Detector | Idea | Known limitation |
|---|---|---|
| `rcma` | Signal when the pooled effect moves by 50% or more. | Undefined for difference measures with a prior effect near zero; returns `not_applicable` rather than a spurious number. |
| `ottawa` | Change in significance (p < 0.04) or effect size (≥ 50%). | Four qualitative signals in the original method are not automated; supplied by the analyst instead. |
| `barrowman` | Ratio of participants contributed by new studies to the number needed to reach significance. | Only applies when the prior meta-analysis was not significant. |
| `sufficiency_changepoint` | Rosenthal's fail-safe N and the stability of the cumulative effect. | Fail-safe N has been discredited since Becker (2005); implemented for fidelity to the published method, not as an endorsement. |
| `simulation` | Simulated power of the next batch of studies. | Flagged none of 80 reviews in the one published comparison of all five methods. |

The two signals above differ because they measure different things:
`rcma` compares the pooled effects, `ottawa` compares the relative risk
reductions. See `?ottawa`.

### Where the implementation departs from its source

The column above lists limitations of the *methods*. Two detectors also depart
from the *procedure* their sources describe, and the results of those two
should not be read as a literal reproduction:

- **`sufficiency_changepoint`** tests stability with a change-point statistic
  (`max_m |Z_m|`) under an order-permutation null, not with the ordinary least
  squares slope of the cumulative series the source specifies. That slope has
  no valid null distribution — the cumulative mean is autocorrelated by
  construction and convergent by the law of large numbers — and fired on 209
  of 300 samples containing no change at all. The published slope is still
  computed and reported in `detail$slope`; it decides nothing. The detector
  carries the substitution in its name for that reason: its sufficiency half
  is the published one, its stability half is not, and a function called
  `sufficiency()` would have implied otherwise.
- **`simulation`** departs in exactly one way, and it cannot be removed: the
  source simulates at the level of **participants** and computes an effect from
  them, while this package never sees participant-level data and simulates the
  effect directly from `yi` and `vi`. Its other three distinctive choices — a
  *t* distribution rather than a normal, one simulated study carrying the
  combined precision of the recent ones rather than one per study, and a strict
  threshold — are not departures at all: each is what Pattanittum et al. (2012),
  Appendix S1, specifies. See `?simulation`.

`vignette("methods")` measures each departure and says where it degrades.

`vignette("methods", package = "staleness")` covers each of these in full,
including the formula, the original source, and the critique in more
detail.

## Backtesting

In the one published comparison that ran all five side by side, on 80
Cochrane reviews, two of them flagged nothing at all: the sufficiency and
stability method
and the simulation method each identified zero out-of-date reviews. The three
that did discriminate — Ottawa (34 reviews), recursive CMA (7) and Barrowman
(7) — agreed at Kappa = 0.14, essentially chance (Pattanittum et al., 2012). `staleness`
exists to make that comparison repeatable, on any body of evidence, with a
design that avoids the trap of defining the evaluation target using the same rule a
detector is scored against: three operational targets are implemented and
described as what they are, contaminated detector-target pairs are marked in
the data itself, and `lead_time()` measures how far ahead of the evidence a detector
actually fires, not just whether it eventually fires at all.

What that repeatability buys is agreement with a stated criterion, measured
over evidence that really accumulated in the order it accumulated. It does not
buy a verdict on which method is correct, because no target here is an
outcome — see the scope section above.

See `vignette("backtesting", package = "staleness")` for the full argument
and a worked backtest.

## License

MIT. See `LICENSE`.

## Citing this work

The findings are described in a preprint, archived on Zenodo under CC BY 4.0:

> Núñez, J. (2026). *Five published signals for updating meta-analyses, applied
> and evaluated: an exploratory study of when they can be used at all.* Zenodo.
> <https://doi.org/10.5281/zenodo.22050352>

**It is a preprint: it has not been peer reviewed.** Every figure in it is
reproduced by the tag `preprint-v1`, not by the default branch.

`10.5281/zenodo.22050352` is that version. `10.5281/zenodo.22050351` is the
concept DOI, which always resolves to the newest version — cite the first when
you mean the text you read, the second when you mean the work.

To cite the software as such — a specific version, or the implementation rather
than the findings:

```r
citation("staleness")
```
