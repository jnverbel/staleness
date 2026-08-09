# Applicability across historical reviews

The published comparison of these five methods (Pattanittum et al. 2012)
applied them to 80 Cochrane reviews **selected for having a non-significant
pooled result**. Two of the five — `barrowman()` and `simulation()` — only
speak when the prior meta-analysis was not significant, so that cohort is the
one place where they can always be asked.

`applicability.R` asks what happens on reviews that were *not* selected that
way. It sweeps every dataset in `metadat`, keeps the ones that can carry a
backtest, and runs all five detectors over each.

The script runs *against* the package, so the package has to be available
first. From a clean clone:

    R CMD INSTALL .
    Rscript inst/applicability/applicability.R

or, while developing, without installing:

    Rscript -e 'pkgload::load_all("."); source("inst/applicability/applicability.R")'

It also needs `metadat`, which is in `Suggests`. It takes a few seconds and
the results are deterministic.

## What it finds

17 of metadat's 110 data frames survive the inclusion criteria. Across their
185 yearly cuts:

| Detector | Can answer at all | Mean sensitivity | Reviews scoring zero |
|---|---|---|---|
| `ottawa` | 17 of 17 | 0.320 | 8 |
| `rcma` | 17 of 17 | 0.309 | 9 |
| `sufficiency` | 17 of 17 | 0.041 | 11 |
| `simulation` | **5 of 17** | 0.167 | 2 of 3 |
| `barrowman` | **4 of 17** | 0.167 | 2 of 3 |

**168 of the 185 cuts (91%) had an already-significant prior meta-analysis.**
In 11 of the 17 reviews, every single cut did. Since `barrowman()` and
`simulation()` require a non-significant prior, they are inapplicable there by
construction — not wrong, but unable to be asked.

That is the result: on historical reviews taken as they come, two of the five
published methods are structurally silent, and the published comparison could
not see this because its cohort was selected to be exactly the case where they
are not.

`sufficiency()` is the opposite failure. It answers in all 17 and detects
almost nothing: a mean sensitivity of 0.041, and zero in 11 of them.

## Inclusion criteria

Applied in this order, and reported by the script as they bite:

1. a per-study year column;
2. a two-group effect measure buildable with `escalc()` — metadat stores the
   same 2×2 table under eight different column conventions, so the mapping is
   explicit rather than guessed;
3. at least 8 studies, so yearly cuts exist;
4. one effect per study-year. Datasets with an explicit `esid`, or repeated
   `(study, year)` pairs, are excluded: `evidence_stream()` treats each row as
   an independent study and would count one study more than once at a cut.

   Repeated *author* names are not nesting. `dat.bcg` has Comstock three
   times, in three different years — three trials, not three effects from one.
5. the backtest must yield at least three uncensored cuts.

Datasets that ship `yi`/`vi` already computed do not record their scale, and
the scale decides which branch `effect_ratio()` takes. Those are resolved one
at a time from metadat's own documentation; the ones that cannot be resolved
are excluded. `dat.hackshaw1998` is a **log odds ratio** — a ratio measure —
and treating it as a difference would be silently wrong.

## What this is not

Sensitivity and specificity here are scored against `truth_shift()`, which is
this package's own definition of what it means for a review to have gone out
of date, not a published outcome. **This compares the five methods against
each other over real evidence. It does not validate any of them against what
actually happened.**

That second claim is the one made by the four cases in
`tests/testthat/test-external-validation.R` — `dat.lau1992`, `dat.li2007`,
`dat.bangertdrowns2004` and `dat.laopaiboon2015` — where the criterion comes
from each source's published result. Those still number four, not seventeen.

One set of parameters (`horizon = 3`, `window = 5`, `min_k = 3`) is used for
every review so the comparison is like for like.

The headline figures do not depend on `horizon`: running the sweep at
`horizon = 6` returns the same 17 reviews, the same 168 of 185 cuts, and the
same 4-of-17 and 5-of-17 coverage. That is expected rather than reassuring —
whether a prior was already significant is a fact about the evidence at a cut,
and `horizon` only governs how far ahead the truth is evaluated — but it is
worth stating, because it means the finding is a property of the reviews and
not of this configuration. No sensitivity analysis to `window` or `min_k` is
claimed.
