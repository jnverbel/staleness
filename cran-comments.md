## Submission

This is a new submission.

## Test environments

* macOS 26.5 (arm64), R 4.6.1 — local
* macOS latest (arm64), R release — GitHub Actions
* Windows Server 2022 (ucrt), R release — GitHub Actions
* Ubuntu 24.04, R devel (4.7.0) — GitHub Actions
* Ubuntu 24.04, R release (4.6.1) — GitHub Actions
* Ubuntu 24.04, R 4.5.3, 4.4.3, 4.3.3 and 4.2.3 — GitHub Actions
* win-builder, R devel (2026-08-08 r90381 ucrt)
* win-builder, R release (4.6.1 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the usual `New submission` from the incoming-feasibility check.
`checking examples ... OK`, `checking tests ... OK` and
`checking PDF version of manual ... OK` on all of the above.

The URLs it flags as 404 all point at the package's own repository, which is
private while this submission is prepared and becomes public on acceptance.

## Maintainer

Javier Núñez <jnverbel@gmail.com>. A personal address rather than an
institutional one, deliberately: it will outlive any affiliation, and CRAN
needs a maintainer address that keeps working.

## Notes for the reviewer

* `Language: en-GB`. The prose is consistently British throughout (21
  occurrences of `-ise`/`-ised` forms, none American). `inst/WORDLIST` holds
  the remaining flagged tokens, which are variable names (`yi`, `vi`, `ni`),
  author surnames from the bibliography, and domain terms (`backtest`,
  `metafor`, `pivotality`).

* `metaanalyses` appears once, in `vignette("methods")`. It is not a typo: it
  is the spelling used in the title of Ioannidis & Lau (2001), *PNAS*
  98(3):831–836, and the citation reproduces it verbatim.

* The package depends on `metafor` for all model fitting and does no fitting
  of its own. `metadat` is used only to supply an example dataset in the
  vignettes and is in `Suggests`.

* `sufficiency()` implements a published method with one declared deviation:
  the stability half is tested with a change-point statistic under a
  permutation null rather than with an OLS slope on the cumulative series.
  The reason is that the published test has no valid null distribution — a
  *t*-test on a cumulative mean is autocorrelated by construction — and fired
  on 209 of 300 samples containing no change. The deviation, the measured
  calibration of the replacement, and the regimes where it loses power are
  documented in `vignette("methods")` and in `?sufficiency`. The published
  slope is still computed and returned.

* No function writes to the user's home filespace, options, or par. The two
  detectors that randomise restore the caller's `.Random.seed` on exit,
  including on error, and leave a session that had no random stream with none.
