# Preprint

**Published on Zenodo, 2026-08-21: <https://doi.org/10.5281/zenodo.22050352>.**
CC BY 4.0, and the basis for a submission to *Research Synthesis Methods* or the
*Journal of Clinical Epidemiology*. `10.5281/zenodo.22050351` is the concept
DOI, which always resolves to the newest version.

The route there is worth recording, because affiliation decided it twice.
medRxiv declined this manuscript on 2026-08-11 without reviewing it: it requires
every author to hold an organizational affiliation "that provides oversight of
research activities", sole authorship as an independent researcher does not
qualify, and no rewrite changes that. SSRN turns out to require the same thing —
"a current, verified institutional affiliation" — and will not open its
submission form without one. arXiv is not a shortcut either: since 2026-01-21
every new author needs a personal endorsement. Zenodo asks for no affiliation,
runs no moderation, and issues the DOI on publication.

`preprint.md` + `preprint.bib`. Render with:

    pandoc preprint.md --citeproc --bibliography=preprint.bib -o preprint.pdf

## Why post a preprint before the journal route

It is free, it takes days rather than months, and it dates the findings. The
Ottawa instability result in particular is the kind of thing another group
could reach independently.

## What this draft deliberately does

Says n = 17, "exploratory" and "nothing is held out" **in the abstract and in
§2**, not in a limitations section at the end. Every objection a methodological
reviewer would raise is stated before the result it qualifies, because that
converts an objection into declared scope.

Reports the 52% conclusion-change rate as a working figure against French's
9%, with four untested explanations, rather than as a finding. If the original
French coding arrives, this section is rewritten around it.

Says which two of the five detectors are not literal reproductions, in the
methods and again in the limitations.

## Every number here is checked against the data, not against memory

The check found two of them wrong before this was committed: 4,519 pairs with
conclusions at both ends was really 4,530 (subtracting all 142 duplicate-
abstract pairs, when not all of them had conclusions at both ends), and "9,862
reviews contribute 6,686 pairs" was misleading — only the 4,132 reviews with
more than one version contribute any.

Re-run the check before any submission:

    python3 -c "..." # see the session notes; it joins the text to corpus/data
