#!/usr/bin/env python3
"""
Screen the version pairs for a likely change of conclusion, and draw a
stratified sample for human adjudication.

This orders the work. It does not decide it. French et al. (2005) define the
outcome as a change that "alter[s] the substance or meaning of a section or
alter[s] the interpretation", judged by two investigators independently, with
minor changes -- "changes in style or wording that do not alter the substance
or meaning" -- explicitly NOT counting. No text statistic can make that
distinction, which is why what follows produces interpretable FLAGS a human
can check rather than a verdict.

Three of the flags read Cochrane's own controlled vocabulary rather than
generic text similarity, which is what makes them defensible: GRADE-informed
conclusions are written to a template, so moving between its tiers is a change
of substance by construction, not a change of style.

Writes corpus/data/screened.jsonl and corpus/data/sample.jsonl.
"""

import json
import math
import os
import random
import re
from collections import Counter
from difflib import SequenceMatcher

HERE = os.path.dirname(__file__)
IN = os.path.join(HERE, "data", "pairs.jsonl")
OUT = os.path.join(HERE, "data", "screened.jsonl")
SAMPLE = os.path.join(HERE, "data", "sample.jsonl")

SAMPLE_PER_STRATUM = 40
SEED = 20260810

# GRADE's informative statements, in tiers. Cochrane asks authors to phrase
# conclusions this way, so a move between tiers is a change in what the review
# claims to know -- not a rewording of it.
CERTAINTY = [
    (4, re.compile(r"\bhigh[- ]certainty\b|\bhigh[- ]quality evidence\b", re.I)),
    (3, re.compile(r"\bmoderate[- ]certainty\b|\bmoderate[- ]quality evidence\b", re.I)),
    (2, re.compile(r"\blow[- ]certainty\b|\blow[- ]quality evidence\b", re.I)),
    (1, re.compile(r"\bvery low[- ]certainty\b|\bvery low[- ]quality evidence\b", re.I)),
]
HEDGE = [
    (3, re.compile(r"\bprobably (?:reduces|increases|improves|makes|results)\b|\blikely (?:to )?(?:reduce|increase)\b", re.I)),
    (2, re.compile(r"\bmay (?:reduce|increase|improve|make|result|have)\b", re.I)),
    (1, re.compile(r"\bthe evidence is very uncertain\b|\bwe are uncertain\b", re.I)),
]
NO_EVIDENCE = re.compile(
    r"\binsufficient evidence\b|\bno (?:randomised |randomized )?(?:evidence|trials|studies)\b"
    r"|\btoo few (?:trials|studies)\b|\bno conclusions? can be drawn\b", re.I)
NO_DIFFERENCE = re.compile(r"\blittle (?:or|to) no difference\b|\bno (?:clear |significant )?difference\b", re.I)


def tier(text, table):
    """Lowest matching tier, since 'very low' also matches 'low'."""
    hits = [t for t, rx in table if rx.search(text or "")]
    return min(hits) if hits else None


def crosses_null(eff):
    if not eff:
        return None
    ratio = eff["measure"] in ("RR", "OR", "HR", "IRR")
    null = 1.0 if ratio else 0.0
    return eff["lo"] <= null <= eff["hi"]


def sign_of(eff):
    if not eff:
        return None
    null = 1.0 if eff["measure"] in ("RR", "OR", "HR", "IRR") else 0.0
    return 0 if eff["est"] == null else (1 if eff["est"] > null else -1)


def screen(r):
    a, b = r["conclusions_from"] or "", r["conclusions_to"] or ""
    ea, eb = r.get("effect_from"), r.get("effect_to")

    flags = {}
    # 1. How much of the text survived. Low similarity is necessary for a major
    #    change and nowhere near sufficient: a full rewrite can leave the
    #    substance untouched, which is precisely French's "minor change".
    flags["similarity"] = round(SequenceMatcher(None, a, b).ratio(), 3)

    # 2. Certainty tier moved. A review that went from low- to high-certainty
    #    is claiming something different about the same effect.
    ca, cb = tier(a, CERTAINTY), tier(b, CERTAINTY)
    flags["certainty_from"], flags["certainty_to"] = ca, cb
    flags["certainty_moved"] = bool(ca and cb and ca != cb)

    # 3. GRADE hedging tier moved: "may reduce" -> "reduces" is a stronger
    #    claim in Cochrane's own template, not a stylistic choice.
    ha, hb = tier(a, HEDGE), tier(b, HEDGE)
    flags["hedge_moved"] = bool(ha and hb and ha != hb)

    # 4. "Insufficient evidence" appearing or disappearing is the clearest
    #    substantive move a conclusion can make.
    na, nb = bool(NO_EVIDENCE.search(a)), bool(NO_EVIDENCE.search(b))
    flags["no_evidence_flipped"] = na != nb

    # 5. And the same for "little or no difference".
    flags["no_difference_flipped"] = bool(NO_DIFFERENCE.search(a)) != bool(NO_DIFFERENCE.search(b))

    # 6-7. What the pooled estimate did. Reported alongside because a
    #      conclusion that changed while the estimate stood still, or the
    #      reverse, is the interesting case for this package specifically.
    flags["significance_flipped"] = (
        None if ea is None or eb is None else crosses_null(ea) != crosses_null(eb))
    flags["direction_flipped"] = (
        None if ea is None or eb is None else sign_of(ea) != sign_of(eb))

    score = 0.0
    score += 2.0 * (1 - flags["similarity"])
    score += 1.0 * flags["certainty_moved"]
    score += 1.0 * flags["hedge_moved"]
    score += 1.5 * flags["no_evidence_flipped"]
    score += 1.0 * flags["no_difference_flipped"]
    score += 1.0 * bool(flags["significance_flipped"])
    score += 1.0 * bool(flags["direction_flipped"])
    flags["score"] = round(score, 3)
    return flags


def main():
    rows = []
    with open(IN) as fh:
        for line in fh:
            r = json.loads(line)
            if not (r["conclusions_from"] and r["conclusions_to"]):
                continue
            # Not an update: the same record indexed twice. See 02-chains.py.
            if r.get("same_abstract"):
                continue
            r["screen"] = screen(r)
            rows.append(r)

    rows.sort(key=lambda r: -r["screen"]["score"])
    with open(OUT, "w") as out:
        for r in rows:
            out.write(json.dumps(r, ensure_ascii=False) + "\n")

    usable = [r for r in rows if r.get("effect_from") and r.get("effect_to")]

    print(f"pares con conclusiones en ambos extremos : {len(rows)}")
    print(f"  de ellos con efecto en ambos           : {len(usable)}")
    print("\nreparto de banderas (sobre los que tienen conclusiones):")
    for k in ("certainty_moved", "hedge_moved", "no_evidence_flipped",
              "no_difference_flipped"):
        n = sum(1 for r in rows if r["screen"][k])
        print(f"  {k:24s} {n:5d}  ({100*n/len(rows):.0f}%)")
    for k in ("significance_flipped", "direction_flipped"):
        d = [r for r in rows if r["screen"][k] is not None]
        n = sum(1 for r in d if r["screen"][k])
        print(f"  {k:24s} {n:5d}  ({100*n/len(d):.0f}% de {len(d)} con efecto)")

    # The GRADE flags only bite on recent pairs, because the vocabulary is
    # recent: 0.0% of pairs ending by 2010 move a certainty tier, 0.5% for
    # 2011-2016, 3.5% from 2017. Reported so nobody reads their low overall
    # rate as evidence that certainty rarely changes.
    print("\nbanderas GRADE por epoca (el vocabulario se estandariza ~2011):")
    for lab, lo, hi in (("<=2010", 0, 2010), ("2011-2016", 2011, 2016), (">=2017", 2017, 9999)):
        sub = [r for r in rows if r["to_date"] and lo <= int(r["to_date"][:4]) <= hi]
        if not sub:
            continue
        cm = sum(1 for r in sub if r["screen"]["certainty_moved"])
        print(f"  {lab:10s} n={len(sub):5d}  certainty_moved {100*cm/len(sub):4.1f}%")

    sims = sorted(r["screen"]["similarity"] for r in rows)
    q = lambda p: sims[int(p * (len(sims) - 1))]
    print(f"\nsimilitud del texto de conclusiones: mediana {q(.5):.2f}, "
          f"decil inferior {q(.1):.2f}, decil superior {q(.9):.2f}")

    # Stratified by score, NOT top-N. Adjudicating only the highest scores
    # would measure nothing about what the screen misses, and the miss rate is
    # what decides whether the screen can be trusted to order the rest.
    rng = random.Random(SEED)
    strata = {"alto": usable[:len(usable)//3],
              "medio": usable[len(usable)//3:2*len(usable)//3],
              "bajo": usable[2*len(usable)//3:]}
    sample = []
    print("\nmuestra estratificada para adjudicacion humana:")
    for name, pool in strata.items():
        take = rng.sample(pool, min(SAMPLE_PER_STRATUM, len(pool)))
        lo = min(r["screen"]["score"] for r in pool)
        hi = max(r["screen"]["score"] for r in pool)
        print(f"  {name:6s} {len(take):3d} de {len(pool):5d}  (score {lo:.2f}-{hi:.2f})")
        for r in take:
            r["stratum"] = name
            sample.append(r)
    rng.shuffle(sample)          # blind the adjudicator to the stratum order
    with open(SAMPLE, "w") as out:
        for r in sample:
            out.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(f"\n{len(sample)} pares en {SAMPLE}")


if __name__ == "__main__":
    main()
