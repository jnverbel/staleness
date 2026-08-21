#!/usr/bin/env python3
"""
Turn the harvested versions into review chains, and report the denominator the
whole study rests on: how many CONSECUTIVE version pairs have an "Authors'
conclusions" section at both ends.

A pair is what carries the outcome. Version 3 alone says nothing; version 3
against version 4 is where a conclusion can have changed. So the count that
matters is pairs, not versions and not reviews, and it is smaller than either.

Writes corpus/data/pairs.jsonl, one line per consecutive pair.
"""

import json
import os
import re
from collections import defaultdict, Counter
from difflib import SequenceMatcher

HERE = os.path.dirname(__file__)
IN = os.path.join(HERE, "data", "versions.jsonl")
OUT = os.path.join(HERE, "data", "pairs.jsonl")

# Same naive parser the spike measured at 46%, kept identical on purpose: the
# figure quoted in the write-up has to be the figure this produces.
EFF = re.compile(
    r"\b(RR|OR|HR|MD|SMD|RD|IRR)\s*[=:]?\s*(-?\d+\.\d+)[,;]?\s*"
    r"(?:95%\s*(?:CI|confidence interval)\s*[:=]?\s*)"
    r"(-?\d+\.\d+)\s*(?:to|,|-|–)\s*(-?\d+\.\d+)", re.I)
K_STUDIES = re.compile(r"(\d+)\s+(?:studies|trials|RCTs)", re.I)
N_PART = re.compile(r"([\d,]{3,})\s+(?:participants|women|men|patients|infants|children)", re.I)


def effect(abstract):
    """First reported effect, plus the text that precedes it.

    The context is not decoration. Comparing an estimate across two versions
    is only meaningful if both quote the SAME outcome, and nothing guaranteed
    that: 42% of pairs turned out to name clearly different ones and 12% did
    not even share a measure. Carrying the preceding words is what lets that
    be checked instead of assumed.
    """
    m = EFF.search(abstract or "")
    if not m:
        return None
    est, lo, hi = float(m.group(2)), float(m.group(3)), float(m.group(4))
    # An interval that does not contain its own estimate is unusable, and the
    # cause turns out not to be the parser: 14 published Cochrane abstracts
    # report one, e.g. "MD -2.76, 95% CI 3.57 to 1.96" with the minus signs
    # dropped, and "RD 0.03, 95% CI -0.01 to -0.07" running backwards. Not
    # repaired -- guessing which digit is wrong would invent data -- and not
    # swapped either, since a sign lost on both bounds is not fixed by
    # reordering them. Refused, and counted.
    if not (lo <= est <= hi):
        return None
    pre = re.sub(r"<[^>]+>", " ", (abstract or "")[max(0, m.start() - 120):m.start()])
    return {"measure": m.group(1).upper(), "est": est, "lo": lo, "hi": hi,
            "context": re.sub(r"\s+", " ", pre).strip()[-110:]}



def counts(abstract):
    """Participant and study counts, for barrowman.

    `barrowman()` sums participants across a snapshot: it needs how many were
    in the prior meta-analysis and how many the new evidence adds. Nothing
    parsed them until now -- N_PART and K_STUDIES were written above and never
    called -- so the detector was unevaluable on this corpus for want of a
    field the pipeline already knew how to look for.

    Every match is kept, not just the first. An abstract that says "3,500
    participants" once is unambiguous; one that says it three times with three
    different numbers is a judgement call, and the caller should see that
    rather than inherit a silent pick. The first is offered as the headline
    because Cochrane abstracts open with the total, but `all` is what says
    whether to trust it.
    """
    text = abstract or ""
    part = [int(m.group(1).replace(",", "")) for m in N_PART.finditer(text)]
    studies = [int(m.group(1)) for m in K_STUDIES.finditer(text)]
    return {"n_participants": part[0] if part else None,
            "n_participants_all": part,
            "k_studies": studies[0] if studies else None}

# Below this, the two estimates are treated as answering different questions.
# A threshold, not a fact: it is a proxy for outcome identity, chosen so that
# the flattering direction needs an argument rather than a default.
CONTEXT_MATCH = 0.80


def comparable(ea, eb):
    """Whether two parsed effects may be compared at all."""
    if not ea or not eb:
        return False
    if ea["measure"] != eb["measure"]:
        return False
    return SequenceMatcher(None, ea["context"], eb["context"]).ratio() >= CONTEXT_MATCH


def main():
    by_review = defaultdict(list)
    with open(IN) as fh:
        for line in fh:
            r = json.loads(line)
            by_review[r["review_id"]].append(r)

    reviews = len(by_review)
    chain_len = Counter()
    pairs = both_concl = both_effect = protocols = same_abs = 0
    comparable_n = 0

    # An unsuffixed DOI was read as version 1, and that was wrong: Europe PMC
    # carries several index records under the same bare DOI for one review, so
    # 551 reviews had multiple "version 1" records. Zipping the sorted list
    # then paired those with each other, producing 681 pairs (9%) whose two
    # ends were the SAME version -- 542 of them with the later end dated
    # earlier, which is how it surfaced. A pair of a version with itself has no
    # outcome to observe and would have entered the analysis silently.
    #
    # Two guards, because either alone leaves a hole: collapse records that
    # share (review, version) keeping the earliest date, and then pair only
    # ends whose versions actually differ.
    collapsed = 0
    with open(OUT, "w") as out:
        for rid, versions in by_review.items():
            best = {}
            for v in versions:
                k = v["version"]
                if k not in best or (v["date"] or "9999") < (best[k]["date"] or "9999"):
                    if k in best:
                        collapsed += 1
                    best[k] = v
                else:
                    collapsed += 1
            versions = sorted(best.values(), key=lambda r: r["version"])
            reals = [v for v in versions if not v["is_protocol"]]
            protocols += len(by_review[rid]) - len(reals) - 0
            chain_len[len(reals)] += 1
            for a, b in zip(reals, reals[1:]):
                if b["version"] <= a["version"]:
                    continue
                pairs += 1
                ca, cb = a.get("conclusions"), b.get("conclusions")
                ea, eb = effect(a.get("abstract")), effect(b.get("abstract"))
                if ca and cb:
                    both_concl += 1
                if ea and eb:
                    both_effect += 1
                comp = comparable(ea, eb)
                comparable_n += comp
                # Two versions whose ENTIRE abstract is byte-identical are not
                # an update that left the conclusions standing -- that case is
                # real and common (662 of them, with Main results visibly
                # rewritten). These 131 are the same record reaching the index
                # twice. Flagged rather than dropped here, so the count stays
                # visible downstream instead of vanishing from a denominator.
                same_abstract = (a.get("abstract") or "") == (b.get("abstract") or "")
                same_abs += same_abstract
                out.write(json.dumps({
                    "review_id": rid,
                    "same_abstract": same_abstract,
                    "from_version": a["version"], "to_version": b["version"],
                    "from_date": a["date"], "to_date": b["date"],
                    "from_doi": a["doi"], "to_doi": b["doi"],
                    "title": b["title"],
                    "conclusions_from": ca, "conclusions_to": cb,
                    "effect_from": ea, "effect_to": eb,
                    "comparable_effects": comp,
                    "counts_from": counts(a.get("abstract")),
                    "counts_to": counts(b.get("abstract")),
                }, ensure_ascii=False) + "\n")

    print(f"versiones cosechadas      : {sum(len(v) for v in by_review.values())}")
    print(f"  registros colapsados    : {collapsed}  (misma revision y version)")
    print(f"revisiones distintas      : {reviews}")
    print(f"\nlongitud de cadena (revisiones por numero de versiones):")
    for n in sorted(chain_len):
        if n <= 8 or chain_len[n] > 5:
            print(f"  {n:2d} version(es): {chain_len[n]:5d}")
    print(f"\nPARES consecutivos        : {pairs}")
    print(f"  con el resumen ENTERO identico (excluibles): {same_abs}")
    print(f"  con conclusiones en AMBOS: {both_concl}  ({100*both_concl/pairs:.0f}%)")
    print(f"  con efecto en AMBOS      : {both_effect}  ({100*both_effect/pairs:.0f}%)")
    print(f"  con efectos COMPARABLES  : {comparable_n}  ({100*comparable_n/pairs:.0f}%)"
          f"  <- denominador de los detectores")
    print(f"\nescrito: {OUT}")


if __name__ == "__main__":
    main()
