#!/usr/bin/env python3
"""
Join the coding back to the sample and estimate the population rate.

The sample is stratified by screen score, so the raw proportion of major
changes in it is NOT an estimate of anything: the high stratum was
deliberately over-represented. The population rate has to be rebuilt by
weighting each stratum by its real size, and that is the only number that
belongs in a sentence about how often conclusions change.

Also reports what the screen bought. A screen that orders work is worth having
only if the low stratum really does contain fewer events than the high one;
if the strata are indistinguishable, the score is noise and the sample should
have been simple random.
"""

import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))


def main():
    codes = {int(k): v for k, v in json.load(open(f"{D}/coding/codes-claude.json")).items()}
    key = {r["n"]: r for r in json.load(open(f"{D}/coding/key.json"))}
    screened = [json.loads(l) for l in open(f"{D}/screened.jsonl")]

    # The population the sample was actually drawn from, copied from
    # 03-screen.py rather than guessed: pairs parsing an effect at BOTH ends,
    # sorted by score, cut into thirds. This line used to filter on
    # comparable_effects instead, which is a different and much smaller set
    # (560 against 1825), and reweighting the strata by ITS size silently
    # relabelled a prevalence estimated for one population as if it belonged
    # to another. The confirmatory run of 2026-08-22 fixed its outcome
    # threshold on the wrong number because of this line; see RESULTS.md.
    usable = [r for r in screened if r.get("effect_from") and r.get("effect_to")]
    usable.sort(key=lambda r: -r["screen"]["score"])
    third = len(usable) // 3
    bounds = {"alto": (0, third), "medio": (third, 2 * third),
              "bajo": (2 * third, len(usable))}
    n_us = len(usable)
    pop = {k: hi - lo for k, (lo, hi) in bounds.items()}

    print(f"población muestreada (pares con efecto en AMBOS extremos): {n_us}")
    print(f"muestra codificada: {len(codes)}\n")

    print(f"{'estrato':8s} {'n':>4s} {'mayor':>6s} {'menor':>6s} {'ninguno':>8s} {'% mayor':>9s}")
    rate = {}
    for s in ("alto", "medio", "bajo"):
        ns = [n for n in codes if key[n]["stratum"] == s]
        c = Counter(codes[n] for n in ns)
        rate[s] = c["major"] / len(ns)
        lo, hi = wilson(c["major"], len(ns))
        print(f"{s:8s} {len(ns):4d} {c['major']:6d} {c['minor']:6d} {c['none']:8d} "
              f"{100*rate[s]:8.0f}%  [{100*lo:.0f}-{100*hi:.0f}]")

    est = sum(rate[s] * pop[s] / n_us for s in pop)
    raw = sum(1 for n in codes if codes[n] == "major") / len(codes)
    print(f"\nproporción CRUDA en la muestra      : {100*raw:.0f}%   <- no estima nada")
    print(f"estimación REPONDERADA por estrato  : {100*est:.0f}%")
    print(f"eventos esperados en los {n_us} pares : ~{round(est*n_us)}")

    # What the screen bought. If alto and bajo are indistinguishable, the score
    # is noise and the stratification was pointless.
    print(f"\nseparación que consigue el cribado: {100*rate['alto']:.0f}% en el alto "
          f"contra {100*rate['bajo']:.0f}% en el bajo")
    if rate["bajo"] > 0:
        print(f"  razón alto/bajo: {rate['alto']/rate['bajo']:.1f}x")

    # And the honest comparison against the figure this was planned around.
    print(f"\nFrench et al. (2005) midió 9% sobre 254 revisiones actualizadas.")
    print(f"Aquí sale {100*est:.0f}%. Las dos no son comparables sin discutir por que;")
    print(f"ver el README.")

    # The set the detectors actually run on is SMALLER and not a random slice
    # of the above: requiring the two effects to be comparable selects pairs
    # whose abstracts repeat the same outcome, and those score low. Its
    # prevalence has to be rebuilt from its own stratum composition, not
    # inherited from the sampled population, and the two differ by a factor
    # of two.
    score_lo = {}
    for name, (lo, hi) in bounds.items():
        score_lo[name] = usable[hi - 1]["screen"]["score"]

    def stratum_of(score):
        if score >= score_lo["alto"]:
            return "alto"
        if score >= score_lo["medio"]:
            return "medio"
        return "bajo"

    detector_set = [r for r in screened if r.get("comparable_effects")]
    comp = Counter(stratum_of(r["screen"]["score"]) for r in detector_set)
    n_det = len(detector_set)
    est_det = sum(rate[s] * comp[s] / n_det for s in comp)

    print(f"\nsubconjunto sobre el que corren los detectores "
          f"(efectos COMPARABLES): {n_det}")
    for s in ("alto", "medio", "bajo"):
        print(f"  {s:6s} {comp[s]:4d}  ({100*comp[s]/n_det:2.0f}%"
              f" frente al {100*pop[s]/n_us:2.0f}% de la población muestreada)")
    print(f"  prevalencia REPONDERADA en este subconjunto: {100*est_det:.0f}%"
          f"   <- NO es el {100*est:.0f}% de arriba")
    print(f"  eventos esperados en los {n_det} pares: ~{round(est_det*n_det)}")


if __name__ == "__main__":
    main()
