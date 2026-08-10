#!/usr/bin/env python3
"""
Tests for the corpus pipeline.

Written after three defects in one afternoon, every one of them silent and
every one caught by eye rather than by anything automatic:

  1. An unsuffixed DOI was read as version 1, so reviews with several bare-DOI
     records paired those with each other. 681 pairs (9%) had the same version
     at both ends and looked perfectly well-formed.
  2. 142 pairs had a byte-identical abstract at both ends -- the same record
     indexed twice, not an update.
  3. The detector denominator was reported as 1,907 when the honest number was
     748, because nothing checked that two versions quoted the SAME outcome.
     A threefold overstatement, in the flattering direction.

The corpus numbers are what would go into a paper, so they need what the R
package already has. Two kinds of test below: unit tests on the parsing, which
always run, and invariants over the produced files, which skip when the data
is absent. The invariants are mostly RANGE checks, because a range check is
what caught defect 1 -- a gap between versions came back as -13 years, and a
negative gap cannot exist.

    python3 -m unittest discover -s corpus -v
"""

import importlib.util
import json
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


chains = load("chains", "02-chains.py")
screen_mod = load("screen", "03-screen.py")
harvest = load("harvest", "01-harvest.py")


def jsonl(name):
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return [json.loads(l) for l in fh]


# --------------------------------------------------------------------------
# Unit tests: the parsing, on fixtures. These always run.
# --------------------------------------------------------------------------

class TestDoiVersion(unittest.TestCase):
    """Defect 1 lived here: the version comes from the DOI, and the bare form
    is not a reliable version 1."""

    def parse(self, doi):
        m = harvest.DOI_RE.search(doi)
        return (m.group(1).upper(), int(m.group(2)) if m.group(2) else 1) if m else None

    def test_suffixed_doi_gives_review_and_version(self):
        self.assertEqual(self.parse("10.1002/14651858.CD005343.pub7"), ("CD005343", 7))
        self.assertEqual(self.parse("10.1002/14651858.CD000259.pub2"), ("CD000259", 2))

    def test_case_is_normalised(self):
        self.assertEqual(self.parse("10.1002/14651858.cd005343.pub6")[0], "CD005343")

    def test_bare_doi_reads_as_version_one(self):
        # Documented rather than trusted: several records can share this form,
        # which is what produced 681 same-version pairs. The pipeline must
        # therefore collapse duplicates, and TestPairInvariants checks it did.
        self.assertEqual(self.parse("10.1002/14651858.CD000510"), ("CD000510", 1))

    def test_non_cochrane_doi_is_rejected(self):
        self.assertIsNone(self.parse("10.1136/bmj.i3507"))


class TestEffectParser(unittest.TestCase):

    def test_extracts_measure_estimate_and_interval(self):
        e = chains.effect("pooled RR 0.72, 95% CI 0.54 to 0.95, 12 studies")
        self.assertEqual(e["measure"], "RR")
        self.assertAlmostEqual(e["est"], 0.72)
        self.assertAlmostEqual(e["lo"], 0.54)
        self.assertAlmostEqual(e["hi"], 0.95)

    def test_accepts_en_dash_as_a_separator(self):
        e = chains.effect("OR 1.30, 95% CI 1.02–1.66")
        self.assertIsNotNone(e)
        self.assertAlmostEqual(e["hi"], 1.66)

    def test_negative_estimates_survive(self):
        e = chains.effect("MD -1.00, 95% CI -1.21 to -0.79")
        self.assertAlmostEqual(e["est"], -1.00)
        self.assertAlmostEqual(e["lo"], -1.21)

    def test_returns_none_when_there_is_no_effect(self):
        self.assertIsNone(chains.effect("We found no eligible trials."))
        self.assertIsNone(chains.effect(""))
        self.assertIsNone(chains.effect(None))

    def test_interval_brackets_the_estimate(self):
        # A parse that puts the estimate outside its own interval has matched
        # the wrong numbers, and every downstream flag would be computed on it.
        for txt in ("RR 0.72, 95% CI 0.54 to 0.95",
                    "MD -1.00, 95% CI -1.21 to -0.79",
                    "HR 0.98, 95% CI 0.82 to 1.18"):
            e = chains.effect(txt)
            self.assertLessEqual(e["lo"], e["est"], txt)
            self.assertLessEqual(e["est"], e["hi"], txt)

    def test_context_is_captured_for_the_outcome_check(self):
        e = chains.effect("Mortality at 28 days was lower: RR 0.72, 95% CI 0.54 to 0.95")
        self.assertIn("Mortality", e["context"])


class TestComparability(unittest.TestCase):
    """Defect 3. Two estimates may only be compared if they answer the same
    question, and this is the guard that says so."""

    def eff(self, measure, context):
        return {"measure": measure, "est": 1.0, "lo": 0.5, "hi": 2.0,
                "context": context}

    def test_same_outcome_and_measure_is_comparable(self):
        a = self.eff("RR", "all-cause mortality at 28 days was lower with treatment:")
        b = self.eff("RR", "all-cause mortality at 28 days was lower with treatment:")
        self.assertTrue(chains.comparable(a, b))

    def test_different_measure_is_never_comparable(self):
        # An RR against an OR is not a borderline case. 12% of pairs did this.
        a = self.eff("RR", "all-cause mortality at 28 days")
        b = self.eff("OR", "all-cause mortality at 28 days")
        self.assertFalse(chains.comparable(a, b))

    def test_different_outcome_is_not_comparable(self):
        # The shape of the case that surfaced the defect: CD007145 went from
        # RR 0.72 to RR 1.14 on what were almost certainly two outcomes.
        a = self.eff("RR", "excessive gestational weight gain was less common")
        b = self.eff("RR", "caesarean delivery occurred more often in the")
        self.assertFalse(chains.comparable(a, b))

    def test_a_missing_effect_is_not_comparable(self):
        a = self.eff("RR", "mortality")
        self.assertFalse(chains.comparable(a, None))
        self.assertFalse(chains.comparable(None, a))
        self.assertFalse(chains.comparable(None, None))


class TestNullAndDirection(unittest.TestCase):
    """A ratio is null at 1 and a difference at 0. Using one rule for both
    would silently mislabel every mean difference in the corpus."""

    def eff(self, m, est, lo, hi):
        return {"measure": m, "est": est, "lo": lo, "hi": hi, "context": ""}

    def test_ratio_measures_use_one_as_the_null(self):
        self.assertTrue(screen_mod.crosses_null(self.eff("RR", 1.0, 0.8, 1.3)))
        self.assertFalse(screen_mod.crosses_null(self.eff("RR", 0.72, 0.54, 0.95)))

    def test_difference_measures_use_zero_as_the_null(self):
        self.assertTrue(screen_mod.crosses_null(self.eff("MD", 0.1, -0.5, 0.7)))
        self.assertFalse(screen_mod.crosses_null(self.eff("MD", -1.0, -1.21, -0.79)))
        # And the case the shared rule would get backwards: an MD interval
        # that excludes zero but contains 1.
        self.assertFalse(screen_mod.crosses_null(self.eff("MD", 0.5, 0.2, 1.4)))

    def test_direction_is_taken_against_the_right_null(self):
        self.assertEqual(screen_mod.sign_of(self.eff("RR", 0.72, 0.5, 0.9)), -1)
        self.assertEqual(screen_mod.sign_of(self.eff("RR", 1.30, 1.1, 1.6)), 1)
        self.assertEqual(screen_mod.sign_of(self.eff("MD", 0.72, 0.5, 0.9)), 1)


class TestCertaintyTiers(unittest.TestCase):
    """"very low-certainty" also matches "low-certainty". Taking the first
    match rather than the lowest would read a downgrade as no change."""

    def test_very_low_is_not_read_as_low(self):
        self.assertEqual(screen_mod.tier("very low-certainty evidence",
                                         screen_mod.CERTAINTY), 1)
        self.assertEqual(screen_mod.tier("low-certainty evidence",
                                         screen_mod.CERTAINTY), 2)

    def test_absent_vocabulary_gives_none(self):
        self.assertIsNone(screen_mod.tier("the trials were small",
                                          screen_mod.CERTAINTY))

    def test_a_tier_move_needs_a_tier_at_both_ends(self):
        r = {"conclusions_from": "high-certainty evidence shows benefit",
             "conclusions_to": "the trials were small", "effect_from": None,
             "effect_to": None}
        self.assertFalse(screen_mod.screen(r)["certainty_moved"])


class TestScoreBounds(unittest.TestCase):

    def test_identical_conclusions_score_zero(self):
        txt = "Antibiotics probably reduce the risk of complications."
        r = {"conclusions_from": txt, "conclusions_to": txt,
             "effect_from": None, "effect_to": None}
        s = screen_mod.screen(r)
        self.assertEqual(s["similarity"], 1.0)
        self.assertEqual(s["score"], 0.0)

    def test_similarity_stays_inside_the_unit_interval(self):
        for a, b in (("", "x"), ("abc", "abc"), ("one thing", "another entirely")):
            s = screen_mod.screen({"conclusions_from": a, "conclusions_to": b,
                                   "effect_from": None, "effect_to": None})
            self.assertGreaterEqual(s["similarity"], 0.0)
            self.assertLessEqual(s["similarity"], 1.0)


# --------------------------------------------------------------------------
# Invariants over the produced files. Skipped when the data is absent, so the
# unit tests above still run on a clean clone.
# --------------------------------------------------------------------------

class TestPairInvariants(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pairs = jsonl("pairs.jsonl")
        if cls.pairs is None:
            raise unittest.SkipTest("corpus/data/pairs.jsonl absent; run 02-chains.py")

    def test_every_pair_moves_forward_a_version(self):
        # Defect 1, made loud. 681 pairs failed this and nothing noticed.
        bad = [r for r in self.pairs if r["to_version"] <= r["from_version"]]
        self.assertEqual(bad, [], f"{len(bad)} pairs do not advance a version")

    def test_almost_no_pair_goes_backwards_in_time(self):
        # The check that actually caught defect 1: a negative gap cannot
        # exist. A handful survive as genuine indexing oddities, and the count
        # is pinned so it cannot grow unnoticed.
        inverted = [r for r in self.pairs
                    if r["from_date"] and r["to_date"] and r["to_date"] < r["from_date"]]
        self.assertLessEqual(len(inverted), 10,
                             f"{len(inverted)} pairs are dated backwards")

    def test_no_pair_repeats_a_review_and_version(self):
        seen = set()
        for r in self.pairs:
            key = (r["review_id"], r["from_version"], r["to_version"])
            self.assertNotIn(key, seen, f"duplicate pair {key}")
            seen.add(key)

    def test_intervals_bracket_their_estimates(self):
        for r in self.pairs:
            for side in ("effect_from", "effect_to"):
                e = r.get(side)
                if e:
                    self.assertLessEqual(e["lo"], e["est"], r["review_id"])
                    self.assertLessEqual(e["est"], e["hi"], r["review_id"])

    def test_comparable_implies_both_effects_and_one_measure(self):
        # Defect 3, made loud. Whatever "comparable" comes to mean, it can
        # never be true where an estimate is missing or the measures differ.
        for r in self.pairs:
            if r.get("comparable_effects"):
                self.assertIsNotNone(r["effect_from"], r["review_id"])
                self.assertIsNotNone(r["effect_to"], r["review_id"])
                self.assertEqual(r["effect_from"]["measure"],
                                 r["effect_to"]["measure"], r["review_id"])

    def test_the_denominators_are_nested(self):
        # comparable <= both effects <= all pairs, and comparable is the one
        # the detector analysis may use. Stated as an invariant so no future
        # edit can quietly report the wider set as the usable one.
        both = [r for r in self.pairs if r["effect_from"] and r["effect_to"]]
        comp = [r for r in self.pairs if r.get("comparable_effects")]
        self.assertLessEqual(len(comp), len(both))
        self.assertLessEqual(len(both), len(self.pairs))
        self.assertGreater(len(comp), 0)

    def test_the_detector_denominator_is_not_silently_inflated(self):
        # The number that was overstated threefold. Pinned to the measured
        # value with room to move, so a parser change that doubles it has to
        # be a deliberate edit to this line.
        comp = sum(1 for r in self.pairs if r.get("comparable_effects"))
        self.assertLess(comp, 1200,
                        f"{comp} comparable pairs: verify the outcome match before believing it")


class TestScreenInvariants(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.rows = jsonl("screened.jsonl")
        cls.sample = jsonl("sample.jsonl")
        if cls.rows is None or cls.sample is None:
            raise unittest.SkipTest("screened/sample absent; run 03-screen.py")

    def test_duplicated_records_are_excluded(self):
        # Defect 2.
        self.assertEqual([r for r in self.rows if r.get("same_abstract")], [])

    def test_every_screened_pair_has_conclusions_at_both_ends(self):
        for r in self.rows:
            self.assertTrue(r["conclusions_from"])
            self.assertTrue(r["conclusions_to"])

    def test_rows_are_ordered_by_score(self):
        scores = [r["screen"]["score"] for r in self.rows]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_the_sample_is_stratified_and_not_top_n(self):
        from collections import Counter
        strata = Counter(r["stratum"] for r in self.sample)
        self.assertEqual(set(strata), {"alto", "medio", "bajo"})
        self.assertEqual(len(set(strata.values())), 1, "strata differ in size")
        # The point of stratifying: the low stratum must really be low, or the
        # sample says nothing about what the screen misses.
        hi = max(r["screen"]["score"] for r in self.sample if r["stratum"] == "bajo")
        lo = min(r["screen"]["score"] for r in self.sample if r["stratum"] == "alto")
        self.assertLess(hi, lo, "strata overlap: the sample cannot estimate misses")

    def test_the_sample_holds_no_duplicates(self):
        keys = [(r["review_id"], r["from_version"], r["to_version"]) for r in self.sample]
        self.assertEqual(len(keys), len(set(keys)))

    def test_the_sample_is_drawn_from_the_screened_set(self):
        known = {(r["review_id"], r["from_version"], r["to_version"]) for r in self.rows}
        for r in self.sample:
            self.assertIn((r["review_id"], r["from_version"], r["to_version"]), known)


if __name__ == "__main__":
    unittest.main(verbosity=2)
