"""Independent record counting, capacity, and compact-controller checks."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
import unittest

from model.observer import capacity as subject


def independent_event(x, seam, coin):
    a, b = seam
    total = x[a] + x[b]
    y = list(x)
    y[a] = total // 2 if coin == 0 else (total + 1) // 2
    y[b] = total - y[a]
    return tuple(y), (coin, x[a] - y[a])


def oracle(binary):
    words = Counter()
    prefixes = {()}
    paths = []
    pairs = combinations(range(12), 2) if binary else combinations_with_replacement(range(12), 2)
    for a, b in pairs:
        x = [0] * 12
        x[a] += 1
        x[b] += 1
        for c1, c2 in product((0, 1), repeat=2):
            y, r = independent_event(x, (0, 1), c1)
            z, s = independent_event(y, (0, 5) if r[1] > 0 else (0, 1), c2)
            words[r, s] += 1
            prefixes.update(((r,), (r, s)))
            paths.append((tuple(x), y, z, r, s))
    return words, prefixes, paths


class CapacityJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = subject.run()
        cls.source, cls.capacity = subject.dependencies()

    def test_reachable_atoms_and_exact_counts(self):
        for binary, name, count in ((True, "binary_66", 11), (False, "full_78", 17)):
            words, prefixes, paths = oracle(binary)
            family = self.result["families"][name]
            reported = {tuple(tuple(r) for r in row["word"]): row["histories"]
                        for row in family["terminal_word_counts"]}
            self.assertEqual(reported, words)
            self.assertEqual(len(words), count)
            self.assertEqual(family["positive_complete_histories"], len(paths))
            self.assertEqual(family["classical_prefix_encoding"]["states"], len(prefixes))

    def test_existing_evaluator_with_full_carrier_and_reachability(self):
        for binary, count in ((True, 11), (False, 17)):
            words, witnesses, _, _, _ = subject.enumerate_source(self.source, binary)
            packet = subject.full_tape_packet(self.source, self.capacity, words, witnesses, "test")
            self.assertEqual(len(packet["observers"]["observer"]), 36)
            self.assertEqual(len(packet["reachability_witnesses"]), count)
            result = self.capacity.evaluate_terminal(packet)
            self.assertEqual(result["exact_zero_error_capacity"], count)
            self.assertEqual(result["capacity_dimension"], 3354)
            self.assertFalse(result["saturation_passed"])
            supports = list(packet["projection_supports"].values())
            self.assertTrue(all(len(s) == 78 and min(s) >= 546 and max(s) < 3354 for s in supports))
            self.assertEqual(len(set().union(*(set(s) for s in supports))), 78 * count)
            broken = deepcopy(packet)
            keys = list(broken["projection_supports"])
            broken["projection_supports"][keys[1]] = broken["projection_supports"][keys[0]]
            self.assertEqual(self.capacity.evaluate_terminal(broken)["status"], "NO_CAPACITY_CARRIER_REPRESENTATION")

    def test_readout_loss_is_exact_confusability(self):
        controls = self.result["families"]["binary_66"]["readout_controls"]
        expected = {"full_tape": 11, "first_letter_only": 4, "second_letter_only": 4,
                    "chosen_seam_only": 2, "forecast_only": 3,
                    "forecast_and_chosen_seam": 4, "complete_erasure": 1}
        self.assertEqual({k: v["M0"] for k, v in controls.items()}, expected)
        self.assertEqual(controls["full_tape"]["confusability_edges"], 0)
        self.assertEqual(controls["complete_erasure"]["confusability_edges"], 55)

    def test_four_bit_controller_executes_every_benchmark_path(self):
        words, witnesses, prefixes, transitions, histories = subject.enumerate_source(self.source, True)
        encoding = subject.encode_prefixes(prefixes, transitions)
        compiled = subject.compile_controller(self.source, histories, encoding)
        self.assertEqual(set(compiled["decode"]), set(range(16)))
        self.assertEqual(encoding["minimum_fixed_width_bits"], 4)
        _, _, independent = oracle(True)
        for initial, first_loads, last_loads, first_record, last_record in independent:
            before = compiled["encode"][()]
            first = compiled["encode"][(first_record,)]
            last = compiled["encode"][(first_record, last_record)]
            self.assertIn((self.source.Q(1, 2), first, first_loads), subject.compact_tick(compiled, before, initial))
            self.assertIn((self.source.Q(1, 2), last, last_loads), subject.compact_tick(compiled, first, first_loads))
            self.assertEqual(subject.compact_tick(compiled, last, last_loads), ((self.source.Q(1), last, last_loads),))
            view = subject.compact_view(compiled, first)
            self.assertEqual(view["phase"], 1)
            self.assertEqual(view["seam"], (0, 5) if first_record[1] > 0 else (0, 1))
            self.assertEqual(view["forecast"], self.source.FORECAST[first_record])
            self.assertEqual(subject.compact_view(compiled, last)["prefix"], (first_record, last_record))
        with self.assertRaises(ValueError):
            subject.compact_tick(compiled, 16, independent[0][0])
        with self.assertRaises(ValueError):
            subject.compact_tick(compiled, 0, (2,) + (0,) * 11)
        absent = next((code, loads) for code in compiled["decode"]
                      for loads, _, _, _, _ in independent if (code, loads) not in compiled["admitted"])
        with self.assertRaises(ValueError):
            subject.compact_tick(compiled, *absent)

    def test_same_atom_count_does_not_make_noisy_capacity_large(self):
        words, _, _ = oracle(True)
        names = [subject.word_name(w) for w in words]
        # Positive common blank support makes every pair confusable, however tiny.
        channels = [{"rows": {name: {name: Fraction(999, 1000), "blank": Fraction(1, 1000)} for name in names}}]
        graph = self.capacity.compound_confusability_graph(names, channels)
        self.assertEqual(len(self.capacity.maximum_independent_set(graph)), 1)

    def test_canonical_receipt(self):
        self.assertEqual((Path(__file__).with_name("results") / "capacity.json").read_text(), subject.canonical(self.result))


if __name__ == "__main__":
    unittest.main()
