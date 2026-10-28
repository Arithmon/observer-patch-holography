"""Exact controls for the source/controller reference distinction."""
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

from . import adaptive_reference as experiment


class AdaptiveReferenceTests(unittest.TestCase):
    def test_retained_receipt_and_mutated_probability(self):
        result = experiment.run()
        experiment.check_receipt(result, experiment.HERE / "adaptive_reference_results.json")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.json"
            changed = json.loads(json.dumps(result))
            changed["second_attempt"]["source_probabilities"][0]["probability"] = "0"
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "differs from exact replay"):
                experiment.check_receipt(result, path)

    def test_adaptive_counterexample(self):
        initial = experiment.initial_joint()
        first = experiment.advance(initial)
        second = experiment.advance(first)
        self.assertEqual(set(experiment.source_marginal(first).values()), {F(1, 66)})
        marginal = experiment.source_marginal(second)
        deviations = [p - F(1, 66) for p in marginal.values()]
        self.assertEqual(sum(marginal.values()), F(1))
        self.assertEqual(sum(d != 0 for d in deviations), 30)
        self.assertEqual(max(map(abs, deviations)), F(1, 960))
        witness = (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
        self.assertEqual(marginal[witness] - F(1, 66), -F(29, 31680))

    def test_fixed_mixtures_preserve_source_reference(self):
        for exploration in (F(1, 16), F(1, 2), F(1)):
            for repair in (F(1, 16), F(1, 2), F(1)):
                joint = experiment.initial_joint()
                for _ in range(3):
                    joint = experiment.advance(joint, exploration=exploration,
                                               repair=repair, adaptive=False)
                self.assertEqual(set(experiment.source_marginal(joint).values()), {F(1, 66)})

    def test_uniform_exploration_removes_feedback_effect(self):
        joint = experiment.initial_joint()
        for _ in range(3):
            joint = experiment.advance(joint, exploration=F(1))
        self.assertEqual(set(experiment.source_marginal(joint).values()), {F(1, 66)})

    def test_identity_holds_do_not_restore_stationarity(self):
        joint = experiment.initial_joint()
        for _ in range(2):
            joint = experiment.advance(joint, repair=F(1, 2))
        self.assertEqual(experiment.summarize(joint)["maximum_absolute_deviation"], "1/3840")

    def test_record_laws_share_actual_diagonal_map(self):
        # This independently checks the classical reduction through the
        # complete recorded C/E maps, rather than assuming the mixture label
        # erases a different source transition law.
        from ..common_source.native_defect_decision import complete_record_maps
        from ..observer import native
        states = native.states()
        maps = complete_record_maps()
        for law in ("C", "E"):
            for seam in experiment.slow_geometry.source().seams():
                for coin in (0, 1):
                    actual = {}
                    for row in maps:
                        if row["law"] == law and tuple(row["seam"]) == seam and row["coin"] == coin:
                            for source, target in row["pairs"]:
                                self.assertNotIn(source, actual)
                                actual[source] = target
                    self.assertEqual(len(actual), 78)
                    for source, state in enumerate(states):
                        out, _ = native.native(state, seam, coin)
                        self.assertEqual(states[actual[source]], out)


if __name__ == "__main__":
    unittest.main()
