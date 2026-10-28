from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

from . import adaptive_reference, stationary_reference


class StationaryReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = stationary_reference.run()
        cls.joint = {(tuple(row["state"]), row["memory"]): F(row["probability"])
                     for row in cls.result["stationary_joint"] if F(row["probability"])}

    def test_exact_stationarity_through_executor(self):
        # Independent replay through native events, rather than the rational
        # matrix residual used by the construction.
        self.assertEqual(adaptive_reference.advance(self.joint), self.joint)

    def test_retained_receipt_and_mutated_support(self):
        stationary_reference.check_receipt(
            self.result, stationary_reference.HERE / "stationary_reference_results.json")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.json"
            changed = json.loads(json.dumps(self.result))
            changed["exact_certificate"]["joint_support_dimension"] = 132
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "differs from exact replay"):
                stationary_reference.check_receipt(self.result, path)

    def test_recurrent_support_and_source_faithfulness(self):
        certificate = self.result["exact_certificate"]
        self.assertEqual(certificate["joint_support_dimension"], 129)
        self.assertEqual(len(certificate["joint_zero_probability_states"]), 3)
        self.assertEqual(certificate["closed_recurrent_components"], 1)
        self.assertEqual(certificate["source_marginal_support_dimension"], 66)
        self.assertEqual(sum(self.joint.values()), F(1))
        self.assertTrue(all(p > 0 for p in self.joint.values()))

    def test_nonuniform_joint_reference_differs_from_blank_product(self):
        source = adaptive_reference.source_marginal(self.joint)
        self.assertGreater(len(set(source.values())), 1)
        self.assertNotEqual(self.joint, adaptive_reference.initial_joint())
        # Resetting the memory bit destroys stationarity: the reference is
        # a joint source/controller state, not just its source marginal.
        blank = {(x, 0): p for x, p in source.items()}
        self.assertNotEqual(adaptive_reference.advance(blank), blank)


if __name__ == "__main__":
    unittest.main()
