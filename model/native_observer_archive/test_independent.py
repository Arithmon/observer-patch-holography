"""Full-domain independent coverage, isolated before producer comparison."""

from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

from . import independent
from .compare_independent import compare, complete_weights, custody, mutations


class IndependentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(independent.__file__).resolve()
        flags = ["-B"] + ([] if __debug__ else ["-O"])
        # A fresh direct-file process cannot inherit production modules from
        # unittest discovery. Its import and file audit guards are exercised.
        run = subprocess.run([sys.executable, *flags, str(path)],
                             capture_output=True, text=True, timeout=240)
        if run.returncode:
            raise RuntimeError(run.stdout + run.stderr)
        cls.oracle = json.loads(run.stdout)
        cls.comparison = compare()

    def test_full_domain_columns_and_public_partitions(self):
        self.assertEqual(self.comparison["prepared_register_columns_compared"], 4096 * 12 * 2 * 5 * 5)
        self.assertEqual(self.oracle["native_dirty_states_both_inverse_orders"], 4096 * 12 * 2 * 2 * 2)
        for row in self.oracle["configurations"]:
            self.assertEqual(row["entries_per_prefix"], [49152] * 5)
            self.assertGreater(row["first_noops"], 0)
        self.assertIn("blocked", self.oracle["isolation"])

    def test_coherence_and_retained_failure_weights(self):
        self.assertGreater(self.oracle["coherence"]["block_density_difference_frobenius"], 0)
        small = self.comparison["small_rail"]
        self.assertGreater(small["completion_weight"], 0)
        self.assertLess(small["completion_weight"], 1)
        self.assertLess(small["unnormalized_block_max_error"], 2e-10)
        weights = list(map(Fraction, [0, 0, 0, 0, 1]))
        with self.assertRaisesRegex(ValueError, "postselection"):
            complete_weights(weights, [Fraction(1, 2), 0, 0, 0, Fraction(1, 2)])

    def test_controls_reach_scientific_predicates(self):
        result = mutations()
        self.assertEqual(result["swapped_address_bit"]["native_target"], 0)
        self.assertEqual(result["swapped_address_bit"]["mutated_target"], 3)
        self.assertIn("dirty inverse", result["dirty_overwrite"]["predicate"])
        self.assertEqual(self.oracle["completed_source_tv"], {"feedback": "1/2", "uphill": "1"})

    def test_input_admission_survives_optimization(self):
        for word, address in ((True, 0), (4096, 0), (0, -1), (0, 12), (0, .5)):
            with self.assertRaises(ValueError):
                independent.prepared(word, address)
        for memory in (True, -1, 2, .5):
            with self.assertRaises(ValueError):
                independent.reversible(0, 0, memory, (0, 0, 0, 0))
        with self.assertRaises(ValueError):
            independent.instructions(0, "unknown")


class PinScopeTests(unittest.TestCase):
    def test_proof_digest_is_snapshot_metadata_and_code_is_bound(self):
        contract, protocol = independent.admit()
        contract["audited_sha256"]["model/native_observer_archive/PROOF.md"] = "0" * 64
        with patch.object(independent, "admit", return_value=(contract, protocol)):
            self.assertEqual(custody(), protocol)
        contract["audited_sha256"]["model/native_memory_clock/native.py"] = "0" * 64
        with patch.object(independent, "admit", return_value=(contract, protocol)):
            with self.assertRaisesRegex(ValueError, "audited source custody"):
                custody()


if __name__ == "__main__":
    unittest.main()
