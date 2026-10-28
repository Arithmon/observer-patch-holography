"""Exact reference attachment and real-operation boundary controls."""
from fractions import Fraction as F
import json
from pathlib import Path
import tempfile
import unittest

import sympy as sp

from . import modular_access as access, observer
from .family import FamilyParameters


class ModularAccessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = access.report()

    def test_actual_reference_ratio_and_modular_orientation(self):
        for witness in self.packet["reference_witnesses"].values():
            high, low = map(F, witness["conditional_positive_weights"])
            self.assertGreater(high, low)
            self.assertGreater(low, 0)
            self.assertEqual(high / low, F(witness["phase_witness"]["ratio_q"]))
            self.assertEqual(witness["phase_witness"]["target_pauli_y_expectation"], "1")
            self.assertEqual(witness["phase_witness"]["trace_distance_lower_bound"], "1/2")

    def test_realness_for_every_actual_branch_on_signed_real_probe(self):
        # Independent execution on a real coherent test input, including a
        # negative off-diagonal coefficient; diagonal-only preservation
        # would not establish this stronger real-operation boundary.
        i, j = self.packet["reference_witnesses"]["C"]["ordered_source_indices_i_j"]
        real_probe = {(i, 0, i, 0): F(1, 2), (j, 0, j, 0): F(1, 2),
                      (i, 0, j, 0): -F(1, 2), (j, 0, i, 0): -F(1, 2)}
        for c_weight in (F(0), F(1, 2), F(1)):
            params = FamilyParameters(c_weight=c_weight, repair_probability=F(3, 4))
            for branch in observer.branches((), params):
                output = observer.apply_branch(real_probe, branch)
                self.assertTrue(all(isinstance(value, F) for value in output.values()))
                mass = observer.trace(output)
                if mass:
                    self.assertTrue(all(isinstance(value / mass, F) for value in output.values()))
                self.assertTrue(all(output.get((a, 0, b, 0), 0) == output.get((b, 0, a, 0), 0)
                                    for a, _, b, _ in output))

    def test_adaptive_continuation_preserves_realness(self):
        i, j = self.packet["reference_witnesses"]["C"]["ordered_source_indices_i_j"]
        probe = {(i, 0, i, 0): F(1, 2), (j, 0, j, 0): F(1, 2),
                 (i, 0, j, 0): F(1, 2), (j, 0, i, 0): F(1, 2)}
        blocks = {(): probe}
        params = FamilyParameters(c_weight=F(1, 2), repair_probability=F(3, 4))
        for _ in range(3):
            blocks = observer.advance_marginal(blocks, params)
        self.assertEqual(observer.trace(observer.summed_density(blocks)), 1)
        self.assertTrue(all(isinstance(p, F) for density in blocks.values() for p in density.values()))

    def test_complex_phase_service_breaks_the_invariant(self):
        half = sp.Rational(1, 2)
        probe = sp.Matrix([[half, half], [half, half]])
        y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
        phase = sp.diag(1, sp.I)
        target = phase * probe * phase.adjoint()
        self.assertEqual(sp.trace(y * target), 1)
        self.assertEqual(sp.trace(y * probe), 0)
        self.assertEqual((target - sp.eye(2) / 2).eigenvals(), {-half: 1, half: 1})

    def test_tracial_pair_has_no_phase_witness(self):
        with self.assertRaises(ValueError):
            access.phase_control(F(1))

    def test_complete_receipt_and_forged_error_bound(self):
        access.check_receipt(self.packet, access.HERE / "modular_access_results.json")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "forged.json"
            changed = json.loads(json.dumps(self.packet))
            changed["reference_witnesses"]["C"]["phase_witness"]["trace_distance_lower_bound"] = "0"
            path.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, "differs from exact replay"):
                access.check_receipt(self.packet, path)


if __name__ == "__main__":
    unittest.main()
