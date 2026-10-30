"""Independent feasibility and exact-matrix adversaries for particle attachments."""
from fractions import Fraction as F
import itertools
import unittest
from unittest.mock import patch

import numpy as np
import sympy as sp
from scipy.optimize import linprog

from model.full_fock_field.field import FieldGeometry, FullFock
from model.particle_comparison import compare as c


class ExactSpectrumTests(unittest.TestCase):
    def test_exact_source_certificate_and_independent_numerical_spectrum(self):
        certificate, frequencies = c.certify_source()
        self.assertEqual(certificate["source_dimension"],16)
        self.assertEqual(len(frequencies),8)
        matrix = FullFock(4).hamiltonian(FieldGeometry((2,1,1))).toarray()
        expected = sorted(float(e.sympy()) for row in c.spectrum() for e in row)
        np.testing.assert_allclose(np.linalg.eigvalsh(matrix),expected,atol=2e-15,rtol=0)

    def test_source_mutation_cannot_reuse_analytic_spectrum(self):
        original = FullFock.hamiltonian
        def changed(space,geometry):
            matrix = original(space,geometry).tolil()
            matrix[0,0] += 1/16
            return matrix.tocsr()
        with patch.object(FullFock,"hamiltonian",changed), self.assertRaises(ArithmeticError):
            c.certify_source()

    def test_number_nonconserving_source_is_rejected(self):
        original = FullFock.hamiltonian
        def changed(space,geometry):
            matrix = original(space,geometry).tolil()
            matrix[0,1] = matrix[1,0] = 1/16
            return matrix.tocsr()
        with patch.object(FullFock,"hamiltonian",changed), self.assertRaises(ArithmeticError):
            c.certify_source()

    def test_no_tolerance_makes_numerical_dust_an_exact_pole(self):
        original = FullFock.hamiltonian
        def changed(space,geometry):
            matrix = original(space,geometry).tolil()
            matrix[0,0] = 1e-17
            return matrix.tocsr()
        with patch.object(FullFock,"hamiltonian",changed), self.assertRaises(ValueError):
            c.certify_source()

    def test_exact_order_and_inverse(self):
        values = [c.Q(F(a,3),F(b,7)) for a in range(-4,5) for b in range(-4,5)]
        for value in values:
            self.assertEqual(value.sign(),int(sp.sign(value.sympy())))
            if value.sign():
                self.assertEqual(value/value,c.Q(1))
        for invalid in (True,.1,float("nan"),"1/2"):
            with self.assertRaises(ValueError):
                c.Q(invalid)


class AssignmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _,cls.frequencies = c.certify_source()
        cls.windows = ((F(206),F(207)),(F(3470),F(3485)))
        cls.rows = c.affine_cases(cls.frequencies,cls.windows)

    @staticmethod
    def independent_linear_feasibility(frequencies,indices,se,windows):
        a,b,d = [float(frequencies[i].sympy()) for i in indices]
        sm,st = np.sign(b-a),np.sign(d-a)
        # Unknowns kappa and offset=-kappa*mu, without eliminating either.
        result = linprog([0,0], A_eq=[[a,1]],b_eq=[se],
            A_ub=[[sm*b,sm],[-sm*b,-sm],[st*d,st],[-st*d,-st]],
            b_ub=[float(windows[0][1]),-float(windows[0][0]),
                  float(windows[1][1]),-float(windows[1][0])],
            bounds=[(0,None),(None,None)],method="highs")
        return result.success

    def test_all_672_exact_cases_against_independent_linear_programs(self):
        self.assertEqual(len(self.rows),672)
        self.assertFalse(any(row["feasible"] for row in self.rows))
        for row in self.rows:
            self.assertEqual(row["feasible"], self.independent_linear_feasibility(
                self.frequencies,row["indices"],row["electron_sign"],self.windows))

    def test_free_chemical_potential_is_not_accidentally_fixed_to_zero(self):
        frequencies = [c.Q(0),c.Q(205),c.Q(3470)]
        rows = c.affine_cases(frequencies,self.windows)
        # kappa=1, mu=-1 yields masses 1,206,3471 exactly.
        selected = next(r for r in rows if r["indices"] == [0,1,2] and r["electron_sign"] == 1)
        self.assertTrue(selected["feasible"])
        self.assertTrue(self.independent_linear_feasibility(frequencies,[0,1,2],1,self.windows))

    def test_common_shift_and_positive_rescaling_preserve_exclusion(self):
        changed = [7*w+c.Q(F(123,5)) for w in self.frequencies]
        self.assertEqual([r["feasible"] for r in self.rows],
                         [r["feasible"] for r in c.affine_cases(changed,self.windows)])

    def test_all_permutations_and_signs_are_disclosed(self):
        expected = {(indices,sign) for indices in itertools.permutations(range(8),3) for sign in (-1,1)}
        self.assertEqual({(tuple(r["indices"]),r["electron_sign"]) for r in self.rows},expected)
        result = c.lepton_no_go(self.frequencies,self.windows)
        self.assertEqual(result["cases_below_tau_window"]+result["cases_above_tau_window"],672)
        self.assertEqual(result["feasible_cases"],0)

    def test_repeated_frequencies_and_bad_windows_fail(self):
        for frequencies,windows in (([c.Q(1)]*3,self.windows),
            (self.frequencies,((1,2),(3,4))),
            (self.frequencies,((206.,207),(3470,3485))),
            (self.frequencies,((3470,3485),(206,207)))):
            with self.assertRaises(ValueError):
                c.affine_cases(frequencies,windows)


class PhysicalScopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = c.build()

    def test_two_candidates_one_calibration_no_absolute_prediction(self):
        pair = self.result["two_candidate_assignment"]
        self.assertEqual(len(pair["prediction"]["candidates"]),2)
        self.assertEqual(pair["independently_predicted_absolute_masses"],0)
        self.assertEqual(pair["normalized_difference_prediction"],"0")
        self.assertTrue(pair["within_reported_conditional_bound"])
        self.assertEqual(pair["candidate_mass_energy_MeV"][0],pair["candidate_mass_energy_MeV"][1])
        self.assertTrue(pair["calibration_errors_fully_correlated"])
        self.assertFalse(pair["prediction"]["physical_rest_mass_interpretation_derived"])

    def test_retrospective_real_data_replay_not_universe_confirmation(self):
        self.assertEqual(self.result["actual_data_replay"]["ACT_posterior_lower_tail_mass"],
                         ["3114/508985"]*2)
        self.assertTrue(self.result["actual_data_replay"]["custody"]["ACT_weighted_projection_recomputed_here"])
        self.assertFalse(self.result["CP11"]["OPH_universe_membership_confirmed"])
        self.assertFalse(self.result["CP11"]["eligible_physical_dictionary"])
        self.assertEqual(self.result["CP11"]["eligible_physical_comparisons_added"],0)
        self.assertFalse(self.result["sealing_after_prior_target_exposure_makes_prospective"])

    def test_both_gap_and_equality_follow_before_calibration(self):
        prediction,gap = c.half_filled_candidates()
        self.assertEqual(gap,c.Q(F(-1,16),F(1,8)))
        self.assertEqual(prediction["ground_multiplicity"],1)
        self.assertEqual(prediction["exact_normalized_mass_difference"],"0")

    def test_calibration_cannot_override_changed_gap_or_equality(self):
        prediction,gap = c.half_filled_candidates()
        observations = c.strict_json.loads((c.ROOT/c.TARGETS).read_text())
        with self.assertRaises(ValueError):
            c.calibrated_pair(prediction,gap+c.Q(1),observations)
        prediction["exact_normalized_mass_difference"] = "1"
        with self.assertRaises(ValueError):
            c.calibrated_pair(prediction,gap,observations)


if __name__ == "__main__":
    unittest.main()
