"""Fixed reference modular time, actual-state dynamics and source attachment."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math
import unittest

import numpy as np
from scipy.linalg import expm

from . import experiment as e, operational as o, scalar_coupling as scalar
from . import thermal_clock as t


def vector(row):
    values = np.asarray(row['unnormalized_terminal_state'])
    return values[:, 0]+1j*values[:, 1]


class ThermalClockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = t.report(16)

    def test_reference_is_actual_H_gibbs_state_without_smoothing(self):
        for d in (2, 4):
            ref = t.reference(d)
            expected = expm(-ref.hamiltonian)
            expected /= np.trace(expected)
            np.testing.assert_allclose(ref.density, expected, atol=2e-14, rtol=0)
            self.assertGreater(min(ref.probabilities), 0)
            self.assertAlmostEqual(sum(ref.probabilities), 1, places=13)
            np.testing.assert_allclose(ref.log_probabilities,
                                       -ref.beta*ref.energies-ref.log_partition, atol=2e-14)
            np.testing.assert_allclose((ref.basis*ref.probabilities)@ref.basis.conj().T,
                                       ref.density, atol=2e-14)

    def test_modular_identity_has_the_correct_sign_on_general_complex_operators(self):
        rng = np.random.default_rng(871)
        for d in (2, 4):
            ref = t.reference(d, beta=.8)
            n = len(ref.density)
            a = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
            s = .37
            time = -ref.beta*s
            expected = expm(1j*time*ref.hamiltonian)@a@expm(-1j*time*ref.hamiltonian)
            np.testing.assert_allclose(t.modular_flow(ref, a, s), expected, atol=2e-12, rtol=0)
            np.testing.assert_allclose(o.modular_flow(ref.density, a, s), expected, atol=2e-11, rtol=0)
            wrong = expm(-1j*time*ref.hamiltonian)@a@expm(1j*time*ref.hamiltonian)
            self.assertGreater(np.linalg.norm(expected-wrong), 1.)

    def test_fixed_reference_is_stationary_while_complete_actual_state_moves(self):
        for row in self.result['levels']:
            result = row['thermal_clock']
            self.assertGreater(result['actual_state_commutator_with_H_hilbert_schmidt'], .1)
            settings = result['settings']
            self.assertTrue(all(s['reference_stationarity_operator_error'] < 2e-13 for s in settings))
            self.assertGreater(settings[-1]['actual_state_displacement_hilbert_schmidt'], .1)
            values = [s['reference_gauge_one_probability'] for s in settings]
            self.assertLess(max(values)-min(values), 2e-13)
            actual = [s['source_clock_bit_marginal'][1] for s in settings]
            self.assertGreater(max(actual)-min(actual), .1)
            self.assertFalse(result['scope']['fixed_reference_is_updated_from_changing_actual_state'])

    def test_actual_source_states_are_evolved_instead_of_replaced_by_reference(self):
        for row in self.result['levels']:
            d, original = row['dimension'], row['actual_history']
            actual = o.terminal_vectors(original)
            pointer, _ = e.pointer(d, 0)
            for setting in row['thermal_clock']['settings']:
                time = setting['model_coordinate_time']
                source = (scalar.source_unitary(d, 16, time=time)[0] if time else np.eye(4*d*d))
                self.assertEqual(len(setting['branches']), 8)
                for branch in setting['branches']:
                    expected = pointer[branch['clock_bit']]@source@actual[tuple(branch['old_record'])]
                    np.testing.assert_allclose(vector(branch), expected, atol=2e-14, rtol=0)
                    self.assertAlmostEqual(branch['probability'], np.linalg.norm(expected)**2, places=13)
                self.assertLess(setting['old_record_marginal_error'], 2e-12)
                self.assertGreater(setting['appended_source_resources_per_realized_history']['native_pulses'], 0)

    def test_declared_reference_record_weights_cancel_and_all_record_projectors_survive(self):
        ref = t.reference(2)
        n = len(ref.density)
        rng = np.random.default_rng(872)
        observables = {record: rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)) for record in t.RECORDS}
        first = t.block_modular_flow(ref, observables, .5)
        second = t.block_modular_flow(ref, observables, .5, (F(4, 10), F(3, 10), F(2, 10), F(1, 10)))
        for record in t.RECORDS:
            np.testing.assert_allclose(first[record], second[record], atol=2e-13)
            np.testing.assert_allclose(first[record], t.modular_flow(ref, observables[record], .5), atol=2e-13)
        for selected in t.RECORDS:
            projector = {record: np.eye(n) if record == selected else np.zeros((n, n)) for record in t.RECORDS}
            after = t.block_modular_flow(ref, projector, -.9)
            for record in t.RECORDS:
                np.testing.assert_allclose(after[record], projector[record], atol=2e-13)

    def test_tracial_reference_cannot_select_nontrivial_H_and_reference_first_keeps_input(self):
        for row in self.result['levels']:
            audit = row['thermal_clock']['reference_selection_audit']
            self.assertGreater(audit['hamiltonian_spectral_spread'], 1)
            self.assertGreater(audit['centered_hamiltonian_hilbert_schmidt_norm'], 1)
            self.assertGreater(audit['gibbs_minus_tracial_density_operator_norm'], .1)
            self.assertGreater(audit['gauge_record_pinching_reference_displacement_operator_norm'], .01)
            self.assertGreater(audit['conditional_gauge_control_reference_displacement_operator_norm'], .01)
            self.assertLess(audit['tracial_reference_modular_flow_error_from_identity'], 2e-13)
            self.assertLess(audit['reference_first_hamiltonian_reconstruction_error'], 2e-13)
            self.assertTrue(audit['constructed_from_supplied_hamiltonian'])
            self.assertFalse(audit['generated_by_existing_repair'])
            self.assertFalse(audit['equivalent_reference_first_choice_removes_physical_input'])

    def test_ideal_split_certificate_and_measured_floating_defect_remain_distinct(self):
        for row in self.result['levels']:
            for setting in row['thermal_clock']['settings']:
                errors = setting['error_budget']
                bound = float(F(errors['ideal_temporal_certificate']['ideal_split_operator_error_upper']))
                self.assertLessEqual(errors['independent_split_versus_exact_operator_error'], bound+2e-12)
                self.assertLessEqual(errors['source_versus_exact_modular_schrodinger_operator_error'],
                                     bound+errors['source_minus_independent_split_operator_error']+2e-12)
                self.assertFalse(errors['measured_execution_defect_is_a_rigorous_roundoff_certificate'])
                self.assertFalse(errors['physical_service_or_SI_clock_error_included'])

    def test_appended_reference_does_not_change_input_history_or_old_local_obstruction(self):
        original = self.result['levels'][0]['actual_history']
        before = deepcopy(original)
        result = t.receipt(original)
        self.assertEqual(original, before)
        self.assertEqual(result['actual_history_branches'], original['branches'])
        local = o.observation_receipt(original)
        self.assertTrue(local['checks']['source_time_is_not_this_stationary_modular_flow'])
        self.assertFalse(result['scope']['earlier_reduced_state_clock_or_negative_result_changed'])
        self.assertFalse(result['scope']['physical_time_or_temperature_calibrated'])

    def test_beta_controls_relative_parameterization_without_fitting_the_state(self):
        for beta in (.5, 1., 2.):
            ref = t.reference(2, beta)
            observable = e.projector(16, 0, 1)
            time = .4
            expected = expm(1j*time*ref.hamiltonian)@observable@expm(-1j*time*ref.hamiltonian)
            np.testing.assert_allclose(t.modular_flow(ref, observable, -time/beta), expected, atol=2e-13)
        self.assertGreater(np.linalg.norm(t.reference(2, .5).density-t.reference(2, 2.).density), .1)

    def test_invalid_reference_inputs_are_rejected_without_a_noise_floor(self):
        for d, beta in ((3, 1), (True, 1), (2, 0), (2, -1), (2, True), (2, math.nan), (2, math.inf), (2, 1e308)):
            with self.subTest(d=d, beta=beta), self.assertRaises(ValueError):
                t.reference(d, beta)
        ref = t.reference(2)
        blocks = {record: np.eye(16) for record in t.RECORDS}
        for weights in ((0, F(1, 3), F(1, 3), F(1, 3)), (1, 1, 1, 1), (1,), (True, 0, 0, 0)):
            with self.subTest(weights=weights), self.assertRaises(ValueError):
                t.block_modular_flow(ref, blocks, .1, weights)
        for parameter in (True, math.nan, math.inf):
            with self.subTest(parameter=parameter), self.assertRaises(ValueError):
                t.modular_flow(ref, np.eye(16), parameter)

    def test_full_report_is_json_safe_and_scoped_to_fixed_finite_unitary_segments(self):
        json.dumps(self.result, allow_nan=False)
        for row in self.result['levels']:
            result = row['thermal_clock']
            self.assertTrue(result['all_pass'])
            self.assertTrue(result['scope']['identity_covers_unitary_segments_not_measurement_feedback_maps'])
            self.assertFalse(result['scope']['reference_prepared_or_substituted_for_actual_state'])
            self.assertFalse(result['scope']['repair_selection_of_thermal_reference_proved'])


if __name__ == '__main__':
    unittest.main()
