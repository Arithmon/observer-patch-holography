"""Complete source instrument, record persistence and finite energy accounting."""
from copy import deepcopy
import json
import math
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import block_diag, expm

from ..observer.source_scalar import compiler as c, field
from . import experiment as sharp
from . import scalar_coupling as scalar
from . import smooth_protocol as smooth


def vector(row):
    state = np.array(row['unnormalized_terminal_state'])
    return state[:, 0]+1j*state[:, 1]


class SmoothProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.histories = {d: smooth.history(d, 16) for d in (2, 4)}

    def test_emitted_pointer_matches_full_analytic_kraus_with_the_circuit_phases(self):
        for d in (2, 4):
            actual, metadata = smooth.pointer(d)
            q = smooth.coordinate(d)
            values, basis = np.linalg.eigh(q)
            ideal = [np.exp(-.25j*math.pi)*(basis*np.cos(values/2-math.pi/4))@basis.T,
                     np.exp(.25j*math.pi)*(basis*np.cos(values/2+math.pi/4))@basis.T]
            for a, b in zip(actual, ideal, strict=True):
                np.testing.assert_allclose(a, b, atol=2e-13, rtol=0)
            self.assertLess(metadata['maximum_full_kraus_operator_error'], 2e-13)
            self.assertGreater(metadata['resources']['native_pulses'], 0)
            self.assertGreater(metadata['resources']['pair_code_pulses'], 0)
            self.assertEqual(metadata['resources']['proper_code_service_events'],
                             2*metadata['resources']['native_pulses'])

    def test_pointer_is_executed_from_retained_source_words(self):
        smooth._pointer.cache_clear()
        with patch.object(c, 'execute', wraps=c.execute) as execute:
            branches, metadata = smooth.pointer(2)
        execute.assert_called_once()
        self.assertEqual(execute.call_args.args[0], metadata['emitted_word'])
        self.assertEqual(execute.call_args.args[1], 5)
        raw = c.execute(metadata['emitted_word'], 5)
        unitary, _ = c.stabilize_executed_unitary(raw)
        for outcome, branch in enumerate(branches):
            np.testing.assert_allclose(branch, unitary[outcome*16:(outcome+1)*16, :16], atol=2e-14)

    def test_sine_effect_and_nonselective_channel_hold_for_coherent_inputs(self):
        rng = np.random.default_rng(708)
        for d in (2, 4):
            branches, _ = smooth.pointer(d)
            q = smooth.coordinate(d)
            minus, plus = expm(-.5j*q), expm(.5j*q)
            state = rng.normal(size=len(q))+1j*rng.normal(size=len(q))
            state /= np.linalg.norm(state)
            rho = np.outer(state, state.conj())
            density = sum(k@rho@k.conj().T for k in branches)
            expected = (minus@rho@minus.conj().T+plus@rho@plus.conj().T)/2
            np.testing.assert_allclose(density, expected, atol=2e-13, rtol=0)
            sinq = (expm(1j*q)-expm(-1j*q))/(2j)
            for b, k in enumerate(branches):
                np.testing.assert_allclose(k.conj().T@k, (np.eye(len(q))+(-1)**b*sinq)/2,
                                           atol=2e-13, rtol=0)

    def test_record_projectors_are_fixed_even_on_inputs_with_record_coherence(self):
        d, size = 2, 16
        branches, _ = smooth.pointer(d)
        second, _ = scalar.source_unitary(d, 16, time=.5)
        controls = [sharp.feedback(d, a, True)[0] for a in (0, 1)]
        continued = [block_diag(*(branch@second@control for control in controls)) for branch in branches]
        for a in (0, 1):
            record = np.diag([int(i//size == a) for i in range(2*size)])
            dual = sum(k.conj().T@record@k for k in continued)
            np.testing.assert_allclose(dual, record, atol=2e-13, rtol=0)

    def test_complete_histories_keep_all_four_states_and_actual_conditional_law(self):
        for result in self.histories.values():
            self.assertEqual({tuple(row['record']) for row in result['branches']},
                             {(0, 0), (0, 1), (1, 0), (1, 1)})
            for row in result['branches']:
                self.assertAlmostEqual(np.linalg.norm(vector(row))**2, row['probability'], places=13)
            self.assertAlmostEqual(sum(result['first_record_marginal']), 1, places=12)
            self.assertLess(result['maximum_old_record_effect_error'], 2e-12)
            self.assertLess(result['maximum_independent_terminal_vector_error'], 2e-12)
            for law in result['conditional_readout_law']:
                p, signal = law['first_record_probability'], law['unnormalized_sine_expectation']
                np.testing.assert_allclose(law['joint_probabilities'], [(p+signal)/2, (p-signal)/2], atol=2e-13)

    def test_all_energy_changes_are_actual_complete_ensemble_expectations(self):
        for d, result in self.histories.items():
            h = scalar.ideal_hamiltonian(d)
            energy = result['energy_ledger']
            total = sum(float(np.vdot(vector(row), h@vector(row)).real) for row in result['branches'])
            self.assertAlmostEqual(energy['stages'][-1]['energy'], total, places=12)
            self.assertEqual(len(energy['stages']), 7)
            for row in energy['stages']:
                self.assertAlmostEqual(row['probability_mass'], 1, places=12)
                self.assertAlmostEqual(row['energy'], sum(b['unnormalized_energy'] for b in row['branches']), places=13)
            self.assertLess(abs(energy['energy_telescope_residual']), 2e-13)
            read = energy['smooth_readout']
            for key in ('source_dual_energy_identity_residual', 'source_minus_analytic_finite_energy_residual',
                        'finite_boundary_identity_residual', 'full_dual_energy_operator_error'):
                self.assertLess(abs(read[key]), 2e-12)
            self.assertGreater(abs(read['analytic_finite_boundary_correction']), 1e-5)
            self.assertFalse(read['finite_or_conditional_injection_assumed_to_equal_one_eighth'])
            for segment in energy['evolution_drift'].values():
                self.assertLess(abs(segment['source_minus_ideal_split_energy_defect']), 2e-12)

    def test_controls_distinguish_interacting_feedback_from_disconnected_gauge_control(self):
        result = smooth.run(16)
        for level in result['levels']:
            self.assertGreater(abs(level['smooth_scalar_feedback_response']), .005)
            self.assertLess(abs(level['disconnected_scalar_feedback_response']), 2e-12)
        self.assertTrue(result['old_sharp_sign_protocol_preserved'])

    def test_sharp_history_is_unchanged_and_has_the_same_pre_detector_record(self):
        old = sharp.history(2, 16)
        before = deepcopy(old)
        new = smooth.history(2, 16)
        self.assertEqual(old, before)
        np.testing.assert_allclose(new['first_record_marginal'], old['first_record_marginal'], atol=2e-13)
        self.assertGreater(abs(new['second_record_marginal'][1]-old['second_record_marginal'][1]), .01)
        self.assertIn('positive site-0', old['meaning']['second_bit'])
        self.assertIn('sin Q0', new['meaning']['second_bit'])

    def test_all_native_current_preparations_retain_complete_source_histories(self):
        probabilities = []
        for current in (-1, 0, 1):
            result = smooth.history(2, 16, native_current=current)
            self.assertEqual(result['native_current'], current)
            self.assertTrue(result['all_pass'])
            probabilities.append(result['second_record_marginal'][0])
        self.assertGreater(max(probabilities)-min(probabilities), .01)

    def test_invalid_protocol_arguments_fail_before_execution(self):
        for kwargs in ({'d': 3}, {'d': True}, {'d': 2, 'steps': 0}, {'d': 2, 'steps': True},
                       {'d': 2, 'feedback': 1}, {'d': 2, 'native_current': True},
                       {'d': 2, 'native_current': 2}, {'d': 2, 'coupling': math.nan}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                smooth.history(**kwargs)

    def test_receipts_are_json_safe_and_do_not_claim_physical_error_certification(self):
        for result in self.histories.values():
            json.dumps(result, allow_nan=False)
            self.assertFalse(result['scope']['old_sharp_sign_protocol_changed'])
            self.assertFalse(result['scope']['physical_service_error_certified'])
            self.assertEqual(result['resources_per_realized_history']['fresh_pointer_blanks'], 2)
            self.assertEqual(result['resources_per_realized_history']['committed_classical_bits'], 2)


if __name__ == '__main__':
    unittest.main()
