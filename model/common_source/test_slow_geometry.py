"""Bind the existing exact response carrier to executed native/quantum records."""
from copy import deepcopy
from fractions import Fraction as F
from collections import Counter
from itertools import product
import math
import unittest
from unittest.mock import patch

import numpy as np

from ..observer import native
from . import native_join
from . import slow_geometry as g
from . import null_tomography as tomography
from .test_native_join import endpoint


def numerical(vector):
    return np.array([float(F(a))+float(F(b))*math.sqrt(5) for a, b in vector])


class SlowGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.joint = native_join.run(d=2, steps=16)
        cls.geometry = g.report(cls.joint)

    def test_all_native_branches_prove_the_actual_uniform_mean(self):
        cal = self.geometry['calibration']
        self.assertEqual(cal['primitive_branches_checked_exactly'], 78*30*2)
        self.assertEqual(cal['native_state_count'], 78)
        mean = [[F(x) for x in row] for row in cal['mean_operator_exact']]
        for state in native.states():
            outcomes = [endpoint(state, seam, coin)[0]
                        for seam in cal['seams'] for coin in (0, 1)]
            independent = [F(sum(row[i] for row in outcomes), 60) for i in range(12)]
            self.assertEqual(independent,
                             [sum(mean[i][j]*state[j] for j in range(12)) for i in range(12)])

    def test_existing_registered_lean_chart_is_checked_not_assumed_equal(self):
        binding = g.registered_binding()
        self.assertEqual(binding['registered_gram_entries_checked_exactly'], 144)
        self.assertEqual(binding['source_to_registered_lean_port'],
                         [8, 10, 1, 3, 7, 11, 0, 4, 6, 9, 2, 5])
        self.assertTrue(binding['all_pass'])
        frame = np.array([numerical(row) for row in g.raw_frame()])
        normalized = frame @ frame.T / ((1+math.sqrt(5))/2+2)
        np.testing.assert_allclose(normalized, g.source().intrinsic_gram(), atol=2e-15)

    def test_three_dimensional_band_is_exact_but_finite_kernels_are_not_rank_three(self):
        cal = self.geometry['calibration']
        self.assertEqual(cal['response_carrier_rank_exact'], 3)
        self.assertEqual(cal['three_seam_raw_determinant_Qsqrt5'], ['4', '0'])
        distances = [row['normalized_kernel_distance_to_gram_frobenius']
                     for row in cal['finite_kernel_diagnostics']]
        self.assertTrue(all(a > b for a, b in zip(distances, distances[1:])))
        for n in (0, 16):
            self.assertEqual(np.linalg.matrix_rank(g.source().normalized_response_kernel(n), tol=1e-12), 11)
        self.assertTrue(all(row['finite_nonconstant_rank_exact'] == 11
                            for row in cal['finite_kernel_diagnostics']))

    def test_relative_readout_uses_actual_committed_currents_and_keeps_quantum_vectors(self):
        attached = self.geometry['joint_readouts']
        self.assertEqual(attached['path_count'], 1056)
        self.assertEqual(attached['quantum_histories'], self.joint['quantum_histories'])
        frame = np.array([numerical(row) for row in g.raw_frame()])
        for old, row in zip(self.joint['paths'], attached['paths'], strict=True):
            self.assertEqual({key: row[key] for key in old}, old)
            initial = np.array(native_join.ALL_STATES[row['initial_source']])
            final = np.array(row['final_native_loads'])
            readout = row['carrier_readout']
            displacement = numerical(readout['total_raw_record_displacement_Qsqrt5'])
            np.testing.assert_allclose(displacement, (final-initial)@frame, atol=2e-15)
            np.testing.assert_allclose(numerical(readout['total_raw_centroid_displacement_Qsqrt5']),
                                       displacement/2, atol=2e-15)
        self.assertEqual(attached['actual_displacement_span_rank'], 2)
        self.assertFalse(attached['uniform_calibration_schedule_used_by_production'])

    def test_new_binary_continuation_executes_a_genuine_third_direction(self):
        new = self.geometry['three_direction_continuation']
        self.assertEqual(new['actual_displacement_span_rank'], 3)
        self.assertEqual(new['three_direction_raw_determinant_Qsqrt5'], ['4', '0'])
        self.assertEqual(new['path_count'], 1056)
        displacements = []
        for old, row in zip(self.joint['paths'], new['paths'], strict=True):
            for name in ('initial_source', 'after_first_source', 'first_native_record', 'quantum_record',
                         'quantum_branch_reference', 'history_probability', 'terminal_density_prefactor'):
                self.assertEqual(row[name], old[name])
            seam = (0, 5) if row['quantum_record'][1] == 0 else (0, 7)
            after = native_join.ALL_STATES[row['after_first_source']]
            final, record = endpoint(after, seam, row['second_native_record'][0])
            self.assertEqual(row['second_seam'], list(seam))
            self.assertEqual(row['second_native_record'], list(record))
            self.assertEqual(row['final_native_loads'], list(final))
            if row['history_probability'] > 0:
                displacements.append(numerical(row['carrier_readout']['total_raw_record_displacement_Qsqrt5']))
        self.assertEqual(np.linalg.matrix_rank(displacements, tol=1e-12), 3)
        self.assertAlmostEqual(sum(row['history_probability'] for row in new['paths']), 1, places=13)

    def test_seam_orientation_and_signed_current_do_not_change_the_record_readout(self):
        for seam in g.source().seams():
            for current in (-1, 0, 1):
                self.assertEqual(g.current_displacement(seam, current),
                                 g.current_displacement(tuple(reversed(seam)), -current))

    def test_hidden_load_coordinates_are_not_an_injective_history_encoding(self):
        zero = (g.source().q5(0),)*3
        antipodal = [state for state in native.states(binary=True) if g.load_readout(state) == zero]
        self.assertEqual(len(antipodal), 6)
        self.assertEqual(len({tuple(row['first_native_record']) for row in self.joint['paths']}), 4)

    def test_detached_path_mutations_are_rejected(self):
        for mutation in ('weight', 'nan', 'reference', 'seam', 'current', 'duplicate', 'missing', 'bool'):
            bad = deepcopy(self.joint)
            row = bad['paths'][0]
            if mutation == 'weight': row['history_probability'] += .01
            elif mutation == 'nan': row['history_probability'] = math.nan
            elif mutation == 'reference': row['quantum_branch_reference']['current'] = 10
            elif mutation == 'seam': row['second_seam'] = [0, 7]
            elif mutation == 'current': row['first_native_record'][1] += 1
            elif mutation == 'duplicate': bad['paths'][1] = deepcopy(row)
            elif mutation == 'missing': bad['paths'].pop()
            elif mutation == 'bool': row['first_native_record'][0] = False
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                g.attach_joint(bad)

    def test_invalid_loads_and_non_native_seams_cannot_be_read(self):
        for loads in ([0]*12, [2]*12, [True, 1]+[0]*10, [-1, 3]+[0]*10, [2]+[0]*10):
            with self.subTest(loads=loads), self.assertRaises(ValueError):
                g.load_readout(loads)
        for seam, current in (((0, 2), 1), ((0, 0), 0), ((0, 1), True), ((0, 1), 2)):
            with self.subTest(seam=seam, current=current), self.assertRaises(ValueError):
                g.current_displacement(seam, current)

    def test_report_is_json_and_explicit_about_physical_scope(self):
        import json
        json.dumps(self.geometry, allow_nan=False)
        self.assertTrue(self.geometry['all_pass'])
        self.assertFalse(self.geometry['scope']['physical_position_identified'])
        self.assertFalse(self.geometry['scope']['physical_momentum_or_ccr_constructed'])
        self.assertFalse(self.geometry['scope']['physical_length_or_time_scale_calibrated'])


class NativeNullTomographyTests(unittest.TestCase):
    """Direct quadratic evaluation and controls for the conditional tensor map."""

    @staticmethod
    def charges(components):
        c = g.source()
        tensor = [[c.q5(0) for _ in range(4)] for _ in range(4)]
        for (i, j), value in zip(tomography.TENSOR_COMPONENTS, components, strict=True):
            tensor[i][j] = tensor[j][i] = value
        result = []
        for u in g.raw_frame():
            k, value = (c.q5(1),)+u, c.q5(0)
            for i in range(4):
                for j in range(4):
                    value = c.q5_add(value, c.q5_mul(k[i], c.q5_mul(tensor[i][j], k[j])))
            result.append(value)
        return result

    def test_exact_rank_kernel_and_native_band_certificate(self):
        import json
        result = tomography.null_tomography_certificate()
        json.dumps(result, allow_nan=False)
        self.assertTrue(result['all_pass'])
        self.assertEqual(result['exact_design_rank'], 9)
        self.assertEqual(result['exact_kernel_dimension'], 1)
        self.assertEqual(result['exact_charge_compatibility_constraints'], 3)
        self.assertEqual(result['antipodal_merged_rank'], 6)
        self.assertEqual(result['antipodal_merged_kernel_dimension'], 4)
        self.assertFalse(result['scope']['physical_null_directions_derived'])
        self.assertFalse(result['scope']['directional_stress_or_geometry_measurements_constructed'])

    def test_all_ten_independent_tensor_components_reconstruct_modulo_metric(self):
        c = g.source()
        inverse_r2 = c.q5_inv(c.q5_add(c.GOLDEN, c.q5(2)))
        for active in range(10):
            for value in (c.q5(1), c.q5(F(2, 7), F(-3, 5))):
                tensor = [value if i == active else c.q5(0) for i in range(10)]
                readings = self.charges(tensor)
                decoded = tomography.decode_null_charges(readings)
                expected = list(tensor)
                expected[0] = c.q5(0)
                for i in (4, 7, 9):
                    expected[i] = c.q5_add(expected[i], c.q5_mul(tensor[0], inverse_r2))
                with self.subTest(active=active, value=value):
                    self.assertEqual(decoded, tuple(expected))
                    self.assertEqual(self.charges(decoded), readings)

    def test_arbitrary_port_readings_are_rejected_not_projected_into_a_tensor(self):
        c = g.source()
        # A single-port response has a nonzero fast-band contribution.
        for active in range(12):
            readings = [c.q5(int(i == active)) for i in range(12)]
            with self.subTest(active=active), self.assertRaisesRegex(ValueError, 'incompatible'):
                tomography.decode_null_charges(readings)
        readings = self.charges([c.q5(i-3, F(i, 11)) for i in range(10)])
        self.assertEqual(self.charges(tomography.decode_null_charges(readings)), readings)
        readings[4] = c.q5_add(readings[4], c.q5(F(1, 1000000)))
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            tomography.decode_null_charges(readings)

    def test_antipodal_merging_erases_three_nonmetric_tensors(self):
        c = g.source()
        frame = g.raw_frame()
        partner = [next(j for j, v in enumerate(frame) if g._add(u, v) == (c.q5(0),)*3)
                   for u in frame]
        for active in (1, 2, 3):
            tensor = [c.q5(int(i == active)) for i in range(10)]
            readings = self.charges(tensor)
            self.assertTrue(any(q != c.q5(0) for q in readings))
            merged = [c.q5_mul(c.q5(F(1, 2)), c.q5_add(readings[i], readings[partner[i]]))
                      for i in range(12)]
            self.assertEqual(merged, [c.q5(0)]*12)
            self.assertEqual(tomography.decode_null_charges(merged), (c.q5(0),)*10)
            self.assertEqual(tomography.decode_null_charges(readings), tuple(tensor))

    def test_changed_frame_or_incorrect_decoder_cannot_pass_certificate(self):
        c = g.source()
        frame = list(g.raw_frame())
        frame[0] = (c.q5(0),)*3
        with patch.object(tomography, 'raw_frame', return_value=tuple(frame)):
            with self.assertRaisesRegex(AssertionError, 'native response Gram'):
                tomography.null_tomography_certificate()
        decoder = tomography._null_tomography_decoder()
        decoder[0] = [c.q5(0)]*12
        with patch.object(tomography, '_null_tomography_decoder', return_value=decoder):
            with self.assertRaisesRegex(AssertionError, 'left inverse'):
                tomography.null_tomography_certificate()

    def test_stability_bounds_are_exact_and_attainable(self):
        c = g.source()
        decoder = tomography._null_tomography_decoder()
        expected = c.q5(1, F(-1, 5))
        self.assertEqual(tomography._max_row_sum(decoder), expected)
        certificate = tomography.null_tomography_certificate()
        self.assertEqual(certificate['compatibility_residual_max_error_per_charge_max_error_Qsqrt5'],
                         ['1/2', '1/2'])
        attained = False
        for row in decoder:
            error = [c.q5(c.q5_sign(x)) for x in row]
            output = g._dot(row, error)
            self.assertLessEqual(c.q5_sign(c.q5_sub(output, expected)), 0)
            attained = attained or output == expected
        self.assertTrue(attained)

    def test_readings_require_exact_finite_field_values(self):
        c = g.source()
        invalid = ([c.q5(0)]*11, [c.q5(0)]*13, [0]*12,
                   [(0.0, 0)]*12, [(False, 0)]*12, [('0', '0')]*12)
        for readings in invalid:
            with self.subTest(readings=readings), self.assertRaises(ValueError):
                tomography.decode_null_charges(readings)


class NativeRecordTensorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = tomography.native_record_tensor_certificate()
        cls.initial = tuple(int(i in (0, 3)) for i in range(12))

    def test_actual_record_kernel_and_exact_withheld_predictions(self):
        self.assertEqual(self.result['native_local_record_branches_checked'], 3960)
        self.assertEqual(self.result['binary_preparations_exactly_compatible'], 6)
        self.assertEqual(self.result['fast_decay_per_uniform_attempt_Qsqrt5'], ['11/12', '-1/60'])
        self.assertEqual(self.result['withheld_error_amplification_Qsqrt5'], ['2', '2'])
        self.assertTrue(self.result['rational_exact_compatibility_requires_antipodal_even_mean'])
        self.assertFalse(self.result['scope']['physical_stress_tensor_or_lorentz_lift_derived'])
        for row in self.result['examples']:
            errors = row['withheld_errors_Qsqrt5']
            self.assertEqual(all(x == ['0', '0'] for x in errors), row['preparation'] == 'antipodal')
            self.assertEqual(sum(F(p) for p in row['endpoint_readout_law']), 1)
            if row['native_attempts'] == 64:
                bound = tuple(F(v) for v in row['uniform_residual_bound_Qsqrt5'])
                self.assertLess(g.source().q5_sign(g.source().q5_sub(bound, g.source().q5(F(1, 10000)))), 0)

    def test_complete_two_event_records_reproduce_joint_endpoint_law(self):
        menu = tuple((seam, coin) for seam in g.source().seams() for coin in (0, 1))
        counts, histories = Counter(), set()
        for word in product(menu, repeat=2):
            record = tomography.execute_record_word(self.initial, word)
            histories.add(record['events'])
            # Independent endpoint oracle: no reuse of the signed-current formula.
            endpoint = list(self.initial)
            for seam, coin in word:
                a, b = seam
                low, high = sum(endpoint[i] for i in seam)//2, (sum(endpoint[i] for i in seam)+1)//2
                endpoint[a], endpoint[b] = (low, high) if coin == 0 else (high, low)
            self.assertEqual(record['readings'], tuple(endpoint))
            for i, tape in enumerate(record['local_records']):
                self.assertTrue(all(i in seam for seam, _ in tape))
                self.assertEqual(tomography.local_occupancy(record['initial_bits'][i], i, tape), endpoint[i])
            counts[record['readings']] += 1
        self.assertEqual(len(histories), 3600)
        law, _ = tomography.native_record_law(self.initial, 2)
        self.assertEqual(law, tuple(F(counts[state], 3600) for state in native.states(binary=True)))

    def test_mean_compatibility_does_not_assert_pathwise_compatibility(self):
        c = g.source()
        law, mean = tomography.native_record_law(self.initial, 1)
        tomography.decode_null_charges([c.q5(x) for x in mean])
        incompatible_mass = F(0)
        for state, p in zip(native.states(binary=True), law, strict=True):
            try:
                tomography.decode_null_charges([c.q5(x) for x in state])
            except ValueError:
                incompatible_mass += p
        self.assertEqual(incompatible_mass, F(1, 6))

    def test_nine_port_fit_cannot_see_a_changed_withheld_reading(self):
        c = g.source()
        _, mean = tomography.native_record_law(self.initial, 16)
        readings = [c.q5(x) for x in mean]
        prediction = tomography.withheld_charge_predictions([readings[i] for i in tomography.FIT_PORTS])
        self.assertEqual(prediction, tuple(readings[i] for i in tomography.WITHHELD_PORTS))
        readings[3] = c.q5_add(readings[3], c.q5(F(1, 1000)))
        self.assertEqual(prediction, tomography.withheld_charge_predictions([readings[i] for i in tomography.FIT_PORTS]))
        self.assertNotEqual(prediction[0], readings[3])
        with self.assertRaises(ValueError):
            tomography.decode_null_charges(readings)

    def test_local_readout_contract_rejects_impossible_or_remote_records(self):
        for initial, port, records in ((None, 0, []), (True, 0, []), (0, True, []),
                                       (0, 0, [((0, 1), 1)]), (1, 0, [((0, 1), -1)]),
                                       (0, 0, [((7, 11), 0)]), (0, 0, [((0, 1), False)])):
            with self.subTest(initial=initial, port=port, records=records), self.assertRaises(ValueError):
                tomography.local_occupancy(initial, port, records)
        self.assertEqual(tomography.local_occupancy(1, 0, [((0, 1), 1)]),
                         tomography.local_occupancy(1, 0, [((1, 0), -1)]))

    def test_sampling_bound_is_for_independent_complete_trials(self):
        self.assertEqual(tomography.independent_history_sampling_bound(3000000, F(1, 100)), F(1, 100))
        self.assertEqual(tomography.independent_history_sampling_bound(1, F(1, 100)), 1)
        for n, error in ((0, 1), (True, 1), (2, F(0)), (2, .01)):
            with self.assertRaises(ValueError):
                tomography.independent_history_sampling_bound(n, error)


if __name__ == '__main__':
    unittest.main()
