"""Bind the existing exact response carrier to executed native/quantum records."""
from copy import deepcopy
from fractions import Fraction as F
import math
import unittest

import numpy as np

from ..observer import native
from . import native_join
from . import slow_geometry as g
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


if __name__ == '__main__':
    unittest.main()
