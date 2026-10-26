"""Source records, packet-blind inverse, and common-clock regression controls."""
from fractions import Fraction as F
import copy
import unittest
from unittest.mock import patch

import numpy as np

from . import calibrated_spectrum_observer as observer


class ObserverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report, cls.arrays = observer.execute()

    def test_nine_running_action_parameters_recovered_from_simulated_records(self):
        rows = self.report['observer_read_running_spectrum']
        self.assertEqual(len(rows), 9)
        for row in rows:
            lower, upper = map(F, row['record_gap_over_common_reference_interval'])
            self.assertLess(lower, upper)
            self.assertTrue(row['calibrated_action_inside_record_interval'])
        self.assertFalse(self.report['scope']['calibration_recovery_is_independent_empirical_confirmation'])

    def test_inverse_replay_has_no_packet_or_source_access(self):
        expected = [row['record_gap_over_common_reference_interval'] for row in self.report['observer_read_running_spectrum']]
        with patch.object(observer.matter, 'parameters', side_effect=AssertionError('packet read')), \
             patch.object(observer.matter, 'mass_coordinates', side_effect=AssertionError('mass read')), \
             patch.object(observer.frame_source, 'source_frame', side_effect=AssertionError('source read')):
            result = observer.decode_records(self.arrays['complete_simulated_record_counts'])
        self.assertEqual([list(map(str, result['mass_over_reference_intervals'][label]))
                          for label in observer.RECORD_LABELS[:-1]], expected)
        np.testing.assert_array_equal(result['analyzer_offsets'], self.arrays['analyzer_offsets'])

    def test_complete_counts_and_timeout_records(self):
        counts = self.arrays['complete_simulated_record_counts']
        self.assertEqual(counts.shape, (21, 2, 32, 512, 4))
        np.testing.assert_array_equal(counts.sum(axis=(2, 3, 4)), np.full((21, 2), observer.TRIALS_PER_SETTING))
        self.assertTrue(np.all(counts[:, :, 30:].sum(axis=(2, 3, 4)) > 0))
        self.assertFalse(self.report['records']['success_postselection_or_probability_renormalization'])

    def test_quantum_columns_and_failed_outputs_retained(self):
        self.assertEqual(self.arrays['field_unitaries'].shape, (21, 9, 4, 4))
        self.assertEqual(self.arrays['success_read_matrices'].shape, (21, 2, 9, 2, 4, 4))
        self.assertEqual(self.arrays['failed_read_matrices'].shape, (21, 2, 9, 2, 4, 4))
        self.assertEqual(self.arrays['clock_reference_outputs'].shape, (21, 4, 4, 4))
        self.assertEqual(self.arrays['failed_clock_reference_outputs'].shape, (21, 4, 4, 4))
        for value in self.arrays.values():
            if value.dtype.kind != 'U':
                self.assertTrue(np.isfinite(value).all())

    def test_simultaneous_conditional_sampling_and_resources(self):
        sampling = self.report['inversion']['simultaneous_sampling_certificate']
        self.assertEqual(sampling['simultaneous_real_statistics'], 420)
        self.assertLess(F(sampling['familywise_error_upper']), F(1, 100))
        self.assertLessEqual(F(sampling['width_exact']), observer.RECORD_RADIUS)
        resources = self.report['source_resources']
        self.assertGreater(resources['total_native_pulses_per_complete_setting_schedule'], 0)
        self.assertEqual(resources['proper_code_load_export_service_events'],
                         2*resources['total_native_pulses_per_complete_setting_schedule'])

    def test_high_winding_pair_keeps_all_column_phase(self):
        angle = (F(10**12)+F(1, 7),)*2
        gate, reduced, winding, error = observer._source_pair(angle)
        self.assertGreater(abs(winding), 10**10)
        self.assertLess(error, observer.NUMERIC_SCORE_ADMISSION)
        self.assertLess(observer._operator_error(gate.fock_unitary, observer._pair_boxes(reduced)), observer.NUMERIC_SCORE_ADMISSION)
        self.assertGreater(observer._operator_error(-gate.fock_unitary, observer._pair_boxes(reduced)), F(1))

    def test_bad_record_schemas_and_absent_visibility_refuse(self):
        counts = self.arrays['complete_simulated_record_counts']
        for bad in (counts.astype(float), np.ma.array(counts), counts[:, :, :31], counts.tolist()):
            with self.assertRaises(ValueError):
                observer.decode_records(bad)
        bad = np.zeros((1, 2, 32, 512, 4), dtype=np.int64)
        bad[:, :, 31, 0, 0] = observer.TRIALS_PER_SETTING
        with self.assertRaisesRegex(ValueError, 'unresolved'):
            observer.decode_records(bad)

    def test_complete_output_checker_rejects_phase_zero_and_schema_forgery(self):
        for kind in ('phase', 'zero', 'extra'):
            arrays = dict(self.arrays)
            if kind == 'extra':
                arrays['uncontrolled_output'] = np.zeros(1)
            else:
                arrays['success_read_matrices'] = arrays['success_read_matrices'].copy()
                if kind == 'phase':
                    arrays['success_read_matrices'][0, 0, 0] *= -1
                else:
                    arrays['success_read_matrices'][0, 0, 0] = 0
            with self.assertRaises(ValueError):
                observer.compare_execution(self.report, arrays, self.report, self.arrays)
        forged = copy.deepcopy(self.report)
        forged['inversion']['hardware_or_finite_shot_coverage_authenticated_here'] = True
        with self.assertRaises(ValueError):
            observer.compare_execution(forged, self.arrays, self.report, self.arrays)


if __name__ == '__main__':
    unittest.main()
