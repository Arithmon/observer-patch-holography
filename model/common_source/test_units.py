"""Exact scale freedom and conditional SI conversion checks."""
from fractions import Fraction as F
import json
import unittest

from . import units


class UnitsTests(unittest.TestCase):
    def test_exact_constraint_ranks_and_nullspace_freedoms(self):
        for names, rank in (([], 0), (['c'], 1), (['c', 'hbar'], 2),
                            (['c', 'hbar', 'time_calibration'], 3)):
            result = units.constraints(names)
            self.assertEqual(result['rank_exact'], rank)
            self.assertEqual(result['independent_scale_freedoms'], 3-rank)
            for vector in result['nullspace_exponent_basis_exact']:
                values = list(map(F, vector))
                self.assertTrue(any(values))
                for row in result['exponent_rows']:
                    self.assertEqual(sum(x*y for x, y in zip(row, values, strict=True)), 0)
        self.assertEqual(units.constraints(['c', 'hbar'])['nullspace_exponent_basis_exact'], [['1', '1', '-1']])

    def test_distinct_time_choices_preserve_exact_c_and_hbar_but_change_scales(self):
        examples = units.report()['conditional_examples']
        for row in examples:
            t = F(row['external_time_scale_seconds_exact'])
            length = F(row['length_scale_metres_exact'])
            energy_pi = F(row['energy_scale_joules_coefficient_of_inverse_pi_exact'])
            self.assertEqual(length/t, F(299792458))
            self.assertEqual(energy_pi*t, F('6.62607015e-34')/2)
            self.assertFalse(row['external_time_scale_selected_by_model'])
        for key, ratio in (('external_time_scale_seconds_exact', F(1000)),
                           ('length_scale_metres_exact', F(1000)),
                           ('energy_scale_joules_coefficient_of_inverse_pi_exact', F(1, 1000))):
            self.assertEqual(F(examples[1][key])/F(examples[0][key]), ratio)

    def test_unit_choice_is_not_a_physical_calibration_or_prediction(self):
        result = units.report()
        json.dumps(result, allow_nan=False)
        self.assertFalse(result['scope']['emitted_model_propagation_speed_derived'])
        self.assertFalse(result['scope']['spatial_readout_linked_to_modular_clock'])
        self.assertFalse(result['scope']['dimensionless_physical_predictions_established_by_unit_choice'])
        self.assertTrue(result['sources']['metre_definition'].startswith('https://www.bipm.org/'))

    def test_invalid_or_nonexact_time_calibrations_are_rejected(self):
        for time in (True, 0, -1, 1e-15, '1e-15'):
            with self.subTest(time=time), self.assertRaises(ValueError):
                units.conditional_si_example(time)
        for names in (['c', 'c'], ['mass'], 'c'):
            with self.subTest(names=names), self.assertRaises(ValueError):
                units.constraints(names)


if __name__ == '__main__':
    unittest.main()
