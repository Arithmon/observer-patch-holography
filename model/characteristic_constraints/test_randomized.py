from fractions import Fraction as F
import unittest
from dataclasses import replace

import numpy as np

from . import bounds, experiment, randomized as r, randomized_source as s


class RandomizedTests(unittest.TestCase):
    def test_linear_coefficients_reproduce_the_complex_gram(self):
        vector = ((F(-3, 2), F(1, 3)), (F(1), F(-1, 7)), (F(2, 3), F(3, 5)))
        means = ((F(1, 4), F(1, 7)), (F(-1, 3), F(2, 9)), (F(1, 5), F(-1, 2)))
        for phase in (F(0), F(4, 25), F(-1, 4), F(25)):
            design = r.witness_design(vector, phase)
            v = np.array([complex(float(a), float(b)) for a, b in vector])
            expected = (v.conj()@bounds.gram(*means, phase)@v).real
            self.assertAlmostEqual(float(design.value(tuple(x for z in means for x in z))), expected, places=13)
            self.assertEqual(sum(design.probabilities), 1)
            self.assertLess(float(design.algebra_error), 1e-25)

    def test_importance_scores_are_unbiased_for_every_state(self):
        rng = np.random.default_rng(18752287)
        for _ in range(30):
            vector = tuple((F(int(x)), F(int(y))) for x, y in rng.integers(-4, 5, size=(3, 2)))
            try:
                design = r.witness_design(vector)
            except ValueError:
                continue
            means = tuple(F(int(x), 100) for x in rng.integers(-100, 101, size=6))
            average = sum((p*((1+m)*design.score(i, 1)+(1-m)*design.score(i, -1))/2
                           for i, (p, m) in enumerate(zip(design.probabilities, means)) if p), F(0))
            self.assertEqual(average, design.value(means))
            self.assertEqual(max(abs(design.score(i, b)-design.constant)
                for i, p in enumerate(design.probabilities) if p for b in (-1, 1)), design.width)

    def test_services_and_fixed_horizon_threshold_fail_closed(self):
        design = r.witness_design(s.VECTOR)
        trials = 100
        total = trials*(design.constant-design.width)
        good = r.decision(design, total, trials, service_errors=(F(0),)*6)
        self.assertTrue(good['rejects_joint_premises'])
        missing = r.decision(design, total, trials, service_errors=(None,)*6)
        self.assertEqual(missing['status'], 'unbounded_services')
        self.assertFalse(missing['rejects_joint_premises'])
        noisy = r.decision(design, total, trials, service_errors=(F(1),)*6)
        self.assertFalse(noisy['rejects_joint_premises'])
        for total_bad in (total-F(1), total+F(1)):
            with self.assertRaises(ValueError):
                r.decision(design, total_bad, trials, service_errors=(F(0),)*6)
        for alpha in (F(0), F(1), .01, True):
            with self.assertRaises(ValueError):
                r.radius(design, trials, alpha)
        with self.assertRaises(ValueError):
            design.score(4, 1)  # zero-probability source setting
        with self.assertRaises(ValueError):
            design.score(0, True)
        self.assertIsNone(r.crossing_budget(design, F(-1), service_errors=(None,)*6))
        self.assertIsNone(r.crossing_budget(design, F(0), service_errors=(F(0),)*6))
        with self.assertRaises(ValueError):
            r.crossing_budget(design, F(-100), service_errors=(F(0),)*6)
        with self.assertRaises(ValueError):
            r.Design(F(0), (F(-1),)*6)
        with self.assertRaises(ValueError):
            replace(design, constant=F(0))
        with self.assertRaises(ValueError):
            replace(design, coefficients=(F(-1),)*6)
        phased = r.witness_design(s.VECTOR, F(1, 4))
        with self.assertRaises(ValueError):
            replace(phased, algebra_error=F(0))

    def test_uniform_bound_and_correlated_reference_execution(self):
        result = s.execute()
        self.assertEqual(len(result['complete_two_trial_law']), 100)
        self.assertEqual(result['maximum_reference_checkpoint_residual'], 0.)
        self.assertLess(result['raw_complete_mass_defect'], 1e-12)
        self.assertLess(result['finite_array_certificate']['conditional_score_defect_upper_float'], 1e-12)
        self.assertEqual(F(result['finite_array_certificate']['completed_randomizer_score_error_upper']), 0)
        self.assertGreater(result['minimum_conditional_second_score'], 0)
        self.assertLess(result['setting_correlated_preparation_negative_control_mean'], -2.)
        rows = result['conditional_second_scores']
        self.assertGreater(max(row['second_score_conditional_mean'] for row in rows)
                           -min(row['second_score_conditional_mean'] for row in rows), .01)
        self.assertFalse(result['preparation_reset_between_trials'])

    def test_operator_certificate_detects_a_wrong_product_probe(self):
        design = r.witness_design(s.VECTOR)
        words = s.source.source_words()
        matrices = {name: s.source.compiler.execute(word, 1) for name, (word, _) in words.items()}
        operators = [s.source.ramsey(words[name][1], phase)[0]
                     for name in ('U', 'V', 'UV') for phase in (0., -np.pi/2)]
        good = s.finite_certificate(design, matrices, operators)
        operators[-1] = tuple(reversed(operators[-1]))
        bad = s.finite_certificate(design, matrices, operators)
        self.assertLess(F(good['conditional_score_defect_upper']), F(1, 10**12))
        self.assertGreater(F(bad['conditional_score_defect_upper']), F(1))
        with self.assertRaises(ValueError):
            s.finite_certificate(design, matrices, operators[:-1])

    def test_portable_certificate_rechecks_ulp_perturbed_arrays_exactly(self):
        design = r.witness_design(s.VECTOR)
        words = s.source.source_words()
        matrices = {name: s.source.compiler.execute(word, 1) for name, (word, _) in words.items()}
        operators = [s.source.ramsey(words[name][1], phase)[0]
                     for name in ('U', 'V', 'UV') for phase in (0., -np.pi/2)]
        roots = [np.sqrt(float(p)) for p in design.probabilities if p]
        baseline = s.execution_certificate(design, matrices, operators, roots)
        diagnostic_key = 'per_execution_floating_diagnostics_not_certificate_inputs'
        stable = {key: value for key, value in baseline.items() if key != diagnostic_key}
        for direction in (np.inf, -np.inf):
            def perturb(array):
                # Simulate low-level evaluation differences in every component.
                real = np.nextafter(np.array(array.real), direction)
                imag = np.nextafter(np.array(array.imag), direction)
                return real+1j*imag
            changed_matrices = {name: perturb(value) for name, value in matrices.items()}
            changed_operators = [tuple(perturb(value) for value in maps) for maps in operators]
            changed_roots = [np.nextafter(value, direction) for value in roots]
            result = s.execution_certificate(design, changed_matrices, changed_operators, changed_roots)
            self.assertEqual({key: value for key, value in result.items() if key != diagnostic_key}, stable)
            self.assertNotEqual(result[diagnostic_key], baseline[diagnostic_key])
            self.assertEqual(F(result['conditional_score_defect_upper']), F(1, 10**13))
        # Errors exceeding a budget are not swallowed by floating replay tolerance.
        changed = [tuple(k.copy() for k in maps) for maps in operators]
        changed[-1][0][0, 0] += 1e-11
        with self.assertRaisesRegex(ArithmeticError, 'exceeds the portable'):
            s.execution_certificate(design, matrices, changed, roots)
        altered_roots = list(roots)
        altered_roots[0] = np.nextafter(altered_roots[0], np.inf)
        with self.assertRaisesRegex(ArithmeticError, 'exceeds the portable'):
            s.execution_certificate(design, matrices, operators, altered_roots)

    def test_public_means_select_future_design_not_a_statistical_test(self):
        axes, interior, _ = experiment.exact_data()
        scored = {str(extent): experiment.score(axes, interior, experiment.predictions(axes, extent))
                  for extent in (4, 5)}
        designs = r.empirical_designs(scored, axes, interior)
        for row in designs.values():
            self.assertLess(F(row['public_plugin_score']), 0)
            self.assertGreater(row['future_zero_service_crossing_trials'], 1000)
            self.assertFalse(row['public_data_statistical_rejection'])
            self.assertFalse(row['power_guaranteed'])
            spec = row['design']
            design = r.witness_design(tuple(tuple(F(x) for x in z) for z in spec['witness_vector']), F(spec['Weyl_phase']))
            trials = row['future_zero_service_crossing_trials']
            gap = -F(row['public_plugin_score'])-design.algebra_error
            self.assertLess(r.radius(design, trials, F(1, 200)), gap)
            self.assertGreaterEqual(r.radius(design, trials-1, F(1, 200)), gap)


if __name__ == '__main__':
    unittest.main()
