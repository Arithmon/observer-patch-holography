from fractions import Fraction as F
import json
import unittest

import numpy as np

from . import bounds, source, experiment


class ConstraintsTests(unittest.TestCase):
    def test_random_unitaries_both_orders_and_reference_states(self):
        rng = np.random.default_rng(4371)
        for dimension in (2, 3, 5):
            for _ in range(12):
                words = [np.linalg.qr(rng.normal(size=(dimension, dimension))
                         +1j*rng.normal(size=(dimension, dimension)))[0] for _ in range(2)]
                state = rng.normal(size=(dimension, 3))+1j*rng.normal(size=(dimension, 3))
                rho = state@state.conj().T
                rho /= np.trace(rho)
                u, v = words
                a, b = (np.trace(rho@word) for word in words)
                radius2 = (1-abs(a)**2)*(1-abs(b)**2)
                for word in (u@v, v@u):
                    self.assertLessEqual(abs(np.trace(rho@word)-a*b)**2, radius2+1e-13)

    def test_taylor_enclosures_at_extremes(self):
        # Binary64 is a diagnostic; exact zero and narrow rational widths are also checked.
        for x in (F(0), F(1, 100), F(-3, 2), F(25), F(-25)):
            for cosine in (False, True):
                lo, hi = bounds.trig_interval(x, cosine)
                self.assertLessEqual(lo, hi)
                self.assertLess(float(hi-lo), 1e-25)
                target = (np.cos if cosine else np.sin)(float(x))
                self.assertAlmostEqual(float(lo), target, places=14)
        self.assertEqual(bounds.trig_interval(F(0)), (F(0), F(0)))
        self.assertEqual(bounds.trig_interval(F(0), True), (F(1), F(1)))
        for bad in (1., F(26), True):
            with self.assertRaises(ValueError):
                bounds.trig_interval(bad)

    def test_certificate_exact_negative_and_positive_controls(self):
        a, b = (F(4, 5), F(0)), (F(4, 5), F(0))
        bad = bounds.witness_floor(a, b, (F(-1), F(0)), F(0))
        self.assertGreater(F(bad['complex_entry_correction_lower_bound']), 0)
        good = bounds.witness_floor(a, b, (F(16, 25), F(0)), F(0))
        self.assertEqual(F(good['complex_entry_correction_lower_bound']), 0)
        with self.assertRaises(ValueError):
            bounds.witness_floor((F(1), F(0)), b, a, F(0))

    def test_axes_only_split_and_all_orders(self):
        axes, interior, _ = experiment.exact_data()
        self.assertEqual((len(axes), len(interior)), (21, 100))
        for extent, expected_count in ((4, 1), (5, 3)):
            prediction = experiment.predictions(axes, extent)
            result = experiment.score(axes, interior, prediction)
            self.assertEqual(result['order_checks'], 200)
            self.assertEqual(result['certified_incompatible_plugin_Grams'], expected_count)
            self.assertEqual({row['product_order_sign'] for row in result['checks']}, {-1, 1})
            # A changed held-out response changes only scoring, never the prediction function.
            changed = dict(interior)
            changed[1, 1] = (F(0), F(0))
            self.assertEqual(experiment.predictions(axes, extent), prediction)
            self.assertNotEqual(experiment.score(axes, changed, prediction), result)
        with self.assertRaises(ValueError):
            experiment.predictions(axes, 3)
        bad = dict(axes)
        bad[1, 0] = (F(11, 10), F(0))
        with self.assertRaises(ValueError):
            experiment.predictions(bad, 4)

    def test_source_compilation_product_order_and_complete_readout(self):
        words = source.source_words()
        matrices = {k: source.compiler.execute(v[0], 1) for k, v in words.items()}
        np.testing.assert_allclose(matrices['UV'], matrices['U']@matrices['V'], atol=2e-14)
        np.testing.assert_allclose(matrices['VU'], matrices['V']@matrices['U'], atol=2e-14)
        self.assertGreater(np.linalg.norm(matrices['UV']-matrices['VU']), .1)
        result = source.execute()
        self.assertEqual(len(result['instrument_runs']), 8)
        self.assertLess(result['maximum_source_Ramsey_mean_residual'], 1e-12)
        self.assertLess(result['maximum_complete_mass_defect'], 1e-12)
        self.assertEqual(result['maximum_reference_checkpoint_residual'], 0.)
        self.assertLess(result['UV_disk_distance'], result['radius'])
        self.assertLess(result['VU_disk_distance'], result['radius'])

    def test_retained_receipt_replays(self):
        with (experiment.FOLDER/'results.json').open() as stream:
            expected = json.load(stream)
        experiment.compare(expected, experiment.run())
        expected['scope']['OPH_specific_confirmation'] = True
        with self.assertRaises(ValueError):
            experiment.compare(expected, experiment.run())


if __name__ == '__main__':
    unittest.main()
