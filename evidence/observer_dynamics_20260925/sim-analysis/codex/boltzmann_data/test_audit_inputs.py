"""Checks that catch calibration, column-order, and parameter-mapping mistakes."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("codex_planck_input_audit", HERE / "audit_inputs.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class PlanckInputAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = audit.make_receipt()

    def test_download_hashes_and_deterministic_receipt(self):
        self.assertEqual(len(self.receipt["inputs"]), 9)
        self.assertEqual(self.receipt, json.loads((HERE / "input_audit.json").read_text()))

    def test_mass_fraction_and_neutrino_density_match_upstream(self):
        p = self.receipt["parameter_mapping"]
        self.assertEqual(p["set_cosmology_equivalent_values"]["YHe"], .2454006)
        self.assertNotEqual(p["set_cosmology_equivalent_values"]["YHe"], p["helium"]["unused_yhe_placeholder"])
        self.assertNotEqual(p["set_cosmology_equivalent_values"]["YHe"], p["helium"]["BBN_nucleon_fraction_not_CAMB_mass_fraction"])
        self.assertEqual(p["neutrino_density_target"], .0006451439)
        self.assertEqual(p["set_cosmology_equivalent_values"]["TCMB"], 2.7255)
        self.assertEqual(p["primordial"]["r"], 0)
        # Printed logA has six decimals, hence up to 5e-7 fractional error in exp.
        self.assertLess(abs(p["primordial"]["As_from_printed_logA"]/p["primordial"]["As_upstream_ini"]-1), 6e-7)

    def test_calibration_direction(self):
        self.assertEqual(audit.coadded_theory(400, 2), 100)
        self.assertAlmostEqual(self.receipt["parameter_mapping"]["raw_CMB_theory_to_coadded_factor"], 1/1.000442**2)
        with self.assertRaises(ValueError):
            audit.coadded_theory(1, 0)

    def test_theory_order_and_polarization_schema(self):
        table = audit.load_table(audit.OBSERVATION / f"{audit.MODEL_STEM}-theory_R3.01.txt", 6)
        self.assertEqual(table[0][1:4], [1016.73, 2.61753, .0308827])
        for spectrum in ["TE", "EE"]:
            p = self.receipt["spectra_schemas"][spectrum]
            self.assertEqual(p["full"]["rows"], 1995)
            self.assertEqual(p["full"]["ell_max"], 1996)
            self.assertEqual(p["binned"]["rows"], 66)
            self.assertEqual(p["binned"]["columns"][-1], "BestFit")
            self.assertFalse(p["binned"]["contains_bin_edges"])

    def test_signed_residuals_use_error_toward_theory(self):
        # Data=10, lower error=2, upper error=4. Predictions above and below.
        data = [[30, 10, 2, 4], [31, 10, 2, 4]]
        result = audit.diagonal_summary(data, [14, 8])
        self.assertEqual(result["sum_squared_standardized_residuals"], 2)
        self.assertEqual(result["mean_signed_standardized_residual"], 0)
        self.assertFalse(result["is_official_likelihood"])
        with self.assertRaises(ValueError):
            audit.diagonal_summary([[30, 1, 0, 1]], [1])

    def test_bin_center_sampling_fails_to_reproduce_official_binning(self):
        # This regression prevents replacing missing window functions silently.
        diagnostic = self.receipt["diagnostics"]["TT"]["bin_center_interpolation_is_not_actual_binning"]
        self.assertGreater(diagnostic["max_abs_difference_from_official_BestFit_in_sigma"], 1.7)
        self.assertLess(diagnostic["max_abs_difference_from_official_BestFit_in_sigma"], 1.71)

    def test_reject_malformed_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad.txt"
            p.write_text("2 1 2\n3 1 2 3\n")
            with self.assertRaises(ValueError):
                audit.load_table(p, 4)
            p.write_text("2 1 2 3\n2 1 2 3\n")
            with self.assertRaises(ValueError):
                audit.load_table(p, 4)
            p.write_text("2 nan 2 3\n")
            with self.assertRaises(ValueError):
                audit.load_table(p, 4)


if __name__ == "__main__":
    unittest.main()
