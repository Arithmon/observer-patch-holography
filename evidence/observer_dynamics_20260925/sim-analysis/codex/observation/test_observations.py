import importlib.util
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("observations", ROOT / "prepare_observations.py")
OBS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OBS)


class ObservationTests(unittest.TestCase):
    def test_coverage_and_units(self):
        data = OBS.prepare()["planck_tt_l2_40.json"]
        self.assertEqual([r["ell"] for r in data["rows"]], list(range(2, 41)))
        self.assertEqual(sum(r["primary_lowell"] for r in data["rows"]), 28)
        first = data["rows"][0]
        self.assertEqual(first["Dl_uK2"], 225.895)
        self.assertAlmostEqual(first["Cl_uK2"], 225.895 * math.pi / 3)
        self.assertNotEqual(first["Dl_error_minus_uK2"], first["Dl_error_plus_uK2"])

    def test_calibration_and_firas_unit_separation(self):
        data = OBS.prepare()
        self.assertEqual(data["planck_tt_l2_40.json"]["bestfit_calPlanck"], 1.000442)
        row = data["firas_monopole.json"]["rows"][0]
        self.assertEqual(row["monopole_MJy_per_sr"], 200.723)
        self.assertEqual(row["residual_kJy_per_sr"], 5)
        self.assertEqual(row["sigma_kJy_per_sr"], 14)
        self.assertAlmostEqual(row["frequency_GHz"], 2.27 * 29.9792458)
        self.assertEqual(len(data["firas_monopole.json"]["rows"]), 43)


if __name__ == "__main__":
    unittest.main()
