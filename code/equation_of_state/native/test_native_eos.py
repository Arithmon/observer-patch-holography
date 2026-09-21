"""Adversarial tests reject mathematically false but rehashed native receipts."""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("native_eos_independent_verifier",
                                              HERE / "verify_native_eos.py")
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


class NativeWorkCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.valid = json.loads((HERE / "native_eos_receipt.json").read_text())

    def reject(self, mutate):
        changed = copy.deepcopy(self.valid)
        mutate(changed)
        payload = {key: value for key, value in changed.items() if key != "payload_sha256"}
        changed["payload_sha256"] = verifier.digest(payload)
        result = verifier.verify(changed)
        self.assertFalse(result["receipt"], result)
        self.assertTrue(result["errors"])

    def test_independent_complete_receipt(self):
        result = verifier.verify(self.valid)
        self.assertTrue(result["receipt"], result)
        self.assertEqual(result["exact_events_replayed"], 240)
        self.assertEqual(result["candidate_work_jets_checked"], 40)
        self.assertIsNone(result["native_thermodynamic_w"])

    def test_full_native_w_promotion(self):
        self.reject(lambda r: r["claim_boundary"].update(native_thermodynamic_w="1/2"))

    def test_quadratic_heat_promotion(self):
        self.reject(lambda r: r["claim_boundary"].update(quadratic_drop_identified_as_heat=True))

    def test_equilibrium_promotion(self):
        self.reject(lambda r: r["claim_boundary"].update(thermodynamic_equilibrium_established=True))

    def test_momentum_no_go_promotion(self):
        self.reject(lambda r: r["claim_boundary"].update(physical_momentum_no_go=True))

    def test_missing_work_assumption(self):
        self.reject(lambda r: r["work_extension"].update(work_law_supplied=False))

    def test_boolean_scope_cannot_be_numeric(self):
        self.reject(lambda r: r["claim_boundary"].update(native_physical_energy_selected=0))

    def test_exact_topology_rejects_nonedge(self):
        self.reject(lambda r: r["seams"].__setitem__(0, [0, 2]))

    def test_omitted_control(self):
        self.reject(lambda r: r["cases"].pop())

    def test_early_stopping(self):
        self.reject(lambda r: r["cases"][0]["events"].pop())

    def test_changed_retained_readback(self):
        self.reject(lambda r: r["cases"][0]["events"][0]["readback_before"].__setitem__(0, "11"))

    def test_changed_update_even_if_load_unchanged(self):
        def mutate(r):
            event = r["cases"][0]["events"][0]
            event["state_after"][0:2] = ["5", "7"]
            event["quadratic_after"] = "37"
            event["quadratic_drop"] = "35"
            event["cumulative_quadratic_drop"] = "35"
        self.reject(mutate)

    def test_lost_dissipation_record(self):
        self.reject(lambda r: r["cases"][0]["events"][0].update(cumulative_quadratic_drop="0"))

    def test_coherent_but_wrong_pressure_jet(self):
        def mutate(r):
            row = r["cases"][0]["reference_work_jets_initial"][2]
            # Preserve p=-dH/dV and p/rho=w, but falsify the supplied alpha derivative.
            row.update(candidate_pressure="12", volume_derivative="-12", candidate_w="1/2")
        self.reject(mutate)

    def test_pressure_sign(self):
        self.reject(lambda r: r["cases"][0]["reference_work_jets_initial"][2]
                    .update(candidate_pressure="-8", candidate_w="-1/3"))

    def test_energy_offset_ratio_rejects_hidden_fixed_w(self):
        self.reject(lambda r: r["energy_offset_control"].update(w_after="1/3"))

    def test_direction_moment_cannot_be_declared_conserved(self):
        self.reject(lambda r: r["linear_direction_readout_control"]
                    .update(candidate_direction_moment_after="-12"))

    def test_missing_scope_key(self):
        self.reject(lambda r: r["claim_boundary"].pop("all_OPH_work_laws_excluded"))

    def test_noncanonical_fraction(self):
        self.reject(lambda r: r["cases"][0]["reference_work_jets_initial"][0]
                    .update(candidate_energy="144/2"))

    def test_invalid_hash_rejected(self):
        changed = copy.deepcopy(self.valid)
        changed["payload_sha256"] = "0" * 64
        self.assertFalse(verifier.verify(changed)["receipt"])

    def test_invalid_source_bytes_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for relative in self.valid["provenance"]["source_sha256"]:
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("changed source\n")
            result = verifier.verify(self.valid, root)
        self.assertFalse(result["receipt"])
        self.assertTrue(any("source_bytes_changed" in error for error in result["errors"]))

    def test_malformed_payloads_fail_closed(self):
        for malformed in (None, [], "payload", {}, {"schema": []}):
            with self.subTest(payload=malformed):
                self.assertFalse(verifier.verify(malformed)["receipt"])


if __name__ == "__main__":
    unittest.main()
