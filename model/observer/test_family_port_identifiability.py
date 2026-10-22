"""Independent record and port-moment controls for the native W5 interface."""
from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import unittest

from model.observer import family_port_identifiability as subject


def direct_event(loads, seam, coin):
    a, b = seam
    low, high = divmod(loads[a] + loads[b], 2)
    endpoints = (low, low + high) if coin == 0 else (low + high, low)
    out = list(loads)
    out[a], out[b] = endpoints
    return tuple(out), (coin, loads[a] - out[a])


def direct_tape_law(loads):
    law = {}
    for first_coin in (0, 1):
        first, r1 = direct_event(loads, (0, 1), first_coin)
        seam = (0, 5) if r1[1] > 0 else (0, 1)
        for second_coin in (0, 1):
            _, r2 = direct_event(first, seam, second_coin)
            key = (r1, r2)
            law[key] = law.get(key, Q(0)) + Q(1, 4)
    return law


class FamilyPortIdentifiabilityTests(unittest.TestCase):
    def test_fixed_tape_cannot_read_distinct_port_quadrupoles(self):
        antipodes = subject.antipodes()
        self.assertEqual({tuple(sorted((i, j))) for i, j in enumerate(antipodes)},
                         {(0, 3), (1, 2), (4, 7), (5, 6), (8, 11), (9, 10)})
        antipodal = subject.binary_load((4, 7))
        nonantipodal = subject.binary_load((4, 8))
        expected_law = {((a, 0), (b, 0)): Q(1, 4)
                        for a in (0, 1) for b in (0, 1)}
        self.assertEqual(direct_tape_law(antipodal), expected_law)
        self.assertEqual(direct_tape_law(nonantipodal), expected_law)
        self.assertEqual(subject.tape_law(antipodal), expected_law)
        self.assertEqual(subject.tape_law(nonantipodal), expected_law)
        self.assertEqual(subject.w5_norm_squared(antipodal, antipodes), Q(5, 3))
        self.assertEqual(subject.w5_norm_squared(nonantipodal, antipodes), Q(2, 3))

    def test_admitted_port_labelled_probe_resolves_witness(self):
        left = subject.binary_load((4, 7))
        right = subject.binary_load((4, 8))
        self.assertEqual(subject.probe_law(left, (7, 8)),
                         {(0, 1): Q(1, 2), (1, 0): Q(1, 2)})
        self.assertEqual(subject.probe_law(right, (7, 8)),
                         {(0, 0): Q(1, 2), (1, -1): Q(1, 2)})

    def test_all_seam_currents_recover_port_loads_and_detect_corruption(self):
        for loads in subject.native.states():
            packet = subject.all_seam_packet(loads)
            self.assertEqual(subject.recover_loads_from_packet(packet), loads)
        packet = subject.all_seam_packet(subject.binary_load((4, 7)))
        packet[(7, 8)] = ((0, 2), packet[(7, 8)][1])
        with self.assertRaises(ValueError):
            subject.recover_loads_from_packet(packet)
        packet = subject.all_seam_packet(subject.binary_load((4, 7)))
        original_contrast = sum(current for _, current in packet[(7, 8)])
        packet[(7, 8)] = ((0, 2), (1, -1))
        self.assertEqual(sum(current for _, current in packet[(7, 8)]), original_contrast)
        with self.assertRaisesRegex(ValueError, "coin-conditioned record"):
            subject.recover_loads_from_packet(packet)
        del packet[(7, 8)]
        with self.assertRaises(ValueError):
            subject.recover_loads_from_packet(packet)

    def test_receipt_replay_and_scope(self):
        result = subject.run()
        self.assertTrue(result["all_pass"])
        self.assertFalse(result["source"]["mass_reference_data_consumed"])
        self.assertEqual(result["source"]["record_law_class_sizes"], [46, 10, 9, 1])
        self.assertEqual(result["port_resolved_readout"]["source_loads_reconstructed_exactly"], 78)
        self.assertEqual(result["witness"]["W5_norm_squared"], ["5/3", "2/3"])
        receipt = Path(__file__).with_name("results") / "family_port_identifiability.json"
        self.assertEqual(receipt.read_text(encoding="utf-8"), subject.canonical(result))


if __name__ == "__main__":
    unittest.main()
