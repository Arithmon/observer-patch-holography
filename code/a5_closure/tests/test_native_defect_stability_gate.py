"""Exact topology and repair-scope controls for the native defect gate."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest
from unittest.mock import patch

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import flux_defect_criterion_certificate as flux  # noqa: E402
import native_defect_stability_gate as gate  # noqa: E402
from model.observer import native  # noqa: E402


class NativeDefectStabilityGateTests(unittest.TestCase):
    def test_closed_one_and_two_puncture_homology(self):
        support, punctures = gate.support_and_punctures()
        self.assertEqual((len(support["edges"]), len(support["face_rows"])),
                         (30, 20))
        self.assertEqual(punctures, (0, 13))
        ranks = gate.homology_ranks(support, punctures)
        for prime in ("F2", "F3"):
            self.assertEqual(ranks[prime], {
                "boundary_one_rank": 11,
                "full_boundary_two_rank": 19,
                "one_puncture_boundary_two_rank": 19,
                "punctured_boundary_two_rank": 18,
                "full_H1_dimension": 0,
                "one_puncture_H1_dimension": 0,
                "two_puncture_H1_dimension": 1,
            })
        refined, child_faces = flux.edge_midpoint_refinement(support)
        child_punctures = (child_faces[punctures[0]][0],
                           child_faces[punctures[1]][0])
        self.assertEqual((refined["vertex_count"], len(refined["edges"]),
                          len(refined["face_rows"])), (42, 120, 80))
        for row in gate.homology_ranks(refined, child_punctures).values():
            self.assertEqual((row["full_H1_dimension"],
                              row["one_puncture_H1_dimension"],
                              row["two_puncture_H1_dimension"]), (0, 0, 1))

    def test_declared_pair_is_flat_off_cores_and_removable_under_link_updates(self):
        support, punctures = gate.support_and_punctures()
        path = flux.dual_path(support, *punctures)
        alternative = flux.dual_path(support, *punctures,
                                     avoid=frozenset({path[0][1]}))
        self.assertNotEqual(path, alternative)
        for charge in range(1, 6):
            first = gate.verify_local_link_path(support, punctures, path, charge)
            second = gate.verify_local_link_path(support, punctures,
                                                 alternative, charge)
            holonomy = flux.face_holonomies(first, support["boundary"], 20)
            self.assertEqual([holonomy[f] for f in punctures],
                             [charge, (-charge) % 6])
            self.assertEqual(sum(bool(value) for i, value in enumerate(holonomy)
                                 if i not in punctures), 0)
            difference = {edge: (first.get(edge, 0) - second.get(edge, 0)) % 6
                          for edge in range(30)}
            gauge = gate.flat_gauge_potential(support, difference)
            self.assertEqual(flux.coboundary_chain(gauge, support["edges"]),
                             difference)

    def test_false_dual_step_and_nonflat_regauge_fail(self):
        support, punctures = gate.support_and_punctures()
        path = flux.dual_path(support, *punctures)
        wrong = list(path)
        wrong[0] = (wrong[0][0], path[-1][1])
        with self.assertRaisesRegex(ArithmeticError, "adjacent seam update"):
            gate.verify_local_link_path(support, punctures, wrong, 1)
        with self.assertRaisesRegex(ValueError, "nonzero face curvature"):
            gate.flat_gauge_potential(support, {path[0][1]: 1})

    def test_native_record_gate_rejects_one_changed_current(self):
        result = gate.run()
        self.assertEqual(result["native_record"]["coin_current_branches_checked"],
                         78 * 30 * 2)
        self.assertFalse(result["native_record"]["current_selected_nontrivial_Z6_defect"])
        self.assertFalse(result["declared_Z6_link_control"]
                         ["single_seam_link_updates_are_native_repair"])
        self.assertFalse(result["one_step_topology_refinement"]
                         ["native_repair_on_child_complex"])
        self.assertEqual(result["one_step_topology_refinement"]
                         ["declared_child_flux_menu"], list(range(6)))
        self.assertFalse(result["one_step_topology_refinement"]
                         ["coarse_link_connection_transport_checked"])
        original = native.native
        loads = native.states()[0]
        seam = gate.support_and_punctures()[0]["edges"][0]

        def changed(source, selected_seam, coin):
            endpoint, record = original(source, selected_seam, coin)
            if source == loads and selected_seam == seam and coin == 0:
                return endpoint, (coin, record[1] + 1)
            return endpoint, record

        with patch.object(native, "native", side_effect=changed):
            with self.assertRaisesRegex(ArithmeticError,
                                        "native current generated a face flux"):
                gate.run()


if __name__ == "__main__":
    unittest.main()
