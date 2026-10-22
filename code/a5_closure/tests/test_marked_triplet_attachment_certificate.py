"""Hostile controls for source marking and strict receipt verification."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import verify_marked_triplet_attachment_independent as verifier


ROOT = Path(__file__).resolve().parents[3]
REFERENCE = ROOT/'code/a5_closure/manifests/marked_triplet_attachment_reference.json'


class MarkedTripletAttachmentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = json.loads(REFERENCE.read_text())

    def rejected(self, packet):
        with self.assertRaises(ValueError):
            verifier.verify(packet)

    def test_reference_replays(self):
        result = verifier.verify(self.packet)
        self.assertTrue(result['verified'])

    def test_source_selected_generator_is_bound_before_character(self):
        bad = deepcopy(self.packet)
        bad['source_marked_order_five_action'] = list(range(12))
        self.rejected(bad)

    def test_marked_index_cannot_be_selected_by_expected_trace(self):
        bad = deepcopy(self.packet)
        bad['marked_action_index'] = (bad['marked_action_index'] + 1) % 60
        self.rejected(bad)

    def test_trace_is_recomputed_for_source_generator(self):
        bad = deepcopy(self.packet)
        bad['marked_traces'][0] = ['999', '0']
        self.rejected(bad)

    def test_boolean_does_not_pass_as_integer_count(self):
        bad = deepcopy(self.packet)
        bad['incidence_recharting_count'] = True
        self.rejected(bad)

    def test_boolean_does_not_pass_as_permutation_entry(self):
        bad = deepcopy(self.packet)
        bad['source_marked_order_five_action'][0] = True
        self.rejected(bad)

    def test_verifier_reconstructs_refreshed_raw_frame_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in verifier.INPUTS:
                target = root/name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ROOT/name, target)
            source = root/'model/common_source/slow_geometry.py'
            contents = source.read_text()
            contents = contents.replace('(minus, phi, zero), (one, phi, zero)',
                                        '(minus, negative_phi, zero), (one, phi, zero)', 1)
            source.write_text(contents)
            bad = deepcopy(self.packet)
            bad['source_sha256'] = {
                name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                for name in verifier.INPUTS
            }
            with self.assertRaises(ValueError):
                verifier.verify(bad, root=root)


if __name__ == '__main__':
    unittest.main()
