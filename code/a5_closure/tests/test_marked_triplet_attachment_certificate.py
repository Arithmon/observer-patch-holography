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
        cls.packet = verifier.load(REFERENCE)

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

    def test_same_character_does_not_allow_substitution_of_the_source_mark(self):
        bad = deepcopy(self.packet)
        old = bad['source_marked_order_five_action']
        inverse = [old.index(i) for i in range(12)]
        self.assertNotEqual(inverse, old)
        index = bad['actions'].index(inverse)
        # The inverse has the same character. A character-only gate would
        # miss this alteration of the independently committed marker.
        self.assertEqual(bad['characters'][index], bad['marked_traces'])
        bad['source_marked_order_five_action'] = inverse
        bad['marked_action_index'] = index
        self.rejected(bad)

    def test_incidence_automorphism_is_not_an_outer_twist(self):
        bad = deepcopy(self.packet)
        bad['outer_twist_example'] = list(range(12))
        self.rejected(bad)

    def test_swapping_response_triplets_does_not_relabel_the_live_spatial_band(self):
        response = verifier.load(ROOT/verifier.RESPONSE)
        for row in response['closed_paths']:
            row['observed_response'].reverse()
        with self.assertRaisesRegex(ValueError, 'upper action'):
            verifier.verify(self.packet, response=response)

    def test_physical_selection_and_universal_nonselection_promotions_fail(self):
        for key, value in (
                ('physical_golden_embedding_selected', True),
                ('scalar_clock_rescaling_changes_representation', True),
                ('scope', 'all possible source laws and every physical embedding')):
            with self.subTest(key=key):
                bad = deepcopy(self.packet)
                bad[key] = value
                self.rejected(bad)

    def test_hom_dimensions_require_original_integers(self):
        for name, correct in self.packet['fixed_action_hom_dimensions'].items():
            for invalid in (bool(correct), float(correct)):
                with self.subTest(name=name, value=invalid):
                    bad = deepcopy(self.packet)
                    bad['fixed_action_hom_dimensions'][name] = invalid
                    self.rejected(bad)

    def test_duplicate_receipt_keys_are_rejected_before_replay(self):
        text = json.dumps(self.packet)
        # A permissive parser accepts the final value and conceals the
        # conflicting earlier physical-promotion claim from the verifier.
        duplicated = text.replace('"physical_golden_embedding_selected": false',
                                  '"physical_golden_embedding_selected": true, '
                                  '"physical_golden_embedding_selected": false')
        self.assertNotEqual(duplicated, text)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'counterfeit.json'
            path.write_text(duplicated)
            with self.assertRaisesRegex(ValueError, 'duplicate JSON key'):
                verifier.load(path)

    def test_duplicate_nested_source_keys_and_nonfinite_values_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'source.json'
            for text in ('{"source":{"generator":0,"generator":1}}',
                         '{"value":NaN}', '{"value":Infinity}', '{"value":-Infinity}'):
                with self.subTest(text=text):
                    path.write_text(text)
                    with self.assertRaises(ValueError):
                        verifier.load(path)

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
