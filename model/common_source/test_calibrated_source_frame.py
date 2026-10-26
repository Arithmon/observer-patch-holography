"""Independent phase, exterior-sector and actual decoded-wire controls."""
from dataclasses import replace
import math

import numpy as np
import pytest

from . import calibrated_source_frame as frame


@pytest.fixture(scope="module")
def executed():
    return frame.source_frame()


def test_complete_eight_column_source_exterior_lift_and_vacuum_phase(executed):
    source, target = executed["CKM_fock_unitary"], executed["CKM_fock_target"]
    np.testing.assert_allclose(source, target, atol=2e-10, rtol=0)
    assert abs(source[0, 0]-1) < 2e-10
    assert abs(source[7, 7]-np.linalg.det(target[np.ix_((1, 2, 4), (1, 2, 4))])) < 2e-10
    assert not executed["diagnostics"]["global_phase_removed_or_fitted"]
    for mask in range(8):
        for other in range(8):
            if mask.bit_count() != other.bit_count():
                assert abs(source[mask, other]) < 2e-10
    with pytest.raises(ValueError):
        source.flags.writeable = True


def test_independent_two_particle_minors_include_fermion_signs(executed):
    source = executed["CKM_fock_unitary"]
    u = np.asarray(frame.matter.ckm_matrix().evalf(20), complex)
    masks = (3, 5, 6)
    subsets = ((0, 1), (0, 2), (1, 2))
    for i, left in enumerate(subsets):
        for j, right in enumerate(subsets):
            assert source[masks[i], masks[j]] == pytest.approx(np.linalg.det(u[np.ix_(left, right)]), abs=2e-10)
    # A bosonic/permutation-only pair lift would miss these signs.
    assert abs(source[3, 5]+u[2, 1].conjugate()) < 2e-10


def test_actual_full90_mapping_has_six_disjoint_CKM_blocks_and_singlet_signs(executed):
    mapping = executed["mapping"]
    assert len(mapping["CKM_blocks"]) == 6
    assert len(mapping["negative_singlet_phases"]) == 18
    addresses = [x for row in mapping["CKM_blocks"] for x in row["CAR_indices"]]
    addresses += [row["CAR_index"] for row in mapping["negative_singlet_phases"]]
    addresses += mapping["identity_CAR_indices"]
    assert sorted(addresses) == list(range(90))
    expected = np.kron(np.eye(2), np.asarray(frame.matter.mass_basis_change().evalf(20), complex))
    np.testing.assert_allclose(executed["one_body_basis_change90"], expected, atol=2e-10, rtol=0)
    assert all(target % 9 == 0 for event in executed["decoded_sourceword"] for target in event["targets"])
    assert executed["resources"]["native_pulses"] > 6*executed["CKM_block_resources"]["native_pulses"]


def test_public_words_and_resources_cannot_poison_a_warm_source_frame(executed):
    expected = executed["resources"]["native_pulses"]
    first = frame.source_frame()
    first["resources"]["native_pulses"] = -1
    first["sourceword"][0]["targets"][0] = -1
    first["mapping"]["CKM_blocks"][0]["CAR_indices"][0] = -1
    again = frame.source_frame()
    assert again["resources"]["native_pulses"] == expected
    assert again["sourceword"][0]["targets"][0] >= 0
    assert again["mapping"]["CKM_blocks"][0]["CAR_indices"][0] >= 0


@pytest.mark.parametrize("corruption", ("complex_phase", "vacuum_phase", "wrong_wire", "bosonic_swap"))
def test_actual_cached_source_words_cannot_hide_phase_or_wiring_corruption(executed, monkeypatch, corruption):
    if corruption == "wrong_wire":
        original = frame.protected.decoded_source_word
        def changed(word, sites):
            result = original(word, sites)
            result[0]["targets"][0] += 1
            return result
        monkeypatch.setattr(frame.protected, "decoded_source_word", changed)
    else:
        original = frame.g.local_gate
        def changed(kind, coefficient=0., duration=0.):
            if corruption == "complex_phase" and kind == "hopping" and abs(complex(coefficient).real) > .1:
                return original(kind, complex(coefficient).conjugate(), duration)
            gate = original(kind, coefficient, duration)
            if corruption == "bosonic_swap" and kind == "fswap":
                word = tuple(gate.word)+(frame.protected.compiler.pulse("J", math.pi, [0, 1]),)
                return replace(gate, word=word)
            if corruption == "vacuum_phase" and kind == "hopping":
                word = tuple(gate.word)+(frame.protected.compiler.pulse("I", .01, [0]),)
                return replace(gate, word=word)
            return gate
        monkeypatch.setattr(frame.g, "local_gate", changed)
    with pytest.raises(ArithmeticError):
        frame.source_frame()


@pytest.mark.parametrize("singlet_sign", (1, -1))
def test_complete_64_column_source_frame_mass_frame_inverse_composition(singlet_sign):
    result = frame.source_defined_down_sector((.13, -.27, .91), singlet_sign=singlet_sign)
    assert result["source_all_input_columns"] == 64
    assert not result["source_global_phase_removed_or_fitted"]
    assert not result["numerical_admission_is_outward_or_physical_certificate"]
    np.testing.assert_allclose(result["weak_unitary"], result["weak_target"], atol=2e-10, rtol=0)
    np.testing.assert_allclose(result["mass_pair_product"], result["mass_target"], atol=2e-10, rtol=0)
    assert max(result["errors"].values()) < 2e-10
    with pytest.raises(ValueError):
        result["weak_unitary"].flags.writeable = True
    # The weak action is not itself the independent diagonal mass product.
    assert np.linalg.norm(result["weak_unitary"]-result["mass_pair_product"]) > .01


def test_64_column_composition_rejects_unreduced_nonfinite_or_boolean_angles():
    for bad in ((True, .2, .3), (5., .2, .3), (float("nan"), .2, .3), [.1, .2, .3]):
        with pytest.raises(ValueError):
            frame.source_defined_down_sector(bad)
    with pytest.raises(ValueError):
        frame.source_defined_down_sector((.1, .2, .3), singlet_sign=False)
