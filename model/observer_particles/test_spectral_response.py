"""Independent spectral identities and countercontrols for visible field weights."""
import numpy as np
import pytest

from .spectral_response import charge_resolved_spectrum, spectrum_summary


def test_dark_modes_and_negative_frequencies_remain_with_unnormalized_weights():
    h = np.diag([-2., 1., 1., 5.])
    q = np.diag([-1., 1., 1., 0.])
    b = np.array([[2, 0], [0, 3], [0, 0], [0, 0]], complex)
    result = charge_resolved_spectrum(h, q, b)
    assert sum(row["cluster_dimension"] for row in result["clusters"]) == 4
    rows = {tuple(row["frequency_interval"]): row for row in result["clusters"]}
    assert rows[(1., 1.)]["cluster_dimension"] == 2
    assert rows[(1., 1.)]["numerical_residue_rank"] == 1
    assert rows[(5., 5.)]["dark_at_numerical_admission"]
    np.testing.assert_allclose(result["residue_sum"], np.diag([4., 9.]), atol=0)
    # At these exact degenerate energies the resolvent has one residue per
    # full eigenspace. This inverse independently checks normalization/sign.
    z = .4+1.7j
    expansion = sum(row["residue"]/(z-row["frequency_interval"][0])
                    for row in result["clusters"])
    np.testing.assert_allclose(expansion, b.conj().T@np.linalg.solve(z*np.eye(4)-h, b))


def test_joint_charge_weights_are_invariant_under_degenerate_basis_changes():
    h = np.diag([1., 1., 2., 2.])
    q = np.diag([-1., 1., -1., 1.])
    b = np.array([[1, 0], [1j, 0], [0, 2], [0, 3j]], complex)
    raw = np.array([[1, 1j, 2, 3], [2j, 1, 1, -1j],
                    [3, 2, 1j, 1], [1, -1, 2j, 3j]], complex)
    u, _ = np.linalg.qr(raw)
    first = charge_resolved_spectrum(h, q, b)
    rotated = charge_resolved_spectrum(u@h@u.conj().T, u@q@u.conj().T, u@b)
    assert len(first["clusters"]) == len(rotated["clusters"]) == 4
    for a, c in zip(first["clusters"], rotated["clusters"], strict=True):
        np.testing.assert_allclose(a["charge_interval"], c["charge_interval"], atol=1e-14)
        np.testing.assert_allclose(a["frequency_interval"], c["frequency_interval"], atol=1e-14)
        np.testing.assert_allclose(a["residue"], c["residue"], atol=2e-14)
    np.testing.assert_allclose(rotated["residue_sum"], b.conj().T@b, atol=2e-14)


def test_unresolved_frequencies_are_not_averaged_into_a_single_pole():
    frequencies = np.array([1., 1.+.75e-8, 1.+1.5e-8])
    result = charge_resolved_spectrum(np.diag(frequencies), np.eye(3), np.eye(3),
                                     energy_resolution=1e-8)
    assert [r["cluster_dimension"] for r in result["clusters"]] == [2, 1]
    np.testing.assert_array_equal(result["clusters"][0]["frequencies"], frequencies[:2])
    summary = spectrum_summary(result)
    assert not summary["scope"]["nearby_frequencies_identified_as_one_exact_pole"]
    assert sum(r["cluster_dimension"] for r in summary["clusters"]) == 3
    assert not summary["scope"]["exact_residue_ranks_certified"]
    assert not any(r["exact_energy_degeneracy_asserted"] for r in summary["clusters"])


def test_nonconserved_charge_and_damping_are_outside_hermitian_contract():
    with pytest.raises(ValueError, match="not conserved"):
        charge_resolved_spectrum(np.diag([1., 2.]), [[0, 1], [1, 0]], np.eye(2))
    with pytest.raises(ValueError, match="Hermitian"):
        charge_resolved_spectrum(np.diag([1.-.2j, 2.-.2j]), np.eye(2), np.eye(2))
    zero = charge_resolved_spectrum(np.diag([1., 2.]), np.eye(2), np.zeros((2, 1)))
    assert all(r["dark_at_numerical_admission"] for r in zero["clusters"])
    assert sum(r["cluster_dimension"] for r in zero["clusters"]) == 2


def test_charged_nambu_symbol_keeps_conjugate_and_spin_multiplicities():
    from model.common_source.chiral_candidate import kinetic

    h = kinetic.bdg_symbol((.2, -.3, .7), np.array([[0, .4], [.4, 0]]), speed=1.)
    charge = kinetic.nambu_charge(np.diag([-1., 1.]))
    # A fixed annihilation coordinate sees less than the entire Nambu space.
    b = np.eye(8, dtype=complex)[:, :1]
    result = charge_resolved_spectrum(h, charge, b)
    assert len(result["clusters"]) == 4
    assert all(r["cluster_dimension"] == 2 for r in result["clusters"])
    assert sum(r["dark_at_numerical_admission"] for r in result["clusters"]) == 2
    np.testing.assert_allclose(result["residue_sum"], [[1]], atol=1e-14)
    assert sum(r["cluster_dimension"] for r in result["clusters"]) == 8
