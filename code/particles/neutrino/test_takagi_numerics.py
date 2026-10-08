"""Independent spectra guard against false massless-neutrino readouts."""

import sys
import json
import subprocess
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from particles.neutrino.build_forward_majorana_matrix import _sorted_takagi
from particles.neutrino.derive_intrinsic_neutrino_exact_eta_map import (
    _build_exact_eta_map,
    _takagi_unitary,
)
from particles.takagi import sorted_takagi


READERS = [sorted_takagi, _sorted_takagi, _takagi_unitary]


@pytest.mark.parametrize("reader", READERS)
@pytest.mark.parametrize("scale", [1e-320, 1e-200, 1.0, 1e200])
def test_diagonal_masses_survive_unit_changes_without_gram_squaring(reader, scale):
    expected = np.array([scale, 2 * scale, 3 * scale])
    # Phases and signs must be removed by congruence, not absolute-value
    # projection of a Hermitian matrix. Every supplied mass is representable.
    matrix = np.diag([expected[0] * 1j, -expected[1], expected[2]])
    masses, unitary = reader(matrix)
    np.testing.assert_allclose(masses / expected, np.ones(3), rtol=2e-15, atol=0)
    scaled = matrix.real / scale + 1j * (matrix.imag / scale)
    np.testing.assert_allclose(unitary.T @ scaled @ unitary, np.diag(masses / scale), atol=2e-14)
    np.testing.assert_allclose(unitary.conj().T @ unitary, np.eye(3), atol=2e-14)


@pytest.mark.parametrize("reader", READERS)
def test_resolved_small_mass_is_not_replaced_by_zero(reader):
    rotation = np.array([[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]])
    expected = np.array([1e-12, 1.0, 3.0])
    matrix = rotation @ np.diag(expected) @ rotation.T
    masses, unitary = reader(matrix)
    # Assembly itself has O(epsilon) absolute error: allow it relative to the
    # small mass, but never the O(sqrt(epsilon)) loss from the normal matrix.
    np.testing.assert_allclose(masses, expected, rtol=2e-4, atol=0)
    np.testing.assert_allclose(unitary.T @ matrix @ unitary, np.diag(masses), atol=2e-14)


@pytest.mark.parametrize("reader", READERS)
def test_complex_congruence_and_null_mode(reader):
    r12 = np.array([[0.6, 0.8, 0], [-0.8, 0.6, 0], [0, 0, 1]])
    r23 = np.array([[1, 0, 0], [0, 0.8, 0.6], [0, -0.6, 0.8]])
    u = r12 @ np.diag(np.exp(1j * np.array([0.2, -0.7, 1.1]))) @ r23
    expected = np.array([0.0, 2.0, 5.0])
    matrix = u.conj() @ np.diag(expected) @ u.conj().T
    masses, actual = reader(matrix)
    np.testing.assert_allclose(masses, expected, atol=3e-15, rtol=2e-15)
    np.testing.assert_allclose(actual.T @ matrix @ actual, np.diag(masses), atol=6e-15)
    np.testing.assert_allclose(actual.conj().T @ actual, np.eye(3), atol=2e-15)


@pytest.mark.parametrize("reader", READERS)
def test_degenerate_block_without_congruence_resolution_is_refused(reader):
    with pytest.raises(ValueError, match="degenerate-block"):
        reader(np.array([[0, 1j, 0], [1j, 0, 0], [0, 0, 3]]))


@pytest.mark.parametrize("reader", READERS)
@pytest.mark.parametrize("bad", [
    [[True, 0], [0, 1]], [[np.nan, 0], [0, 1]], [[np.inf, 0], [0, 1]],
    [[1, 0, 0], [0, 1, 0]],
    [[1e-200, 1e-200], [0, 1e-200]],
])
def test_invalid_or_relatively_nonsymmetric_input_is_refused(reader, bad):
    with pytest.raises(ValueError):
        reader(bad)


def test_intrinsic_cubic_cancellation_cannot_erase_two_massive_modes():
    a, rho = 1.0 + 1e-8, 1.0
    result = _build_exact_eta_map(a, rho, 0.0, np.zeros(3))
    # At zero cycle and eta, M=(a-rho)I+rho*11^T. Its eigenvalues are
    # exact without diagonalizing either M or M†M.
    expected = np.array([a - rho, a - rho, a + 2 * rho])
    np.testing.assert_allclose(result.masses, expected, rtol=5e-8, atol=0)
    np.testing.assert_allclose(result.masses_squared, expected**2, rtol=1e-7, atol=0)
    assert np.all(result.masses_squared > 0)
    np.testing.assert_allclose(result.u_takagi.T @ result.majorana @ result.u_takagi,
                               np.diag(result.masses), atol=2e-14)


@pytest.mark.parametrize("mode", ["canonical_selector", "real_seed"])
def test_forward_serialization_preserves_complex_phases_at_small_mass_scale(tmp_path, mode):
    root = Path(__file__).resolve().parents[2]
    script = root / "particles/neutrino/build_forward_majorana_matrix.py"
    anchor = json.loads((root / "particles/runs/neutrino/neutrino_scale_anchor.json").read_text())
    results = []
    for scale in (1.0, 1e-20):
        anchor["anchors"]["m_star_gev"] = scale
        anchor_path, output = tmp_path / "anchor.json", tmp_path / "majorana.json"
        anchor_path.write_text(json.dumps(anchor))
        subprocess.run([sys.executable, str(script), "--scale-anchor", str(anchor_path),
                        "--mode", mode, "--out", str(output)], check=True, capture_output=True)
        payload = json.loads(output.read_text())
        matrix = np.array(payload["majorana_matrix_real"]) + 1j * np.array(payload["majorana_matrix_imag"])
        unitary = np.array(payload["U_nu_real"]) + 1j * np.array(payload["U_nu_imag"])
        masses = np.array(payload["masses_sorted_gev"])
        np.testing.assert_allclose(unitary.T @ (matrix / scale) @ unitary,
                                   np.diag(masses / scale), atol=3e-14)
        assert (payload["eigenvalues_raw_gev"] is None) == (mode == "canonical_selector")
        assert payload["public_surface_candidate_allowed"] is False
        results.append(masses / scale)
    np.testing.assert_allclose(results[0], results[1], rtol=2e-15)


def test_real_shortcut_rejects_a_nonsymmetric_matrix_instead_of_using_one_triangle(tmp_path):
    root = Path(__file__).resolve().parents[2]
    family = json.loads((root / "particles/runs/neutrino/family_response_tensor.json").read_text())
    family["C_nu_hat_real"][0][1] += 0.25
    family_path, output = tmp_path / "family.json", tmp_path / "majorana.json"
    family_path.write_text(json.dumps(family))
    result = subprocess.run([sys.executable, str(root / "particles/neutrino/build_forward_majorana_matrix.py"),
                             "--mode", "real_seed", "--family", str(family_path), "--out", str(output)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "Majorana matrix must be complex symmetric" in result.stderr
    assert not output.exists()


def _run_splittings(tmp_path, masses):
    root = Path(__file__).resolve().parents[2]
    majorana, output = tmp_path / "majorana.json", tmp_path / "gaps.json"
    majorana.write_text(json.dumps({"masses_sorted_gev": masses}))
    result = subprocess.run([sys.executable, str(root / "particles/neutrino/build_forward_splittings.py"),
                             "--majorana", str(majorana), "--out", str(output)],
                            capture_output=True, text=True)
    return result, output


@pytest.mark.parametrize("scale", [1.0, 1e-20, 1e-150, 1e150])
def test_forward_gap_ratio_is_invariant_under_resolved_unit_changes(tmp_path, scale):
    result, output = _run_splittings(tmp_path, [scale, 2 * scale, 3 * scale])
    assert result.returncode == 0, result.stderr
    payload = json.loads(output.read_text())
    for key, expected in {"s1_minus_s0": 3, "s2_minus_s0": 8, "s2_minus_s1": 5}.items():
        assert payload["ascending_mass_sq_gaps_gev2"][key] / (scale * scale) == pytest.approx(expected, rel=2e-15)
    assert payload["ascending_gap_ratio_s10_over_s20"] == pytest.approx(3 / 8, rel=2e-15)
    assert payload["public_surface_candidate_allowed"] is False


def test_forward_gap_does_not_subtract_separately_rounded_large_squares(tmp_path):
    lower = 1e154
    upper = float(np.nextafter(lower, np.inf))
    result, output = _run_splittings(tmp_path, [0, lower, upper])
    assert result.returncode == 0, result.stderr
    payload = json.loads(output.read_text())
    # Independently factored difference avoids the two nearly equal squares.
    expected = (upper - lower) * (upper + lower)
    assert payload["ascending_mass_sq_gaps_gev2"]["s2_minus_s1"] == pytest.approx(expected, rel=2e-15)


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_forward_gaps_refuse_underflow_or_overflow_instead_of_false_zeros(tmp_path, scale):
    result, output = _run_splittings(tmp_path, [scale, 2 * scale, 3 * scale])
    assert result.returncode != 0
    assert "squared mass gap exceeds finite binary64 range" in result.stderr
    assert not output.exists()


def test_equal_masses_keep_exact_zero_gaps_and_undefined_ratio(tmp_path):
    result, output = _run_splittings(tmp_path, [2, 2, 2])
    assert result.returncode == 0, result.stderr
    payload = json.loads(output.read_text())
    assert set(payload["ascending_mass_sq_gaps_gev2"].values()) == {0}
    assert payload["ascending_gap_ratio_s10_over_s20"] is None
