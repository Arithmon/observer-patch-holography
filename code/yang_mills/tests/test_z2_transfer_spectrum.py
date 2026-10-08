"""Original-input controls for the finite Wilson logarithm and Perron state.

The free Hamiltonian is known before forming a transfer matrix.  Its small
transfer eigenvalues and its small leading eigengap stress different numerical
operations, so checking the Hamiltonian gap alone does not certify this receipt.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from pathlib import Path
import sys

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

import z2_finite_transfer_receipt as z2  # noqa: E402


def _original_input_dual_coupling(beta: float) -> float:
    """Compute log(coth(beta)) from the supplied float in scalar decimal math."""
    with localcontext() as ctx:
        ctx.prec = 80
        q = (-2 * Decimal.from_float(beta)).exp()
        return float(((1 + q) / (1 - q)).ln())


@pytest.mark.parametrize("beta_t", [0.01, 0.02, 0.1, 0.5, 10.0, 18.0])
def test_free_wilson_full_generator_and_stationary_diagnostics(
    beta_t: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    orbits = z2.Z2GaugeOrbits(2)
    n = orbits.n_orbits
    rate = _original_input_dual_coupling(beta_t)
    identity = np.eye(n)
    # On the 2 x 2 torus the shortest nonempty closed electric loop has two
    # links.  The character spectrum therefore gives gap(H) = 2 * rate,
    # independently of the producer's transfer eigendecomposition.
    expected_h = rate / 2 * sum(
        identity - identity[orbits.flip[link]]
        for link in range(orbits.n_links)
    )
    observed = {}
    original_doob_transform = z2.doob_transform

    def record_doob_transform(H, omega, e0):
        observed["H"] = H.copy() - e0 * identity
        observed["omega"] = omega.copy()
        result = original_doob_transform(H, omega, e0)
        observed["L"] = result.copy()
        return result

    monkeypatch.setattr(z2, "doob_transform", record_doob_transform)
    result = z2.evaluate(orbits, "wilson", beta_s=0.0, beta_t=beta_t)

    # Normalize by the original-input rate before comparing: an absolute
    # 1e-9 check would accept a wholly wrong Hamiltonian when beta_t is large.
    np.testing.assert_allclose(observed["H"] / rate, expected_h / rate,
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(observed["L"] / rate, expected_h / rate,
                               rtol=2e-12, atol=2e-12)
    np.testing.assert_allclose(observed["omega"] ** 2, np.full(n, 1 / n),
                               rtol=2e-12, atol=0)
    assert result["doob_generator_rows_sum_zero"]
    assert result["doob_generator_offdiagonal_nonpositive"]
    fit = result["constant_rate_fit"]
    assert fit["relative_frobenius_residual"] < 2e-12
    np.testing.assert_allclose(np.asarray(fit["rates"]) / rate, 1,
                               rtol=2e-12, atol=0)
    fiber = result["fiber_dependent_rates"]
    assert fiber["rate_min"] / rate == pytest.approx(1, rel=2e-12)
    assert fiber["rate_max"] / rate == pytest.approx(1, rel=2e-12)
    assert fiber["spread_max_over_min"] == pytest.approx(1, rel=2e-12)
    assert fiber["offdiagonal_mass_outside_single_flip"] / rate < 2e-10
    assert result["dobrushin"]["eta_star"] < 2e-12
    assert result["dobrushin"]["dobrushin_condition_holds"]
    assert result["dobrushin"]["unit_rate_floor_c_star_times_1_minus_eta"] == (
        pytest.approx(1, rel=2e-12)
    )
    assert result["spectral"]["gap_H"] / rate == pytest.approx(2, rel=2e-12)
    assert result["spectral"]["gap_unit_rate_heat_bath"] == pytest.approx(2, rel=2e-12)
    assert result["spectral"]["pi_min"] == pytest.approx(1 / n, rel=2e-12)
    assert result["spectral"]["pi_max"] == pytest.approx(1 / n, rel=2e-12)


@pytest.mark.parametrize("beta_s,beta_t", [(0.01, 0.01), (0.001, 10.0), (0.0, 18.0)])
def test_generic_log_refuses_unresolved_spectrum_or_perron_state(
    beta_s: float, beta_t: float
) -> None:
    # The first input loses the bottom of the transfer spectrum.  The others
    # have positive, well-scaled eigenvalues but an unresolved leading gap.
    # Either failure must prevent all downstream physical classifications.
    transfer = z2.Z2GaugeOrbits(2).wilson_transfer(beta_s, beta_t)
    with pytest.raises((ValueError, RuntimeError),
                       match=r"(?i)(resol|precision|condition|spectral|perron)"):
        z2.symmetric_log_hamiltonian(transfer)


@pytest.mark.parametrize("scale", [1e-200, 1.0, 1e200])
@pytest.mark.parametrize("second_eigenvalue", [0.5, 1 - 1e-5])
def test_generic_log_keeps_resolved_repeated_and_close_eigenvalues(
    scale: float, second_eigenvalue: float
) -> None:
    # An exact dyadic orthogonal basis isolates the repeated-eigenvalue and
    # energy-origin cases from geometry, orbit enumeration and the producer.
    basis = np.array([[1, 1, 1, 1], [1, -1, 1, -1],
                      [1, 1, -1, -1], [1, -1, -1, 1]], dtype=float) / 2
    eigenvalues = np.array([1.0, second_eigenvalue, 0.5, 0.25])
    transfer = scale * ((basis * eigenvalues) @ basis.T)
    expected = (basis * -np.log(eigenvalues)) @ basis.T
    H, omega, maximum = z2.symmetric_log_hamiltonian(transfer)
    np.testing.assert_allclose(H, expected, rtol=2e-11, atol=2e-14)
    np.testing.assert_allclose(omega ** 2, np.full(4, 0.25),
                               rtol=2e-10, atol=0)
    assert maximum / scale == pytest.approx(1, rel=2e-14)
    assert z2.spectral_gap(H) == pytest.approx(-np.log(second_eigenvalue), rel=2e-10)
