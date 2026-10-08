"""Original-input controls for the finite Wilson logarithm and Perron state.

The free Hamiltonian is known before forming a transfer matrix.  Its small
transfer eigenvalues and its small leading eigengap stress different numerical
operations, so checking the Hamiltonian gap alone does not certify this receipt.
"""

from __future__ import annotations

from decimal import Decimal, localcontext
from pathlib import Path
import sys

import mpmath as mp
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


def _original_input_wilson_control(beta_s: float, beta_t: float, dps: int):
    """Build the 2 x 2 transfer directly from integer incidence and scalar exp.

    No producer geometry, matrix, eigenvalues, support decision or stationary
    distribution enters this oracle.  Repeating the precision is an accuracy
    check, not an interval certificate.
    """
    def link(x, y, direction):
        return 2 * ((x % 2) * 2 + (y % 2)) + direction

    stars = [
        sum(1 << j for j in {link(x, y, 0), link(x, y, 1),
                            link(x - 1, y, 0), link(x, y - 1, 1)})
        for x in range(2) for y in range(2)
    ]
    gauge = {0}
    for star in stars:
        gauge |= {element ^ star for element in gauge}
    representatives = sorted({min(config ^ g for g in gauge) for config in range(256)})
    plaquettes = [
        sum(1 << j for j in {link(x, y, 0), link(x + 1, y, 1),
                            link(x, y + 1, 0), link(x, y, 1)})
        for x in range(2) for y in range(2)
    ]
    potential = [sum((-1) ** ((r & p).bit_count()) for p in plaquettes)
                 for r in representatives]
    with mp.workdps(dps):
        bs, bt = mp.mpf(beta_s), mp.mpf(beta_t)
        weights = [mp.exp(bs * p / 2) for p in potential]
        kinetic = [mp.exp(bt * (8 - 2 * distance)) for distance in range(9)]
        transfer = mp.matrix([
            [weights[i] * weights[j] * mp.fsum(
                kinetic[(a ^ b ^ g).bit_count()] for g in gauge
            ) for j, b in enumerate(representatives)]
            for i, a in enumerate(representatives)
        ])
        eigenvalues, eigenvectors = mp.eigsy(transfer)
        energies = [-mp.log(value / eigenvalues[31]) for value in eigenvalues]
        H = eigenvectors * mp.diag(energies) * eigenvectors.T
        omega = eigenvectors[:, 31]
        if sum(omega) < 0:
            omega = -omega
        return (np.array(H.tolist(), dtype=float),
                np.array(list(omega), dtype=float),
                float(eigenvalues[31]), float(energies[30]))


@pytest.mark.parametrize("beta_s,beta_t", [(0.01, 0.01), (0.1, 0.1), (0.3, 0.8),
                                          (-0.3, 0.8), (0.001, 4.0)])
def test_source_wilson_hamiltonian_against_original_input_controls(
    beta_s: float, beta_t: float
) -> None:
    low = _original_input_wilson_control(beta_s, beta_t, 60)
    high = _original_input_wilson_control(beta_s, beta_t, 90)
    for first, second in zip(low, high):
        np.testing.assert_allclose(first, second, rtol=1e-14, atol=1e-45)
    expected_H, expected_omega, expected_maximum, expected_gap = high
    H, omega, maximum = z2.wilson_hamiltonian(z2.Z2GaugeOrbits(2), beta_s, beta_t)
    h_scale = np.max(np.abs(expected_H))
    np.testing.assert_allclose(H / h_scale, expected_H / h_scale,
                               rtol=2e-11, atol=2e-12)
    np.testing.assert_allclose(omega ** 2, expected_omega ** 2, rtol=2e-11, atol=0)
    assert maximum == pytest.approx(expected_maximum, rel=2e-12)
    result = z2.evaluate(z2.Z2GaugeOrbits(2), "wilson", beta_s=beta_s, beta_t=beta_t)
    assert result["spectral"]["gap_H"] == pytest.approx(expected_gap, rel=2e-11)
    assert result["spectral"]["pi_min"] == pytest.approx(min(expected_omega ** 2), rel=2e-11)
    assert result["spectral"]["pi_max"] == pytest.approx(max(expected_omega ** 2), rel=2e-11)


def test_l3_source_spectrum_near_free_limit_resolves_small_transfer_eigenvalues(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # A nonzero spatial coupling exercises the source-factor path.  The free
    # limit has an analytic whole-H control even though tanh(.1)^18 ~ 1e-18:
    # forming a dense transfer Gram matrix first cannot resolve its bottom.
    orbits = z2.Z2GaugeOrbits(3)
    rate = _original_input_dual_coupling(0.1)
    identity = np.eye(orbits.n_orbits)
    expected_H = rate / 2 * sum(
        identity - identity[orbits.flip[link]] for link in range(orbits.n_links)
    )
    original_svd = z2.dgejsv
    source_calls = []

    def record_source_svd(*args, **kwargs):
        source_calls.append(args[0].shape)
        return original_svd(*args, **kwargs)

    monkeypatch.setattr(z2, "dgejsv", record_source_svd)
    H, omega, _ = z2.wilson_hamiltonian(orbits, 1e-15, 0.1)
    # The full-H comparison is the numerical control; this transparent spy
    # additionally checks that a tiny nonzero coupling was not silently set to zero.
    assert source_calls == [(orbits.n_orbits, orbits.n_orbits)]
    np.testing.assert_allclose(H / rate, expected_H / rate,
                               rtol=2e-10, atol=2e-10)
    np.testing.assert_allclose(omega ** 2, np.full(orbits.n_orbits, 1 / orbits.n_orbits),
                               rtol=2e-11, atol=0)


@pytest.mark.parametrize("public_evaluate", [False, True])
def test_source_refuses_unresolved_perron_probabilities(public_evaluate: bool) -> None:
    # At this original-input point the first factor prototype had excellent
    # absolute pi agreement but 4.69% relative population error.  Source factor
    # accuracy alone cannot justify the pointwise divisions in the Doob law.
    orbits = z2.Z2GaugeOrbits(2)
    with pytest.raises(RuntimeError, match=r"(?i)perron.*(separation|unresolved|support)"):
        if public_evaluate:
            z2.evaluate(orbits, "wilson", beta_s=0.001, beta_t=10.0)
        else:
            z2.wilson_hamiltonian(orbits, 0.001, 10.0)


@pytest.mark.parametrize("beta_s,beta_t", [(3.0, 3.0), (3.5, 2.5)])
def test_small_support_doob_diagnostics_are_accurate_or_explicitly_unresolved(
    beta_s: float, beta_t: float,
) -> None:
    # This point has a well-separated transfer Perron root, but pi_min ~2e-20.
    # Tiny matrix-log errors are amplified in the pointwise Doob division.
    # A false conservation flag must not be presented as physical evidence.
    H, omega, _, _ = _original_input_wilson_control(beta_s, beta_t, 90)
    expected_doob = H * omega[None, :] / omega[:, None]
    np.testing.assert_allclose(expected_doob.sum(axis=1), 0, atol=5e-14)
    orbits = z2.Z2GaugeOrbits(2)
    mask = np.eye(orbits.n_orbits, dtype=bool)
    for flip in orbits.flip:
        mask[np.arange(orbits.n_orbits), flip] = True
    expected_outside_mass = np.abs(expected_doob[~mask]).sum()
    try:
        result = z2.evaluate(orbits, "wilson", beta_s=beta_s, beta_t=beta_t)
    except RuntimeError as error:
        assert any(word in str(error).lower() for word in ("unresolved", "precision", "resolution"))
        return
    assert result["doob_generator_rows_sum_zero"]
    assert result["fiber_dependent_rates"]["offdiagonal_mass_outside_single_flip"] == (
        pytest.approx(expected_outside_mass, rel=1e-7, abs=0)
    )


@pytest.mark.parametrize("fault", ["info", "rank", "nonzero_rank", "denormal", "singular", "normalization"])
def test_source_does_not_report_an_incomplete_svd_as_physical_evidence(
    fault: str, monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Fault injection checks the fail-closed interface, not SVD accuracy.
    # The unmodified real SVD and the same supported source point are checked
    # against original-input high precision in the tests above.
    original_svd = z2.dgejsv

    def incomplete_svd(*args, **kwargs):
        singular, left, right, work, rank, info = original_svd(*args, **kwargs)
        singular, work, rank = singular.copy(), work.copy(), rank.copy()
        if fault == "info":
            info = 1
        elif fault == "rank":
            rank[0] -= 1
        elif fault == "nonzero_rank":
            rank[1] -= 1
        elif fault == "denormal":
            rank[2] = 1
        elif fault == "singular":
            singular[-1] = np.nan
        else:
            work[0] = np.nan
        return singular, left, right, work, rank, info

    monkeypatch.setattr(z2, "dgejsv", incomplete_svd)
    with pytest.raises(RuntimeError, match=r"(?i)(unresolved|singular|precision)"):
        z2.evaluate(z2.Z2GaugeOrbits(2), "wilson", beta_s=0.1, beta_t=0.1)
