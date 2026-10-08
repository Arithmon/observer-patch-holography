"""Original-input Kogut--Susskind controls across Perron and response scales.

The exact integer L=2 geometry and scalar high-precision Hamiltonian are built
without producer masks, orbit labels, ground states, or probabilities.  In
particular, influence differences are evaluated before converting populations
to float64.  Repeated-precision agreement is not an interval certificate.
"""

from __future__ import annotations

import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import z2_finite_transfer_receipt as z2  # noqa: E402


def _original_input_control(lam: float, dps: int) -> dict[str, float]:
    def link(x: int, y: int, direction: int) -> int:
        return 2 * ((x % 2) * 2 + y % 2) + direction

    stars = [
        sum(1 << i for i in {link(x, y, 0), link(x, y, 1),
                            link(x - 1, y, 0), link(x, y - 1, 1)})
        for x in range(2) for y in range(2)
    ]
    gauge = {0}
    for star in stars:
        gauge |= {element ^ star for element in gauge}
    representatives = sorted({min(c ^ g for g in gauge) for c in range(256)})
    labels = {r ^ g: i for i, r in enumerate(representatives) for g in gauge}
    plaquettes = [
        sum(1 << i for i in {link(x, y, 0), link(x + 1, y, 1),
                            link(x, y + 1, 0), link(x, y, 1)})
        for x in range(2) for y in range(2)
    ]
    potential = [sum((-1) ** ((r & p).bit_count()) for p in plaquettes)
                 for r in representatives]
    with mp.workdps(dps):
        coupling = mp.mpf(lam)  # preserve the original supplied float exactly
        h = mp.zeros(32)
        for i, r in enumerate(representatives):
            h[i, i] = -potential[i]
            for link_index in range(8):
                h[i, labels[r ^ (1 << link_index)]] -= coupling
        eigenvalues, eigenvectors = mp.eigsy(h)
        omega = list(eigenvectors[:, 0])
        if sum(omega) < 0:
            omega = [-value for value in omega]
        assert all(value > 0 for value in omega)
        pi = [value * value for value in omega]
        kernels = [
            [pi[i] / (pi[i] + pi[labels[r ^ (1 << l)]])
             for i, r in enumerate(representatives)]
            for l in range(8)
        ]
        influences = [
            mp.fsum(
                max(abs(kernels[l][i] - kernels[l][labels[r ^ (1 << u)]])
                    for i, r in enumerate(representatives))
                for u in range(8) if u != l
            )
            for l in range(8)
        ]
        rates = [
            coupling * (omega[i] / omega[labels[r ^ (1 << l)]]
                        + omega[labels[r ^ (1 << l)]] / omega[i])
            for l in range(8) for i, r in enumerate(representatives)
        ]
        return {
            "gap_H": float(eigenvalues[1] - eigenvalues[0]),
            "ground_energy": float(eigenvalues[0]),
            "pi_min": float(min(pi)), "pi_max": float(max(pi)),
            "rate_min": float(min(rates)), "rate_max": float(max(rates)),
            "eta_star": float(max(influences)),
        }


def _compare_reported_components(result: dict, expected: dict[str, float]) -> None:
    actual = {
        "gap_H": result["spectral"]["gap_H"],
        "ground_energy": result["ground_energy"],
        "pi_min": result["spectral"]["pi_min"],
        "pi_max": result["spectral"]["pi_max"],
        "rate_min": result["fiber_dependent_rates"]["rate_min"],
        "rate_max": result["fiber_dependent_rates"]["rate_max"],
        "eta_star": result["dobrushin"]["eta_star"],
    }
    for key, value in expected.items():
        assert actual[key] == pytest.approx(value, rel=1e-7, abs=0), key
    assert result["doob_generator_rows_sum_zero"]
    assert result["doob_generator_offdiagonal_nonpositive"]


@pytest.mark.parametrize("lam", [0.05, 0.5, 1.0, 2.0, 1e4, 5e-5, 1e10, 1e16])
def test_ks_reported_components_resolve_original_input_or_refuse(lam: float) -> None:
    low = _original_input_control(lam, 90)
    high = _original_input_control(lam, 120)
    for key in high:
        np.testing.assert_allclose(low[key], high[key], rtol=1e-14, atol=0,
                                   err_msg=key)
    try:
        result = z2.evaluate(z2.Z2GaugeOrbits(2), "kogut_susskind", lam=lam)
    except RuntimeError as error:
        # The default grid and independent resolved controls must compute.
        # Stress cases may decline a diagnostic whose precision is lost.
        if lam in (0.05, 0.5, 1.0, 2.0, 1e4):
            raise
        assert any(word in str(error).lower()
                   for word in ("unresolved", "precision", "resolution"))
        return
    _compare_reported_components(result, high)
