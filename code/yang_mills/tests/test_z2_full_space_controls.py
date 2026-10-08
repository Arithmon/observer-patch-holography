"""Compare the orbit receipt with an independent 256-configuration L=2 model.

Geometry is built from signed edge endpoints, gauge transformations from site
signs, and heat-bath probabilities from the physical +1/-1 link states.  None
of the producer's masks, orbit members, flips, or probabilities construct the
controls.  These ordinary-coupling checks concern finite mathematical
conventions, not the resolution of extreme floating-point parameters.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import z2_finite_transfer_receipt as z2  # noqa: E402


@pytest.fixture(scope="module")
def full_space() -> dict:
    sites = [(x, y) for x in range(2) for y in range(2)]
    edges = [
        (a, ((a[0] + (d == 0)) % 2, (a[1] + (d == 1)) % 2), d)
        for a in sites
        for d in (0, 1)
    ]
    n = len(edges)
    size = 2**n
    spins = np.array(
        [[1 - 2 * ((c // 2**i) % 2) for i in range(n)] for c in range(size)]
    )
    boundaries = []
    for x, y in sites:
        boundary = [
            ((x, y), 0), (((x + 1) % 2, y), 1),
            ((x, (y + 1) % 2), 0), ((x, y), 1),
        ]
        boundaries.append([
            next(i for i, (a, _, d) in enumerate(edges) if (a, d) == key)
            for key in boundary
        ])
    potential = np.array([
        sum(math.prod(int(s[i]) for i in p) for p in boundaries)
        for s in spins
    ])
    orbit_sets = []
    remaining = set(range(size))
    while remaining:
        c = min(remaining)
        members = set()
        for g in range(2**len(sites)):
            signs = {a: 1 - 2 * ((g // 2**i) % 2) for i, a in enumerate(sites)}
            transformed = [
                spins[c, i] * signs[a] * signs[b]
                for i, (a, b, _) in enumerate(edges)
            ]
            members.add(sum(2**i for i, s in enumerate(transformed) if s == -1))
        orbit_sets.append(members)
        remaining -= members
    restriction = np.zeros((size, len(orbit_sets)))
    for o, members in enumerate(orbit_sets):
        restriction[list(members), o] = 1 / math.sqrt(len(members))
    return {
        "n_links": n, "size": size, "spins": spins,
        "potential": potential, "orbit_sets": orbit_sets,
        "restriction": restriction,
    }


def test_orbit_basis_matches_site_sign_action(full_space: dict) -> None:
    orbits = z2.Z2GaugeOrbits(2)
    members = full_space["orbit_sets"]
    q = full_space["restriction"]
    assert len(members) == 32
    assert all(len(orbit) == 8 for orbit in members)
    assert [min(orbit) for orbit in members] == list(orbits.reps)
    assert all(orbit == set(orbits.members[o]) for o, orbit in enumerate(members))
    np.testing.assert_allclose(q.T @ q, np.eye(32), atol=1e-15)
    np.testing.assert_array_equal(
        full_space["potential"][[min(orbit) for orbit in members]],
        orbits.plaquette_sum,
    )
    # Distinct single-link moves must remain distinct after the quotient;
    # otherwise extracting each collar rate from one matrix entry overcounts.
    for c in range(full_space["size"]):
        destinations = [
            next(o for o, orbit in enumerate(members) if c ^ 2**l in orbit)
            for l in range(full_space["n_links"])
        ]
        assert len(set(destinations)) == full_space["n_links"]


@pytest.mark.parametrize("transfer, params", [
    ("wilson", {"beta_s": 0.1, "beta_t": 0.1}),
    ("wilson", {"beta_s": 0.5, "beta_t": 0.5}),
    ("wilson", {"beta_s": 0.3, "beta_t": 0.8}),
    ("kogut_susskind", {"lam": 0.5}),
    ("kogut_susskind", {"lam": 1.0}),
    ("kogut_susskind", {"lam": 2.0}),
])
def test_full_configuration_operator_and_consumers(
    full_space: dict, transfer: str, params: dict[str, float]
) -> None:
    n, size = full_space["n_links"], full_space["size"]
    q, p = full_space["restriction"], full_space["potential"]
    spins = full_space["spins"]
    orbits = z2.Z2GaugeOrbits(2)
    if transfer == "wilson":
        weight = np.exp(params["beta_s"] * p / 2)
        t = (weight[:, None] * np.exp(params["beta_t"] * (spins @ spins.T))
             * weight[None, :])
        eigenvalues, eigenvectors = np.linalg.eigh(t)
        assert np.all(eigenvalues > 0)
        h = (eigenvectors * -np.log(eigenvalues / eigenvalues[-1])) @ eigenvectors.T
        omega, e0 = eigenvectors[:, -1], 0.0
        orbit_t = orbits.wilson_transfer(**params)
        np.testing.assert_allclose(q.T @ t @ q, orbit_t, atol=1e-10, rtol=1e-13)
        orbit_h, orbit_omega, _ = z2.symmetric_log_hamiltonian(orbit_t)
        orbit_e0 = 0.0
    else:
        h = -np.diag(p.astype(float))
        for c in range(size):
            for l in range(n):
                h[c, c ^ 2**l] = -params["lam"]
        eigenvalues, eigenvectors = np.linalg.eigh(h)
        omega, e0 = eigenvectors[:, 0], eigenvalues[0]
        orbit_h = orbits.kogut_susskind(**params)
        orbit_omega, orbit_e0 = z2.ground_state(orbit_h)
    if omega.sum() < 0:
        omega = -omega
    assert np.all(omega > 0)
    pi = omega**2
    orbit_pi = orbit_omega**2
    orbit_pi /= orbit_pi.sum()
    np.testing.assert_allclose(q.T @ h @ q, orbit_h, atol=1e-9, rtol=1e-10)
    np.testing.assert_allclose(q.T @ omega, orbit_omega, atol=1e-12)
    gen = (h - e0 * np.eye(size)) * omega[None, :] / omega[:, None]
    orbit_gen = z2.doob_transform(orbit_h, orbit_omega, orbit_e0)
    np.testing.assert_allclose(q.T @ gen @ q, orbit_gen, atol=1e-9, rtol=1e-10)
    np.testing.assert_allclose(gen.sum(axis=1), 0, atol=1e-10)
    np.testing.assert_allclose(orbit_gen.sum(axis=1), 0, atol=1e-10)

    projectors, plus_probabilities = [], []
    for l in range(n):
        expectation = np.zeros((size, size))
        plus_probability = np.zeros(size)
        for c in range(size):
            cplus, cminus = c & ~(2**l), c | 2**l
            prob = pi[cplus] / (pi[cplus] + pi[cminus])
            expectation[c, cplus], expectation[c, cminus] = prob, 1 - prob
            plus_probability[c] = prob
        projectors.append(expectation)
        plus_probabilities.append(plus_probability)
    for full_e, orbit_e in zip(projectors, z2.heat_bath_projectors(orbits, orbit_pi)):
        np.testing.assert_allclose(q.T @ full_e @ q, orbit_e, atol=1e-12)
        np.testing.assert_allclose(full_e @ full_e, full_e, atol=1e-14)
        np.testing.assert_allclose(pi @ full_e, pi, atol=1e-14)
    influence = np.zeros((n, n))
    for l in range(n):
        for u in range(n):
            if l != u:
                influence[l, u] = max(
                    abs(plus_probabilities[l][c] - plus_probabilities[l][c ^ 2**u])
                    for c in range(size)
                )
    dob = z2.dobrushin_influence(orbits, orbit_pi)
    np.testing.assert_allclose(influence.sum(axis=1), dob["row_sums"], atol=1e-12)
    receipt = z2.evaluate(orbits, transfer, **params)
    assert receipt["doob_generator_rows_sum_zero"]
    assert receipt["dobrushin"]["eta_star"] == pytest.approx(
        float(influence.sum(axis=1).max()), abs=1e-12
    )
    quotient_hb = q.T @ sum(np.eye(size) - e for e in projectors) @ q
    symmetric_hb = quotient_hb * orbit_omega[:, None] / orbit_omega[None, :]
    expected_hb_gap = np.linalg.eigvalsh((symmetric_hb + symmetric_hb.T) / 2)[1]
    expected_h_gap = np.linalg.eigvalsh(q.T @ h @ q)[1] - e0
    assert receipt["spectral"]["gap_H"] == pytest.approx(expected_h_gap, abs=1e-9)
    assert receipt["spectral"]["gap_unit_rate_heat_bath"] == pytest.approx(
        expected_hb_gap, abs=1e-12
    )
    if transfer == "kogut_susskind":
        rates = np.array([
            [params["lam"] * (omega[c] / omega[c ^ 2**l]
                              + omega[c ^ 2**l] / omega[c]) for c in range(size)]
            for l in range(n)
        ])
        reconstruction = sum(
            rate[:, None] * (np.eye(size) - e) for rate, e in zip(rates, projectors)
        )
        np.testing.assert_allclose(gen, reconstruction, atol=1e-11, rtol=1e-12)
        reported = receipt["fiber_dependent_rates"]
        assert reported["rate_min"] == pytest.approx(float(rates.min()), abs=1e-12)
        assert reported["rate_max"] == pytest.approx(float(rates.max()), abs=1e-12)
        assert reported["offdiagonal_mass_outside_single_flip"] < 1e-12
        assert receipt["doob_generator_offdiagonal_nonpositive"]
