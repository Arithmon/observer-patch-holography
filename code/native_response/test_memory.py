"""Full-shell counting, unreduced private overlaps and kinetic controls."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from native_response import memory
from native_response.geometry import bare_contrast, load_response
from native_response.source import edges, laplacian, source_generators
from native_response.test_geometry import conjugated_kraus


@pytest.fixture(scope="module", params=(1 / 16, 1 / 8))
def data(request):
    return memory.closure_data(request.param)


def shell():
    return np.array([[int(p in occupied) for p in range(12)]
                     for occupied in combinations(range(12), 2)])


def test_reference_covariance_and_all_overlap_counts_are_full_shell_counts():
    states = shell()
    centered = states - states.mean(axis=0)
    np.testing.assert_allclose(centered.T @ centered / 66,
                               (5 / 33) * (np.eye(12) - np.ones((12, 12)) / 12),
                               rtol=0, atol=1e-15)
    for e, f in combinations(edges(), 2):
        shared = set(e) & set(f)
        if not shared:
            continue
        p = shared.pop()
        both_odd = (states[:, e[0]] != states[:, e[1]]) & (
            states[:, f[0]] != states[:, f[1]])
        assert Counter(states[both_odd, p]) == {1: 9, 0: 1}
    for a, b in edges():
        assert np.count_nonzero(states[:, a] != states[:, b]) == 20


def test_two_edge_pair_orbits_exhaust_all_120_overlaps():
    packet = json.loads((Path(__file__).resolve().parents[1]
                        / "source_selection_model/response.json").read_text())
    adjacent = {frozenset((e, f)) for e, f in combinations(edges(), 2) if set(e) & set(f)}
    orbits = []
    for representative in (((0, 2), (0, 4)), ((0, 2), (0, 8))):
        orbit = set()
        for row in packet["closed_paths"]:
            action = row["port_action"]
            orbit.add(frozenset(tuple(sorted((action[a], action[b])))
                                for a, b in representative))
        assert len(orbit) == 60
        assert all(sum(e in pair for pair in orbit) == 4 for e in edges())
        orbits.append(orbit)
    assert orbits[0].isdisjoint(orbits[1])
    assert orbits[0] | orbits[1] == adjacent


def kraus_residual(edge, theta):
    z = bare_contrast()
    q = sum(k.conj().T @ z @ k for k in conjugated_kraus(edge, theta))
    beta = np.trace(z @ q).real / 72
    return q - beta * z


def embed_private_pair(operator, edge, ports):
    """Explicit 216-dimensional embedding, with no normalized partial trace."""
    coordinates = np.array(list(np.ndindex(6, 6, 6)))
    first, second = (ports.index(p) for p in edge)
    spectator = ({0, 1, 2} - {first, second}).pop()
    local = 6 * coordinates[:, first] + coordinates[:, second]
    return (operator[np.ix_(local, local)]
            * (coordinates[:, spectator, None] == coordinates[None, :, spectator]))


def test_local_overlap_reduction_matches_unreduced_three_private_trace(data):
    seams = edges()
    # Shared port first/first, second/first, and second/second: retain signs.
    pairs = (((0, 2), (0, 4)), ((0, 2), (2, 5)), ((2, 5), (3, 5)))
    states = shell()
    for e, f in pairs:
        assert e in seams and f in seams
        left, right = (kraus_residual(edge, data["theta"]) for edge in (e, f))
        ports = sorted(set(e) | set(f))
        # Count actual initial full configurations, rather than insert9 and1.
        patterns = Counter(tuple(x[ports]) for x in states
                           if x[e[0]] != x[e[1]] and x[f[0]] != x[f[1]])
        direct = 0.0
        for pattern, count in patterns.items():
            e_sign = 1 - pattern[ports.index(e[0])]
            f_sign = 1 - pattern[ports.index(f[0])]
            eb = left[36 * e_sign:36 * (e_sign + 1), 36 * e_sign:36 * (e_sign + 1)]
            fb = right[36 * f_sign:36 * (f_sign + 1), 36 * f_sign:36 * (f_sign + 1)]
            E = embed_private_pair(eb, e, ports)
            F = embed_private_pair(fb, f, ports)
            direct += count * np.einsum("ij,ji->", E, F).real / (66 * 216)
        assert data["edge_gram"][seams.index(e), seams.index(f)] == pytest.approx(
            direct, rel=2e-11, abs=1e-15)
    residual = kraus_residual(seams[0], data["theta"])
    odd_count = sum(x[seams[0][0]] != x[seams[0][1]] for x in states)
    direct_diagonal = odd_count * np.vdot(residual, residual).real / (66 * 72)
    assert data["edge_gram"][0, 0] == pytest.approx(direct_diagonal, rel=2e-12)
    # Same-edge load coherences contribute: deleting them must change the norm.
    assert np.linalg.norm(residual[:36, 36:]) > 1e-3


def test_all_seams_and_bands_have_the_certified_structure(data):
    gram = data["edge_gram"]
    for row, edge in zip(gram, edges()):
        values = [abs(x) for x, other in zip(row, edges()) if edge != other and set(edge) & set(other)]
        assert len(values) == 8
        distinct = sorted(values)
        np.testing.assert_allclose(distinct[:4], distinct[0], rtol=1e-10, atol=1e-15)
        np.testing.assert_allclose(distinct[4:], distinct[4], rtol=1e-10, atol=1e-15)
        assert all(x == 0 for x, other in zip(row, edges()) if not set(edge) & set(other))
    margin = np.min(np.diag(gram) - (np.abs(gram).sum(axis=1) - np.diag(gram)))
    assert margin > (0.0037 if data["theta"] == 1 / 8 else 0.00018)
    boundary = memory.incidence()
    np.testing.assert_array_equal(boundary.T @ boundary, laplacian())
    np.testing.assert_allclose(data["normalized_leakage"],
                               boundary.T @ gram @ boundary / (3600 * (5 / 33)),
                               rtol=1e-14, atol=1e-18)
    np.testing.assert_allclose(data["normalized_leakage"] @ np.ones(12), 0, atol=1e-18)
    spectrum = np.linalg.eigvalsh(data["normalized_leakage"])
    assert np.count_nonzero(spectrum > 1e-7) == 11
    # The equilibrium rate and the excess are checked separately, not just a ratio.
    expected = {1 / 8: (0.9759867045598407, 4.17822334397240e-5, 0.0206705310245),
                1 / 16: (0.9982364126230368, 3.29877975132457e-6, 0.00156003481520)}
    alpha, excess, ratio = expected[data["theta"]]
    assert data["alpha"] == pytest.approx(alpha, rel=0, abs=2e-13)
    np.testing.assert_allclose(spectrum[1:4], excess, rtol=2e-10, atol=0)
    rate = data["alpha"] * (5 - np.sqrt(5)) / 60
    np.testing.assert_allclose(spectrum[1:4] / rate**2, ratio, rtol=2e-10, atol=0)
    # The positive operator lower bound applies to arbitrary superpositions.
    lower = (11 * margin / 6000) * laplacian()
    assert np.linalg.eigvalsh(data["normalized_leakage"] - lower).min() > -1e-17


def test_single_port_correction_from_unprojected_event_moments(data):
    """Obtain c2 before projecting, then compare the corrected local statistic."""
    states = shell()
    incident = [edge for edge in edges() if 0 in edge]
    assert len(incident) == 5 and all(edge[0] == 0 for edge in incident)
    z = bare_contrast()
    # A f=(1/60) sum_e (Z_e-Phi_e(Z_e)) over the five incident edges.
    # These are not the projected residuals Phi_e(Z_e)-beta Z_e.
    defects = [z - sum(k.conj().T @ z @ k
                       for k in conjugated_kraus(edge, data["theta"]))
               for edge in incident]
    gram = np.zeros((5, 5))
    first = 0.0
    for i, (edge, defect) in enumerate(zip(incident, defects)):
        odd = states[:, edge[0]] != states[:, edge[1]]
        for state in states[odd]:
            sign = 1 - state[0]
            block = defect[36 * sign:36 * (sign + 1),
                           36 * sign:36 * (sign + 1)]
            first += (state[0] - 1 / 6) * np.trace(block).real / (66 * 36 * 60)
        gram[i, i] = np.count_nonzero(odd) * np.vdot(defect, defect).real / (66 * 72)
        for j in range(i):
            other = incident[j]
            ports = sorted(set(edge) | set(other))
            patterns = Counter(tuple(x[ports]) for x in states
                               if x[edge[0]] != x[edge[1]]
                               and x[other[0]] != x[other[1]])
            overlap = 0.0
            for pattern, count in patterns.items():
                sign = 1 - pattern[ports.index(0)]
                block = slice(36 * sign, 36 * (sign + 1))
                # Full three-private-factor trace; do not use the producer's
                # endpoint partial traces or its projected edge Gram.
                left = embed_private_pair(defect[block, block], edge, ports)
                right = embed_private_pair(defects[j][block, block], other, ports)
                overlap += count * np.einsum("ij,ji->", left, right).real / (66 * 216)
            gram[i, j] = gram[j, i] = overlap
    variance = 5 / 36
    r = first / variance
    c1 = 1 - r
    c2 = 1 - 2 * r + gram.sum() / (3600 * variance)
    corrected = c2 - c1**2 - (1 - c1)**2 / 10
    expected = data["leakage_matrix"][0, 0] / variance
    assert r == pytest.approx(data["alpha"] / 11, rel=3e-14, abs=0)
    assert corrected == pytest.approx(expected, rel=0, abs=4e-16)
    # Uncorrected scalar variance contains the known mixture of carrier bands.
    assert c2 - c1**2 - expected == pytest.approx(data["alpha"]**2 / 1210,
                                               rel=0, abs=4e-16)
    # Lower endpoints already established by test_certificate.py.
    lower = {1 / 16: Fraction(186661676, 10**12),
             1 / 8: Fraction(3778430107, 10**12)}[data["theta"]]
    assert corrected > float(lower / 100) > 0


def test_native_local_mode_mixture_is_not_hidden_memory():
    """Exact 66-state events give the false positive and its corrected zero."""
    states = shell()
    indices = {tuple(state): i for i, state in enumerate(states)}
    counts = np.zeros((66, 66), dtype=int)
    for i, state in enumerate(states):
        for a, b in edges():
            swapped = state.copy()
            swapped[a], swapped[b] = state[b], state[a]
            counts[i, i] += 1
            counts[i, indices[tuple(swapped)]] += 1
    np.testing.assert_array_equal(counts.sum(axis=1), 60)
    np.testing.assert_array_equal(counts, counts.T)
    f = 6 * states[:, 0] - 1  # Integer multiple of N_0-1/6.
    next_f = counts @ f
    norm = int(f @ f)
    c1 = Fraction(int(f @ next_f), 60 * norm)
    c2 = Fraction(int(next_f @ next_f), 3600 * norm)
    assert c1 == Fraction(10, 11)
    assert c2 == Fraction(91, 110)
    assert c2 - c1**2 == Fraction(1, 1210)
    assert c2 - c1**2 - (1 - c1)**2 / 10 == 0


@pytest.mark.parametrize("keep", ("g", "h"))
def test_full_network_one_term_control_has_zero_memory(monkeypatch, keep):
    ds = source_generators()

    def control(edge, theta):
        a, b = edge
        private = (ds[b] - ds[a] if keep == "g"
                   else ds[a] + ds[b] - ds.sum(axis=0) / 6)
        private = np.kron(private, np.eye(6)) + np.kron(np.eye(6), private)
        pauli = (np.array([[0, -1j], [1j, 0]]) if keep == "g"
                 else np.array([[0, 1], [1, 0]]))
        W = expm(theta * np.kron(pauli, private))
        pullback = W.conj().T @ bare_contrast() @ W
        R = (pullback[:36, :36] + pullback[36:, 36:]) / 2
        return {"alpha": 1 - np.vdot(R, R).real / 36,
                "output": W @ np.kron(np.eye(2), R) @ W.conj().T}

    monkeypatch.setattr(memory, "load_response", control)
    result = memory.closure_data(1 / 8)
    np.testing.assert_allclose(result["normalized_leakage"], 0, atol=1e-28, rtol=0)


def test_native_zero_coupling_and_public_band_output():
    native = memory.closure_data(0)
    assert native["alpha"] == 1
    np.testing.assert_array_equal(native["normalized_leakage"], np.zeros((12, 12)))
    bands = memory.band_diagnostics(1 / 8)
    assert [b["multiplicity"] for b in bands] == [3, 5, 3]
    expected = (4.17822334397240e-5, 6.65553092475258e-5, 6.91432819728189e-5)
    for band, value in zip(bands, expected):
        np.testing.assert_allclose(band["second_moment_excess"], value, rtol=2e-10, atol=0)
        assert all(x > 0 for x in band["clock_invariant_excess"])
