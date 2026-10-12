"""Independent source, local-channel, and equilibrium geometry controls."""
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import block_diag, expm
import sympy as sp

from native_response.geometry import (
    alpha_lower_bound, bare_contrast, checked_theta, edge_dirichlet,
    exchange_amplitude, exchange_generator, load_response, odd_unitary,
)
from native_response.source import edges, laplacian, source_generators, vertices


ROOT = Path(__file__).resolve().parents[2]
ATOL = 2e-12


def decoded_blocks(raw):
    def scalar(value):
        a, b, c, d = map(Fraction, value)
        return complex(float(a) + float(b) * np.sqrt(5),
                       float(c) + float(d) * np.sqrt(5))
    return block_diag(*[np.array([[scalar(x) for x in row] for row in block])
                        for block in raw])


def conjugated_kraus(edge, theta):
    """Construct all four qubit depolarizer branches, not a partial trace."""
    unitary = odd_unitary(edge, theta)
    branches = []
    for i in range(2):
        for j in range(2):
            local = np.zeros((2, 2))
            local[i, j] = 1 / np.sqrt(2)
            branch = np.kron(local, np.eye(36))
            branches.append(unitary @ branch @ unitary.conj().T)
    return branches


def test_source_formula_matches_authoritative_response_and_proper_covariance():
    packet = json.loads((ROOT / "code/source_selection_model/response.json").read_text())
    source = source_generators()
    recorded = np.array([decoded_blocks(row) for row in packet["generators"]])
    np.testing.assert_allclose(source, recorded, atol=ATOL, rtol=0)
    np.testing.assert_allclose(source.conj().transpose(0, 2, 1), -source,
                               atol=ATOL, rtol=0)
    np.testing.assert_allclose(source.sum(axis=0),
                               block_diag(1j * np.eye(3), np.zeros((3, 3))),
                               atol=ATOL, rtol=0)
    assert len(packet["closed_paths"]) == 60
    actions, orbit = set(), set()
    seam_set = set(edges())
    for row in packet["closed_paths"]:
        action = tuple(row["port_action"])
        transport = decoded_blocks(row["observed_response"])
        actions.add(action)
        orbit.add(tuple(sorted((action[0], action[2]))))
        assert {tuple(sorted((action[a], action[b]))) for a, b in edges()} == seam_set
        np.testing.assert_allclose(transport @ transport.conj().T, np.eye(6),
                                   atol=ATOL, rtol=0)
        for offset in (0, 3):
            assert abs(np.linalg.det(transport[offset:offset + 3,
                                              offset:offset + 3]) - 1) < ATOL
        for port in range(12):
            np.testing.assert_allclose(transport @ source[port] @ transport.conj().T,
                                       source[action[port]], atol=ATOL, rtol=0)
    assert len(actions) == 60
    assert orbit == seam_set


def test_exact_source_commutator_coefficient():
    """Independent symbolic source formula fixes the theta^4 coefficient."""
    root5 = sp.sqrt(5)
    phi, conjugate = (1 + root5) / 2, (1 - root5) / 2

    def hat(w):
        x, y, z = w
        return sp.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])

    def response(sign):
        v = sp.Matrix([0, sign, phi])
        w = sp.Matrix([0, sign, conjugate])
        quadrupole = v * v.T / 2 - (v.dot(v) / 6 - sp.Rational(1, 12)) * sp.eye(3)
        return sp.diag(hat(v / 2) + sp.I * quadrupole, hat(w / 2))

    first, second = response(1), response(-1)
    center = sp.diag(sp.I * sp.eye(3), sp.zeros(3))
    g, h = second - first, first + second - center / 6
    q = sp.simplify(sp.I * (g * h - h * g))
    assert sp.simplify(sp.trace(q)) == 0
    assert sp.simplify(sp.trace(q * q) - (49 + 19 * root5)) == 0
    independent = np.array((h - sp.I * g).evalf(), dtype=complex)
    expected = np.kron(independent, np.eye(6)) + np.kron(np.eye(6), independent)
    np.testing.assert_allclose(exchange_amplitude((0, 2)), expected, atol=ATOL, rtol=0)


@pytest.mark.parametrize("theta", (0, 1 / 64, 1 / 16, 1 / 8, -1 / 8, 1 / 2, 1))
def test_local_output_matches_four_kraus_branches_and_private_blocks(theta):
    response = load_response((0, 2), theta)
    branches = conjugated_kraus((0, 2), theta)
    z = bare_contrast()
    direct = sum(branch.conj().T @ z @ branch for branch in branches)
    np.testing.assert_allclose(direct, response["output"], atol=ATOL, rtol=0)
    np.testing.assert_allclose(sum(branch.conj().T @ branch for branch in branches),
                               np.eye(72), atol=ATOL, rtol=0)
    retained = response["retained_private"]
    np.testing.assert_allclose(retained, retained.conj().T, atol=ATOL, rtol=0)
    for first in range(2):
        for second in range(2):
            indices = [6 * i + j for i in range(3 * first, 3 * first + 3)
                       for j in range(3 * second, 3 * second + 3)]
            block = retained[np.ix_(indices, indices)]
            assert abs(np.trace(block)) < ATOL
            assert np.linalg.eigvalsh(block).max() <= 1 + ATOL
            assert np.linalg.eigvalsh(block).min() >= -1 - ATOL
            assert np.trace(block @ block).real <= 8 + ATOL
    alpha = response["alpha"]
    assert alpha_lower_bound(theta) - ATOL <= alpha <= 1 + ATOL
    output = response["output"]
    assert abs(np.trace(z @ output).real / 72 - (1 - alpha)) < ATOL
    assert abs(np.trace(output @ output).real / 72 - (1 - alpha)) < ATOL


def test_source_dirichlet_geometry_and_rank_three_gram():
    frame = vertices()
    carrier = laplacian()
    assert len(edges()) == 30 and np.all(np.diag(carrier) == 5)
    assert (0, 2) in edges() and (0, 1) not in edges()
    np.testing.assert_allclose(carrier @ frame, (5 - np.sqrt(5)) * frame,
                               atol=ATOL, rtol=0)
    result = edge_dirichlet(1 / 8)
    alpha = load_response((0, 2), 1 / 8)["alpha"]
    np.testing.assert_allclose(result, alpha * carrier / 396, atol=ATOL, rtol=0)
    np.testing.assert_allclose(60 * result / np.trace(result), carrier,
                               atol=ATOL, rtol=0)
    eigenvalues, eigenvectors = np.linalg.eigh(result)
    assert np.count_nonzero(np.isclose(eigenvalues, eigenvalues[1], rtol=0,
                                      atol=1e-12)) == 3
    slow_projector = eigenvectors[:, 1:4] @ eigenvectors[:, 1:4].T
    gram = frame @ frame.T / (frame[0] @ frame[0])
    np.testing.assert_allclose(4 * slow_projector, gram, atol=ATOL, rtol=0)


def test_one_step_paired_load_measurement_has_twice_the_dirichlet_form():
    """One isolated tracial protocol is valid; repeated dephasing changes it."""
    theta = 1 / 8
    branches = conjugated_kraus((0, 2), theta)
    # Joint probabilities Tr[P_out K P_in (I72/72) P_in Kdagger].
    flip_probability = sum(np.linalg.norm(branch[:36, 36:])**2
                           + np.linalg.norm(branch[36:, :36])**2
                           for branch in branches) / 72
    alpha = load_response((0, 2), theta)["alpha"]
    assert abs(flip_probability - alpha / 2) < ATOL
    # Total-two shell: odd endpoint probability20/66; uniform30-seam attempts.
    increment_form = (20 / 66) * flip_probability * laplacian() / 30
    np.testing.assert_allclose(increment_form, 2 * edge_dirichlet(theta),
                               atol=ATOL, rtol=0)
    # A conditioned tracial load preparation develops load coherence.
    initial = np.diag([1.] * 36 + [0.] * 36) / 36
    evolved = sum(branch @ initial @ branch.conj().T for branch in branches)
    assert np.linalg.norm(evolved[:36, 36:]) > 1e-4


@pytest.mark.parametrize("keep", ("g", "h"))
def test_one_product_controls_leave_alpha_one(keep):
    ds = source_generators()
    g, h = ds[2] - ds[0], ds[0] + ds[2] - ds.sum(axis=0) / 6
    private = g if keep == "g" else h
    private = np.kron(private, np.eye(6)) + np.kron(np.eye(6), private)
    load = (np.array([[0, -1j], [1j, 0]]) if keep == "g"
            else np.array([[0, 1], [1, 0]]))
    unitary = expm(np.kron(load, private) / 8)
    transformed = unitary.conj().T @ bare_contrast() @ unitary
    retained = (transformed[:36, :36] + transformed[36:, 36:]) / 2
    np.testing.assert_allclose(retained, np.zeros((36, 36)), atol=ATOL, rtol=0)


def test_endpoint_orientation_and_adjoint_are_not_silent_conventions():
    theta = 1 / 8
    forward = load_response((0, 2), theta)
    reverse = load_response((2, 0), theta)
    swap = np.kron(np.array([[0, 1], [1, 0]]), np.eye(36))
    np.testing.assert_allclose(reverse["output"], -swap @ forward["output"] @ swap,
                               atol=ATOL, rtol=0)
    assert abs(forward["alpha"] - reverse["alpha"]) < ATOL
    unitary = odd_unitary((0, 2), theta)
    np.testing.assert_allclose(exchange_generator((0, 2)).conj().T,
                               -exchange_generator((0, 2)), atol=ATOL, rtol=0)
    wrong_pullback = unitary.T @ bare_contrast() @ unitary.conj()
    assert np.linalg.norm(wrong_pullback - forward["pullback"]) > 1


@pytest.mark.parametrize("theta", (True, np.bool_(False), None, [], 1j,
                                   float("nan"), float("inf"),
                                   np.ma.array(.1, mask=True)))
def test_invalid_theta_rejected(theta):
    with pytest.raises(ValueError, match="finite real"):
        checked_theta(theta)


@pytest.mark.parametrize("edge", ((0, 1), (0, 0), (0., 2), (True, 2),
                                  (0, 12), None, np.ma.array([0, 2], mask=True)))
def test_invalid_source_edge_rejected(edge):
    with pytest.raises(ValueError, match="source|ports|seam"):
        exchange_amplitude(edge)


@pytest.mark.parametrize("theta", (1.0001, -1.0001, 1e308))
def test_numerical_diagnostic_domain_is_separate_from_all_angle_bound(theta):
    with pytest.raises(ValueError, match="abs\\(theta\\) <= 1"):
        odd_unitary((0, 2), theta)
    assert alpha_lower_bound(theta) == 1 / 9
