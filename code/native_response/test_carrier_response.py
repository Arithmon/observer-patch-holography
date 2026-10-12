"""Small independent controls for the analytic collective-source extension.

The arbitrary-factor statement is a proof, not a giant matrix calculation.
Two private factors expose lost phases, wrong Cayley lifts, and the distinction
between an inner response and a chart permutation of central records.
"""
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as sp

from native_response.carrier import decode, encode, refresh_address
from native_response.geometry import exchange_generator, odd_unitary
from native_response.source import source_generators


ROOT = Path(__file__).resolve().parents[2]


def exact_blocks(raw):
    def scalar(row):
        a, b, c, d = map(sp.Rational, row)
        return a + b * sp.sqrt(5) + sp.I * (c + d * sp.sqrt(5))
    return sp.diag(*[sp.Matrix([[scalar(x) for x in row] for row in block])
                     for block in raw])


def real_column(matrix):
    values = list(matrix)
    return sp.Matrix([sp.re(x).expand() for x in values]
                     + [sp.im(x).expand() for x in values])


def order(action):
    current = list(range(12))
    for n in range(1, 61):
        current = [action[p] for p in current]
        if current == list(range(12)):
            return n
    raise AssertionError("not a source permutation")


@pytest.fixture(scope="module")
def packet():
    return json.loads((ROOT / "code/source_selection_model/response.json").read_text())


def collective(matrix):
    identity = np.eye(len(matrix))
    return np.kron(matrix, identity) + np.kron(identity, matrix)


def test_exact_source_span_is_the_full_closed_twelve_dimensional_algebra(packet):
    generators = [exact_blocks(row) for row in packet["generators"]]
    # This rank calculation uses SymPy's independent number-field arithmetic,
    # not the source author's elimination or recorded bracket coefficients.
    columns = sp.Matrix.hstack(*(real_column(x) for x in generators))
    assert columns.to_DM(extension=True).rank() == 12
    for x, actual in zip(generators, source_generators(), strict=True):
        assert (x + x.conjugate().T).applyfunc(sp.expand) == sp.zeros(6)
        assert x[:3, 3:] == sp.zeros(3)
        assert x[3:, :3] == sp.zeros(3)
        assert x[3:, 3:].applyfunc(sp.im) == sp.zeros(3)
        np.testing.assert_allclose(actual, np.array(x, dtype=complex),
                                   atol=2e-15, rtol=0)
    # Every generator is in u(3) + so(3), a closed space of real dimension 12;
    # rank 12 therefore proves equality and closure without a cutoff.


def test_collective_response_preserves_brackets_and_has_a_left_inverse():
    source = source_generators()
    lifted = np.array([collective(x) for x in source])
    # Tr_2 Delta(X)/6 = X + Tr(X) I/6. Removing its scalar part recovers X.
    for x, y in zip(source, lifted, strict=True):
        reduced = np.trace(y.reshape(6, 6, 6, 6), axis1=1, axis2=3) / 6
        recovered = reduced - np.trace(reduced) * np.eye(6) / 12
        np.testing.assert_allclose(recovered, x, atol=2e-15, rtol=0)
    for a, b in ((0, 2), (0, 4), (1, 9)):
        bracket = source[a] @ source[b] - source[b] @ source[a]
        np.testing.assert_allclose(lifted[a] @ lifted[b] - lifted[b] @ lifted[a],
                                   collective(bracket), atol=2e-14, rtol=0)
    x = source[0] + 2 * source[4]
    np.testing.assert_allclose(expm(collective(x) / 7),
                               np.kron(expm(x / 7), expm(x / 7)),
                               atol=2e-15, rtol=0)


def test_full_effects_detect_the_relative_phase_on_one_private_factor():
    # The response centre is a relative phase, not a scalar on C^6.
    centre = source_generators().sum(axis=0)
    effect = np.zeros((6, 6), complex)
    effect[0, 3] = 1
    probe = np.kron(effect, np.eye(6))
    observed = collective(centre) @ probe - probe @ collective(centre)
    np.testing.assert_allclose(observed, 1j * probe, atol=4e-15, rtol=0)
    # If all effects were block diagonal, precisely this response is silent.
    diagonal = np.diag([1, 2, 3, 4, 5, 6])
    np.testing.assert_allclose(centre @ diagonal - diagonal @ centre,
                               0, atol=4e-15, rtol=0)


def test_executed_paths_tensor_lift_but_collective_cayley_is_not_the_lift(packet):
    path = next(row for row in packet["closed_paths"]
                if order(row["port_action"]) == 3)
    unitary = np.array(exact_blocks(path["observed_response"]), complex)
    identity = np.eye(6)
    cayley_control = (unitary - identity) @ np.linalg.inv(unitary + identity)
    lifted_control = collective(cayley_control)
    wrong = ((np.eye(36) + lifted_control)
             @ np.linalg.inv(np.eye(36) - lifted_control))
    right = np.kron(unitary, unitary)
    assert np.linalg.norm(wrong - right) > 1
    for p, x in enumerate(source_generators()):
        np.testing.assert_allclose(right @ collective(x) @ right.conj().T,
                                   collective(source_generators()[path["port_action"][p]]),
                                   atol=8e-15, rtol=0)


def test_inner_response_fixes_centre_while_chart_relabels_it(packet):
    path = next(row for row in packet["closed_paths"]
                if order(row["port_action"]) == 5)
    action = path["port_action"]
    unitary = np.array(exact_blocks(path["observed_response"]), complex)

    def inner(blocks):
        return unitary @ blocks @ unitary.conj().T

    def chart(blocks):
        result = np.empty_like(blocks)
        result[action] = inner(blocks)
        return result

    p = next(p for p in range(12) if action[p] != p)
    central = np.zeros((12, 6, 6), complex)
    central[p] = np.eye(6)
    np.testing.assert_allclose(inner(central), central, atol=4e-16, rtol=0)
    expected = np.zeros_like(central)
    expected[action[p]] = np.eye(6)
    np.testing.assert_allclose(chart(central), expected, atol=4e-16, rtol=0)
    state = np.zeros_like(central)
    vector = np.array([1, 0, 0, 1j, 0, 0]) / np.sqrt(2)
    state[p] = np.outer(vector, vector.conj())
    np.testing.assert_allclose(decode(chart(state)), inner(decode(state)), atol=2e-16)
    np.testing.assert_allclose(chart(encode(decode(state))),
                               encode(inner(decode(state))), atol=2e-16)
    np.testing.assert_allclose(chart(refresh_address(state)),
                               refresh_address(chart(state)), atol=2e-16)
    # Private-factor reindexing also commutes with every uniform response.
    for x in source_generators():
        lifted = collective(x).reshape(6, 6, 6, 6)
        np.testing.assert_array_equal(lifted.transpose(1, 0, 3, 2), lifted)


def test_marked_slow_and_lower_triplets_have_no_intertwiner(packet):
    paths = packet["closed_paths"]
    five = next(row for row in paths if order(row["port_action"]) == 5)
    u = exact_blocks(five["observed_response"])
    slow, lower = u[:3, :3], u[3:, 3:]
    traces = {sp.expand(sp.trace(slow)), sp.expand(sp.trace(lower))}
    assert traces == {(1 + sp.sqrt(5)) / 2, (1 - sp.sqrt(5)) / 2}
    # Order five alone shares a fixed axis. All order-three paths remove
    # that spurious intertwiner; this is an exact kernel, not just a trace test.
    equations = []
    for path in [five] + [row for row in paths if order(row["port_action"]) == 3]:
        matrix = exact_blocks(path["observed_response"])
        a, b = matrix[:3, :3], matrix[3:, 3:]
        equations.append(sp.kronecker_product(sp.eye(3), a)
                         - sp.kronecker_product(b.T, sp.eye(3)))
    assert sp.Matrix.vstack(*equations).to_DM(extension=True).rank() == 9


def test_complete_chart_implementer_has_endogenous_centralizer_decomposition(packet):
    # A genuine C2 subgroup exchanges both endpoints of this seam. Its odd
    # load sector and two private factors suffice, without a 6^12 matrix.
    path = next(row for row in packet["closed_paths"]
                if order(row["port_action"]) == 2
                and row["port_action"][0] == 2 and row["port_action"][2] == 0)
    unitary = np.array(exact_blocks(path["observed_response"]), complex)
    load_swap = np.array([[0, 1], [1, 0]])
    private_swap = np.eye(36).reshape(6, 6, 6, 6).transpose(1, 0, 2, 3).reshape(36, 36)
    endogenous = np.kron(np.eye(2), np.kron(unitary, unitary))
    centralizer = np.kron(load_swap, private_swap)
    chart = np.kron(load_swap, private_swap @ np.kron(unitary, unitary))
    np.testing.assert_allclose(chart, endogenous @ centralizer, atol=1e-15, rtol=0)
    np.testing.assert_allclose(chart @ chart, np.eye(72), atol=2e-15, rtol=0)
    source = source_generators()
    for p, x in enumerate(source):
        lifted = np.kron(np.eye(2), collective(x))
        np.testing.assert_array_equal(centralizer @ lifted, lifted @ centralizer)
        expected = np.kron(np.eye(2), collective(source[path["port_action"][p]]))
        np.testing.assert_allclose(chart @ lifted @ chart.conj().T,
                                   expected, atol=8e-15, rtol=0)
    # This centralizer is response-trivial, not silent on every observable.
    contrast = np.kron(np.diag([1, -1]), np.eye(36))
    np.testing.assert_array_equal(centralizer @ contrast @ centralizer.conj().T,
                                  -contrast)
    # Internal conjugation reverses g; the load swap restores the canonical
    # seam orientation. Dropping that permutation fails for the actual T.
    generator = exchange_generator((0, 2))
    np.testing.assert_allclose(chart @ generator @ chart.conj().T,
                               generator, atol=8e-15, rtol=0)
    assert np.linalg.norm(endogenous @ generator @ endogenous.conj().T - generator) > 1
    repair_unitary = odd_unitary((0, 2), 1 / 8)

    def repair(value):
        pulled = repair_unitary.conj().T @ value @ repair_unitary
        retained = (pulled[:36, :36] + pulled[36:, 36:]) / 2
        return repair_unitary @ np.kron(np.eye(2), retained) @ repair_unitary.conj().T

    vector = np.arange(1, 73) + 1j * np.arange(72, 0, -1)
    vector /= np.linalg.norm(vector)
    state = np.outer(vector, vector.conj())
    np.testing.assert_allclose(repair(chart @ state @ chart.conj().T),
                               chart @ repair(state) @ chart.conj().T,
                               atol=2e-16, rtol=0)


def test_independent_address_controls_would_enlarge_the_response():
    # Two address blocks already double the rank. Such controls are excluded
    # from the witness grammar, rather than silently treated as twelve probes.
    source = source_generators().reshape(12, -1)
    real = np.concatenate((source.real, source.imag), axis=1).T
    independent = np.kron(np.eye(2), real)
    assert np.linalg.matrix_rank(real) == 12
    assert np.linalg.matrix_rank(independent) == 24
