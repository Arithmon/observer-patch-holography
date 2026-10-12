"""Actual record branches, exact shell controls and a coherent N4 witness."""
from fractions import Fraction
from itertools import combinations

import numpy as np
import pytest

from native_response.geometry import load_response, odd_unitary
from native_response.port_readout import (
    BINARY_LOADS, correlation_from_equal_ports, correlation_from_same_port,
    corrected_memory_from_equal_ports, corrected_memory_from_same_port,
    expand_port_record, forget_port_record, luders_branches, measure_occupied_port,
    occupied_port_effects, occupied_port_kraus, sector_indices,
    unread_port_readout,
)
from native_response.source import edges


def independent_loads():
    return tuple(tuple(int(p in pair) for p in range(12))
                 for pair in combinations(range(12), 2))


def independent_selectors(d=1):
    """Construct rectangular branch maps from binary occupation strings."""
    rows = independent_loads()
    result = []
    for port in range(12):
        selected = [d * i + j for i, row in enumerate(rows)
                    if row[port] for j in range(d)]
        matrix = np.zeros((11 * d, 66 * d))
        for output, source in enumerate(selected):
            matrix[output, source] = 1 / np.sqrt(2)
        result.append(matrix)
    return result


def coherent_state(d=2):
    vector = np.zeros(66 * d, complex)
    for pair, private, amplitude in (((0, 1), 0, 1), ((0, 2), 1, 1j),
                                    ((3, 4), 0, -1), ((3, 4), 1, -1j)):
        vector[d * BINARY_LOADS.index(pair) + private] = amplitude / 2
    return np.outer(vector, vector.conj())


def test_effects_are_a_povm_and_output_sectors_are_not_native_projections():
    effects = occupied_port_effects(2)
    np.testing.assert_array_equal(sum(effects), np.eye(132))
    for effect, kraus, port in zip(effects, occupied_port_kraus(2), range(12)):
        assert kraus.shape == (22, 132)
        assert len(sector_indices(port, 2)) == 22
        # Irrational 1/sqrt(2) is rounded; the dyadic effect and direct
        # instrument branch formulas are separately checked exactly below.
        np.testing.assert_allclose(kraus.conj().T @ kraus, effect, atol=2e-16, rtol=0)
        np.testing.assert_array_equal(effect @ effect, effect / 2)
        assert not np.array_equal(effect @ effect, effect)
    # Two port effects overlap on the load {0,1}; their central output
    # projections instead live in different direct-sum blocks.
    assert np.trace(effects[0] @ effects[1]) == Fraction(1, 2)
    assert sum(len(sector_indices(p)) for p in range(12)) == 132
    assert len(BINARY_LOADS) == 66


def test_instrument_matches_independent_kraus_branches_and_keeps_private_phase():
    rho = coherent_state()
    actual = measure_occupied_port(rho)
    kraus = independent_selectors(2)
    for block, branch in zip(actual, kraus):
        np.testing.assert_allclose(block, branch @ rho @ branch.conj().T,
                                   atol=3e-17, rtol=0)
        assert np.linalg.eigvalsh(block).min() >= -1e-16
    assert sum(np.trace(block) for block in actual) == 1
    loads = independent_loads()
    overlap = np.array([[sum(a * b for a, b in zip(x, y)) / 2
                         for y in loads] for x in loads])
    expected = rho * np.kron(overlap, np.ones((2, 2)))
    np.testing.assert_array_equal(forget_port_record(actual), expected)
    i = 2 * BINARY_LOADS.index((3, 4))
    assert expected[i, i + 1] == rho[i, i + 1] != 0  # private coherence
    assert np.linalg.norm(expected - rho) > 0.1  # native load coherence changes


def test_trace_is_preserved_but_fixed_algebra_has_66_load_blocks():
    rho = np.eye(132) / 132
    recorded = measure_occupied_port(rho)
    np.testing.assert_array_equal(forget_port_record(recorded), rho)
    for block in recorded:
        np.testing.assert_array_equal(block, np.eye(22) / 264)
        assert np.trace(block) == pytest.approx(1 / 12, abs=2e-17)
    # Exact support of every fixed matrix unit: x=y, arbitrary private indices.
    coefficients = [[Fraction(len(set(x) & set(y)), 2)
                     for y in BINARY_LOADS] for x in BINARY_LOADS]
    assert all((value == 1) == (i == j) for i, row in enumerate(coefficients)
               for j, value in enumerate(row))
    counts = {value: sum(row.count(value) for row in coefficients)
              for value in (Fraction(0), Fraction(1, 2), Fraction(1))}
    assert counts == {Fraction(0): 2970, Fraction(1, 2): 1320, Fraction(1): 66}
    private = np.array([[1, 1j], [-1j, 2]])
    fixed = np.kron(np.diag(np.arange(66)), private)
    np.testing.assert_array_equal(unread_port_readout(fixed), fixed)
    unit = np.zeros((66, 66))
    unit[0, 1] = 1  # {0,1} and {0,2} share exactly one port.
    once = unread_port_readout(unit)
    twice = unread_port_readout(once)
    np.testing.assert_array_equal(once, unit / 2)
    np.testing.assert_array_equal(twice, unit / 4)


def test_expanded_record_keeps_its_label_and_continues_in_the_address_carrier():
    from native_response.carrier import decode, tick
    from native_response.test_carrier import native_binary

    loads, channels, _ = native_binary()  # independent actual (coin,current) Kraus maps
    left, right = BINARY_LOADS.index((0, 2)), BINARY_LOADS.index((0, 4))
    rho = np.zeros((66, 66), complex)
    rho[left, left] = rho[right, right] = 0.5
    rho[left, right], rho[right, left] = -0.5j, 0.5j
    incoming = np.zeros((12, 66, 66), complex)
    incoming[9] = rho
    compressed = measure_occupied_port(decode(incoming))
    expanded = expand_port_record(compressed)
    assert expanded.shape == incoming.shape
    for port, occupation in enumerate(loads.T):
        projector = np.diag(occupation)
        np.testing.assert_array_equal(expanded[port], projector @ rho @ projector / 2)
    weights = np.trace(expanded, axis1=1, axis2=2)
    assert weights[0] == 0.5 and weights[2] == weights[4] == 0.25
    assert weights[9] == 0  # the requested old address is not another retained record
    changed_address = np.zeros_like(incoming)
    changed_address[7] = rho
    np.testing.assert_array_equal(
        expand_port_record(measure_occupied_port(decode(changed_address))), expanded)
    native_after_read = decode(expanded)
    np.testing.assert_array_equal(native_after_read, forget_port_record(compressed))
    assert np.linalg.norm(native_after_read - rho) > 0.1
    expected = sum(channels.apply(edge, native_after_read) for edge in edges()) / 30
    continued = tick(expanded, channels)
    np.testing.assert_allclose(decode(continued), expected, atol=2e-15, rtol=0)
    # Tick refreshes the address as declared. Keeping the initial tag for a
    # later endpoint comparison requires the separately supplied record log.
    np.testing.assert_allclose(np.trace(continued, axis1=1, axis2=2),
                               np.full(12, 1 / 12), atol=2e-16, rtol=0)


def test_one_central_tag_does_not_store_a_full_pair_or_an_endpoint_history():
    # Uniform full-shell preparation: each observed port leaves eleven
    # equally likely load pairs. Best one-tag exact-pair decoding succeeds1/11.
    joint = [[Fraction(int(port in pair), 132) for pair in BINARY_LOADS]
             for port in range(12)]
    assert all(sum(row) == Fraction(1, 12) for row in joint)
    assert all(sum(value > 0 for value in row) == 11 for row in joint)
    assert sum(max(row) for row in joint) == Fraction(1, 11)
    # Full pair and two unrestricted endpoint labels require more than one
    # twelve-valued record. This instrument supplies only its current tag.
    assert len(BINARY_LOADS) == 66 > 12
    assert len(tuple((p, q) for p in range(12) for q in range(12))) == 144 > 12


def test_disjoint_superpositions_are_lost_even_with_all_classical_outcomes():
    left, right = BINARY_LOADS.index((0, 1)), BINARY_LOADS.index((2, 3))
    plus = np.zeros((66, 66))
    plus[left, left] = plus[right, right] = Fraction(1, 2)
    plus[left, right] = plus[right, left] = Fraction(1, 2)
    minus = plus.copy()
    minus[left, right] = minus[right, left] = -Fraction(1, 2)
    outputs = [measure_occupied_port(state) for state in (plus, minus)]
    for a, b in zip(*outputs):
        np.testing.assert_array_equal(a, b)
    mixture = forget_port_record(outputs[0])
    distance = lambda a, b: np.abs(np.linalg.eigvalsh(a - b)).sum() / 2
    assert distance(plus, minus) == 1
    assert distance(plus, mixture) == distance(minus, mixture) == Fraction(1, 2)
    # Each branch's zero effect on at least one input vector is the support
    # argument for ANY instrument. Output rotations cannot restore that term.
    for p, branch in enumerate(independent_selectors()):
        assert not (np.any(branch[:, left]) and np.any(branch[:, right]))
        phases = np.diag([1j ** (j + p) for j in range(11)])
        changed = phases @ branch[::-1]
        np.testing.assert_array_equal(changed @ plus @ changed.conj().T,
                                      changed @ minus @ changed.conj().T)
    # The 1/2 bound is minimax over this pair, not a claimed global optimum.


def test_coherent_dilation_is_reversible_until_the_pointer_is_classical():
    rho = coherent_state()
    branches = independent_selectors(2)
    isometry = np.vstack(branches)
    coherent = isometry @ rho @ isometry.conj().T
    np.testing.assert_allclose(isometry.conj().T @ coherent @ isometry, rho,
                               atol=1e-16, rtol=0)
    classical = np.zeros_like(coherent)
    for p, expected in enumerate(measure_occupied_port(rho)):
        sector = slice(22 * p, 22 * (p + 1))
        classical[sector, sector] = coherent[sector, sector]
        np.testing.assert_allclose(classical[sector, sector], expected,
                                   atol=3e-17, rtol=0)
    assert np.linalg.norm(classical - coherent) > 0.1
    # A sector projector is not central in the coherent encoded algebra.
    pointer = np.diag([1] * 22 + [0] * (264 - 22))
    assert np.linalg.norm(pointer @ coherent - coherent @ pointer) > 0.1


def native_counts():
    """Integer lazy-swap transitions from actual fair seam events."""
    loads = independent_loads()
    index = {row: i for i, row in enumerate(loads)}
    counts = np.zeros((66, 66), dtype=np.int64)
    for i, row in enumerate(loads):
        for a, b in edges():
            changed = list(row)
            changed[a], changed[b] = changed[b], changed[a]
            counts[i, i] += 1
            counts[index[tuple(changed)], i] += 1
    return counts


def test_endpoint_record_law_and_corrected_native_null_from_actual_branches():
    counts = native_counts()
    projectors = np.array(independent_loads(), dtype=np.int64)
    rho = np.eye(66) / 66
    branches = measure_occupied_port(rho)
    effects = occupied_port_effects()
    exact = []
    for attempts in range(3):
        transition = np.linalg.matrix_power(counts, attempts)
        numerator = projectors.T @ transition @ projectors
        denominator = 4 * 66 * 60**attempts
        q = np.array([[Fraction(int(x), denominator) for x in row]
                      for row in numerator], dtype=object)
        exact.append(q[0, 0])
        for p, branch in enumerate(branches):
            # Begin with the ACTUAL selective state on the original shell.
            state = np.zeros((66, 66), complex)
            rows = [i for i, load in enumerate(independent_loads()) if load[p]]
            state[np.ix_(rows, rows)] = branch
            for _ in range(attempts):
                state = np.diag(counts @ state.diagonal() / 60)
            for r in range(12):
                actual = np.trace(effects[r] @ state).real
                assert actual == pytest.approx(float(q[p, r]), abs=2e-17)
        assert all(sum(row) == Fraction(1, 12) for row in q)
        assert all(q[p, p] == q[0, 0] for p in range(12))
        assert correlation_from_equal_ports(sum(q[p, p] for p in range(12))) == (
            correlation_from_same_port(q[0, 0]))
    assert exact == [Fraction(1, 24), Fraction(61, 1584), Fraction(113, 3168)]
    c1, c2 = map(correlation_from_same_port, exact[1:])
    assert (c1, c2) == (Fraction(10, 11), Fraction(91, 110))
    assert c2 - c1**2 == Fraction(1, 1210)
    assert corrected_memory_from_same_port(*exact[1:]) == 0
    assert corrected_memory_from_equal_ports(*(12 * q for q in exact[1:])) == 0
    assert correlation_from_same_port(exact[0]) == 1
    # Reusing a stored label produces 1/12 at zero separation; two fresh
    # occupied-port draws instead produce 1/24. They are different protocols.
    with pytest.raises(ValueError, match="endpoint range"):
        correlation_from_same_port(Fraction(1, 12))


def test_port_transition_noise_must_be_subtracted_to_recover_geometry():
    occupations = np.array(independent_loads(), dtype=np.int64)
    counts = native_counts()
    q0 = occupations.T @ occupations
    q1 = occupations.T @ counts @ occupations
    # Use integer numerators at common denominator 4*66*60.
    difference = 60 * q0 - q1
    laplacian = np.zeros((12, 12), dtype=np.int64)
    for a, b in edges():
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    np.testing.assert_array_equal(difference, 10 * laplacian)
    # Raw label increments have nonzero noise even when no attempt occurs.
    for i in range(12):
        for j in range(12):
            increment = 2 * (Fraction(int(i == j), 12) - Fraction(int(q0[i, j]), 264))
            assert increment == (Fraction(int(i == j), 11) - Fraction(1, 132))


def test_stationarity_is_a_required_protocol_premise_not_inferred_from_one_number():
    # A reset to load {0,1} is CPTP but not tracial-state preserving. It gives
    # q_00=1/24, a scalar also allowed in the stationary model, while destroying
    # every centered initial/final correlation. The decoder cannot infer the
    # actual protocol from this one number; final marginals expose this control.
    initial_tag_probability = Fraction(1, 12)
    final_tag_probability = Fraction(1, 2)
    q = initial_tag_probability * final_tag_probability
    assert q == Fraction(1, 24)
    true_covariance = q - initial_tag_probability * final_tag_probability
    assert true_covariance == 0
    assert correlation_from_same_port(q) == 1
    assert final_tag_probability != initial_tag_probability


@pytest.mark.parametrize("theta", (0, 1 / 16, 1 / 8))
def test_actual_n4_odd_instrument_has_endpoint_law_and_midpoint_disturbance(theta):
    # Restrict the full shell to {a,s} and {b,s}: the spectator s is occupied
    # in both patterns. The actual E_p POVM has effects N_a/2, N_b/2, I/2.
    d = 36
    na = np.diag([1] * d + [0] * d)
    nb = np.eye(72) - na
    effects = [na / 2, nb / 2, np.eye(72) / 2]
    readout = [na / np.sqrt(2), nb / np.sqrt(2), np.eye(72) / np.sqrt(2)]
    unitary = odd_unitary((0, 2), theta)
    channel_kraus = []
    for i in range(2):
        for j in range(2):
            matrix_unit = np.zeros((2, 2))
            matrix_unit[i, j] = 1
            channel_kraus.append(unitary @ np.kron(matrix_unit, np.eye(d))
                                  @ unitary.conj().T / np.sqrt(2))
    apply = lambda state: sum(k @ state @ k.conj().T for k in channel_kraus)
    unread = lambda state: sum(k @ state @ k.conj().T for k in readout)
    rho = np.eye(72) / 72
    for n in (0, 1, 2):
        for initial, ep in zip(readout, effects):
            selective = initial @ rho @ initial.conj().T
            for _ in range(n):
                selective = apply(selective)
            for eq in effects:
                pulled = eq.copy()
                for _ in range(n):
                    pulled = sum(k.conj().T @ pulled @ k for k in channel_kraus)
                actual = np.trace(eq @ selective).real
                expected = np.trace(ep @ pulled).real / 72
                assert actual == pytest.approx(expected, abs=2e-16)
    selective = readout[0] @ rho @ readout[0].conj().T
    uninterrupted = np.trace(effects[0] @ apply(apply(selective))).real
    intervened = np.trace(effects[0] @ apply(unread(apply(selective)))).real
    q = load_response((0, 2), theta)["output"]
    off_norm = 2 * np.vdot(q[:d, d:], q[:d, d:]).real / 72
    assert uninterrupted - intervened == pytest.approx(off_norm / 32, abs=2e-16)
    if theta:
        assert uninterrupted - intervened > 1e-6
    else:
        assert abs(uninterrupted - intervened) < 2e-16
        np.testing.assert_allclose(q, 0, atol=2e-16, rtol=0)


@pytest.mark.parametrize("value", (True, False, 1.0, -1, 0, None, "2", np.ma.masked))
def test_invalid_private_dimensions_are_rejected(value):
    with pytest.raises(ValueError, match="positive integer"):
        occupied_port_effects(value)


@pytest.mark.parametrize("value", (True, 0.0, -1, 12, None, "0", np.ma.masked))
def test_invalid_port_labels_are_rejected(value):
    with pytest.raises(ValueError, match="port must"):
        sector_indices(value)


@pytest.mark.parametrize("value", (True, None, 1j, [], float("nan"), float("inf"),
                                   np.ma.array(0.01, mask=True), Fraction(-1, 100),
                                   Fraction(1, 24) + Fraction(1, 10**100)))
def test_invalid_endpoint_probabilities_are_rejected(value):
    with pytest.raises(ValueError, match="probability"):
        correlation_from_same_port(value)


@pytest.mark.parametrize("value", (True, None, float("nan"), Fraction(-1, 100),
                                   Fraction(1, 2) + Fraction(1, 10**100),
                                   np.ma.array(0.1, mask=False)))
def test_invalid_equal_tag_probabilities_are_rejected(value):
    with pytest.raises(ValueError, match="probability"):
        correlation_from_equal_ports(value)


def test_invalid_matrix_states_and_record_blocks_are_rejected():
    with pytest.raises(ValueError, match="66"):
        measure_occupied_port(np.eye(12) / 12)
    with pytest.raises(ValueError, match="trace one"):
        measure_occupied_port(np.eye(66))
    bad = np.eye(66) / 66
    bad[0, 0], bad[1, 1] = -0.1, 2 / 66 + 0.1
    with pytest.raises(ValueError, match="positive"):
        measure_occupied_port(bad)
    with pytest.raises(ValueError, match="masked"):
        measure_occupied_port(np.ma.array(np.eye(66) / 66, mask=False))
    for embedding in (forget_port_record, expand_port_record):
        with pytest.raises(ValueError, match="twelve"):
            embedding([np.eye(11)] * 11)
        with pytest.raises(ValueError, match="11"):
            embedding([np.eye(12)] * 12)
    # This legitimate nontracial preparation is outside the decoder's
    # tracial endpoint model: its same-port identity joint probability is1/4.
    pure = np.zeros((66, 66))
    pure[0, 0] = 1
    blocks = measure_occupied_port(pure)
    assert np.trace(blocks[0]) == Fraction(1, 2)
    with pytest.raises(ValueError, match="endpoint range"):
        correlation_from_same_port(Fraction(1, 4))


def test_measurement_refuses_erasing_a_subnormal_population():
    # This diagonal population survives the shared density validator. Its
    # two half-weight branches cannot be represented in binary64 separately.
    tiny = np.nextafter(0., 1.)
    state = np.zeros((66, 66))
    state[0, 0], state[1, 1] = 1., tiny
    from quantum_information.states import density_matrix
    assert density_matrix(state)[1, 1] == tiny
    with np.errstate(under="raise"):
        with pytest.raises(ValueError, match="branch underflows"):
            measure_occupied_port(state)


@pytest.mark.parametrize("imaginary", [False, True])
def test_unread_channel_refuses_erasing_a_fixed_private_operator(imaginary):
    # Within one load sector, arbitrary private matrix units are exact fixed
    # points of the forgotten-record channel. An unresolved numerical split
    # must not turn a nonzero fixed operator into an apparently exact zero.
    value = np.zeros((132, 132), complex)
    tiny = np.nextafter(0., 1.)
    if imaginary:
        value[0, 1] = 1j * tiny
    else:
        value[0, 1] = tiny
    with pytest.raises(ValueError, match="branch underflows"):
        unread_port_readout(value)
    # The next larger population splits exactly, so the guard does not
    # replace all subnormal entries with an arbitrary numerical cutoff.
    value *= 2
    np.testing.assert_array_equal(unread_port_readout(value), value)
