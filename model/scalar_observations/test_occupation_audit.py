"""Independent Fock-matrix and residual-quadrature audits of the cutoff bridge.

The reference below forms oscillator powers before any occupation projection.
It does not use the certificate's ladder recurrence or residual Gram formula.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, prod, sqrt

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.hermite import hermgauss, hermval
from scipy.linalg import cosm, expm, sinm

from . import occupation_exact as exact
from . import occupation_readout as readout
from . import occupation_adjoint as adjoint
from . import occupation_free as free
from . import occupation
from . import protocol
from .lattice import GalerkinParameters, SevenModeLattice
from ..observer import native


def _normalized(vector):
    return {key: complex(*value) * sqrt(prod(factorial(n) for n in key[:7]))
            / vector.denominator for key, value in vector.coefficients.items()}


def _padded_full_h(occupation, coupling):
    """One full H action, reconstructed using padded normalized matrices."""
    size = max(occupation) + 7
    annihilation = np.diag(np.sqrt(np.arange(1, size)), 1)
    q = (annihilation + annihilation.T) / 2
    p = 1j * (annihilation.T - annihilation)
    q2, p2, q4 = q @ q, p @ p, np.linalg.matrix_power(q, 4)
    # Independent torus quadrature fixes every quartic coefficient.
    points = np.asarray(list(product(range(5), repeat=3))) * (2*np.pi/5)
    modes = np.column_stack([np.ones(len(points)), *[
        f(points[:, axis]) * sqrt(2) for axis in range(3)
        for f in (np.cos, np.sin)]])
    out = {}

    def add(factors, multiplier):
        choices = [[(row, op[row, occupation[mode]])
                    for row in range(size) if abs(op[row, occupation[mode]]) > 1e-13]
                   for mode, op in factors]
        for entries in product(*choices):
            key = list(occupation)
            value = multiplier
            for (mode, _), (row, coefficient) in zip(factors, entries, strict=True):
                key[mode] = row
                value *= coefficient
            out[tuple(key)] = out.get(tuple(key), 0j) + value

    for mode in range(7):
        add([(mode, p2)], .5)
        add([(mode, q2)], .5 if mode == 0 else 1.)
        add([(mode, q4)], coupling*np.mean(modes[:, mode]**4)/24)
    for left, right in combinations(range(7), 2):
        add([(left, q2), (right, q2)],
            coupling*np.mean(modes[:, left]**2*modes[:, right]**2)/4)
    # Centering affects only a common Hamiltonian phase.
    out[occupation] = out.get(occupation, 0j) - (F(41, 8) + F(49, 128)*coupling)
    return {key: value for key, value in out.items() if abs(value) > 2e-13}


@pytest.mark.parametrize("coupling", [0, 1])
def test_full_h_boundary_action_matches_independent_padded_quadrature(coupling):
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    generator = exact.Generator(dimensions, "H", coupling)
    for occupation in ((0,)*7, (1, 3, 1, 0, 1, 0, 1)):
        monomial = exact.Vector({occupation: (1, 0)})
        input_norm = sqrt(prod(factorial(n) for n in occupation))
        expected = {key: input_norm*value
                    for key, value in _padded_full_h(occupation, coupling).items()}
        actual = {**_normalized(generator.apply(monomial)),
                  **_normalized(generator.apply(monomial, exterior=True))}
        for key in set(actual) | set(expected):
            assert actual.get(key, 0) == pytest.approx(expected.get(key, 0), abs=2e-12)
        inside = set(generator.apply(monomial).coefficients)
        outside = set(generator.apply(monomial, exterior=True).coefficients)
        assert inside.isdisjoint(outside)
        assert all(all(n < d for n, d in zip(key, dimensions, strict=True)) for key in inside)
        assert all(any(n >= d for n, d in zip(key, dimensions, strict=True)) for key in outside)


@pytest.mark.parametrize("kind,time", [("H", F(1, 20)), ("Q", F(-1, 3))])
def test_residual_integral_matches_direct_differentiation_and_quadrature(kind, time):
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    basis = [tuple([0, n, 0, 0, 0, 0, 0]) for n in range(4)]
    columns = []
    for occupation in basis:
        if kind == "H":
            columns.append(_padded_full_h(occupation, 1))
        else:
            n = occupation[1]
            column = {tuple([0, n+1, 0, 0, 0, 0, 0]): sqrt(n+1)/2}
            if n:
                column[tuple([0, n-1, 0, 0, 0, 0, 0])] = sqrt(n)/2
            columns.append(column)
    rows = sorted(set(basis).union(*(set(column) for column in columns)))
    full = np.asarray([[column.get(row, 0j) for column in columns] for row in rows])
    embedding = np.asarray([[int(row == column) for column in basis] for row in rows])
    finite = embedding.T @ full
    # Different trailing records remain orthogonal throughout the flow.
    vector = exact.Vector({basis[0]+(0,): (3, 1), basis[3]+(0,): (-1, 2),
                           basis[1]+(1,): (1, -2)}, 7)
    v = np.zeros((4, 2), complex)
    for key, value in _normalized(vector).items():
        v[key[1], key[7]] = value
    degree = 3
    powers = [v]
    for j in range(1, degree+1):
        powers.append((-1j*float(time)/j)*(finite @ powers[-1]))
    nodes, weights = leggauss(8)
    integrated = 0.
    for u, weight in zip((nodes+1)/2, weights/2, strict=True):
        value = sum(term*u**j for j, term in enumerate(powers))
        derivative = sum(j*powers[j]*u**(j-1)/float(time) for j in range(1, degree+1))
        residual = embedding @ derivative + 1j*full @ value
        integrated += weight * float(np.linalg.norm(residual)**2) * float(time)**2
    _, receipt = exact.polynomial_step(exact.Generator(dimensions, kind), vector,
                                       time, degree=degree, bits=40)
    assert float(receipt["full_integrated_squared_residual_times_duration"]) == pytest.approx(
        integrated, rel=3e-12, abs=2e-15)


def test_coherent_weyl_recombination_preserves_old_records_and_correct_read_phase():
    a = (0, 0, 0, 0, 0, 0, 0, 0)
    b = (0, 3, 0, 0, 0, 0, 0, 1)
    positive = exact.Vector({a: (2, 3), b: (-1, 2)}, 7)
    negative = exact.Vector({a: (1, -2), b: (3, 1)}, 11)
    actual = _normalized(exact._read_combine(positive, negative))
    plus, minus = _normalized(positive), _normalized(negative)
    for key in (a, b):
        # Actual Kraus maps include another common 1/sqrt(2), tracked by weight.
        assert actual[key+(0,)]/sqrt(2) == pytest.approx(
            ((1+1j)*plus[key]+(1-1j)*minus[key])/(2*sqrt(2)))
        assert actual[key+(1,)]/sqrt(2) == pytest.approx(
            ((1-1j)*plus[key]+(1+1j)*minus[key])/(2*sqrt(2)))
    assert exact.norm_squared(exact._read_combine(positive, negative)) == (
        exact.norm_squared(positive)+exact.norm_squared(negative))


def test_endpoint_rounding_charges_factorial_weight_and_deleted_coefficients():
    vector = exact.Vector({(0,)*7: (10, 3), (0, 4, 0, 0, 0, 0, 0): (1, -1)}, 97)
    rounded, bound = exact.quantize(vector, bits=3)
    before, after = _normalized(vector), _normalized(rounded)
    actual = sqrt(sum(abs(after.get(key, 0)-value)**2 for key, value in before.items()))
    assert actual <= float(bound)+2e-16
    assert actual > .05
    assert (0, 4, 0, 0, 0, 0, 0) not in rounded.coefficients
    assert exact.norm_squared(rounded) != 1


@pytest.mark.parametrize("coupling,sign", [(0, -1), (0, 1), (1, -1), (1, 1)])
def test_complete_certificate_contains_independent_finite_history(coupling, sign):
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    certificate = exact.certificate(dimensions, coupling, sign, degree=4, substeps=2, bits=32)
    lattice = SevenModeLattice(dimensions, GalerkinParameters(coupling=coupling))
    reference = protocol.history(lattice, sign=sign, source_execution=False)
    error = float(certificate["finite_vector_error_upper"])
    assert error < .001
    actual_half_l1 = 0.
    for row, observed in zip(certificate["branches"], reference, strict=True):
        assert row["bits"] == observed["bits"]
        lo, hi = map(float, row["finite_probability_interval"])
        assert lo-2e-13 <= observed["probability"] <= hi+2e-13
        actual_half_l1 += abs(observed["probability"]-float(row["raw_polynomial_probability"]))/2
    assert actual_half_l1 <= error+error*error/2
    assert abs(float(certificate["raw_polynomial_mass"])-1) <= 2*error+error*error
    assert certificate["raw_state_squared_weight"] == F(1, 4)
    assert certificate["occupied_polynomial_coefficients"] <= 16


def test_positive_native_current_has_exact_weight_for_a_named_record_event():
    counts = {j: 0 for j in (-1, 0, 1)}
    for state in native.states(binary=True):
        for coin in (0, 1):
            current = native.native(state, native.E, coin)[1][1]
            counts[current] += 1
    assert F(counts[1], sum(counts.values())) == F(5, 66)
    # The final native repair retains both field bits and first current, so
    # the event is observable on the final complete record alphabet.
    probabilities = [.13, .21, .29, .37]
    law = [{"bits": list(bits), "probability": probability}
           for bits, probability in zip(product((0, 1), repeat=2), probabilities, strict=True)]
    records = protocol.native_records([law, law])
    event = sum(row["probability"] for row in records
                if row["records"][1] > 0 and row["records"][2:4] == [1, 1])
    assert event == pytest.approx(float(F(5, 66))*.37, abs=1e-15)


def test_unrestricted_weyl_elements_match_independent_hermite_quadrature():
    nodes, weights = hermgauss(80)
    weights /= sqrt(np.pi)
    phase = np.exp(1j*nodes/sqrt(2))
    polynomials = [hermval(nodes, [0]*n+[1])/sqrt(2**n) for n in range(6)]
    for m, n in product(range(6), repeat=2):
        expected = np.dot(weights, polynomials[m]*polynomials[n]*phase)
        real, imaginary = readout.weyl_element(m, n)
        actual = complex(float(real), float(imaginary))*np.exp(-1/8)
        assert actual == pytest.approx(expected, abs=4e-12)
    lower, upper = readout.exp_minus_eighth(12)
    # Compare at an order whose interval remains visible at float precision.
    assert float(lower) <= np.exp(-1/8) <= float(upper)


def _readout_quadrature(vector, weight):
    """Integrate |sum c_n H_n(x)/sqrt(2^n)|^2 against Gaussian effects."""
    nodes, weights = hermgauss(80)
    weights /= sqrt(np.pi)
    groups = {}
    for key, (real, imaginary) in vector.coefficients.items():
        spectator = key[:1]+key[2:7]
        group = (key[7], spectator)
        groups.setdefault(group, np.zeros(len(nodes), complex))
        groups[group] += ((real+1j*imaginary)/vector.denominator
                         *hermval(nodes, [0]*key[1]+[1])/sqrt(2**key[1]))
    out = {(first, last): 0. for first in (0, 1) for last in (0, 1)}
    for (first, spectator), field in groups.items():
        density = float(weight)*prod(factorial(n) for n in spectator)*abs(field)**2
        for last in (0, 1):
            effect = (1+(1 if last else -1)*np.sin(nodes/sqrt(2)))/2
            out[first, last] += np.dot(weights, density*effect)
    return out


def test_exact_terminal_effect_preserves_complex_coherence_and_record_weights():
    vector = exact.Vector({
        (0, 0, 0, 0, 0, 0, 0, 0): (3, 2),
        (0, 3, 0, 0, 0, 0, 0, 0): (-2, 1),
        (0, 5, 0, 0, 0, 0, 0, 0): (1, -1),
        (2, 0, 0, 0, 0, 0, 0, 0): (1, 2),
        (2, 1, 0, 0, 0, 0, 0, 0): (-2, -1),
        (0, 0, 0, 0, 0, 0, 0, 1): (2, -3),
        (0, 1, 0, 0, 0, 0, 0, 1): (1, 2),
    }, 32)
    weight = F(1, 2)
    expected = _readout_quadrature(vector, weight)
    receipt = readout.terminal_probabilities(vector, weight, F(0))
    for row in receipt["branches"]:
        lo, hi = row["polynomial_effect_probability_interval"]
        value = expected[tuple(row["bits"])]
        assert float(lo)-3e-14 <= value <= float(hi)+3e-14
    assert float(receipt["raw_polynomial_mass"]) == pytest.approx(sum(expected.values()), abs=3e-14)
    assert receipt["raw_polynomial_mass"] != 1
    assert receipt["branches_renormalized"] is False
    assert receipt["terminal_quantum_output_certified"] is False


def test_terminal_effect_intervals_cover_a_normalized_true_flagged_state():
    # True amplitudes are (2+i)|0,0>+(1+2i)|0,1>, with physical weight1/10.
    # Every true bit-pair probability is1/4. Add explicit dyadic perturbations
    # in the original component, measured mode and a factorial-weighted spectator.
    vector = exact.Vector({
        (0, 0, 0, 0, 0, 0, 0, 0): (33, 16),
        (0, 0, 0, 0, 0, 0, 0, 1): (16, 32),
        (0, 1, 0, 0, 0, 0, 0, 0): (1, -2),
        (2, 0, 0, 0, 0, 0, 0, 1): (1, 0),
    }, 16)
    error = exact.sqrt_upper(F(1+5+2, 2560))
    receipt = readout.terminal_probabilities(vector, F(1, 10), error)
    assert receipt["raw_polynomial_mass"] != 1
    for row in receipt["branches"]:
        lo, hi = row["unrestricted_probability_interval"]
        assert lo <= F(1, 4) <= hi


@pytest.mark.parametrize("bit", [0, 1])
def test_reverse_read_recombination_requires_sqrt_two_error_factor(bit):
    # Four independent Weyl/pointer error columns, two existing pointer rows.
    # The forward fresh-pointer2x2 norm-one estimate cannot apply to this map.
    empty = exact.Vector({})
    matrix = np.zeros((2, 4), complex)
    for sign_index, old in product(range(2), repeat=2):
        pointer = [1, 1]
        pointer[bit] = old
        key = (0,)*7+tuple(pointer)
        column = exact.Vector({key: (1, 0)})
        positive, negative = (column, empty) if sign_index == 0 else (empty, column)
        transformed = _normalized(adjoint._reverse_read_combine(positive, negative, bit))
        for output, value in transformed.items():
            assert output[7+1-bit] == 1
            matrix[output[7+bit], 2*sign_index+old] = value
    np.testing.assert_allclose(np.linalg.svd(matrix, compute_uv=False), sqrt(2), atol=2e-15)
    # At Q=0 the two Weyl values coincide; the physical inverse is R(pi/4)^T.
    vector = np.asarray([2+1j, -1+3j])
    expected = np.asarray([[1, 1], [-1, 1]]) @ vector/sqrt(2)
    np.testing.assert_allclose(matrix @ np.tile(vector, 2)/sqrt(2), expected, atol=2e-15)


def test_full_reverse_instrument_matches_dense_unitary_without_discarding_bits():
    dimensions = (2, 4, 2, 2, 2, 2, 2)
    model = SevenModeLattice(dimensions)
    indices = [2*n for n in range(4)]
    h = model.hamiltonian[indices][:, indices].toarray()
    h -= float(F(705, 128))*np.eye(4)
    q = model.operators[1]["q"]
    angle = (q+np.pi/2*np.eye(4))/2
    k0, k1 = cosm(angle), sinm(angle)
    vector = exact.Vector({(0, n, 0, 0, 0, 0, 0, first, last):
                           (1+n-2*first, 1+last-n)
                           for n, first, last in product(range(4), range(2), range(2))}, 16)
    weight = F(1, 4)
    expected = np.zeros((4, 2, 2), complex)
    for key, value in _normalized(vector).items():
        expected[key[1], key[7], key[8]] = sqrt(float(weight))*value

    def inverse_read(state, bit):
        result = np.zeros_like(state)
        for other in (0, 1):
            left = state[:, 0, other] if bit == 0 else state[:, other, 0]
            right = state[:, 1, other] if bit == 0 else state[:, other, 1]
            zero, one = k0 @ left+k1 @ right, -k1 @ left+k0 @ right
            if bit == 0:
                result[:, 0, other], result[:, 1, other] = zero, one
            else:
                result[:, other, 0], result[:, other, 1] = zero, one
        return result

    interval, sign = F(1, 8), -1
    backwards = expm(1j*float(interval)*h)
    expected = inverse_read(expected, 1)
    expected = np.einsum("ij,jab->iab", backwards, expected)
    expected[:, 1, :] = expm(-1j*q/4) @ expected[:, 1, :]
    expected = inverse_read(expected, 0)
    expected = np.einsum("ij,jab->iab", backwards, expected)
    expected = np.einsum("ij,jab->iab", expm(-1j*sign*q/4), expected)
    result = adjoint.reverse_history(vector, weight, dimensions, 1, sign, interval,
                                     degree=8, substeps=4, bits=80)
    actual = np.zeros_like(expected)
    for key, value in _normalized(result["state"]).items():
        actual[key[1], key[7], key[8]] = sqrt(float(result["raw_state_squared_weight"]))*value
    np.testing.assert_allclose(actual, expected, atol=3e-12, rtol=0)
    assert result["raw_state_squared_weight"] == F(1, 16)
    assert np.linalg.norm(actual[:, 1, :]) > .1
    assert np.linalg.norm(actual[:, :, 1]) > .1
    prefixes = [stage["reverse_prefix_vector_error_upper"] for stage in result["stages"]]
    assert prefixes == sorted(prefixes)


def test_adjoint_event_bound_contains_independent_free_single_mode_observation():
    # In the free action, all six spectators factor out of every record
    # probability. A much larger, independent single-mode calculation supplies
    # the numerical reference without reusing any polynomial certificate code.
    size = 48
    annihilation = np.diag(np.sqrt(np.arange(1, size)), 1)
    q = (annihilation+annihilation.T)/2
    p = 1j*(annihilation.T-annihilation)
    h = (p @ p+2*q @ q)/2
    u = expm(-1j*h/8)
    angle = (q+np.pi/2*np.eye(size))/2
    kraus = cosm(angle), sinm(angle)
    kick = expm(1j*q/4)
    state = u @ kick[:, 0]
    probability = 0.
    for first in (0, 1):
        branch = kraus[first] @ state
        if first:
            branch = kick @ branch
        final = kraus[1-first] @ u @ branch
        probability += float(np.vdot(final, final).real)
    forward = exact.certificate((2, 4, 2, 2, 2, 2, 2), coupling=0,
                                degree=6, substeps=4, bits=64, return_terminal_state=True)
    certificate = adjoint.adjoint_certificate(forward)
    lo, hi = certificate["unrestricted_event_probability_interval"]
    assert float(lo) <= probability <= float(hi)
    assert certificate["event_probability_error_upper"] < certificate["ordinary_probability_error_upper"]
    assert certificate["forward_integrated_full_defect_and_jumps_upper"] >= (
        certificate["forward_endpoint_vector_error_upper"])


@pytest.mark.parametrize("dimension", [4, 12])
def test_reduced_free_generator_matches_full_padded_oscillator(dimension):
    size = dimension+6
    annihilation = np.diag(np.sqrt(np.arange(1, size)), 1)
    q = (annihilation+annihilation.T)/2
    p = 1j*(annihilation.T-annihilation)
    h = p @ p/2+q @ q-.75*np.eye(size)
    generator = free.ObservedFreeHamiltonian(dimension)
    for n in range(dimension):
        key = (0, n, 0, 0, 0, 0, 0, 1)
        vector = exact.Vector({key: (1, -1)})
        actual = {**_normalized(generator.apply(vector)),
                  **_normalized(generator.apply(vector, exterior=True))}
        for destination in range(size):
            expected = (1-1j)*sqrt(factorial(n))*h[destination, n]
            observed = actual.get((0, destination, 0, 0, 0, 0, 0, 1), 0)
            assert observed == pytest.approx(expected, abs=2e-10)
    with pytest.raises(ValueError, match="spectator"):
        generator.apply(exact.Vector({(2, 0, 0, 0, 0, 0, 0): (1, 0)}))


@pytest.mark.parametrize("sign", [-1, 1])
def test_factorized_free_records_match_seven_mode_prefix_and_unrestricted_effect(sign):
    # Four spectator levels permit real spectator evolution; freezing the
    # spectator state is not how the reference below obtains its probabilities.
    dimensions = (4,)*7
    model = SevenModeLattice(dimensions, GalerkinParameters(coupling=0))
    q = model.operators[1]["q"]
    kick = expm(1j*q/4)
    state = model.apply_local(model.vacuum_state(), [1], expm(1j*sign*q/4))
    state = model.evolve(state, .125)
    angle = (q+np.pi/2*np.eye(4))/2
    maps = cosm(angle), sinm(angle)
    branches = np.column_stack([model.apply_local(state, [1], k) for k in maps])
    branches[:, 1] = model.apply_local(branches[:, 1], [1], kick)
    branches = model.evolve(branches, .125)
    # Integrate the unrestricted final effects in the finite Hermite basis.
    nodes, weights = hermgauss(80)
    weights /= sqrt(np.pi)
    modes = np.asarray([hermval(nodes, [0]*n+[1])/sqrt(2**n*factorial(n)) for n in range(4)])
    effects = [np.einsum("mx,nx,x->mn", modes, modes,
                        weights*(1+sgn*np.sin(nodes/sqrt(2)))/2) for sgn in (-1, 1)]
    expected = {}
    for first in (0, 1):
        tensor = branches[:, first].reshape(dimensions, order="F")
        matrix = np.moveaxis(tensor, 1, 0).reshape(4, -1)
        density = matrix @ matrix.conj().T
        for last in (0, 1):
            expected[first, last] = float(np.trace(effects[last] @ density).real)
    spectator_blank_indices = [4*n for n in range(4)]
    spectator_excitations = (float(np.linalg.norm(branches)**2)
                            -float(np.linalg.norm(branches[spectator_blank_indices])**2))
    assert spectator_excitations > .01
    certificate = free.free_prefinal_certificate(4, sign=sign, degree=8, substeps=4,
                                                 bits=80, return_state=True)
    # This comparison uses finite-prefix error because the seven-mode
    # reference has the same observed finite prefix and unrestricted final effect.
    finite = readout.terminal_probabilities(certificate["state"],
                                           certificate["raw_state_squared_weight"],
                                           certificate["finite_vector_error_upper"])
    assert certificate["finite_vector_error_upper"] < F(1, 10**8)
    for row in finite["branches"]:
        lower, upper = row["unrestricted_probability_interval"]
        assert float(lower)-5e-13 <= expected[tuple(row["bits"])] <= float(upper)+5e-13
    assert certificate["scope"]["unobserved_modes_traced_out_exactly"] == 6
    assert certificate["scope"]["seven_mode_quantum_output_state_certified"] is False


@pytest.mark.parametrize("sign", [-1, 1])
def test_factorized_free_certificate_covers_large_independent_observed_mode(sign):
    size = 48
    annihilation = np.diag(np.sqrt(np.arange(1, size)), 1)
    q = (annihilation+annihilation.T)/2
    p = 1j*(annihilation.T-annihilation)
    u = expm(-1j*(p @ p/2+q @ q)/8)
    kick = expm(1j*q/4)
    state = u @ expm(1j*sign*q/4)[:, 0]
    angle = (q+np.pi/2*np.eye(size))/2
    maps = cosm(angle), sinm(angle)
    expected = {}
    for first in (0, 1):
        branch = maps[first] @ state
        if first:
            branch = kick @ branch
        for last in (0, 1):
            final = maps[last] @ u @ branch
            expected[first, last] = float(np.vdot(final, final).real)
    certificate = free.free_prefinal_certificate(12, sign=sign)
    assert certificate["unrestricted_vector_error_upper"] < F(1, 10**4)
    for row in certificate["terminal_effect_certificate"]["branches"]:
        lower, upper = row["unrestricted_probability_interval"]
        assert float(lower) <= expected[tuple(row["bits"])] <= float(upper)


@pytest.mark.parametrize("arguments", [{"substeps": 0}, {"substeps": True},
                                       {"degree": 0}, {"bits": 0}])
def test_no_zero_work_helper_can_issue_a_zero_error_flow_certificate(arguments):
    generator = exact.Generator((2,)*7, "Q")
    with pytest.raises(ValueError):
        exact.flow(generator, exact.Vector({(0,)*7: (1, 0)}), F(1, 4), **arguments)


@pytest.mark.parametrize("time", [.125, True, -.5])
def test_exact_flow_rejects_inexact_or_boolean_time(time):
    generator = exact.Generator((2,)*7, "Q")
    vector = exact.Vector({(0,)*7: (1, 0)})
    with pytest.raises(ValueError):
        exact.flow(generator, vector, time)
    with pytest.raises(ValueError):
        exact.polynomial_step(generator, vector, time)


@pytest.mark.parametrize("mass", [F(4, 5), F(6, 5)])
def test_event_bound_charges_raw_mass_defect_in_addition_to_half_l1(mass):
    # The normalized true distribution is certain on01. Its raw comparison
    # has mass0.8 or1.2 there. Half-L1 alone is0.1 and cannot enclose the
    # event's0.2 discrepancy; the exact mass correction supplies the other0.1.
    proof = {"branches": [{"bits": list(bits),
                            "polynomial_effect_probability_interval":
                            [mass, mass] if bits == (0, 1) else [F(0), F(0)]}
                           for bits in product((0, 1), repeat=2)],
             "raw_polynomial_mass": mass,
             "raw_probability_half_l1_error_upper": F(1, 10)}
    event = occupation.terminal_event(proof)
    assert event["raw_mass_correction"] == F(1, 10)
    assert event["event_error_upper"] == F(1, 5)
    lower, upper = event["unrestricted_event_probability_interval"]
    assert lower <= 1 <= upper
    assert event["branches_renormalized"] is False
