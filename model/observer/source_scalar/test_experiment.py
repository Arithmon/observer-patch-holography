"""Independent physical targets and complete source/field instrument checks."""
from decimal import Decimal, localcontext
from fractions import Fraction
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import field


def ladder(d):
    """Independent occupation matrix, explicitly populated by its action."""
    result = np.zeros((d, d))
    for occupied in range(1, d):
        result[occupied-1, occupied] = math.sqrt(occupied)
    return result


def physical_targets(d):
    # Form the physical squares BEFORE cutting down to the working code.
    a = ladder(d+2)
    q = (a+a.T)/math.sqrt(2)
    p = (a-a.T)/(1j*math.sqrt(2))
    q2, p2 = (q@q)[:d, :d], (p@p)[:d, :d]
    q, p = q[:d, :d], p[:d, :d]
    onsite = (p2+(41/32)*q2)/2
    return q, p, q2, p2, onsite


def physical_vacuum(d):
    """Build the normal-mode Fock series by creation-operator application."""
    create = ladder(d).T
    difference_creation = (np.kron(np.eye(d), create)-np.kron(create, np.eye(d)))/math.sqrt(2)
    blank = np.eye(d*d)[0]
    state = np.zeros(d*d)
    occupation = blank.copy()
    for n in range(d-1):
        if n % 2 == 0:
            pair = n//2
            amplitude = (-1/9)**pair*math.sqrt(math.comb(2*pair, pair))/2**pair
            state += amplitude*occupation
        occupation = difference_creation@occupation/math.sqrt(n+1)
    return state/np.linalg.norm(state)


def complex_vector(encoded):
    array = np.asarray(encoded)
    if array.ndim == 2 and array.shape[-1] == 2:
        return array[:, 0]+1j*array[:, 1]
    return array.astype(complex)


@pytest.mark.parametrize("d", (2, 4, 8, 16))
def test_projected_action_against_larger_physical_occupation_space(d):
    q, p, qq, pp, onsite = physical_targets(d)
    actual_q, actual_p = field.quadratures(d)
    actual_qq, actual_pp = field.projected_squares(d)
    np.testing.assert_allclose(actual_q, q, atol=1e-14)
    np.testing.assert_allclose(actual_p, p, atol=1e-14)
    np.testing.assert_allclose(actual_qq, qq, atol=1e-14)
    np.testing.assert_allclose(actual_pp, pp, atol=1e-14)
    expected = np.kron(np.eye(d), onsite)+np.kron(onsite, np.eye(d))-(9/32)*np.kron(q, q)
    np.testing.assert_allclose(field.physical_hamiltonian(d), expected, atol=2e-14)
    boundary = np.zeros((d, d))
    boundary[-1, -1] = d/2
    np.testing.assert_allclose(qq-q@q, boundary, atol=1e-14)
    np.testing.assert_allclose(pp-p@p, boundary, atol=1e-14)


def test_cutoff_countercontrol_rejects_naive_shear_and_missing_boundary():
    q, p, qq, pp, _ = physical_targets(2)
    np.testing.assert_allclose(q@q, np.eye(2)/2, atol=1e-14)
    np.testing.assert_allclose(p@p, np.eye(2)/2, atol=1e-14)
    np.testing.assert_allclose((qq+pp)/2, np.diag([.5, 1.5]), atol=1e-14)
    assert np.linalg.norm(qq@q-q@qq) > .9
    # The top-level commutator has the opposite sign from canonical CCR.
    np.testing.assert_allclose(q@p-p@q, 1j*np.diag([1., -1.]), atol=1e-14)
    initial = np.array([1., 0.])
    shear = expm(-.37j*(q@q))
    np.testing.assert_allclose(shear.conj().T@p@shear, p, atol=1e-14)
    assert np.linalg.norm((shear.conj().T@p@shear-p+.74*q)@initial) > .5


@pytest.mark.parametrize("d", (4, 8, 16))
def test_vacuum_target_and_rational_bound_against_independent_fock_series(d):
    state, metadata = field.vacuum_vector(d)
    np.testing.assert_allclose(state, physical_vacuum(d), atol=2e-15)
    assert abs(np.linalg.norm(state)-1) < 2e-15
    # This is an entangled site state, not independent onsite squeezing.
    singular = np.linalg.svd(state.reshape(d, d), compute_uv=False)
    assert singular[1] > .04
    for n1 in range(d):
        for n0 in range(d):
            if n0+n1 > d-2 or (n0+n1) % 2:
                assert state[n0+d*n1] == 0
    with localcontext() as context:
        context.prec = 80
        total = sum((Fraction(math.comb(2*n, n), 4**n*81**n) for n in range(d//2)), Fraction())
        retained = 4*Decimal(5).sqrt()/9*Decimal(total.numerator)/Decimal(total.denominator)
        tail = 1-retained
        declared_tail = Fraction(metadata["tail_probability_upper"])
        declared_distance = Fraction(metadata["infinite_vacuum_trace_distance_upper"])
        assert Decimal(declared_tail.numerator)/Decimal(declared_tail.denominator) >= tail
        assert declared_distance**2 >= declared_tail
        assert abs(Decimal(metadata["retained_probability"])-retained) < Decimal("2e-15")
        threshold = Fraction(1, 60) if d == 4 else Fraction(1, 6000)
        assert declared_distance < threshold


def test_infinite_reference_from_two_mode_covariance_and_response():
    stiffness = np.array([[41/32, -9/32], [-9/32, 41/32]])
    eigenvalues, rotation = np.linalg.eigh(stiffness)
    frequencies = np.sqrt(eigenvalues)
    covariance = rotation@np.diag(1/(2*frequencies))@rotation.T
    for time in (.25, .5, 1., 1.25):
        response = rotation@np.diag(np.sin(time*frequencies)/frequencies)@rotation.T
        for current in (-1, 0, 1):
            mean = current*response[1, 0]
            expected = (1+math.exp(-covariance[1, 1]/2)*math.sin(mean))/2
            result = field.ideal_prediction(time=time, current=current)
            assert abs(result["plus_probability"]-expected) < 1e-15
            assert abs(result["variance"]-covariance[1, 1]) < 1e-15


@pytest.mark.parametrize("d", (4, 8))
def test_exact_nested_commutator_certificate_against_orthonormal_matrices(d):
    from .bounds import nested_commutator_bounds, split_error_bound
    q, _, _, _, onsite = physical_targets(d)
    a = np.kron(np.eye(d), onsite)+np.kron(onsite, np.eye(d))
    b = -(9/32)*np.kron(q, q)
    commutator = a@b-b@a
    squared = [float(np.linalg.norm(matrix, "fro")**2) for matrix in
               (a@commutator-commutator@a, b@commutator-commutator@b)]
    exact = nested_commutator_bounds(d)
    np.testing.assert_allclose([float(x) for x in exact], squared, rtol=2e-14, atol=1e-12)
    bound = split_error_bound(d, 64)
    finer = split_error_bound(d, 128)
    assert Fraction(bound["unitary_operator_error_upper"]) == 4*Fraction(finer["unitary_operator_error_upper"])
    for value, key in zip(exact, ("onsite_nested_commutator_norm_upper", "bond_nested_commutator_norm_upper")):
        rounded = Fraction(bound[key])
        assert rounded**2 >= value
        assert (rounded-Fraction(1, 10**15))**2 < value
    half = expm(-1j*a/128)
    step = half@expm(-1j*b/64)@half
    error = np.linalg.norm(np.linalg.matrix_power(step, 64)-expm(-1j*(a+b)), 2)
    assert 0 < error < float(Fraction(bound["unitary_operator_error_upper"]))


def test_rational_split_bound_rejects_inexact_or_invalid_clock_inputs():
    from .bounds import split_error_bound
    for steps, time in ((True, Fraction(1)), (0, Fraction(1)), (4, 1.), (4, Fraction(-1))):
        with pytest.raises(ValueError):
            split_error_bound(4, steps, time)


@pytest.fixture(scope="module")
def source_histories():
    from .experiment import source_history
    return {current: source_history(4, 4, current) for current in (-1, 0, 1)}


def independent_terminal_branches(d, steps, current):
    q, _, _, _, onsite = physical_targets(d)
    dt = 1/steps
    half = expm(-.5j*dt*onsite)
    half_both = np.kron(half, half)
    bond = expm(1j*dt*(9/32)*np.kron(q, q))
    split = half_both@bond@half_both
    state = np.kron(np.eye(d), expm(1j*current*q))@physical_vacuum(d)
    state = np.linalg.matrix_power(split, steps)@state
    minus, plus = expm(-.5j*q), expm(.5j*q)
    branches = {}
    for outcome in (-1, 1):
        kraus = np.kron((minus-1j*outcome*plus)/2, np.eye(d))
        branches[outcome] = kraus@state
    return branches


def test_source_preparation_and_all_pointer_branches_against_direct_matrices(source_histories):
    target = physical_vacuum(4)
    for current, result in source_histories.items():
        prepared = complex_vector(result["prepared_state"])
        assert abs(np.vdot(prepared, prepared)-1) < 2e-10
        np.testing.assert_allclose(np.outer(prepared, prepared.conj()), np.outer(target, target), atol=2e-10)
        branches = {row["outcome"]: row for row in result["branches"]}
        assert set(branches) == {-1, 1}
        assert abs(sum(row["probability"] for row in branches.values())-1) < 2e-10
        for outcome, oracle in independent_terminal_branches(4, 4, current).items():
            row = branches[outcome]
            state = complex_vector(row["state"])
            probability = float(np.vdot(oracle, oracle).real)
            assert abs(row["probability"]-probability) < 2e-10
            assert abs(np.vdot(state, state)-1) < 2e-10
            # Complete conditional density, not just the detector count.
            np.testing.assert_allclose(row["probability"]*np.outer(state, state.conj()),
                                       np.outer(oracle, oracle.conj()), atol=2e-10)


def test_actual_detector_parity_and_nonzero_response(source_histories):
    def plus(current):
        return next(row["probability"] for row in source_histories[current]["branches"] if row["outcome"] == 1)
    assert abs(plus(0)-.5) < 2e-10
    assert abs(plus(1)+plus(-1)-1) < 2e-10
    assert plus(1) > .51


def test_repeated_source_routines_resist_tiny_backend_norm_bias(monkeypatch):
    """Reproduce a harmless per-routine drift that becomes visible on repetition."""
    from . import compiler, experiment
    execute = compiler.execute

    def biased_execute(*args, **kwargs):
        return (1-2e-12)*execute(*args, **kwargs)

    # Uncorrected, the four onsite uses plus bond on every step accumulate
    # probability loss far beyond the unchanged 2e-10 instrument check.
    assert 1-(1-2e-12)**(2*5*128) > 2e-9
    monkeypatch.setattr(compiler, "execute", biased_execute)
    experiment._executed.cache_clear()
    try:
        result = experiment.source_history(4, 128, 1)
        assert abs(result["post_read_total_probability"]-1) < 2e-10
        for row in result["executed_routine_roundoff_diagnostics"].values():
            assert row["projection_change"] > 1e-12
        oracle = independent_terminal_branches(4, 128, 1)
        for row in result["branches"]:
            state = complex_vector(row["state"])
            expected = oracle[row["outcome"]]
            np.testing.assert_allclose(row["probability"]*np.outer(state, state.conj()),
                                       np.outer(expected, expected.conj()), atol=2e-10)
    finally:
        experiment._executed.cache_clear()


def test_preparation_and_pointer_phase_cannot_be_omitted(source_histories):
    q, _, _, _, onsite = physical_targets(4)
    half = np.kron(expm(-.125j*onsite), expm(-.125j*onsite))
    step = half@expm(1j*(9/128)*np.kron(q, q))@half
    propagation = np.linalg.matrix_power(step, 4)@np.kron(np.eye(4), expm(1j*q))
    correct_read = np.kron((expm(-.5j*q)-1j*expm(.5j*q))/2, np.eye(4))
    wrong_read = np.kron((expm(-.5j*q)+expm(.5j*q))/2, np.eye(4))
    no_preparation = float(np.linalg.norm(correct_read@propagation@np.eye(16)[0])**2)
    no_pointer_phase = float(np.linalg.norm(wrong_read@propagation@physical_vacuum(4))**2)
    actual = next(row["probability"] for row in source_histories[1]["branches"] if row["outcome"] == 1)
    assert abs(actual-no_preparation) > 1e-4
    assert abs(actual-no_pointer_phase) > .1


@pytest.mark.parametrize("args", ((2, 4, 1), (4, True, 1), (4, 0, 1),
                                  (4, 4, True), (4, 4, 2), (4, 4, 1.)))
def test_source_history_rejects_wrong_code_clock_or_current(args):
    from .experiment import source_history
    with pytest.raises(ValueError):
        source_history(*args)


@pytest.mark.parametrize("dimension", (True, 0, 3, 4., 32))
def test_invalid_onsite_dimension_rejected(dimension):
    with pytest.raises(ValueError):
        field.vacuum_vector(dimension)
