"""Independent exact and spectral challenges to the local-clock field join."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as sp

from ..observer import native
from . import core as c


@pytest.fixture(scope='module')
def matrices():
    data = c.source()
    a = np.array(data.remote_counts, dtype=float)/data.denominator
    b = np.array(data.incident_counts, dtype=float)/data.denominator
    d = np.eye(66)-a
    values, vectors = np.linalg.eigh(d)
    inverse_root = (vectors/np.sqrt(values))@vectors.T
    transfer = inverse_root@b@inverse_root
    exact = c.completed()
    completed = np.array([[float(F(value, exact.denominator)) for value in row] for row in exact.counts])
    return a, b, d, transfer, completed


def test_native_split_and_triangle_are_exact_positive_source_operators():
    data = c.source()
    assert data.region == (0, 1) and data.triangle == (0, 1, 5)
    assert len(data.incident_edges) == 9 and len(data.remote_edges) == 21
    assert data.tick_probability == F(3, 10)
    triangle = sp.Matrix(data.triangle_counts)
    assert (triangle-3*sp.eye(66))*(triangle-6*sp.eye(66)) == sp.zeros(66)
    facts = c.spectral_certificate()
    assert facts['native_attempt_transfer_lower'] == F(2, 3)
    assert facts['D_upper'] == F(19, 30)
    assert facts['effective_transfer_lower'] == F(3, 38)
    for value in (data, c.completed(), c.field_law()):
        assert not hasattr(value, '__dict__')
        with pytest.raises(AttributeError):
            value.source_digest = 'mutated'
    with pytest.raises(TypeError):
        data.remote_counts[0][0] = 0


def test_field_controller_needs_only_disclosed_local_loads_and_actual_current():
    data = c.source()
    cases = 0
    for initial in data.states:
        local = tuple(initial[p] for p in data.region)
        for seam in data.seams:
            for coin in (0, 1):
                final, (_, current) = native.native(initial, seam, coin)
                observed = c.update_field(local, seam, current)
                assert observed == tuple(final[p] for p in data.region)
                if not set(seam).intersection(data.region):
                    assert observed == local
                cases += 1
    assert cases == 3960
    # This boundary current changes the disclosed field. Discarding it
    # gives a different field even though its other endpoint is remote.
    assert c.update_field((1, 0), (0, 5), 1) == (0, 0)


def test_exact_completed_kernel_retains_remote_evolution_in_the_right_order():
    data, kernel = c.source(), c.completed()
    a = sp.Matrix(data.remote_counts)/data.denominator
    b = sp.Matrix(data.incident_counts)/data.denominator
    k = sp.Matrix(kernel.counts)/kernel.denominator
    assert k*(sp.eye(66)-a) == b
    assert (sp.eye(66)-a)*k != b
    assert k*sp.ones(66, 1) == sp.ones(66, 1)
    assert sp.ones(1, 66)*k == sp.ones(1, 66)
    assert k != k.T
    assert all(value >= 0 for value in k)


def test_completed_field_laws_preserve_every_record_and_full_native_population():
    base, changed = c.field_law(), c.field_law(feedback=True)
    expected = tuple((initial, word) for initial in product((0, 1), repeat=2)
                     for word in product((0, 1), repeat=4))
    assert base.records == changed.records == expected
    assert len(base.records) == 64
    for law in (base, changed):
        assert all(len(state) == 66 and all(type(value) is int and value >= 0 for value in state)
                   for state in law.counts)
        assert sum(map(sum, law.counts)) == law.denominator
        probabilities = law.probabilities()
        assert sum(probabilities.values()) == 1
        marginal = {initial: sum(p for (label, _), p in probabilities.items() if label == initial)
                    for initial in product((0, 1), repeat=2)}
        assert marginal == {(0, 0): F(15, 22), (0, 1): F(5, 33),
                            (1, 0): F(5, 33), (1, 1): F(1, 66)}
    assert c.total_variation(base, changed) == F(5, 99)
    assert c.field_statistics(base)['joint_cumulant'] != 0
    assert c.field_statistics(base)['means'] == (F(1, 6),)*4
    assert c.field_statistics(changed)['expected_native_attempts'] == F(27, 2)
    assert c.field_statistics(changed)['expected_incident_ticks'] == F(25, 6)
    assert c.read_clock((1, 0, 1, 0), feedback=True) == (1, 3, 4, 5)
    assert c.read_clock((0, 1, 0, 1), feedback=True) == (1, 2, 3, 4)
    probabilities = base.probabilities()
    probabilities.clear()
    assert len(c.field_law().probabilities()) == 64


def _independent_cylinders(kernel, ports, feedback):
    data = c.source()
    intervention = np.array(data.internal_counts, dtype=float)/2
    branches = {(label, ()): np.array([int(tuple(state[p] for p in data.region) == label)/66
                                     for state in data.states])
                for label in product((0, 1), repeat=2)}
    for cursor, port in enumerate(ports):
        future = {}
        for (initial, word), state in branches.items():
            evolved = kernel@state
            for bit in (0, 1):
                retained = evolved*np.array([int(s[port] == bit) for s in data.states])
                if feedback and cursor == 0 and bit:
                    retained = intervention@retained
                future[initial, (*word, bit)] = retained
        branches = future
    return branches


@pytest.mark.parametrize('feedback', (False, True))
def test_full_adaptive_cylinder_dictionary_matches_the_symmetric_effective_model(matrices, feedback):
    _, _, _, transfer, native_kernel = matrices
    exact = c.field_law(feedback=feedback)
    native = _independent_cylinders(native_kernel, exact.read_ports, feedback)
    effective = _independent_cylinders(transfer, exact.read_ports, feedback)
    for record, expected in exact.probabilities().items():
        assert abs(native[record].sum()-float(expected)) < 4e-15
        assert abs(effective[record].sum()-float(expected)) < 4e-15
    # The equality concerns observable cylinder laws. It is not equality of
    # hidden native and reconstructed terminal configuration distributions.
    assert max(np.linalg.norm(native[key]-effective[key]) for key in native) > 1e-4


def test_effective_transfer_is_positive_stochastic_with_a_bounded_logarithm(matrices):
    _, _, _, transfer, _ = matrices
    eigenvalues = np.linalg.eigvalsh(transfer)
    assert eigenvalues.min() >= float(c.spectral_certificate()['effective_transfer_lower'])-1e-14
    assert eigenvalues.max() <= 1+1e-14
    assert transfer.min() >= -2e-15
    np.testing.assert_allclose(transfer, transfer.T, atol=2e-15, rtol=0)
    np.testing.assert_allclose(transfer.sum(axis=0), 1, atol=5e-15, rtol=0)
    degree = 128
    polynomial = sum((1-eigenvalues)**j/j for j in range(1, degree+1))
    assert max(abs(-np.log(eigenvalues)-polynomial)) <= float(c.logarithm_remainder(degree))+1e-13
    assert c.logarithm_remainder(256) < F(1, 10**9)


def test_outside_field_has_no_automatic_same_dictionary(matrices):
    _, _, d, transfer, native_kernel = matrices
    facts = c.dictionary_certificate()
    assert facts['outside_port'] == 2 and facts['outside_commutator_exact'] != 0
    outside = np.diag([state[2] for state in c.source().states])
    assert np.linalg.norm(d@outside-outside@d) > .1
    native = _independent_cylinders(native_kernel, (2, 2, 2), False)
    effective = _independent_cylinders(transfer, (2, 2, 2), False)
    assert sum(abs(native[key].sum()-effective[key].sum()) for key in native)/2 > 1e-5


def test_identity_holding_cancels_but_real_remote_evolution_does_not():
    data, kernel = c.source(), c.completed()
    a = sp.Matrix(data.remote_counts)/data.denominator
    b = sp.Matrix(data.incident_counts)/data.denominator
    activity = sp.Rational(1, 2)
    held = activity*b*(sp.eye(66)-((1-activity)*sp.eye(66)+activity*a)).inv(method='DM')
    assert held == sp.Matrix(kernel.counts)/kernel.denominator
    assert c.holding_control()['complete_kernel_unchanged']
    assert c.scalar_remote_control()['four_read_record_total_variation'] > F(3, 1000)


def test_finite_native_stopping_keeps_its_complete_timeout_state():
    data = c.source()
    a = sp.Matrix(data.remote_counts)/data.denominator
    b = sp.Matrix(data.incident_counts)/data.denominator
    cap = 3
    success = sum((b*a**n for n in range(cap)), sp.zeros(66))
    timeout = a**cap
    assert sp.ones(1, 66)*(success+timeout) == sp.ones(1, 66)
    assert sp.ones(1, 66)*timeout == (1-data.tick_probability)**cap*sp.ones(1, 66)
    # A selected exact initial state challenges the factor order, which
    # stationary scalar two-point functions alone need not distinguish.
    wrong = sum((a**n*b for n in range(cap)), sp.zeros(66))
    assert success != wrong


def test_actual_change_clock_has_an_exact_negative_mode_instead():
    control = c.movement_clock_control()
    assert control['completed_adjoint_eigenvalue'] == '-1'
    matrix = sp.Matrix(c.movement.event_block('C'))
    sign = sp.Matrix([1-2*state[c.movement.CLOCK_PORT] for state in c.quantum.STATES])
    assert matrix.T*sign == -sign


@pytest.mark.parametrize('ports', ([0, 1], (), (True,), (2,), (0,)*7))
def test_read_contract_fails_closed_even_after_the_exact_cache_is_warm(ports):
    c.field_law()
    with pytest.raises(ValueError):
        c.field_law(ports)


@pytest.mark.parametrize('feedback', (0, 1, None, 'yes'))
def test_feedback_contract_is_not_a_truthy_flag(feedback):
    with pytest.raises(ValueError):
        c.field_law(feedback=feedback)


@pytest.mark.parametrize('region', ((True, 1), (0, 0), (0, 12), [0, 1], (0,)))
def test_public_field_counter_rejects_ambiguous_regions(region):
    with pytest.raises(ValueError):
        c.update_field((1, 0), (0, 1), 1, region=region)
