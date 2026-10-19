"""Independent primitive, algebra, timeout and modular-isometry controls."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as c
from . import passive as p


def z(real=0, imag=0):
    return F(real), F(imag)


def dense(value, dimension=66):
    answer = np.zeros((dimension, dimension), dtype=complex)
    for (i, j), (real, imag) in value.items():
        answer[i, j] = complex(real, imag)
    return answer


@pytest.mark.parametrize('law', ['C', 'E'])
def test_independent_primitive_runtime_matches_complex_map_and_adjoint(law):
    # Independent existing primitive Kraus arrays, including cross-degree
    # coherence and a non-Hermitian matrix unit. No event-column oracle.
    x = {(0, 0): z(F(2, 7)), (1, 1): z(F(5, 7)),
         (0, 1): z(F(1, 9), F(2, 11)), (1, 0): z(F(1, 9), F(-2, 11)),
         (17, 38): z(F(-3, 19), F(1, 23))}
    instrument = c.marked(law).instrument
    for current in (0, 1, None):
        branches = (instrument.kraus[current] if current is not None
                    else (*instrument.kraus[0], *instrument.kraus[1]))
        actual = sum(k@dense(x)@k.conj().T for k in branches)
        adjoint = sum(k.conj().T@dense(x)@k for k in branches)
        np.testing.assert_allclose(dense(p.native_map(law, x, current=current)),
                                   actual, rtol=0, atol=2e-15)
        np.testing.assert_allclose(dense(p.native_map(law, x, current=current, adjoint=True)),
                                   adjoint, rtol=0, atol=2e-15)
    assert p.native_map(law, p.identity()) == p.identity()
    assert p.native_map(law, p.identity(), adjoint=True) == p.identity()


@pytest.mark.parametrize('law', ['C', 'E'])
def test_exact_hilbert_schmidt_duality_and_untouched_reference(law):
    # Four reference blocks test the extension of the same primitive, rather
    # than replacing the reference by a classical mixture.
    i, j = 0, 1
    blocks = {(0, 0): {(i, i): z(F(1, 2))},
              (1, 1): {(j, j): z(F(1, 2))},
              (0, 1): {(i, j): z(F(1, 2))},
              (1, 0): {(j, i): z(F(1, 2))}}
    state = np.zeros((132, 132), dtype=complex)
    expected = np.zeros_like(state)
    for (a, b), block in blocks.items():
        state[a*66:(a+1)*66, b*66:(b+1)*66] = dense(block)
        expected[a*66:(a+1)*66, b*66:(b+1)*66] = dense(p.native_map(law, block))
    actual = sum(np.kron(np.eye(2), k)@state@np.kron(np.eye(2), k.conj().T)
                 for group in c.marked(law).instrument.kraus for k in group)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-15)
    a = p._matrix({(0, 1): z(1, 2), (17, 17): z(F(2, 7))}, 66)
    b = p._matrix({(0, 1): z(3, -1), (17, 17): z(F(3, 11))}, 66)
    image = p._matrix(p.native_map(law, p._sparse(b)), 66)
    adjoint = p._matrix(p.native_map(law, p._sparse(a), adjoint=True), 66)
    assert ((a.H*image).trace()-(adjoint.H*b).trace()).expand() == 0


def test_modular_invariant_algebra_test_distinguishes_rotated_commutative_algebra():
    unit = p.identity(2)
    x = {(0, 1): z(1), (1, 0): z(1)}
    diagonal = {(0, 0): z(1)}
    reference = F(1, 3), F(2, 3)
    rotated = p.inspect_algebra((unit, x), reference, lambda a: a)
    assert not rotated['modular_invariant_for_every_parameter']
    centralizer = p.inspect_algebra((unit, diagonal), reference, lambda a: a)
    assert centralizer['modular_invariant_for_every_parameter']
    full = p.inspect_algebra((unit, diagonal, x, {(0, 1): z(0, 1), (1, 0): z(0, -1)}),
                             reference, lambda a: a)
    assert full['modular_eigenvalue_ratios'] == ['1/2', '1', '2']
    assert full['modular_invariant_for_every_parameter']


def test_isometry_alone_is_not_mistaken_for_modular_sufficiency():
    # For the tracial reference every modular action is identity, whereas
    # conjugation by the swap is a nonidentity Hilbert-Schmidt isometry.
    unit = p.identity(2)
    a = {(0, 0): z(1)}
    swap = lambda value: {(1-i, 1-j): parts for (i, j), parts in value.items()}
    result = p.inspect_algebra((unit, a), (F(1, 2), F(1, 2)), swap)
    assert swap(a) != a
    assert not result['strict_loss_obstructs_every_unitary_modular_parameter']
    assert not result['passing_isometry_check_is_a_sufficiency_claim']
    assert not result['CPTP_admissibility_checked_here']


@pytest.mark.parametrize('law', ['C', 'E'])
def test_actual_centralizer_loses_norm_despite_commuting_with_modular_reference(law):
    degree, reference = p.palm_data()
    projector = {(i, i): z(1) for i, value in enumerate(degree) if value == 8}
    result = p.inspect_algebra((p.identity(), projector), reference,
                              lambda a: p.native_map(law, a, adjoint=True))
    assert result['modular_invariant_for_every_parameter']
    assert result['modular_eigenvalue_ratios'] == ['1']
    assert result['strict_loss_obstructs_every_unitary_modular_parameter']
    assert F(result['basis_checks'][1]['image_HS_norm_squared']) < F(
        result['basis_checks'][1]['centered_HS_norm_squared'])
    scalar = p.inspect_algebra((p.identity(),), reference,
                              lambda a: p.native_map(law, a, adjoint=True))
    assert not scalar['strict_loss_obstructs_every_unitary_modular_parameter']


@pytest.mark.parametrize('cap', [0, 1, 2, 8])
def test_complete_palm_event_record_separates_modular_control_with_same_timeouts(cap):
    degree, reference = p.palm_data()
    c_law = p.event_history(cap, 'C')
    assert c_law == p.event_history(cap, 'E')
    target = p.event_history(cap, modular_control=True)
    assert sum(c_law.values()) == sum(target.values()) == 1
    assert c_law['failed_herald',] == F(28, 33)
    for first, mass, q in ((8, F(2, 5), F(13, 15)), (10, F(3, 5), F(5, 6))):
        assert c_law['timeout', first, cap] == F(5, 33)*mass*q**cap
        assert c_law['timeout', first, cap] == target['timeout', first, cap]
        initial = {(i, i): z(F(5, 33)*reference[i]) for i in range(66) if degree[i] == first}
        for step in range(1, cap+1):
            # Independent all-matrix-unit stopped-channel implementation.
            event = c.apply('C', initial, first_at=step)
            for last in (8, 10):
                assert c_law['complete', first, step, last] == sum(
                    value[0] for (i, j), value in event.items() if i == j and degree[i] == last)
            assert sum(c_law['complete', first, step, last] for last in (8, 10)) == sum(
                target['complete', first, step, last] for last in (8, 10))
    assert sum(abs(c_law[key]-target[key]) for key in c_law)/2 == (
        2-F(13, 15)**cap-F(5, 6)**cap)/33


def test_stationary_event_reference_is_not_stationary_native_attempt_reference():
    _, reference = p.palm_data()
    sigma = {(i, i): z(value) for i, value in enumerate(reference)}
    assert c.apply('C', sigma) == sigma
    moved = p.native_map('C', sigma)
    assert sum(abs(moved[i, i][0]-reference[i]) for i in range(66))/2 == F(1, 150)
    assert sum(abs(value-F(1, 66)) for value in reference)/2 == F(3, 55)
    assert p.orbit_gap(F(65, 66), 221) >= F(32467, 33000)
    assert p.orbit_gap(F(3, 55), 221) >= F(589, 11000)
    assert p.orbit_gap(F(65, 66), 0) == p.observer_HS_gap(0) == 0
    assert p.observer_HS_gap(1) == F(921, 20000)
    assert p.observer_HS_gap(10, F(1, 2)) < p.observer_HS_gap(10)


@pytest.mark.parametrize('basis,reference', [
    (({(0, 0): z(1)},), (F(1, 2), F(1, 2))),  # no ambient identity
    ((p.identity(2), p.identity(2)), (F(1, 2), F(1, 2))),
    ((p.identity(2), {(0, 1): z(1)}), (F(1, 2), F(1, 2))),  # no adjoint
    ((p.identity(3), {(0, 0): z(1), (1, 1): z(2)}), (F(1, 3),)*3),  # no products
    ((p.identity(2),), (F(1), F(0))),
    ((p.identity(2),), (F(1, 2), F(1, 3))),
])
def test_algebra_admission_rejects_invalid_domains(basis, reference):
    with pytest.raises(ValueError):
        p.inspect_algebra(basis, reference, lambda a: a)


def test_map_and_operator_inputs_are_strict_and_copy_safe():
    with pytest.raises(ValueError, match='unital'):
        p.inspect_algebra((p.identity(2),), (F(1, 2), F(1, 2)), lambda a: {})
    value = F(1, 3)
    original = {(0, 0): (value, F(0))}
    copied = p.operator(original)
    value._numerator = 2
    assert copied[0, 0][0] == F(1, 3)
    malformed = F(1, 2)
    malformed._numerator = 2
    for value in (True, 0.5, malformed):
        with pytest.raises(ValueError):
            p.operator({(0, 0): (value, F(0))})
        with pytest.raises(ValueError):
            p.orbit_gap(value, 1)
    for count in (-1, True, 1.0):
        with pytest.raises(ValueError):
            p.event_history(count)
        with pytest.raises(ValueError):
            p.observer_HS_gap(count)
    for coverage in (0, F(3, 2), True):
        with pytest.raises(ValueError):
            p.orbit_gap(F(1, 2), 2, coverage)
    with pytest.raises(ValueError):
        p.native_map('C', {}, current=True)
    with pytest.raises(ValueError):
        p.native_map('other', {})
