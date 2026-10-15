"""Same-history algebra attachments and controls against false regional autonomy."""
import numpy as np
import pytest

from . import algebra as a


X = np.array([[0., 1.], [1., 0.]])
Z = np.diag([1., -1.])
P0, P1 = np.diag([1., 0.]), np.diag([0., 1.])


def mixed_view():
    rho = np.array([[.3, .1j, 0, 0], [-.1j, .2, 0, 0],
                    [0, 0, .2, .08], [0, 0, .08, .3]], complex)
    return a.CQView(((0,), (1,), (2,)), (.3*rho, .7*np.eye(4)/4, np.zeros((4, 4))), 4)


def test_restriction_keeps_actual_record_weights_and_observable_pairing():
    view = mixed_view()
    probe = np.array([[.3, .5j], [.4-.2j, -.7]])
    for region in ((0,), (1,), (), (1, 0)):
        local = a.restrict(view, (2, 2), region)
        assert local.probabilities == pytest.approx(view.probabilities)
        operator = probe if local.dimension == 2 else np.eye(local.dimension)
        assert a.expectation(local, operator) == pytest.approx(
            a.expectation(view, a.embed_operator(operator, (2, 2), region)), abs=2e-15)
    assert a.public_expectation(view, {(0,): 2., (1,): -1., (2,): 800.}) == pytest.approx(-.1)
    assert not view.blocks[0].flags.writeable
    with pytest.raises(ValueError):
        view.blocks[0].setflags(write=True)


def test_supplied_factor_order_uses_actual_tensor_indices_and_disjoint_regions_commute():
    factors = (2, 3, 2)
    vector = np.zeros(12)
    vector[np.ravel_multi_index((1, 2, 0), factors)] = 1
    rho = np.outer(vector, vector)
    np.testing.assert_array_equal(a.partial_trace(rho, factors, (0,)), P1)
    np.testing.assert_array_equal(a.partial_trace(rho, factors, (2,)), P0)
    left = a.embed_operator(X, factors, (0,))
    right = a.embed_operator(Z, factors, (2,))
    np.testing.assert_array_equal(left@right, right@left)
    np.testing.assert_array_equal(a.embed_operator(np.eye(4), factors, (2, 0)), np.eye(12))


def test_conditional_expectation_is_cp_bimodule_projection_but_not_homomorphism():
    # Independent Kraus implementation of E(X)=Tr_B(X) tensor I_B / 2.
    rng = np.random.default_rng(61015)
    raw = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    positive = raw@raw.conj().T
    kraus = [np.kron(np.eye(2), np.outer(np.eye(2)[i], np.eye(2)[j]))/np.sqrt(2)
             for i in range(2) for j in range(2)]
    enlarged = sum(np.kron(k, np.eye(2))@positive@np.kron(k.conj().T, np.eye(2)) for k in kraus)
    assert np.linalg.eigvalsh(enlarged)[0] > -1e-12
    matrix = raw[:4, :4]
    direct = sum(k@matrix@k.conj().T for k in kraus)
    expected = a.conditional_expectation(matrix, (2, 2), (0,))
    np.testing.assert_allclose(direct, expected, atol=1e-14)
    np.testing.assert_allclose(a.conditional_expectation(expected, (2, 2), (0,)), expected)
    left, right = np.kron(X, np.eye(2)), np.kron(Z, np.eye(2))
    np.testing.assert_allclose(a.conditional_expectation(left@matrix@right, (2, 2), (0,)),
                               left@expected@right, atol=1e-14)
    erased = np.kron(np.eye(2), X)
    np.testing.assert_array_equal(a.conditional_expectation(erased, (2, 2), (0,)), np.zeros((4, 4)))
    np.testing.assert_array_equal(a.conditional_expectation(erased@erased, (2, 2), (0,)), np.eye(4))


def test_direct_sum_gns_uses_actual_mixed_blocks_and_retains_zero_record():
    view = mixed_view()
    gns = a.gns(view)
    operators = {record: np.eye(4)*(i+1)+np.kron(X, Z)
                 for i, record in enumerate(view.records)}
    assert gns.components[2] is None
    assert gns.dimension == 32
    assert gns.expectation(operators) == pytest.approx(a.expectation(view, operators), abs=2e-15)
    square = {record: operator.conj().T@operator for record, operator in operators.items()}
    assert gns.inner(operators, operators) == pytest.approx(a.expectation(view, square), abs=1e-14)
    for component in gns.components:
        if component is not None:
            with pytest.raises(ValueError):
                component.root_columns.setflags(write=True)


def test_incompatible_independent_numerical_rank_cuts_cannot_claim_an_isometry():
    # The small eigenvalues are individually below the common GNS rank cut,
    # but their regional sum lies above it. Blindly returning J would give
    # ||J*J-I||=1, despite the mathematically exact GNS inclusion existing.
    density = np.diag([1-1.2e-12, 0., .6e-12, .6e-12])
    view = a.CQView(((),), (density,), 4)
    with pytest.raises(ValueError, match='support cuts'):
        a.gns_embedding(view, (2, 2), (0,))


@pytest.mark.parametrize('pure', (False, True))
def test_regional_gns_is_an_isometry_into_the_same_common_state(pure):
    view = mixed_view() if not pure else a.CQView(((),), (np.diag([1., 0., 0., 0.]),), 4)
    for region in ((0,), (1,)):
        inclusion, local, common = a.gns_embedding(view, (2, 2), region)
        np.testing.assert_allclose(inclusion.conj().T@inclusion, np.eye(local.dimension), atol=2e-14)
        np.testing.assert_allclose(inclusion@local.lift(np.eye(2)), common.lift(np.eye(4)), atol=2e-14)
        np.testing.assert_allclose(inclusion@local.lift(X+1j*Z),
                                   common.lift(a.embed_operator(X+1j*Z, (2, 2), region)), atol=2e-14)
        if pure:
            # A rank-one regional state has two GNS dimensions, not a support corner of size one.
            assert local.dimension == 2


def test_modular_flow_uses_actual_faithful_blocks_fixes_records_and_selects_no_clock():
    view = a.CQView(((0,), (1,)), (.2*np.diag([.7, .3]), .8*np.diag([.4, .6])), 2)
    first = a.modular_flow(view, X, .3)
    composed = a.modular_flow(view, first, -.8)
    direct = a.modular_flow(view, X, -.5)
    for record in view.records:
        np.testing.assert_allclose(composed[record], direct[record], atol=1e-14)
    assert a.expectation(view, first) == pytest.approx(a.expectation(view, X), abs=1e-14)
    public = {(0,): np.eye(2), (1,): -np.eye(2)}
    for record, matrix in a.modular_flow(view, public, 412.).items():
        np.testing.assert_allclose(matrix, public[record], atol=1e-14)
    with pytest.raises(ValueError, match='faithful'):
        a.modular_flow(mixed_view(), np.eye(4), .2)
    with pytest.raises(ValueError, match='faithful'):
        a.modular_flow(a.CQView(((),), (P0,), 2), X, .2)


def test_equal_local_states_can_have_opposite_interacting_future_observations():
    # Bell plus/minus have identical marginals. Inverse Bell preparation
    # converts their correlation phase into an exactly distinguishable first bit.
    plus = np.array([1., 0., 0., 1.])/np.sqrt(2)
    minus = np.array([1., 0., 0., -1.])/np.sqrt(2)
    views = [a.CQView(((),), (np.outer(v, v),), 4) for v in (plus, minus)]
    locals_ = [a.restrict(view, (2, 2), (0,)) for view in views]
    np.testing.assert_allclose(locals_[0].blocks[0], locals_[1].blocks[0])
    cnot = np.eye(4)[[0, 1, 3, 2]]
    hadamard = np.array([[1., 1.], [1., -1.]])/np.sqrt(2)
    unitary = np.kron(hadamard, np.eye(2))@cnot
    instrument = ((np.kron(P0, np.eye(2))@unitary,), (np.kron(P1, np.eye(2))@unitary,))
    laws = [a.forecast(view, instrument) for view in views]
    assert laws[0] == pytest.approx({(0,): 1., (1,): 0.}, abs=1e-14)
    assert laws[1] == pytest.approx({(0,): 0., (1,): 1.}, abs=1e-14)
    assert a.forecast(locals_[0], ((P0,), (P1,))) == pytest.approx(
        a.forecast(locals_[1], ((P0,), (P1,))))
    with pytest.raises(ValueError, match='2 by 2'):
        a.forecast(locals_[0], instrument)


def test_forecast_retains_each_old_record_and_rejects_incomplete_or_nonfinite_maps():
    view = mixed_view()
    instrument = ((np.kron(P0, np.eye(2)),), (np.kron(P1, np.eye(2)),))
    law = a.forecast(view, instrument)
    for record, mass in view.probabilities.items():
        assert sum(law[(*record, bit)] for bit in (0, 1)) == pytest.approx(mass)
    assert len(law) == 6
    with pytest.raises(ValueError, match='complete'):
        a.forecast(view, instrument[:1])
    with pytest.raises(ValueError, match='finite'):
        a.forecast(view, ((np.full((4, 4), np.nan),),))
    with pytest.raises(ValueError, match='finite'):
        a.forecast(view, ((np.eye(4)*1e308,),))


def test_tiny_retained_branch_has_its_own_gns_without_complex_reciprocal_overflow():
    mass = np.nextafter(0., 1.)
    view = a.CQView(((0,), (1,)), (P0*mass, P1), 2)
    assert view.probabilities[(0,)] == mass
    gns = a.gns(view)
    assert gns.components[0].rank == 1
    assert gns.dimension == 4
    assert np.isfinite(gns.lift(X)).all()


@pytest.mark.parametrize('factors,region', (((2, 2), (0, 0)), ((True, 4), (0,)),
                                          ((2, 2), (True,)), ((4,), (1,)), ((4,), [0])))
def test_ambiguous_factorizations_are_rejected(factors, region):
    with pytest.raises(ValueError):
        a.restrict(mixed_view(), factors, region)


@pytest.mark.parametrize('value', (np.diag([1.1, -.1]), np.eye(2), [[1., 1.], [0., 0.]],
                                   np.ones((2, 2), dtype=bool), [[True, 0.], [0., 0.]],
                                   np.ma.array(np.eye(2)/2)))
def test_invalid_states_are_rejected_without_normalization_or_kernel_filling(value):
    with pytest.raises(ValueError):
        a.CQView(((),), (value,), 2)


def test_same_runtime_history_reference_restriction_and_continuation_forecast():
    from .runtime import Event, Instrument, Program, execute
    # An untouched reference is entangled with the work; it remains outside access.
    bell = np.array([1., 0., 0., 1.])/np.sqrt(2)
    read = Instrument(kraus=((P0,), (P1,)), visible=True)
    program = Program('algebra-test', 2, np.eye(2)/2,
                      events=(Event('read', instruments=(read,)),))
    initial = execute(program, initial=np.outer(bell, bell), stop_after=0)
    predicted = a.forecast(a.observer_view(initial), read.kraus)
    actual = execute(program, checkpoint=initial)
    view = a.observer_view(actual)
    assert predicted == pytest.approx(view.probabilities, abs=1e-14)
    assert view.dimension == 2
    assert actual.reference_dimension == 2
    assert a.public_expectation(view, {(0,): 1., (1,): -1.}) == pytest.approx(0.)
    np.testing.assert_allclose(view.blocks[0], P0/2, atol=1e-14)
    np.testing.assert_allclose(view.blocks[1], P1/2, atol=1e-14)


def test_regional_observers_share_the_actual_source_compiled_interacting_history():
    from fractions import Fraction
    from ..common_source import open_environment, protocol
    from . import common
    program = protocol.continuation_program(2, 2)
    env = open_environment.Environment((open_environment.PhaseNoise(1, 0, Fraction(1, 64), 4),))
    history = common.execute(program, env)
    complete = a.observer_view(history)
    # Existing source basis is scalar1 tensor scalar0 tensor gauge/matter/link.
    local = a.observer_view(history, (2, 2, 4), (1,))
    assert len(complete.records) == len(local.records) == 16
    assert local.probabilities == pytest.approx(complete.probabilities, abs=2e-14)
    probe = X+.3*Z
    assert a.expectation(local, probe) == pytest.approx(
        a.expectation(complete, a.embed_operator(probe, (2, 2, 4), (1,))), abs=2e-14)
    inclusion, regional_gns, common_gns = a.gns_embedding(complete, (2, 2, 4), (1,))
    np.testing.assert_allclose(inclusion.conj().T@inclusion,
                               np.eye(regional_gns.dimension), atol=3e-10)
    np.testing.assert_allclose(inclusion@regional_gns.lift(probe),
                               common_gns.lift(a.embed_operator(probe, (2, 2, 4), (1,))), atol=3e-11)
