"""Challenge finite all-input work bounds and their exact execution ledger."""
from dataclasses import replace
from fractions import Fraction as F
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.linalg import expm

from ..common_source import protocol, scalar_coupling
from ..primitive_observer import runtime
from . import work_bounds as w


@pytest.fixture(scope='module')
def step():
    return w.certify_step(F(1, 4), 128)


@pytest.fixture(scope='module', params=(-1, 1))
def suffix(request):
    return w.certify_suffix(request.param, 128)


def _nuclear(array):
    return float(np.linalg.svd(array, compute_uv=False).sum())


def test_step_has_small_all_input_error_against_independent_continuous_matrix(step):
    assert tuple(record for record, _ in step.maps) == ((),)
    actual = step.maps[0][1]
    ideal = expm(-.25j*scalar_coupling.ideal_hamiltonian(2, .25))
    assert np.linalg.norm(actual-ideal, 2) <= float(step.stinespring_error_upper)+1e-13
    assert F(0) < step.half_diamond_upper < F(1, 10**6)
    assert step.half_diamond_upper == step.stinespring_error_upper*(2+step.stinespring_error_upper)/2
    assert step.channel_norm_upper == (1+step.stinespring_error_upper)**2
    assert step.includes_current_kick is False
    assert not actual.flags.writeable


def test_complete_suffix_matches_independent_all_input_reference_instrument(suffix):
    ideal = protocol.execute(suffix.program, backend='continuous')
    targets = {branch.record: branch.kraus for branch in ideal.branches}
    assert tuple(record for record, _ in suffix.maps) == ideal.records
    stacked = np.vstack([matrix-targets[record] for record, matrix in suffix.maps])
    assert np.linalg.norm(stacked, 2) <= float(suffix.stinespring_error_upper)+2e-13
    assert F(0) < suffix.half_diamond_upper < F(1, 50000)
    assert suffix.polynomial_remainder_upper < F(1, 10**16)
    rng = np.random.default_rng(160391)
    state = rng.normal(size=(16, 3))+1j*rng.normal(size=(16, 3))
    state /= np.linalg.norm(state)
    distance = 0.
    for record, matrix in suffix.maps:
        actual, target = (matrix@state).ravel(), (targets[record]@state).ravel()
        distance += _nuclear(np.outer(actual, actual.conj())-np.outer(target, target.conj()))/2
    assert 0 < distance <= float(suffix.half_diamond_upper)+2e-12
    assert suffix.includes_current_kick is True


def test_common_runtime_uses_the_same_complete_operators_on_correlated_input(suffix):
    # The reference carries native/logical information independently of the
    # finite work register. This is an execution check, not a rounded-state
    # certificate; later applications have their explicit exact residual API.
    state = np.zeros((16, 2), complex)
    state[0, 0], state[15, 1], state[6, 0] = .5, .5j, np.sqrt(.5)
    initial = np.outer(state.ravel(), state.ravel().conj())
    actual = runtime.execute(suffix.compiled_program, initial=initial)
    for branch, (record, matrix) in zip(actual.branches, suffix.maps, strict=True):
        assert branch.record == record
        expected = (matrix@state).ravel()
        np.testing.assert_allclose(branch.density, np.outer(expected, expected.conj()), atol=2e-13, rtol=0)


def test_wrong_complete_postmeasurement_action_does_not_pass_certificate(monkeypatch):
    original = w.common.compile_program
    def wrong(program):
        compiled = original(program)
        event = compiled.events[-1]
        instrument = event.instruments[0]
        phase = np.diag([1 if i//4 % 2 == 0 else -1 for i in range(16)])
        modified = runtime.Instrument(kraus=tuple((phase@group[0],) for group in instrument.kraus),
                                      visible=instrument.visible, resources=instrument.resources,
                                      binding=instrument.binding)
        return runtime.Program(compiled.name, 16, compiled.initial_density,
                               (*compiled.events[:-1], replace(event, instruments=(modified,))),
                               compiled.preparation_resources, compiled.binding)
    monkeypatch.setattr(w.common, 'compile_program', wrong)
    result = w.certify_program(protocol.Program(operations=(protocol.Operation('gauge_read'),)))
    assert result.half_diamond_upper > 1


def test_integer_lanes_preserve_uint64_and_do_not_share_cached_mutable_arrays():
    array = np.zeros((16, 16), dtype=np.uint64)
    array[0, 0] = 2**63+1
    real, imag, denominator = w.exact_matrix(array)
    assert real[0, 0] == 2**63+1 and denominator == 1
    assert not np.any(imag)
    real[0, 0] = 0
    imag[1, 1] = 12
    again = w.exact_matrix(array)
    assert again[0][0, 0] == 2**63+1 and not np.any(again[1])


def test_exact_work_boundaries_reject_subclass_that_forges_snapshot_bytes():
    class HiddenBytes(np.ndarray):
        def tobytes(self, *args, **kwargs):
            return bytes(self.nbytes)

    actual = np.eye(16)
    hidden = actual.view(HiddenBytes)
    assert w.nuclear_upper(actual) == 16
    calls = (
        lambda: w.exact_matrix(hidden),
        lambda: w.nuclear_upper(hidden),
        lambda: w.certify_preparation(hidden),
        lambda: w.exact_product_residual(hidden, actual, actual),
        lambda: w.exact_product_residual(actual, hidden, actual),
        lambda: w.exact_product_residual(actual, actual, hidden),
        lambda: w.exact_density_residual(hidden, (actual,), actual),
        lambda: w.exact_density_residual(actual, (hidden,), actual),
        lambda: w.exact_density_residual(actual, (actual,), hidden),
    )
    for call in calls:
        with pytest.raises(ValueError):
            call()


def test_exact_product_residual_includes_imaginary_and_underflow_errors():
    left = np.eye(16, dtype=complex)*(1+1j)/4
    right = np.eye(16, dtype=complex)*(1-1j)/8
    product = np.eye(16, dtype=complex)/16
    assert w.exact_product_residual(product, left, right) == 0
    product[2, 7] += 1j/32
    assert w.exact_product_residual(product, left, right) == F(1, 32)
    tiny = np.eye(16)*1e-200
    error = w.exact_product_residual(np.zeros((16, 16)), tiny, tiny)
    assert error > 0
    assert error >= 4*F(1e-200)**2


def test_exact_density_residual_checks_original_nonhermitian_arrays_and_factor_order():
    first = np.eye(16, dtype=complex)
    first[0, 1] = 1j/4
    last = np.eye(16, dtype=complex)
    last[1, 0] = 1/2
    initial = np.zeros((16, 16), complex)
    initial[0, 0], initial[1, 1] = .75, .25
    expected = last@first@initial@first.conj().T@last.conj().T
    assert w.exact_density_residual(expected, (first, last), initial) == 0
    assert w.exact_density_residual(expected, (last, first), initial) > 0
    bad = expected.copy()
    bad[5, 2] += 1j/64
    assert w.exact_density_residual(bad, (first, last), initial) == F(1, 32)
    assert w.exact_density_residual(bad, (first, last), initial) >= F(1, 128)


def test_finite_preparation_bound_includes_raw_nonpositive_mass_and_phase_errors(step):
    assert w.certify_preparation(step.compiled_program.initial_density) < F(1, 10**12)
    exact = np.zeros((16, 16), complex)
    exact[0, 0] = 1
    assert w.certify_preparation(exact) == 0
    bad = exact.copy()
    bad[5, 5], bad[3, 7] = -1/16, 1j/8
    bound = w.certify_preparation(bad)
    assert bound >= F(3, 32)
    assert float(bound) >= _nuclear(bad-exact)/2


@pytest.mark.parametrize('array', [np.zeros((15, 15)), np.zeros(16), [[1]],
    np.ones((16, 16), dtype=bool), np.full((16, 16), float('nan')),
    np.ma.array(np.eye(16), mask=False)])
def test_exact_array_boundary_rejects_ambiguous_or_nonfinite_inputs(array):
    with pytest.raises(ValueError):
        w.exact_matrix(array)


@pytest.mark.parametrize('time', [True, .25, 1, F(1, 3)])
def test_step_requires_declared_exact_time(time):
    with pytest.raises(ValueError):
        w.certify_step(time)


@pytest.mark.parametrize('current', [True, -1., 0, 2])
def test_suffix_requires_signed_current(current):
    with pytest.raises(ValueError):
        w.certify_suffix(current)


def test_finite_certificate_rejects_other_cutoffs_sharp_and_hidden_inputs():
    with pytest.raises(ValueError):
        w.certify_program(protocol.Program(dimension=4))
    with pytest.raises(ValueError):
        w.certify_program(protocol.Program(coupling=F(0)))
    with pytest.raises(ValueError):
        w.certify_program(protocol.Program(operations=(protocol.Operation('sharp_read'),)))
    with pytest.raises(ValueError):
        w.certify_program(protocol.Program(native_current=1), include_kick=False)
    with pytest.raises(ValueError):
        w.exact_density_residual(np.eye(16), [np.eye(16)], np.eye(16))


def test_exact_finite_generator_matches_independent_tensor_hamiltonian_on_every_column():
    # 128 H at occupation dimension two has integer entries. In particular,
    # projecting the infinite quadratic term is not squaring a truncated Q.
    identity = np.eye(2, dtype=int)
    flip = np.array([[0, 1], [1, 0]])
    number = np.diag([0, 1])
    gauge_flip = np.kron(identity, flip)+np.kron(flip, identity)
    electric = np.diag([0, 1, 1, 0])
    occupation = np.kron(identity, number)+np.kron(number, identity)
    expected = (128*np.kron(np.eye(4, dtype=int), gauge_flip+electric)
                +146*np.kron(occupation+np.eye(4, dtype=int), np.eye(4, dtype=int))
                -18*np.kron(np.kron(flip, flip), np.eye(4, dtype=int))
                +32*np.kron(np.kron(identity, flip), np.diag([1, -1, 1, -1])))
    actual = np.zeros((16, 16), dtype=int)
    for column in range(16):
        terms = w.sh._compress(w.ih._hamiltonian_numerator(w._basis(column).coefficients), 2)
        for (n0, n1, gauge), (real, imag) in terms.items():
            assert imag == 0
            actual[4*(n0+2*n1)+gauge, column] = real
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(actual, actual.T)
    np.testing.assert_allclose(actual/128, scalar_coupling.ideal_hamiltonian(2, .25), atol=3e-16, rtol=0)


def test_work_accountant_retains_gain_in_repeated_steps_suffix_and_preparation(monkeypatch):
    from . import experiment as e

    # Stress the algebra with exact dyadic CP maps that amplify trace. These
    # are deliberately outside the execution completeness gate. The oracle
    # is their analytic output mass, not a restatement of the recurrence.
    step_gain, suffix_gain, mass_error = F(5, 4), F(9, 8), F(1, 16)
    source = float(step_gain)*np.eye(16)
    zero = np.zeros((16, 16))
    zero[0, 0] = 1
    original = (1+2*float(mass_error))*zero
    projectors = tuple(np.diag([int((i % 2) == bit) for i in range(16)]) for bit in (0, 1))
    maps = tuple(((bit,), float(suffix_gain)*projector) for bit, projector in enumerate(projectors))

    def history(n, with_suffix):
        density = float(step_gain**(2*n))*original
        branches = (tuple(SimpleNamespace(record=record, density=k@density@k.T) for record, k in maps)
                    if with_suffix else (SimpleNamespace(record=(), density=density),))
        return SimpleNamespace(branches=branches, initial_density=original,
                               reference_dimension=1, preparation_mode='declared')

    prefixes = tuple(history(n, False) for n in range(4))
    prefixes[0].program = SimpleNamespace(events=(SimpleNamespace(
        instruments=(SimpleNamespace(kraus=((source,),)),)),))
    work = SimpleNamespace(horizon=3, delta=F(1, 4), steps=128, prefixes=prefixes,
        successes=tuple((n, j, history(n, True)) for n in range(1, 4) for j in (-1, 1)))
    monkeypatch.setattr(e.c, 'validate_work', lambda value: value)
    monkeypatch.setattr(e, 'work_binding', lambda value: 'synthetic-gain-control')
    monkeypatch.setattr(e.wb, 'certify_step', lambda *args: SimpleNamespace(
        maps=(((), source),), stinespring_error_upper=step_gain-1))
    monkeypatch.setattr(e.wb, 'certify_suffix', lambda *args: SimpleNamespace(
        maps=maps, stinespring_error_upper=suffix_gain-1, compiled_program=None))
    monkeypatch.setattr(e.wb, 'certify_preparation', lambda initial: mass_error)
    monkeypatch.setattr(e.runtime, 'execute', lambda program: history(0, True))
    monkeypatch.setattr(e.report, '_gate', lambda *args: None)
    bounds = e.work_certificates(work)

    for n in range(1, 4):
        exact_prefix_error = (step_gain**(2*n)*(1+2*mass_error)-1)/2
        exact_suffix_error = (suffix_gain**2*step_gain**(2*n)*(1+2*mass_error)-1)/2
        assert bounds['errors'][n, 0] >= exact_prefix_error
        assert bounds['errors'][n, -1] >= exact_suffix_error
        assert bounds['errors'][n, 1] >= exact_suffix_error
        assert bounds['errors'][n, 1]-exact_suffix_error < F(1, 10**14)
        assert bounds['arithmetic'][n, 1] == 0
    assert bounds['prefix_operator_errors'][3] > 3*(step_gain-1)


def test_joint_accountant_requires_native_work_cross_term_and_timeout(monkeypatch):
    from . import experiment as e

    # Exact rank-one positive factors saturate the claimed tensor bound.
    # Their excess masses must not be removed by normalization or clipping.
    records = {(1, -1): F(1, 4), (1, 1): F(1, 4), (1, 0): F(1, 2)}
    native_error, work_error = F(1, 16), F(1, 8)
    checks = {'distance_upper': len(records)*native_error,
              'per_record': {key: {'ideal_mass': mass, 'distance_upper': native_error}
                             for key, mass in records.items()}}
    monkeypatch.setattr(e.clock, 'compare_outputs', lambda *args, **kwargs: checks)
    monkeypatch.setattr(e.c, 'validate_work', lambda value: value)
    monkeypatch.setattr(e, 'work_binding', lambda value: 'synthetic-tensor-control')
    monkeypatch.setattr(e.clock, 'statistics', lambda *args, **kwargs:
                        {'timeout_probability': records[1, 0]})
    monkeypatch.setattr(e.report, '_gate', lambda *args: None)
    bounds = {'errors': {key: work_error for key in records}, 'work_binding': 'synthetic-tensor-control'}
    result = e.joint_certificate(SimpleNamespace(horizon=1), bounds, 'C', {}, cap=1)
    actual_distance = sum(((mass+2*native_error)*(1+2*work_error)-mass)/2
                          for mass in records.values())
    assert result['finite_error'] == actual_distance
    assert result['finite_error'] > checks['distance_upper']+work_error
    assert result['eventual_error'] == actual_distance+records[1, 0]


def test_probability_interval_charges_full_nuclear_error_for_unnormalized_arrays():
    from . import experiment as e

    # A rank-one excess 2e has half nuclear error e but changes this event
    # probability by 2e. The normalized-state Helstrom factor is unavailable.
    error, exact_probability = F(1, 16), F(1, 2)
    computed_probability = exact_probability+2*error
    lower, upper = map(F, e.interval(computed_probability, error))
    assert lower <= exact_probability <= upper
    assert computed_probability-error > exact_probability
