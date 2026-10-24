"""Source-routed effects, complete failure laws and identification controls."""
from fractions import Fraction as F
from itertools import combinations
from unittest.mock import patch

import numpy as np
import pytest

from . import tomography as t


def test_exact_design_and_every_route():
    proof = t.certificate()
    assert proof['exact_design']['real_hermitian_rank'] == 4356
    assert proof['exact_routes_checked'] == 2145
    assert proof['settings'] == 4291
    assert F(proof['complete_record_half_l1_upper']) < F(1, 10**7)
    assert proof['native_attempts_per_trial'] == 2520
    assert proof['all_failures_retained']


@pytest.mark.parametrize('indices', [(True, 2), (0, False), (0., 1), (-1, 3), (0, 66), (2, 2), (3, 2)])
def test_pair_domain(indices):
    with pytest.raises(ValueError):
        t.word(*indices, 'mixer')


def test_incorrect_route_is_not_certified():
    with patch.object(t, 'route', return_value=()):
        with pytest.raises(ArithmeticError, match='isolate'):
            t.certificate()


def test_wrong_source_phase_is_not_certified():
    original = t.compiled.target
    with patch.object(t.compiled, 'target', side_effect=lambda kind: original('mixer')):
        with pytest.raises(ArithmeticError, match='quadrature'):
            t.design_certificate()


def _density_factor():
    rng = np.random.default_rng(10049)
    factor = rng.normal(size=(66, 3))+1j*rng.normal(size=(66, 3))
    factor /= np.linalg.norm(factor)
    return factor


def test_reconstruction_from_independent_effects_and_missing_phase_control():
    factor = _density_factor()
    rho = factor@factor.conj().T
    x, y = {}, {}
    for a, b in combinations(range(66), 2):
        # Explicit two-dimensional Pauli effects independent of decoder.
        block = rho[np.ix_((a, b), (a, b))]
        x[a, b] = float(np.trace(np.array([[0, 1], [1, 0]])@block).real)
        y[a, b] = float(np.trace(np.array([[0, -1j], [1j, 0]])@block).real)
    result = t.reconstruct(list(map(float, rho.diagonal().real)), x, y)
    np.testing.assert_allclose(result, rho, atol=1e-17, rtol=0)
    conjugate = rho.conj()
    assert np.linalg.norm(rho-conjugate) > .1
    np.testing.assert_array_equal(conjugate.real, rho.real)
    # Omitting the phase analyzer leaves these conjugate inputs indistinguishable.
    with pytest.raises(ValueError):
        t.reconstruct(list(map(float, rho.diagonal().real)), x, {})


@pytest.mark.parametrize('a,b', [(0, 1), (7, 42), (64, 65)])
@pytest.mark.parametrize('kind', ['mixer', 'y_analyzer'])
def test_whole_source_word_effect_and_reference_outputs(a, b, kind):
    rng = np.random.default_rng(a*66+b)
    factor = rng.normal(size=(66, 2, 2))+1j*rng.normal(size=(66, 2, 2))
    factor /= np.linalg.norm(factor)
    actual = t.records(factor, a, b, kind)
    target = t.records(factor, a, b, kind, raw=False)
    rho = np.einsum('irk,jrk->ij', factor, factor.conj())
    score = np.trace(target[str(a)]-target[str(b)]).real
    expected = 2*(rho[a, b].real if kind == 'mixer' else -rho[a, b].imag)
    assert abs(score-expected) < 1e-15
    reference = np.einsum('irk,isk->rs', factor, factor.conj())
    np.testing.assert_allclose(sum(target.values()), reference, atol=5e-16, rtol=0)
    distance = sum(np.linalg.svd(actual[k]-target[k], compute_uv=False).sum() for k in actual)/2
    assert distance < 1e-7
    assert np.trace(actual['timeout']).real > 0
    assert set(actual) == {str(i) for i in range(66)} | {'timeout', 'export_abort'}


def test_padding_abort_cannot_be_removed_by_renormalization():
    factor = _density_factor()[:, None, :]
    # This route mixes some native spectator pairs with padding. Their mass
    # is retained as a real terminal event even in the ideal model.
    found = False
    for a, b in [(0, 1), (7, 42), (64, 65)]:
        outcomes = t.records(factor, a, b, 'mixer', raw=False)
        aborted = np.trace(outcomes['export_abort']).real
        if aborted > .01:
            found = True
            assert sum(np.trace(outcomes[str(i)]).real for i in range(66)) < .99
    assert found


def test_inverse_routing_is_charged_and_not_assumed_an_adjoint():
    original = t.compiler.execute
    def corrupt(tape, qubits, state=None):
        out = original(tape, qubits, state)
        if qubits == 1 and tape[0]['kind'] == 'X' and tape[0]['angle'] < 0:
            out = out.copy(); out[0, 0] += .01
        return out
    with patch.object(t.compiler, 'execute', side_effect=corrupt):
        with pytest.raises(ArithmeticError, match='primitive'):
            t.certificate()


def test_sampling_budget_charges_bias_and_every_setting():
    result = t.sampling_bound(10**14, F(1, 10000), F(1, 10**8))
    assert result['simultaneous_failure_upper'] == F(4291, 10**6)
    assert result['copies_used'] == 4291*10**14
    ideal = t.sampling_bound(10**14, F(1, 10000), F())
    assert result['operator_error_upper'] > ideal['operator_error_upper']
    assert t.sampling_bound(1, F(1, 10000), F())['simultaneous_failure_upper'] == 1


@pytest.mark.parametrize('args', [(True,F(1),F()), (1,0.1,F()), (1,F(),F()), (1,F(1),F(-1))])
def test_bad_sampling_domain(args):
    with pytest.raises(ValueError):
        t.sampling_bound(*args)


def test_modular_bound_requires_certified_faithfulness_and_scales_with_parameter():
    assert t.modular_bound(F(1,1000),F(1,100),F(1)) == F(1,9)
    assert t.modular_bound(F(1,1000),F(1,100),F(-2)) == F(2,9)
    with pytest.raises(ValueError):
        t.modular_bound(F(1,100),F(1,100),F(1))


def test_linear_estimate_is_not_silently_projected_to_a_density():
    keys = tuple(combinations(range(66), 2))
    diagonal = [0.0]*66; diagonal[0] = 1.1
    result = t.reconstruct(diagonal, dict.fromkeys(keys, 0.), dict.fromkeys(keys, 0.))
    assert np.trace(result) == 1.1
