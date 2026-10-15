"""Exact nonidentifiability, useful probes and complete native execution."""
from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

from tools import strict_json
from ..observer.scalar_response.experiment import compare_receipt
from ..observer.source_scalar import compiler
from . import response, source
from .experiment import run


def test_four_level_null_directions_are_exact_and_not_a_precision_issue():
    report = response.complete_design()
    assert report['before']['rank'] == 3
    assert report['before']['nullspace_exact'] == [
        ['1', '0', '0', '0', '0', '0'],
        ['0', '0', '0', '0', '1', '0'],
        ['0', '0', '0', '0', '0', '1']]
    assert response.exact_blind_compression() == {
        'binom(N,4)': ['0']*4, 'binom(N,5)': ['0']*4}


def test_whole_calibration_evolution_is_unchanged_with_noncommuting_background():
    rng = np.random.default_rng(1107)
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    h = raw+raw.conj().T
    correction = np.diag([math.comb(n, 4)+2*math.comb(n, 5) for n in range(4)])
    assert np.array_equal(h, h+correction)
    for time in (0., .1, 2.3):
        assert np.array_equal(expm(-1j*time*h), expm(-1j*time*(h+correction)))


def test_next_occupations_close_only_the_declared_family():
    report = response.complete_design()
    assert report['added_occupations'] == [4, 5]
    assert report['after']['rank'] == 5
    assert report['after']['nullspace_exact'] == [['1', '0', '0', '0', '0', '0']]
    wider = response.complete_design(candidates=(4, 5), order=6)
    assert wider['remaining_non_gauge_directions'] == 1
    assert wider['after']['nullspace_exact'][-1] == ['0']*6+['1']


def test_exact_inverse_retains_conditioning_cost():
    inverse = sp.Matrix(response.finite_difference_inverse())
    matrix = response.response_matrix(range(1, 6))[:, 1:]
    assert inverse*matrix == sp.eye(5)
    assert [sum(abs(value) for value in inverse.row(i)) for i in range(5)] == [1, 3, 7, 15, 31]


@pytest.mark.parametrize('extra,time', [(0., source.TIME), (1., source.TIME), (-.7, .2), (2.3, .9)])
def test_native_full_evolution_and_complete_read_law_agree_with_reference(extra, time):
    result = source.execute(extra, time)
    np.testing.assert_allclose(result['complete_read_probabilities'], source.independent_reference(extra, time), atol=2e-13, rtol=0)
    assert result['full_phase_evolution_target_residual'] < 2e-13
    assert result['complete_probability_mass'] == pytest.approx(1, abs=2e-13)
    assert len(result['complete_read_probabilities']) == 8
    assert result['stage_pulses'] == {'prepare': 1, 'evolve': 38, 'read': 2}
    word = source.words(extra, time)['evolve']
    np.testing.assert_allclose(compiler.execute(word, 3), expm(-1j*time*np.diag(source.energies(extra))), atol=2e-13, rtol=0)


def test_phase_sensitive_read_is_necessary_and_preserves_every_record():
    coherent = [source.execute(extra) for extra in (0., 1.)]
    assert coherent[1]['outcome_zero_probability']-coherent[0]['outcome_zero_probability'] == pytest.approx(1/math.sqrt(2), abs=1e-13)
    for extra in (0., 1.):
        result = source.execute(extra, commit_before_read=True)
        joint = np.array(result['joint_old_occupation_then_read_code'])
        assert joint.shape == (8, 8)
        assert joint.sum() == pytest.approx(1, abs=1e-13)
        assert result['outcome_zero_probability'] == pytest.approx(.5, abs=1e-13)
        np.testing.assert_allclose(joint.sum(axis=0), result['complete_read_probabilities'], atol=1e-13, rtol=0)


def test_ramsey_slope_identity_for_arbitrary_noncommuting_hermitian_background():
    rng = np.random.default_rng(541)
    raw = rng.normal(size=(8, 8))+1j*rng.normal(size=(8, 8))
    h = raw+raw.conj().T
    for n in range(1, 8):
        psi = np.zeros(8, complex)
        read = np.zeros(8, complex)
        psi[[0, n]] = 1/math.sqrt(2)
        read[0], read[n] = 1/math.sqrt(2), -1j/math.sqrt(2)
        rho = np.outer(psi, psi.conj())
        effect = np.outer(read, read.conj())
        slope = np.trace(effect@(-1j*(h@rho-rho@h))).real
        assert 2*slope == pytest.approx((h[n, n]-h[0, 0]).real, abs=2e-13)


def test_fixed_time_witness_has_ideal_sampling_guarantee_and_explicit_bias_limit():
    result = response.bernoulli_design(source.TIME, 1.)
    assert result['worst_case_one_hypothesis_error_upper'] == pytest.approx(math.exp(-19/4))
    assert result['worst_case_one_hypothesis_error_upper'] < .01
    robust = response.bernoulli_design(source.TIME, 1., trials=20, probability_error=.01)
    assert robust['worst_case_one_hypothesis_error_upper'] < .01
    assert response.bernoulli_design(source.TIME, 1., probability_error=.5)['worst_case_one_hypothesis_error_upper'] == 1.
    assert response.derivative_noise_radius(1/128, 200, .01) > 20
    # A repeated phase can alias: rank at time zero is no global finite-time
    # uniqueness theorem. The design keeps this negative control explicit.
    assert response.bernoulli_design(2*math.pi, 1.)['ideal_binary_tv'] < 1e-14


@pytest.mark.parametrize('bad', [True, float('nan'), float('inf'), '1', None])
def test_invalid_candidate_parameters_fail_closed(bad):
    with pytest.raises(ValueError):
        source.execute(bad)


@pytest.mark.parametrize('call', [
    lambda: response.response_matrix([1, 1]),
    lambda: response.response_matrix([0]),
    lambda: response.response_matrix([1.2]),
    lambda: response.bernoulli_design(-1, 1),
    lambda: response.bernoulli_design(1, 1, probability_error=-.1),
    lambda: response.bernoulli_design(1, 1, trials=True),
    lambda: response.derivative_noise_radius(0, 200, .01),
    lambda: response.derivative_noise_radius(.1, 200, 1),
    lambda: source.execute(1, commit_before_read=1),
])
def test_invalid_design_inputs_fail_closed(call):
    with pytest.raises(ValueError):
        call()


def test_retained_receipt_replays_with_exact_provenance_and_design():
    expected = strict_json.load(Path(__file__).with_name('results.json'))
    actual = run()
    compare_receipt(expected, actual)
    changed = deepcopy(actual)
    changed['response_design']['before']['rank'] += 1
    with pytest.raises(ValueError):
        compare_receipt(expected, changed)
    changed = deepcopy(actual)
    changed['input_sha256'][next(iter(changed['input_sha256']))] = '0'*64
    with pytest.raises(ValueError):
        compare_receipt(expected, changed)


def test_producer_rejects_wrong_global_phase_even_when_output_laws_agree(monkeypatch):
    original = compiler.compile_diagonal
    def wrong_phase(phases):
        return original(phases)+[compiler.pulse('I', .001, [0])]
    monkeypatch.setattr(compiler, 'compile_diagonal', wrong_phase)
    # The unconditional measurement cannot see a global phase, but the
    # declared complete source-word target check must still reject it.
    row = source.execute(0.)
    np.testing.assert_allclose(row['complete_read_probabilities'],
                               source.independent_reference(0.), atol=2e-12, rtol=0)
    assert row['full_phase_evolution_target_residual'] > .0009
    with pytest.raises(ArithmeticError, match='Intervention design check failed'):
        run()


@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -1.])
def test_producer_rejects_invalid_full_phase_diagnostics(monkeypatch, invalid):
    original = source.execute
    def malformed(*args, **kwargs):
        row = original(*args, **kwargs)
        row['full_phase_evolution_target_residual'] = invalid
        return row
    monkeypatch.setattr(source, 'execute', malformed)
    with pytest.raises(ArithmeticError, match='Intervention design check failed'):
        run()


@pytest.mark.parametrize('candidate', [0., 1.])
@pytest.mark.parametrize('outcome', range(source.DIMENSION))
@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), -float('inf')])
def test_producer_rejects_nonfinite_independent_reference_laws(monkeypatch, candidate, outcome, invalid):
    original = source.independent_reference

    def malformed(extra, *args, **kwargs):
        law = original(extra, *args, **kwargs)
        if extra == candidate:
            law[outcome] = invalid
        return law

    monkeypatch.setattr(source, 'independent_reference', malformed)
    # A later JSON encoder rejection is not a substitute for a failed
    # scientific acceptance gate. In particular, max([finite, NaN]) can
    # silently retain its finite first element.
    with pytest.raises(ArithmeticError, match='Intervention design check failed'):
        run()
