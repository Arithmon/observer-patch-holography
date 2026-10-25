"""Exact native quadrature effects, source words and continuing poststates."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import phase_readout as p

REFERENCE = (F(2, 225), F(1, 36), F(2, 45))
OTHER_REFERENCE = (F(1, 100), F(1, 50), F(3, 20))
WEIGHTS = (F(11), F(37), F(59))


def _density(weights=REFERENCE):
    return np.diag([float(weights[sum(x[i] for i in p.local.CORE)])
                    for x in p.local.native.states(binary=True)]).astype(complex)


def _moments(output):
    return tuple(sum(p.score(record)**power*np.trace(rho).real for record, rho in output.items())
                 for power in (1, 2))


def test_exact_native_effects_and_gaussian_integer_quadrature_sign():
    proof = p.exact_effects()
    assert proof['exact_gaussian_integer_pauli_conjugations_checked']
    assert len(proof['rows']) == 4
    assert all(row['matrix_units'] == 4356 and row['selected_spectator_pairs'] == 9
               and row['diagonal_input_closes_after_native_read'] for row in proof['rows'])
    assert proof['native_current_and_coin_only']
    assert not proof['hidden_spectator_selection_or_projective_occupation_read']


@pytest.mark.parametrize('index,law', [(0, 'C'), (0, 'E'), (1, 'C'), (1, 'E')])
def test_actual_two_native_records_match_both_quadrature_formulas(index, law):
    for reference in (REFERENCE, OTHER_REFERENCE):
        for parameter in (F(-1), F(-3, 4), F(), F(3, 4), F(1)):
            for quadrature in p.QUADRATURES:
                out = p.execute(_density(reference), WEIGHTS, law, index=index,
                                quadrature=quadrature, s=parameter, raw=False)
                mean, second = _moments(out)
                lo, hi = p.response_enclosure(WEIGHTS, reference, parameter, quadrature=quadrature)
                assert float(lo)-2e-14 <= mean <= float(hi)+2e-14
                assert second == pytest.approx(float(9*(reference[0]+reference[1])), abs=2e-14)
                assert set(out) == set(p.record_alphabet(law, index))
                assert sum(np.trace(a).real for a in out.values()) == pytest.approx(1, abs=2e-14)


def test_quadrature_orientation_and_the_real_analyzer_countercontrol():
    observed = {}
    for parameter in (F(-3, 4), F(3, 4)):
        for quadrature in p.QUADRATURES:
            observed[parameter, quadrature] = _moments(p.execute(_density(), WEIGHTS, 'C',
                quadrature=quadrature, s=parameter, raw=False))[0]
    assert observed[F(-3, 4), 'cos'] == pytest.approx(observed[F(3, 4), 'cos'])
    assert observed[F(-3, 4), 'sin'] == pytest.approx(-observed[F(3, 4), 'sin'])
    assert observed[F(-3, 4), 'sin'] < -.01 < 0 < observed[F(3, 4), 'sin']


def test_core_tag_uses_disturbed_first_output_not_an_independent_reset():
    rho = _density()
    first, _ = p.cap._read_stage(rho, 0, 'C',
        p.stage(WEIGHTS, F(3, 4), quadrature='sin', raw=False))
    wrong = 0.
    for record, state in first.items():
        for branch in p.interface.fixed_maps('C', p.local.CORE):
            key = ('read', *record, branch.coin, branch.label, branch.current)
            wrong += p.score(key)*np.trace(state).real*np.trace(branch.matrix@rho@branch.matrix.conj().T).real
    right = _moments(p.execute(rho, WEIGHTS, 'C', quadrature='sin', s=F(3, 4), raw=False))[0]
    assert abs(right-wrong) > 1e-3


@pytest.mark.parametrize('index,law', [(0, 'C'), (1, 'E')])
def test_ideal_population_kernel_is_the_complete_diagonal_quantum_instrument(index, law):
    initial = np.arange(1, 67, dtype=float)
    initial /= initial.sum()
    for quadrature in p.QUADRATURES:
        records, kernel = p.population_instrument(WEIGHTS, F(-3, 4), index=index, quadrature=quadrature, law=law)
        out = p.execute(np.diag(initial).astype(complex), WEIGHTS, law,
                        index=index, quadrature=quadrature, s=F(-3, 4), raw=False)
        assert kernel.shape == (66*len(records), 66)
        assert np.all(kernel >= 0)
        np.testing.assert_allclose(kernel.sum(axis=0), 1, atol=2e-14)
        for i, record in enumerate(records):
            if record[0] == 'abort':
                assert not np.any(kernel[66*i:66*(i+1)])
                assert not np.any(out[record])
            else:
                np.testing.assert_allclose(out[record], np.diag(kernel[66*i:66*(i+1)]@initial), atol=2e-14)


def test_complete_instrument_preserves_untouched_entangled_reference():
    rng = np.random.default_rng(193)
    vector = rng.normal(size=132)+1j*rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    original = np.einsum('aras->rs', density.reshape(66, 2, 66, 2))
    for law in ('C', 'E'):
        out = p.execute(density, WEIGHTS, law, index=1, quadrature='sin', s=F(-3, 4), raw=False)
        marginal = sum(a if record[0] == 'abort' else np.einsum('aras->rs', a.reshape(66, 2, 66, 2))
                       for record, a in out.items())
        np.testing.assert_allclose(marginal, original, atol=2e-14)
        assert out['abort', 1].shape == (2, 2)
        assert all(np.linalg.eigvalsh(a).min() >= -1e-14 for a in out.values())


def test_preparation_keeps_unread_core_and_all_failures(monkeypatch):
    original = p.tomography.population_counts
    calls = []
    def recorded():
        calls.append(True)
        return original()
    monkeypatch.setattr(p.tomography, 'population_counts', recorded)
    out = p.prepared(WEIGHTS, 'C', initial=17, quadrature='sin', s=F(-3, 4), raw=False)
    assert calls
    assert set(out) == set(p.record_alphabet('C', prepared=True))
    assert np.trace(out['preparation', 'herald_failure']).real > .8
    assert sum(np.trace(a).real for a in out.values()) == pytest.approx(1, abs=5e-10)


def test_statistical_budget_charges_source_gain_and_finite_preparation():
    b = p.complete_trial_budget()
    preparation, q, maximum = p.tomography._record_contract()
    error, gain = p.source_bounds()
    assert b['complete_trial_halftrace_upper'] == gain*preparation+q*error
    assert b['score_mean_bias_upper'] == 4*b['complete_trial_halftrace_upper']
    assert b['score_second_moment_upper'] == 9*q/20+8*b['complete_trial_halftrace_upper']
    assert b['native_attempts_per_trial_upper'] == maximum+3
    assert b['raw_source_certification_required']
    assert gain > 1 and b['score_mean_bias_upper'] < F(1, 10**6)


@pytest.mark.parametrize('record', [
    ('read', True, 0, 0, 0, 0, 0), ('read', 0, 0, 1), ('abort', False),
    ('preparation', 'missing'), ('read', 0, 0, 3, 0, 0, 0), [],
])
def test_bad_records_are_not_silently_zero_scores(record):
    with pytest.raises(ValueError):
        p.score(record)


def test_zero_scores_include_all_failure_and_unmatched_native_branches():
    assert p.score(('numerical_abort',)) == 0
    assert p.score(('abort', 1)) == 0
    for flag in p.PREPARATION_FLAGS:
        assert p.score(('preparation', flag)) == 0
    assert p.score(('read', 0, 1, 1, 0, 0, 0)) == 2
    assert p.score(('read', 1, -1, -1, 1, -1, -1)) == -2
    assert p.score(('read', 0, 1, 1, 1, 0, 0)) == 0
    assert p.score(('read', 0, 1, 1, 0, 1, 1)) == 0


@pytest.mark.parametrize('index', [0, 1])
def test_raw_quadrature_word_certifies_every_original_source_address(index):
    raw = p.quadrature_rotation(index)
    assert p.certify_quadrature(index, raw) <= p.ROTATION_GATE
    assert len(p.quadrature_word(index)) == len(p.cap.rotation_word(index, 1))+20
    assert not raw.flags.writeable
    assert np.linalg.norm(raw-p.ideal_quadrature(index), 'fro') < float(p.ROTATION_GATE)
    modified = raw.copy()
    native_codes = set(p.full.support()[0])
    padding = next(i for i in range(512) if i not in native_codes)
    modified[padding, padding] += 1e-4
    assert p.certify_quadrature(index, modified) > p.ROTATION_GATE
    assert p.certify_quadrature(index, raw.conj()) > p.ROTATION_GATE
    with pytest.raises(ValueError):
        p.certify_quadrature(index, raw.view(type('Hostile', (np.ndarray,), {})))


def test_raw_source_instrument_and_complete_word_budget():
    raw = p.execute(_density(), WEIGHTS, 'C', quadrature='sin', s=F(-3, 4))
    ideal = p.execute(_density(), WEIGHTS, 'C', quadrature='sin', s=F(-3, 4), raw=False)
    budget = p.certificate(WEIGHTS, F(-3, 4), quadrature='sin')
    assert p.history.distance(raw, ideal) <= float(budget['source_history_halfdiamond_upper'])
    assert budget['source_pulses'] == (len(p.cap.rotation_word(0, -1))
        +len(p.common.diagonal_word(F(-3, 4), p.learned.ratios(WEIGHTS)))+len(p.quadrature_word(0)))
    assert budget['proper_code_service_events'] == 2*budget['source_pulses']
    assert budget['native_current_reads'] == 2 and budget['native_exports'] == 1
    assert not budget['intermediate_pulses_proved_spatially_local']



@pytest.mark.parametrize('mutation', ['phase', 'output', 'amplitude'])
def test_exact_effect_theorem_rejects_changed_native_arrays(monkeypatch, mutation):
    from dataclasses import replace
    original = p.interface.fixed_maps
    def corrupted(law, seam=None):
        branches = original(law, seam)
        first, *rest = branches
        matrix = first.matrix.copy()
        before, after = first.mapping[0]
        if mutation == 'phase':
            matrix[after, before] *= 1j
        elif mutation == 'output':
            matrix[after, before] = 0
            matrix[(after+1) % 66, before] = 1/np.sqrt(2)
        else:
            matrix[after, before] /= 2
        return (replace(first, matrix=matrix), *rest)
    monkeypatch.setattr(p.interface, 'fixed_maps', corrupted)
    with pytest.raises(ArithmeticError, match='partial injection'):
        p.exact_effects()


def test_exact_effect_theorem_rejects_a_missing_native_record(monkeypatch):
    original = p.interface.fixed_maps
    monkeypatch.setattr(p.interface, 'fixed_maps', lambda law, seam=None: original(law, seam)[:-1])
    with pytest.raises(ArithmeticError, match='complete branch'):
        p.exact_effects()
