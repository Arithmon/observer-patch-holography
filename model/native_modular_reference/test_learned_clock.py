"""Independent native-space comparisons for record-inferred source clocks."""
from fractions import Fraction as F
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from . import learned_clock as l

WEIGHTS = (F(2, 225), F(1, 36), F(2, 45))


def test_ratios_do_not_use_known_reference_or_success_probability(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('Compilation must depend on the inferred weights only')
    monkeypatch.setattr(l.local, 'reference', forbidden)
    other = (F(1, 12), F(3, 20), F(1, 18))
    ratios = l.ratios(other)
    codes = l.full.support()[0]
    for code, state in zip(codes, l.local.native.states(binary=True), strict=True):
        assert ratios[code] == other[state[0]+state[1]]/other[0]
    assert all(ratios[i] == 1 for i in set(range(512))-set(codes))
    assert ratios == l.ratios(tuple(F(7, 13)*x for x in other))


def test_inference_log_bound_covers_every_box_corner_and_ignores_scale():
    box = tuple((w*F(99, 100), w*F(101, 100)) for w in WEIGHTS)
    bound = l.inference_bound(WEIGHTS, box, F(2))
    for mask in range(8):
        truth = tuple(box[q][(mask >> q) & 1] for q in range(3))
        delta = 2*np.log(np.array(list(map(float, WEIGHTS)))/np.array(list(map(float, truth))))
        delta -= (delta.max()+delta.min())/2
        assert np.max(np.abs(np.exp(1j*delta)-1)) <= float(bound)+1e-15
    scale = F(5, 33)
    assert bound == l.inference_bound(tuple(scale*x for x in WEIGHTS),
                                     tuple(tuple(scale*x for x in b) for b in box), F(2))
    assert l.inference_bound(WEIGHTS, box, F(0)) == 0


@pytest.mark.parametrize('bad', [(F(0), F(1), F(1)), (True, F(1), F(1)),
                                 [F(1)]*3, (F(1),)*2, (float('nan'), F(1), F(1))])
def test_inferred_weight_admission_fails_closed(bad):
    with pytest.raises(ValueError):
        l.ratios(bad)


@pytest.mark.parametrize('box', [((F(0), F(1)),)*3, ((F(2), F(1)),)*3,
                                ((0., 1.),)*3, ((F(1), F(2)),)*3,
                                [(F(1), F(2))]*3])
def test_nonfaithful_or_invalid_confidence_region_cannot_be_compiled(box):
    with pytest.raises(ValueError):
        l.inference_bound(WEIGHTS, box)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_source_inferred_action_on_entangled_reference_and_complete_feedback(law):
    rng = np.random.default_rng(4127)
    vector = rng.normal(size=132)+1j*rng.normal(size=132)
    vector /= np.linalg.norm(vector)
    density = np.outer(vector, vector.conj())
    # Perturb all three weights so the known-reference source word is wrong.
    learned = tuple(w*d for w, d in zip(WEIGHTS, (F(101, 100), F(99, 100), F(103, 100))))
    raw = l.execute(density, learned, law)
    ideal = l.execute(density, learned, law, raw=False)
    assert l.history.distance(raw, ideal) < float(l.certificate(learned)['source_history_halfdiamond_upper'])
    assert abs(sum(l.history.probabilities(raw).values())-1) < 1e-9
    assert set(raw) == set(ideal)
    assert all(a.shape == ((2, 2) if key[0] == 'abort' else (132, 132)) for key, a in raw.items())
    truth = l.cap.execute(density, law, raw=False)
    box = tuple((min(a, b), max(a, b)) for a, b in zip(learned, WEIGHTS))
    assert l.history.distance(ideal, truth) <= float(2*l.inference_bound(learned, box))
    assert l.history.distance(ideal, l.execute(density, learned, law, raw=False, feedback=False)) > 1e-4
    # The output states agree with a directly exponentiated native generator,
    # independently of the 512-address gauge and padding extension.
    native_weights = [float(learned[x[0]+x[1]]) for x in l.local.native.states(binary=True)]
    native_u = np.diag(np.exp(-1j*np.log(native_weights)))
    codes = l.full.support()[0]
    encoded = np.diag(np.exp(-1j*np.log(list(map(float, l.ratios(learned))))))
    gauged = encoded[np.ix_(codes, codes)]*complex(float(learned[0]))**(-1j)
    np.testing.assert_allclose(gauged, native_u, rtol=0, atol=1e-15)


def test_source_certificate_and_arrays_cannot_be_poisoned():
    cert = l.certificate(WEIGHTS)
    pulses = cert['source_pulses']
    cert['source_pulses'] = 0
    cert['source_history_halfdiamond_upper'] = F(0)
    assert l.certificate(WEIGHTS)['source_pulses'] == pulses > 0
    with pytest.raises(ValueError):
        l.clock(WEIGHTS)[0, 0] = 7


def test_original_full_array_certificate_rejects_wrong_padding():
    actual = l.clock(WEIGHTS).copy()
    pad = next(iter(set(range(512))-set(l.full.support()[0])))
    actual[pad, pad] = .5
    error = l.common.certify_diagonal(F(1), l.ratios(WEIGHTS), actual)
    assert error > F(1, 3)


def _retained_counts():
    from . import reference_tomography as readout
    return l.strict_json.load(readout.OUTPUT)['continuing_execution']['execution']['counts']


def test_fixed_admission_cap_is_not_the_observed_error_plus_failure():
    counts = _retained_counts()
    decision = l.admission(counts)
    assert decision['status'] == 'accept'
    assert decision['observed_confidence_event_bound'] < l.ACCEPTANCE_CAP
    assert decision['complete_protocol_halftrace_upper'] == (
        l.ACCEPTANCE_CAP+F(3, 512)*(l.source_bounds()[1]+1)/2)
    counts['success_1'] -= 10000
    counts['success_0'] += 10000
    changed = l.admission(counts)
    assert changed['inference']['status'] == 'faithful'
    assert changed['reason'] == 'confidence_box_too_wide'
    assert changed['complete_protocol_halftrace_upper'] == decision['complete_protocol_halftrace_upper']
    assert l.admission({k: 2*v for k, v in counts.items()})['reason'] == 'fixed_trial_budget_required'


@pytest.mark.parametrize('failure', ['unresolved', 'wide', 'compilation'])
def test_rejected_training_preserves_carrier_and_does_not_start_clock(monkeypatch, failure):
    counts = _retained_counts()
    if failure == 'unresolved':
        counts = {k: l.TRAINING_TRIALS if k == 'herald_failure' else 0 for k in counts}
    if failure == 'wide':
        counts['success_1'] -= 10000
        counts['success_0'] += 10000
    def failed_certificate(*args, **kwargs):
        if failure != 'compilation':
            raise AssertionError('rejected records must not attempt a clock')
        raise ArithmeticError('original source array missed its gate')
    def forbidden(*args, **kwargs):
        raise AssertionError('rejected records must not start future preparation')
    monkeypatch.setattr(l, 'certificate', failed_certificate)
    monkeypatch.setattr(l, 'prepared', forbidden)
    decision, outputs = l.prepared_from_records(counts, 'C', initial=37)
    assert decision['status'] == 'stop'
    assert len(outputs) == 1
    np.testing.assert_array_equal(next(iter(outputs.values())), np.diag([int(i == 37) for i in range(66)]))


@pytest.mark.parametrize('mutation', ['probability_sum', 'distance', 'boolean', 'count'])
def test_receipt_rejects_inconsistent_complete_history_diagnostics(mutation):
    actual = l.strict_json.load(Path(l.__file__).with_name('learned_clock_results.json'))
    changed = deepcopy(actual)
    if mutation == 'probability_sum':
        for row in changed['histories']:
            row['record_probabilities'] = {k: v+1.9e-10 for k, v in row['record_probabilities'].items()}
    elif mutation == 'distance':
        changed['histories'][0]['complete_output_halftrace_diagnostic'] = .5
    elif mutation == 'boolean':
        changed['training_counts']['warmup_timeout'] = False
    else:
        changed['training_counts']['success_1'] += 1
    with pytest.raises(ValueError):
        l.compare_receipt(changed, actual)
