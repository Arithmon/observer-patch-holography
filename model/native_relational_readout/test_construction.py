"""Actual current-controlled stopping, complete tails and coherent controls."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as r


@pytest.fixture(scope='module')
def sources():
    return {law: r.source(law) for law in ('C', 'E')}


def test_native_moving_maps_coincide_but_zero_current_coherence_differs(sources):
    assert sources['C'].moving_maps == sources['E'].moving_maps
    assert sources['C'].no_current_counts != sources['E'].no_current_counts
    for marked in sources.values():
        assert set(marked.no_current_counts[i][i] for i in range(66)) == {50, 52}
        assert len(marked.moving_maps) == 60
        for mapping in marked.moving_maps:
            assert len(set(dict(mapping).values())) == len(mapping)


@pytest.mark.parametrize('law,coherence,plus', [
    ('C', F(401, 1365), F(883, 1365)),
    ('E', F(1667, 8855), F(5261, 8855))])
def test_exact_infinite_event_readout_and_independent_phase_degree_formula(sources, law, coherence, plus):
    marked = sources[law]
    state = r.exact_stopped_state(marked)
    assert r.trace(state) == 1
    assert r.phase_probabilities(state) == (plus, 1-plus)
    assert r.phase_probabilities(r.exact_stopped_state(marked, sign=-1)) == (1-plus, plus)
    denominator = 10 if law == 'C' else 18
    independently = (2*F(3, denominator+3)+4*F(4, denominator+4)+4*F(5, denominator+5))/10
    assert coherence == independently == 2*plus-1
    audit = r.exact_audit(marked)
    assert audit['phase_sector_moving_degree_histogram'] == {'3': 2, '4': 4, '5': 4}
    assert audit['positive_global_idle_rescaling_preserves_complete_stopped_quantum_instrument']


@pytest.mark.parametrize('law', ['C', 'E'])
def test_resolvent_matches_independent_finite_native_CP_composition(sources, law):
    marked = sources[law]
    initial, _, _, _ = r.quantum.phase_apparatus()
    survival = initial
    stopped = np.zeros_like(initial)
    for n in range(1, 9):
        current = sum(k@survival@k.conj().T for k in marked.instrument.kraus[1])
        stopped += current
        survival = sum(k@survival@k.conj().T for k in marked.instrument.kraus[0])
        np.testing.assert_allclose(current, r.dense(r.exact_stopped_state(marked, first_at=n)), atol=2e-14, rtol=0)
        np.testing.assert_allclose(stopped, r.dense(r.exact_stopped_state(marked, through=n)), atol=2e-14, rtol=0)
        np.testing.assert_allclose(survival, r.dense(r.exact_tail(marked, n)), atol=2e-14, rtol=0)
        assert abs(np.trace(stopped+survival)-1) < 3e-14
        assert r.trace(r.exact_tail(marked, n)) <= F(13, 15)**n


def test_same_waiting_clock_does_not_imply_same_coherent_measurement(sources):
    attached = r.clock_attachment(sources)
    assert attached['prepared_class_law'] == ['2/5', '1/2', '1/10']
    assert attached['mean_first_current_attempts'] == '33/5'
    c, e = attached['joint_law']
    assert c['no_current_after_sixteen_attempts'] == e['no_current_after_sixteen_attempts']
    for left, right in zip(c['joint_first_current_phase_law'], e['joint_first_current_phase_law'], strict=True):
        assert sum(map(F, left['joint_phase_probabilities'])) == sum(map(F, right['joint_phase_probabilities']))
    assert c['joint_first_current_phase_law'][0]['joint_phase_probabilities'] == ['67/600', '1/24']
    assert e['joint_first_current_phase_law'][0] == c['joint_first_current_phase_law'][0]
    assert c['joint_first_current_phase_law'][1] != e['joint_first_current_phase_law'][1]


def test_exact_event_sums_are_not_normalized_after_discarding_unresolved_trials(sources):
    for marked in sources.values():
        finite = r.exact_stopped_state(marked, through=6)
        tail = r.exact_tail(marked, 6)
        assert 0 < r.trace(finite) < 1
        assert r.trace(finite)+r.trace(tail) == 1
        assert sum(r.phase_probabilities(finite)+r.phase_probabilities(tail)) == 1
        normalized_success = r.phase_probabilities(finite)[0]/r.trace(finite)
        assert normalized_success != r.phase_probabilities(r.exact_stopped_state(marked))[0]


def test_complete_latched_runtime_retains_impossible_flags_tail_and_checkpoint():
    result = r.run()
    assert result['all_pass']
    assert result['C_minus_E_plus_probability'] == '3644/69069'
    for row in result['finite_executions']:
        assert len(row['complete_flag_and_phase_records']) == 128
        assert sum(F(item['exact_probability']) for item in row['complete_flag_and_phase_records']) == 1
        assert row['impossible_flag_words_mass_diagnostic'] == 0
        assert row['checkpoint_complete_cq_distance_diagnostic'] == 0
        assert row['actual_branch_native_costs_checked']
        assert row['complete_terminal_density_residual_diagnostic'] < 3e-14
        assert row['complete_probability_residual_diagnostic'] < 3e-14
        assert F(row['unresolved_tail_mass']) <= F(row['uniform_all_input_tail_bound'])
    assert not result['scope']['physical_halt_device_or_clock_calibration_derived']
    assert not result['scope']['laboratory_measurement_identified_or_empirical_agreement_claimed']


def test_mistaking_no_current_for_identity_changes_the_measurement(sources):
    # Native zero-current events attenuate phase. Merely normalizing the
    # first-attempt moving outcome is a different, postselected protocol.
    for marked in sources.values():
        one = r.exact_stopped_state(marked, first_at=1)
        conditioned = r.phase_probabilities(one)[0]/r.trace(one)
        assert conditioned == F(67, 92)
        assert conditioned != r.phase_probabilities(r.exact_stopped_state(marked))[0]


@pytest.mark.parametrize('law', ['C', 'E'])
def test_two_attempt_latch_preserves_complete_reference_correlations(sources, law):
    marked = sources[law]
    p = r.program(marked, 2)
    initial = r.quantum.reference_preparation()
    history = r.runtime.execute(p, initial=initial)
    def apply(bit, density):
        return sum((np.kron(k, np.eye(2))@density@np.kron(k.conj().T, np.eye(2))
                    for k in marked.instrument.kraus[bit]), np.zeros_like(density))
    zero = apply(0, initial)
    prefixes = {(0, 0): apply(0, zero), (0, 1): apply(1, zero),
                (1, 0): np.zeros_like(initial), (1, 1): apply(1, initial)}
    read = p.events[-1].instruments[0]
    for branch in history.branches:
        flags, phase = branch.record[:2], branch.record[2]
        k = np.kron(read.kraus[phase][0], np.eye(2))
        expected = k@prefixes[flags]@k.conj().T
        np.testing.assert_allclose(branch.density, expected, atol=2e-14, rtol=0)
    total = sum(branch.density for branch in history.branches).reshape(66, 2, 66, 2)
    np.testing.assert_allclose(np.trace(total, axis1=0, axis2=2), np.eye(2)/2, atol=2e-14, rtol=0)


def test_corrupted_native_record_fails_before_relational_readout(monkeypatch):
    original = r.native.native
    def wrong(state, seam, coin):
        output, record = original(state, seam, coin)
        return output, (1-coin, record[1])
    monkeypatch.setattr(r.native, 'native', wrong)
    with pytest.raises(ArithmeticError, match='grammar'):
        r.source('C')


@pytest.mark.parametrize('count', [True, -1, 0, 7, 1., None])
def test_runtime_horizon_requires_a_bounded_integer(sources, count):
    with pytest.raises(ValueError):
        r.program(sources['C'], count)


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1., .01])
def test_invalid_numerical_diagnostics_cannot_pass(value):
    with pytest.raises(ArithmeticError):
        r.bounded(value, 'mutant')
