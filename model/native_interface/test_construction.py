"""Native herald/coin readout, coherent-service boundary and full outputs."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as c


@pytest.fixture(scope='module')
def result():
    return c.run()


def dephase_work(density, reference=1):
    tensor = density.reshape(66, reference, 66, reference)
    return (tensor*np.eye(66)[:, None, :, None]).reshape(density.shape)


def test_exact_native_branch_audit_and_readout_effects(result):
    for law, count in (('C', 120), ('E', 180)):
        row = result['native_reachability_and_record_audit'][law]
        assert row['native_branch_addresses_checked'] == 3960
        assert row['recorded_partial_injection_maps_checked'] == count
        assert row['exact_herald_probability_per_sign'] == '5/66'
        assert row['exact_total_useful_herald_probability'] == '5/33'
        assert row['exact_failure_probability'] == '28/33'
        values = list(map(F, row['orientation_pointer_plus_effect_diagonal']))
        assert values.count(F(1)) == values.count(F(0)) == 10
        assert values.count(F(1, 2)) == 46


def test_shared_mixer_pairing_is_the_actual_native_current_support():
    pairs = set(c.phase_pairs())
    assert pairs == set(c.quantum.phase_apparatus()[3])
    for law in ('C', 'E'):
        branch, = (branch for branch in c.fixed_maps(law) if branch.current == -1)
        assert pairs == {(output, initial) for initial, output in branch.mapping}
    b = c.mixer()
    np.testing.assert_allclose(b@b, np.eye(66), atol=3e-16, rtol=0)
    # The new service violates the one-output-per-column native grammar.
    assert np.count_nonzero(np.count_nonzero(b, axis=0) == 2) == 20


def test_all_preparation_failures_and_pointer_records_are_executed(result):
    assert len(result['complete_heralded_measurements']) == 6
    for row in result['complete_heralded_measurements']:
        probabilities = list(map(F, row['exact_joint_probabilities']))
        assert probabilities[0] == F(28, 33) and probabilities[1] == 0
        assert sum(probabilities) == 1
        assert sum(probabilities[2:]) == F(5, 33)
        assert probabilities[2:4] == probabilities[4:6]
        assert row['runtime_joint_probabilities'] == pytest.approx(list(map(float, probabilities)), abs=2e-14)
        assert row['conditional_native_preparation_vs_previous_phase_service_residual'] < 2e-14
        assert row['complete_terminal_density_residual_diagnostic'] < 2e-14
        assert row['checkpoint_complete_cq_distance_diagnostic'] == 0
    at_four = {row['candidate']: row for row in result['complete_heralded_measurements']
               if row['unread_native_attempts_after_herald'] == 4}
    assert at_four['C']['conditional_prepared_phase_match_probabilities'][0] == '1921/2592'
    assert at_four['E']['conditional_prepared_phase_match_probabilities'][0] == '12401/20000'
    joint = sum(F(at_four['C']['exact_joint_probabilities'][i])-F(at_four['E']['exact_joint_probabilities'][i])
                for i in (2, 4))
    assert joint == F(12259, 668250)
    assert result['four_attempt_discrimination']['unconditional_complete_record_total_variation'] == str(joint)


def test_matching_effects_do_not_select_the_old_terminal_instrument(result):
    for row in result['native_phase_pointer']:
        assert row['complete_effect_matrix_residual_diagnostic'] < 2e-14
        assert row['certain_plus_state_native_vs_old_pointer_trace_distance'] == pytest.approx(.5, abs=2e-14)
        assert not row['same_complete_postmeasurement_instrument_as_old_square_root_pointer']
        assert not row['certain_state_invariance_required_for_Lueders_selection_satisfied']
        assert not row['whole_carrier_effects_are_projective']


def test_forgetting_the_native_coin_weakens_the_pointer():
    # A fair assignment on j=0 loses the source coin correlation. It gives
    # the unsharp effect I/2+S/4, not the claimed phase effect I/2+S/2.
    rho, s, _, _ = c.quantum.phase_apparatus()
    b = c.mixer()
    effect = sum((1 if branch.current > 0 else 0 if branch.current < 0 else .5)
                 *b.conj().T@branch.matrix.conj().T@branch.matrix@b
                 for branch in c.fixed_maps('C'))
    np.testing.assert_allclose(effect, np.eye(66)/2+s/4, atol=3e-16, rtol=0)
    assert np.trace(effect@rho).real == pytest.approx(.75, abs=2e-14)


@pytest.mark.parametrize('law', ['C', 'E'])
def test_native_adaptive_record_laws_ignore_work_coherence_but_not_full_outputs(law):
    def instrument(seam):
        maps = c.fixed_maps(law, seam)
        return c.runtime.Instrument(kraus=tuple(tuple(b.matrix for b in maps if b.current == current)
                                                for current in (-1, 0, 1)), visible=True)
    seams = c.quantum.mixing.seams()
    first = instrument(seams[0])
    p = c.runtime.Program('native_only_adaptive_boundary', 66, np.eye(66)/66,
        (c.runtime.Event('first_current', instruments=(first,)),
         c.runtime.Event('adaptive_current', instruments=tuple(instrument(seams[i]) for i in (0, 1, 2)), control=0)))
    initial = c.reference_input()
    coherent = c.runtime.execute(p, initial=initial)
    diagonal = c.runtime.execute(p, initial=dephase_work(initial, 2))
    for left, right in zip(coherent.branches, diagonal.branches, strict=True):
        assert left.probability == pytest.approx(right.probability, abs=2e-14)
        np.testing.assert_allclose(dephase_work(left.density, 2), right.density, atol=2e-14, rtol=0)
        # Even the surviving reference state conditional on a native record
        # depends only on the input diagonal work/reference blocks.
        np.testing.assert_allclose(np.trace(left.density.reshape(66, 2, 66, 2), axis1=0, axis2=2),
                                   np.trace(right.density.reshape(66, 2, 66, 2), axis1=0, axis2=2), atol=2e-14, rtol=0)
    assert c.runtime.complete_distance(coherent, diagonal) > .1


@pytest.mark.parametrize('law', ['C', 'E'])
def test_entanglement_survives_successful_herald_and_complete_reference_audit(result, law):
    p = c.program(law, 1)
    prefix = c.runtime.execute(p, initial=c.reference_input(), stop_after=1)
    success = next(branch for branch in prefix.branches if branch.record == (1,))
    assert success.probability == pytest.approx(.25, abs=2e-14)
    normalized = success.density/success.probability
    partial_transpose = normalized.reshape(66, 2, 66, 2).transpose(0, 3, 2, 1).reshape(132, 132)
    assert np.linalg.eigvalsh(partial_transpose)[0] == pytest.approx(-.5, abs=2e-14)
    row = next(row for row in result['reference_checks'] if row['candidate'] == law)
    assert row['complete_reference_output_residual_diagnostic'] < 2e-14
    assert row['untouched_reference_marginal_residual_diagnostic'] < 2e-14
    assert not row['stationary_herald_probability_asserted_for_entangled_override']


def test_real_mixer_still_leaves_opposite_imaginary_phases_unresolved(result):
    for row in result['remaining_complex_phase_boundary']:
        assert row['opposite_imaginary_phase_state_trace_distance'] == pytest.approx(1., abs=2e-14)
        assert np.asarray(row['complete_native_readout_probabilities']) == pytest.approx(np.full((2, 2), .5), abs=2e-14)
        assert not row['full_complex_tomography_supplied_by_one_real_mixer']


def test_native_record_corruption_cannot_redefine_herald_or_pointer(monkeypatch):
    original = c.native.native
    def changed(state, seam, coin):
        output, record = original(state, seam, coin)
        return output, (1-coin, record[1])
    monkeypatch.setattr(c.native, 'native', changed)
    with pytest.raises(ArithmeticError, match='grammar'):
        c.fixed_maps('C')


@pytest.mark.parametrize('count', [None, True, 0, -1, 65, 1.])
def test_unsupported_evolution_count_is_rejected(count):
    with pytest.raises(ValueError):
        c.program('C', count)


@pytest.mark.parametrize('law', ['X', 'c', None, 1, True])
def test_unknown_native_law_is_rejected(law):
    with pytest.raises(ValueError):
        c.fixed_maps(law)


@pytest.mark.parametrize('value', [float('inf'), float('nan'), -1., .01])
def test_failed_numerical_contract_cannot_pass(value):
    with pytest.raises(ArithmeticError):
        c.bounded(value, 'mutant')


def test_portable_source_paths_and_scope(result):
    assert all('\\' not in path and not path.startswith('/') for path in result['sources'])
    scope = result['scope']
    assert scope['native_preparation_herald_and_binary_effects_constructed']
    assert scope['one_shared_balanced_mixer_service_supplied']
    for key in ('mixer_derived_as_an_executable_native_operation',
                'independent_target_density_or_POVM_received_by_runtime',
                'full_complex_observer_interface_or_physical_clock_derived',
                'new_empirical_measurement_or_validation'):
        assert not scope[key]
