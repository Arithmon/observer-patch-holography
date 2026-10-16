"""Independent reference, complete source-map and retained-outcome controls."""
from collections import Counter
from fractions import Fraction as F
import math
import json

import numpy as np
import pytest

from . import construction as c


def test_palm_reference_is_faithful_and_pair_census_is_not_uniform():
    degree, palm, pairs, mask = c.reference()
    assert Counter(degree) == {8: 30, 10: 36}
    assert sum(degree) == 600
    assert all(p == F(d, 600) > 0 for d, p in zip(degree, palm, strict=True))
    ratios = Counter(palm[a]/palm[b] for a, b in pairs)
    assert ratios == {F(4, 5): 2, F(5, 4): 2, F(1): 6}
    assert mask == (10, 18, 25, 33, 41, 42, 57, 58)
    assert sum(palm[i] for pair in pairs for i in pair) == F(23, 75)


@pytest.mark.parametrize('setting,expected', [(0, F(1)), (1, F(4, 5)), (2, F(3, 5)), (4, F(1))])
def test_exact_pair_interference_law_and_original_binary_bound(setting, expected):
    degree, _, pairs, mask = c.reference()
    # Independent state-vector calculation; no receipt or returned formula.
    match = sum(abs(((1j**setting if degree[a] == 8 else 1)
                    +(1j**setting if degree[b] == 8 else 1))/2)**2 for a, b in pairs)/10
    assert abs(match-float(expected)) < 1e-15
    assert c.probabilities(setting) == expected
    matrix, target = c.local_matrix(setting), c.ideal_local(setting)
    bound = c.certify(setting)
    squared = sum((F(float((x-y).real))**2+F(float((x-y).imag))**2
                   for x, y in zip(matrix.flat, target.flat, strict=True)), F())
    # Difference subtraction happens in binary above, so check an independent
    # triangle margin; the production bound converts each operand first.
    assert math.sqrt(float(squared)) <= float(bound['operator_error_upper'])+1e-16
    assert 0 <= bound['operator_error_upper'] < F(1, 10**9)
    assert bound['half_diamond_upper'] == bound['operator_error_upper']*(2+bound['operator_error_upper'])/2


def test_modular_target_identification_is_only_the_active_corner():
    degree, _, pairs, _ = c.reference()
    codes = [c.compiled.circuit.encode(state) for state in c.quantum.STATES]
    active = {i for pair in pairs for i in pair}
    target = np.kron(np.eye(4), c.ideal_local(1))
    for i in active:
        assert target[codes[i], codes[i]] == (1j if degree[i] == 8 else 1)
    omitted = [i for i in range(66) if i not in active and degree[i] == 8]
    assert len(omitted) == 22
    assert all(target[codes[i], codes[i]] == 1 for i in omitted)
    assert all(target[codes[i], codes[i]] != 1j for i in omitted)


def test_all_512_source_columns_are_the_declared_spectator_lift():
    full = c.full_matrix(1)
    for spectator in range(4):
        np.testing.assert_array_equal(full[128*spectator:128*(spectator+1),
                                          128*spectator:128*(spectator+1)], c.local_matrix(1))
    assert np.count_nonzero(full[:128, 128:]) == 0
    # A padding error cannot hide behind correctness on twenty active columns.
    damaged = c.local_matrix(1).copy(); damaged[127, 126] = 2**-16
    assert c.certify(1, damaged)['operator_error_upper'] >= F(1, 2**16)


def test_unitary_inverse_phase_and_common_phase_fail_complete_operator_target():
    correct = c.ideal_local(1)
    for wrong in (correct.conj(), -correct):
        np.testing.assert_array_equal(wrong.conj().T@wrong, np.eye(128))
        assert c.certify(1, wrong)['operator_error_upper'] > 1
    # The cosine readout alone aliases inverse phase and whole turns.
    assert c.probabilities(0) == c.probabilities(4)


def test_three_gate_bound_accounts_for_amplified_raw_maps():
    # Scalar multiples of identity saturate this telescope exactly, so an
    # unweighted sum fails even in the simplest complete channel example.
    eps_b, eps_d = F(1, 4), F(1, 3)
    def bound(eps):
        return {'half_diamond_upper': eps*(2+eps)/2, 'channel_norm_upper': (1+eps)**2}
    b, d = bound(eps_b), bound(eps_d)
    assert c.composite_bound(d, b) == (((1+eps_b)**2*(1+eps_d))**2-1)/2
    assert c.composite_bound(d, b) > 2*b['half_diamond_upper']+d['half_diamond_upper']


def test_complete_native_extensions_keep_padding_and_zero_outcomes():
    herald, pointer = c.interface.herald_and_pointer('C')
    code = np.zeros(512); code[511] = 1
    for instrument, expected in ((c._lift_native(herald, herald=True), (1., 0., 0.)),
                                 (c._lift_native(pointer), (.5, .5))):
        masses = tuple(sum(float(np.vdot(k@code, k@code).real) for k in group)
                       for group in instrument.kraus)
        np.testing.assert_allclose(masses, expected, atol=2e-16, rtol=0)


def test_full_runtime_preserves_failed_preparation_records_resources_and_resume():
    result = c.execute_setting('C', 1, checkpoint=True)
    assert result['exact_joint_probabilities'] == ['28/33', '0', '2/33', '1/66', '2/33', '1/66']
    assert result['checkpoint_complete_distance_diagnostic'] < 1e-12
    assert result['full_512_terminal_density_residual_diagnostic'] < 1e-12
    for branch in result['realized_branch_resources']:
        if branch['record'][0] == 0:
            assert 'compiled_modular_pulses' not in branch['resources']
        else:
            assert branch['resources']['compiled_mixer_pulses'] == 68


def test_zero_probability_control_retains_raw_leakage_and_all_hidden_factors():
    result = c.execute_setting('C', 0)
    assert result['exact_joint_probabilities'] == ['28/33', '0', '5/66', '0', '5/66', '0']
    assert len(result['executed_joint_probabilities']) == 6
    assert all(p >= 0 for p in result['executed_joint_probabilities'])
    # Nominally zero success/read-one branches retain source defects instead
    # of being silently replaced by their exact ideal zero probabilities.
    assert result['executed_joint_probabilities'][3] > 0
    assert result['executed_joint_probabilities'][5] > 0


def test_factor_input_is_original_binary_mixture_and_preparation_costs_survive():
    embedding = c.compiled.embedding()
    initial = embedding@embedding.T/66
    p = c.runtime.Program('factor_identity', 512, initial, (),
                          preparation_resources=(('prepared_inputs', 3),))
    branch, = c.factor_execute(p)
    np.testing.assert_array_equal(branch.density, initial)
    assert branch.resources == (('prepared_inputs', 3),)
    codes = np.flatnonzero(np.diag(initial))
    for coherent in (False, True):
        changed = initial.copy()
        if coherent:
            changed[codes[0], codes[1]] = changed[codes[1], codes[0]] = 1e-4
        else:
            changed[codes[0], codes[0]] += 1e-4
            changed[codes[1], codes[1]] -= 1e-4
        wrong = c.runtime.Program('wrong_factor_input', 512, changed, ())
        with pytest.raises(ValueError, match='exactly the declared binary uniform native input'):
            c.factor_execute(wrong)


def test_uniform_reference_and_native_evolution_are_distinct_controls():
    uniform = c.program('C', 2, uniform_reference=True)
    assert json.loads(uniform.binding)['preparation_state'] == 'I66/66, not Palm'
    assert json.loads(uniform.binding)['reference_state'] == 'I66/66'
    assert json.loads(uniform.events[2].instruments[1].binding)['quarter_turns'] == 0
    assert c.quantum.RATES['C'] != c.quantum.RATES['E']
    assert all((1+rate)/2 != c.probabilities(1) for rate in c.quantum.RATES.values())


@pytest.mark.parametrize('bad', [True, -1, 3, 1., '1', None])
def test_phase_settings_reject_alias_types_and_undeclared_values(bad):
    for operation in (c.word, c.local_matrix, c.full_matrix, c.ideal_local, c.certify, c.probabilities):
        with pytest.raises(ValueError):
            operation(bad)


def test_hostile_subclass_is_rejected_before_metadata_and_cache_is_immutable():
    class Hostile(np.ndarray):
        def __getattribute__(self, name):
            if name in ('shape', 'dtype', 'ndim', 'copy'):
                raise AssertionError('untrusted dispatch')
            return super().__getattribute__(name)
    with pytest.raises(ValueError, match='base ndarray'):
        c.certify(1, np.eye(128).view(Hostile))
    with pytest.raises(ValueError):
        c.local_matrix(1)[0, 0] = 9
    with pytest.raises(ValueError):
        c.certify(1, np.full((128, 128), np.nan))
