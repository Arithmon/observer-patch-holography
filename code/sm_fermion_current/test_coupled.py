"""Exact-action derivatives, omitted forces and coherent coupled-read attacks."""
import copy
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_coupled as v


@pytest.fixture(scope='module')
def receipt():
    return v.load()


@pytest.fixture(scope='module')
def parent():
    return v.old.load(v.ROOT / v.PARENT_PATH)


def test_whole_coupled_history(receipt):
    result = v.verify(receipt)
    assert result['exact_rational_replay'] and result['same_declared_hybrid_action']
    assert result['continuation_events'] == 1020
    assert result['links'] == 144 and result['plaquettes'] == 108
    assert result['mean_Gauss'] and not result['operator_Gauss']
    assert result['occupied_family'] == 0 and result['vacuum_families'] == [1, 2]
    assert result['arbitrary_complex_Yukawa_coefficients'] == 27
    lo, hi = map(v.F, result['feedback_current_difference_interval'])
    assert 0 < lo < hi
    assert not result['quadratic_Maxwell_kinetic'] and not result['global_Z6_quotient_selected']


def test_fixed_law_derivatives_and_nonabelian_expected_reduction():
    assert v.fixed_electric_law_check()
    assert all(row == [0] * 11 for row in v.omitted_first_variation_check().values())


def test_pair_expectation_zero_does_not_annihilate_quantum_Higgs_source():
    # Basis vacuum, first occupied, second occupied, both occupied. A=cc+h.c.
    # has zero expectation in each fixed-number basis ket but nonzero variance
    # on vacuum/both-filled states. Thus a zero mean Higgs force is not a
    # quantum truncation with H identically zero as an operator.
    import sympy as s
    source = s.zeros(4)
    source[0, 3] = source[3, 0] = 1
    number = s.diag(0, 1, 1, 2)
    pair = s.zeros(4)
    pair[0, 3] = 1
    assert number * pair - pair * number == -2 * pair
    assert all(source[i, i] == 0 for i in range(4))
    assert (source * source)[0, 0] == 1


def test_unequal_internal_occupations_break_omitted_gauge_force():
    import sympy as s
    T = s.diag(s.Rational(1, 2), -s.Rational(1, 2))
    assert s.trace(T * s.eye(2)) == 0
    assert s.trace(T * s.diag(1, 0)) == s.Rational(1, 2)
    # Equal internal populations still have local nonabelian charge fluctuations.
    assert 2 * v.F(1, 64) * v.F(63, 64) == v.F(63, 2048) > 0


@pytest.mark.parametrize('flag', list(v.SCOPE))
def test_scope_rewrites_fail(receipt, flag):
    packet = copy.deepcopy(receipt)
    packet['scope'][flag] = not packet['scope'][flag]
    with pytest.raises(ValueError, match='scope'):
        v.validate_custody(packet)


@pytest.mark.parametrize('field', ['gauge_group', 'electric_kinetic', 'kinetic_lambda', 'electric_duration',
                                  'kinetic_parameter_boundary', 'Wilson_loop', 'magnetic_kick', 'schedule'])
def test_common_action_cannot_be_silently_changed(receipt, field):
    packet = copy.deepcopy(receipt)
    packet['law'][field] = 'quadratic Maxwell or retuned step'
    with pytest.raises(ValueError, match='action'):
        v.validate_custody(packet)


@pytest.mark.parametrize('field', ['families', 'occupied_family', 'vacuum_families', 'Yukawa_coefficients',
                                  'anomalous_covariance', 'nonabelian_links', 'Higgs_and_conjugate_momentum'])
def test_omitted_force_hypotheses_are_required(receipt, field):
    packet = copy.deepcopy(receipt)
    packet['mean_field_restriction'][field] = 'unrestricted nonzero value'
    with pytest.raises(ValueError, match='restriction'):
        v.validate_custody(packet)


def reseal(run):
    previous = run['prefix']['final_event_hash']
    for event in run['events']:
        event['previous_hash'] = previous
        event['event_hash'] = v.old.hashed({k: x for k, x in event.items() if k != 'event_hash'})
        previous = event['event_hash']
    run['final_event_hash'] = previous


@pytest.mark.parametrize('case', ['drift_sign', 'disabled_drift', 'missing_electric_read', 'stale_E_writer',
                                  'Wilson_sign', 'missing_loop_link', 'nonunitary_return', 'return_current_sign',
                                  'stale_return_link', 'stale_return_E', 'fake_read_parent', 'omit_plaquette',
                                  'Gauss_readout', 'energy_interval', 'feedback_witness', 'changed_prefix',
                                  'fake_control_dependency'])
def test_resealed_wrong_dynamics_and_consumption_fail(receipt, parent, case):
    run = copy.deepcopy(receipt['runs'][3 if case == 'fake_control_dependency' else 1])
    drift = next(e for e in run['events'] if e['metadata'] == {'edge': 14} and 'drift' in e['operation'])
    matter = next(e for e in run['events'] if e['operation'] == 'coupled_matter_return')
    if case == 'drift_sign':
        drift['writes'][0]['value'][1] = str(-v.F(drift['writes'][0]['value'][1]))
    elif case == 'disabled_drift':
        drift['writes'][0]['value'] = ['1', '0']
    elif case == 'missing_electric_read':
        drift['reads'] = [r for r in drift['reads'] if not r['port'].startswith('electric/')]
        drift['parents'] = sorted({r['writer'] for r in drift['reads']})
    elif case == 'stale_E_writer':
        row = next(r for r in drift['reads'] if r['port'] == 'electric/14')
        row['writer'], row['version'], row['value'] = 334, 0, '0'
        drift['parents'] = sorted({r['writer'] for r in drift['reads']})
    elif case in ('Wilson_sign', 'missing_loop_link'):
        kick = next(e for e in run['events'] if e['operation'] == 'Wilson_kick' and
                    any(v.F(w['value']) for w in e['writes']))
        if case == 'Wilson_sign':
            row = next(w for w in kick['writes'] if v.F(w['value']))
            row['value'] = str(-v.F(row['value']))
        else:
            kick['reads'] = [r for r in kick['reads'] if r['port'] != next(r['port'] for r in kick['reads'] if r['port'].startswith('link/'))]
            kick['parents'] = sorted({r['writer'] for r in kick['reads']})
    elif case == 'nonunitary_return':
        next(w for w in matter['writes'] if w['port'].startswith('orbital/e_c'))['value'][0] = '2'
    elif case == 'return_current_sign':
        row = next(w for w in matter['writes'] if w['port'] == 'factor_current/14')
        row['value'] = str(-v.F(row['value']))
    elif case in ('stale_return_link', 'stale_return_E'):
        key = 'link/14' if case == 'stale_return_link' else 'electric/14'
        row = next(r for r in matter['reads'] if r['port'] == key)
        former = next(e for e in reversed(run['prefix']['events']) if any(w['port'] == key for w in e['writes']))
        value = next(w for w in former['writes'] if w['port'] == key)
        row.update(writer=former['id'], version=value['version'], value=value['value'])
        matter['parents'] = sorted({r['writer'] for r in matter['reads']})
    elif case == 'fake_read_parent':
        matter['parents'].append(0)
    elif case == 'omit_plaquette':
        run['events'].pop(next(i for i, e in enumerate(run['events']) if e['operation'] == 'Wilson_kick'))
    elif case == 'Gauss_readout':
        run['checkpoints'][-1]['readout']['gauss_residual'][21] = '1'
    elif case == 'energy_interval':
        run['checkpoints'][0]['readout']['magnetic_energy_interval'] = ['0', '0']
    elif case == 'feedback_witness':
        run['feedback_witness']['electric_drift_event'] -= 1
    elif case == 'changed_prefix':
        run['prefix']['checkpoints'][-1]['readout']['charge'][21] = '0'
    else:
        drift['reads'].append({'port': 'electric/14', 'version': 1, 'writer': 480, 'value': '0'})
        drift['parents'] = sorted({r['writer'] for r in drift['reads']})
    reseal(run)
    with pytest.raises(ValueError):
        v.replay(run, receipt['carrier'], parent)


@pytest.mark.parametrize('case', ['hash', 'bytes', 'missing', 'parent', 'extra'])
def test_source_and_parent_custody(receipt, case):
    packet = copy.deepcopy(receipt)
    key = next(iter(packet['source_pins']))
    if case == 'hash':
        packet['source_pins'][key]['sha256'] = '0' * 64
    elif case == 'bytes':
        packet['source_pins'][key]['bytes'] += 1
    elif case == 'missing':
        del packet['source_pins'][key]
    elif case == 'parent':
        packet['parent']['sha256'] = '0' * 64
    else:
        packet['invented'] = True
    with pytest.raises(ValueError):
        v.validate_custody(packet)
