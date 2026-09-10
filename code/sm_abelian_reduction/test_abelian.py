"""Actual replay, scope controls and rehashed operation/custody attacks."""
import copy
import json
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_abelian as verifier
import pilot


@pytest.fixture(scope='module')
def receipt(): return verifier.load()


@pytest.fixture(scope='module')
def replayed(receipt): return verifier.verify(receipt)


def test_fresh_independent_replay(replayed):
    assert replayed['verified'] and replayed['mathematical_replay']
    assert replayed['sites'] == 64 and replayed['links'] == 144
    assert replayed['exact_first_electric_kick'] == '1/25000'
    assert replayed['neighbor_intensity_response'] > 1e-6
    assert replayed['nonvacuous_off_support_sites'] == 57
    assert not replayed['continuous_history_error_enclosure']


def test_custody_is_not_mathematical_replay(receipt):
    bad = copy.deepcopy(receipt); bad['exact_first_edge_control']['current'] = '0'
    assert verifier.validate_custody(bad)['mathematical_replay'] is False
    with pytest.raises(ValueError, match='exact current witness'):
        verifier.verify(bad)


@pytest.mark.parametrize('field', list(verifier.SCOPE))
def test_scope_rewrite_fails(receipt, field):
    bad = copy.deepcopy(receipt); bad['scope'][field] = not bad['scope'][field]
    with pytest.raises(ValueError, match='scope'): verifier.validate_custody(bad)


@pytest.mark.parametrize('field', list(verifier.UNITS))
def test_units_and_clock_rewrite_fails(receipt, field):
    bad = copy.deepcopy(receipt); bad['units'][field] = 'measured seconds'
    with pytest.raises(ValueError, match='units'): verifier.validate_custody(bad)


@pytest.mark.parametrize('field', ['kappa_1', 'w_standard', 'gauss', 'current', 'schedule', 'log_chart', 'boundary'])
def test_action_law_rewrite_fails(receipt, field):
    bad = copy.deepcopy(receipt); bad['law'][field] = 'uncertified alternative'
    with pytest.raises(ValueError, match='law'): verifier.validate_custody(bad)


@pytest.mark.parametrize('case', ['parent_hash', 'parent_bytes', 'parent_bool', 'source_hash', 'source_missing', 'parent_commit', 'extra_root'])
def test_source_custody_fails(receipt, case):
    bad = copy.deepcopy(receipt)
    if case == 'parent_commit': bad['parent_commit'] = '0'*40
    elif case == 'extra_root': bad['invented'] = True
    else:
        field = 'parent_pins' if case.startswith('parent') else 'source_pins'; key = next(iter(bad[field]))
        if case == 'source_missing': del bad[field][key]
        elif case == 'parent_bool': bad[field][key]['bytes'] = True
        elif case == 'parent_bytes': bad[field][key]['bytes'] += 1
        else: bad[field][key]['sha256'] = '0'*64
    with pytest.raises(ValueError): verifier.validate_custody(bad)


@pytest.mark.parametrize('case', ['basis', 'metric', 'current', 'yukawa', 'omitted', 'false_negative'])
def test_symbolic_reduction_fails(receipt, case):
    row = copy.deepcopy(receipt['reduction'])
    if case == 'basis': row['code_q'] = ['1/2', '1/2']
    if case == 'metric': row['standard_K'] = ['k2', 'k1']
    if case == 'current': row['weak_and_hypercharge_currents'][0] = '1'
    if case == 'yukawa': row['unrestricted_complex_yukawa_slots'] = 0
    if case == 'omitted': row['omitted_equations'].remove('lower_Higgs')
    if case == 'false_negative': row['hypercharge_only_orthogonal_gauge_defect_at_k1_k2_1'] = '0'
    with pytest.raises(ValueError, match='reduction'): verifier.reduction_check(row)


@pytest.mark.parametrize('case', ['source_record', 'address', 'mass', 'conductance', 'electric', 'curl_sign', 'magnetic', 'bool_index'])
def test_source_geometry_fails(receipt, case):
    geo = copy.deepcopy(receipt['geometry'])
    if case == 'source_record': geo['sites'][21]['source_record'][0] += 1
    if case == 'address': geo['sites'][21]['address_Qphi'][0][0] += 1
    if case == 'mass': geo['sites'][21]['mass'] *= 2
    if case == 'conductance': geo['edges'][0]['conductance'] *= 2
    if case == 'electric': geo['edges'][0]['electric_mass'] *= 2
    if case == 'curl_sign': geo['plaquettes'][0]['boundary'][0][1] *= -1
    if case == 'magnetic': geo['plaquettes'][0]['magnetic_weight'] *= 2
    if case == 'bool_index': geo['sites'][0]['index'][0] = True
    with verifier.mp.workdps(60), pytest.raises(ValueError, match='geometry'):
        verifier.geometry_check(geo)


def reseal(run):
    previous = '0'*64
    for ev in run['events']:
        ev['previous_hash'] = previous
        ev['event_hash'] = verifier.hashed({k: v for k, v in ev.items() if k != 'event_hash'})
        previous = ev['event_hash']
    run['final_event_hash'] = previous


@pytest.mark.parametrize('case', ['wrong_current', 'missing_feedback', 'stale_writer', 'stale_version', 'missing_parent',
                                  'missing_read', 'wrong_cohort', 'extra_write', 'omitted_event', 'bool_id', 'copied_checkpoint', 'stale_real_same_port_writer'])
def test_resealed_event_forgeries_fail(receipt, case):
    run = copy.deepcopy(receipt['runs'][1])
    ev = next(e for e in run['events'] if e['operation'] == 'edge_kick' and e['metadata']['edge'] == 54)
    # Edge54 is the first outgoing edge at source21; each mutation is nontrivial.
    if case in ('wrong_current', 'missing_feedback'):
        row = next(w for w in ev['writes'] if w['port'].startswith('P/'))
        assert row['value'] != 0
        row['value'] = -row['value'] if case == 'wrong_current' else 0.0
    if case == 'stale_writer': ev['reads'][0]['writer'] = 0
    if case == 'stale_version': ev['reads'][0]['version'] += 1
    if case == 'missing_parent': ev['parents'].pop()
    if case == 'missing_read': ev['reads'].pop()
    if case == 'wrong_cohort': run['cohort'] = 'gauge_copy'
    if case == 'extra_write': ev['writes'].append({'port': 'address/21', 'version': 1, 'value': 0.0})
    if case == 'omitted_event': run['events'].pop(ev['id'])
    if case == 'bool_id': run['events'][0]['id'] = False
    if case == 'copied_checkpoint': run['checkpoints'][0]['values']['psi/21'] = [0.2, 0.0]
    if case == 'stale_real_same_port_writer':
        ev = next(e for e in run['events'] if e['operation'] == 'edge_kick' and e['metadata'] == {'edge': 54, 'step': 2, 'part': 'first'})
        row = next(r for r in ev['reads'] if r['port'] == 'a/54')
        row['writer'] = 64+54; row['version'] = 0
        ev['parents'] = sorted({r['writer'] for r in ev['reads']})
    reseal(run)
    with verifier.mp.workdps(60):
        geo = verifier.geometry_check(receipt['geometry'])
        with pytest.raises(ValueError): verifier.replay(run, geo)


@pytest.mark.parametrize('text', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":-Infinity}', '{"x":1e999}'],
                         ids=['duplicate', 'nan', 'inf', 'negative_inf', 'overflow'])
def test_strict_json(tmp_path, text):
    p = tmp_path/'bad.json'; p.write_text(text)
    with pytest.raises(ValueError): verifier.load(p)


def test_oversize_rejected(tmp_path):
    p = tmp_path/'bad.json'
    with p.open('wb') as f: f.truncate(18000001)
    with pytest.raises(ValueError, match='size'): verifier.load(p)


def test_disjoint_actual_edge_maps_commute():
    # Derive two disjoint actual potential kicks from the producer's operations.
    geo = pilot.geometry(); run = pilot.execute(geo, 'phase_intervention')
    events = [e for e in run['events'] if e['operation'] == 'edge_kick']
    a = next(e for e in events if e['metadata']['edge'] == 54)
    b = next(e for e in events if e['metadata']['edge'] == 0)
    assert not ({r['port'] for r in a['reads']} & {r['port'] for r in b['reads']})
    def apply(order):
        state = {r['port']: (verifier.mp.mpc(*r['value']) if isinstance(r['value'], list) else verifier.mp.mpf(r['value'])) for event in order for r in event['reads']}
        for event in order:
            e = event['metadata']['edge']; row = geo['edges'][e]; i, j = row['ends']
            qi, qj, pi, pj, akey, P = f'psi/{i}', f'psi/{j}', f'pi/{i}', f'pi/{j}', f'a/{e}', f'P/{e}'
            U = verifier.mp.exp(verifier.mp.j*state[akey]); c = verifier.mp.mpf(row['conductance'])
            fi = c*(U*state[qj]-state[qi]); fj = c*(verifier.mp.conj(U)*state[qi]-state[qj])
            J = -2*c*verifier.mp.im(verifier.mp.conj(state[qi])*U*state[qj])
            state[pi] += fi/400; state[pj] += fj/400; state[P] += J/400
        return state
    with verifier.mp.workdps(60): assert apply([a, b]) == apply([b, a])
    # Their actual state transforms, not their audit hashes, commute.
    assert a['event_hash'] != b['event_hash']


def test_diagnostic_does_not_claim_independent_replay(receipt):
    bad = copy.deepcopy(receipt)
    bad['producer_only_diagnostic']['reverse_potential_order_final_max_difference'] = 0.123
    assert verifier.validate_custody(bad)['accepted_custody']
    bad['producer_only_diagnostic']['independently_replayed'] = True
    with pytest.raises(ValueError): verifier.validate_custody(bad)
