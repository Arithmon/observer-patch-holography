"""Independent CAR/midpoint checks and resealed scientific mutations."""
import copy
import itertools
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_current as v


def ladder(state, mode, creation):
    """Exterior-algebra creation/annihilation on ordered occupation bits."""
    occupied = bool(state & (1 << mode))
    if occupied == creation:
        return None, 0
    sign = (-1) ** ((state & ((1 << mode) - 1)).bit_count())
    return state ^ (1 << mode), sign


def composed(state, first, second):
    middle, a = ladder(state, *second)
    if not a:
        return {}
    final, b = ladder(middle, *first)
    return {final: a * b} if b else {}


def test_exterior_algebra_has_car_and_pauli_exclusion():
    # Exhaustive four-mode Fock check, including all fermionic permutation signs.
    # The actual rank15 state uses the same exterior-algebra construction, not
    # independent commuting oscillator coordinates or an unnormalized wave.
    for state in range(16):
        for i, j in itertools.product(range(4), repeat=2):
            for dag_j in (False, True):
                ab = composed(state, (i, False), (j, dag_j))
                ba = composed(state, (j, dag_j), (i, False))
                actual = {k: ab.get(k, 0) + ba.get(k, 0) for k in ab.keys() | ba.keys()}
                actual = {k: n for k, n in actual.items() if n}
                assert actual == ({state: 1} if dag_j and i == j else {})
            assert composed(state, (i, True), (i, True)) == {}


@pytest.mark.parametrize('axis', range(3))
@pytest.mark.parametrize('charge', [1, -4, 2, -3, 6])
def test_midpoint_unitarity_and_current_continuity(axis, charge):
    # A dense complex two-spin witness with a nontrivial compact connection.
    U = (v.F(3, 5), v.F(4, 5))
    transport = v.power(U, charge)
    old_u = ((v.F(1, 3), v.F(2, 7)), (v.F(-3, 11), v.F(1, 5)))
    old_v = ((v.F(2, 9), v.F(-1, 4)), (v.F(4, 13), v.F(1, 8)))
    # Solve the four complex linear equations by independent SymPy elimination.
    import sympy as s
    z = lambda a: s.Rational(a[0].numerator, a[0].denominator) + s.I * s.Rational(a[1].numerator, a[1].denominator)
    sigmas = (s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]),
              s.diag(1, -1))
    K = s.I * sigmas[axis] * z(transport)
    h = s.zeros(4)
    h[:2, 2:] = K
    h[2:, :2] = K.conjugate().T
    old = s.Matrix([z(a) for a in old_u + old_v])
    new = (s.eye(4) + s.I * h / 4).inv() * (s.eye(4) - s.I * h / 4) * old
    pair = lambda a: (v.F(str(s.re(s.expand(a)))), v.F(str(s.im(s.expand(a)))))
    new_u, new_v = tuple(pair(a) for a in new[:2]), tuple(pair(a) for a in new[2:])
    imag = v.midpoint_residual(old_u, old_v, new_u, new_v, axis, transport)
    current = -2 * charge * imag
    change = charge * (sum(map(v.norm, new_u)) - sum(map(v.norm, old_u)))
    assert change == -v.DT * current
    assert current != 0
    mutated = (v.add(new_u[0], (v.F(1, 1000), v.F(0))), new_u[1])
    with pytest.raises(ValueError, match='midpoint'):
        v.midpoint_residual(old_u, old_v, mutated, new_v, axis, transport)


@pytest.mark.parametrize('bad', ['0.5', '1/2/3', '2/4', '+1', '01', '1/0', True, 1])
def test_noncanonical_rational_rejected(bad):
    with pytest.raises(ValueError):
        v.fraction(bad)


@pytest.mark.parametrize('text', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":Infinity}', '{"x":0.5}'])
def test_strict_json(tmp_path, text):
    path = tmp_path / 'bad.json'
    path.write_text(text)
    with pytest.raises(ValueError):
        v.load(path)


@pytest.fixture(scope='module')
def receipt():
    return v.load()


def test_full_producer_independent_replay(receipt):
    answer = v.verify(receipt)
    assert answer['exact_rational_replay']
    assert answer['occupied_fermions'] == answer['endpoint_Slater_projector_rank'] == 15
    assert answer['one_particle_modes'] == 1920
    assert answer['events'] == 1833 and answer['local_factors'] == 432
    assert answer['expectation_Gauss'] and not answer['operator_Gauss']
    assert answer['initial_local_charge_variance'] == '945/512'


@pytest.mark.parametrize('field', list(v.SCOPE))
def test_scope_upgrade_rejected(receipt, field):
    changed = copy.deepcopy(receipt)
    changed['scope'][field] = not changed['scope'][field]
    with pytest.raises(ValueError, match='scope'):
        v.validate_custody(changed)


@pytest.mark.parametrize('field', list(v.LAW))
def test_supplied_action_cannot_be_rewritten(receipt, field):
    changed = copy.deepcopy(receipt)
    changed['law'][field] = 'different or physically selected law'
    with pytest.raises(ValueError, match='law'):
        v.validate_custody(changed)


@pytest.mark.parametrize('case', ['missing', 'charge', 'multiplicity', 'chirality', 'spin_multiplicity'])
def test_incomplete_or_wrong_generation_rejected(receipt, case):
    rows = copy.deepcopy(receipt['multiplets'])
    if case == 'missing':
        rows.pop()
    elif case == 'charge':
        rows[1]['integer_charge'] = 4
    elif case == 'multiplicity':
        rows[0]['multiplicity'] = 3
    elif case == 'chirality':
        rows[1]['color_dimension_signed'] = 3
    else:
        rows[0]['multiplicity'] *= 2
    with pytest.raises(ValueError, match='generation'):
        v.multiplet_check(rows)


@pytest.mark.parametrize('case', ['source', 'orientation', 'axis', 'missing_link', 'bool_index'])
def test_source_geometry_mutations_rejected(receipt, case):
    changed = copy.deepcopy(receipt['geometry'])
    if case == 'source':
        changed['sites'][21]['source_record'][0] += 1
    elif case == 'orientation':
        changed['edges'][14]['ends'].reverse()
    elif case == 'axis':
        changed['edges'][14]['axis'] = 2
    elif case == 'missing_link':
        changed['edges'].pop()
    else:
        changed['sites'][0]['index'][0] = True
    with pytest.raises(ValueError, match='geometry'):
        v.geometry_check(changed)


@pytest.mark.parametrize('case', ['digest', 'bytes', 'missing', 'bool_size', 'extra'])
def test_exact_source_content_custody(receipt, case):
    changed = copy.deepcopy(receipt)
    key = next(iter(changed['source_pins']))
    if case == 'digest':
        changed['source_pins'][key]['sha256'] = '0' * 64
    elif case == 'bytes':
        changed['source_pins'][key]['bytes'] += 1
    elif case == 'missing':
        del changed['source_pins'][key]
    elif case == 'bool_size':
        changed['source_pins'][key]['bytes'] = True
    else:
        changed['invented'] = True
    with pytest.raises(ValueError):
        v.validate_custody(changed)


def reseal(run):
    previous = '0' * 64
    for event in run['events']:
        event['previous_hash'] = previous
        event['event_hash'] = v.hashed({k: x for k, x in event.items() if k != 'event_hash'})
        previous = event['event_hash']
    run['final_event_hash'] = previous


@pytest.mark.parametrize('case', ['current_sign', 'missing_feedback', 'nonunitary', 'nonunitary_link',
                                  'wrong_spin', 'stale_writer', 'stale_version', 'missing_read',
                                  'missing_parent', 'coherent_stale_read', 'omit_factor', 'extra_write',
                                  'bool_id', 'readout_gauss', 'readout_projector', 'zero_variance',
                                  'gauge_without_links', 'gauge_wrong_charge', 'corrupt_checkpoint_join'])
def test_resealed_mathematical_and_read_from_mutations(receipt, case):
    run = copy.deepcopy(receipt['runs'][2 if case.startswith('gauge_') else 1])
    event = next(e for e in run['events'] if e['operation'] == 'edge_factor' and e['metadata']['edge'] == 14)
    by_port = {r['port']: r for r in event['writes']}
    if case == 'current_sign':
        row = by_port['factor_current/14']
        assert v.F(row['value']) != 0
        row['value'] = str(-v.F(row['value']))
    elif case == 'missing_feedback':
        row = by_port['electric/14']
        assert v.F(row['value']) != 0
        row['value'] = '0'
    elif case == 'nonunitary':
        by_port['orbital/e_c/21/0']['value'][0] = '2'
    elif case == 'wrong_spin':
        by_port['orbital/e_c/21/0']['value'], by_port['orbital/e_c/21/1']['value'] = (
            by_port['orbital/e_c/21/1']['value'], by_port['orbital/e_c/21/0']['value'])
    elif case == 'nonunitary_link':
        link = next(e for e in run['events'] if e['operation'] == 'prepare_link')
        next(w for w in link['writes'] if w['port'].startswith('link/'))['value'] = ['2', '0']
    elif case in ('stale_writer', 'stale_version', 'missing_read', 'missing_parent'):
        if case == 'stale_writer':
            event['reads'][0]['writer'] = 0
        elif case == 'stale_version':
            event['reads'][0]['version'] += 1
        elif case == 'missing_read':
            event['reads'].pop()
        else:
            event['parents'].pop()
    elif case == 'coherent_stale_read':
        read = next(r for r in event['reads'] if r['port'] == 'orbital/e_c/21/0')
        former = next(e for e in run['events'] if e['operation'] == 'prepare_orbital' and e['metadata'] == {'multiplet': 'e_c', 'site': 21})
        value = next(w for w in former['writes'] if w['port'] == read['port'])
        read.update(writer=former['id'], version=value['version'], value=value['value'])
        event['parents'] = sorted({r['writer'] for r in event['reads']})
    elif case == 'omit_factor':
        run['events'].pop(event['id'])
    elif case == 'extra_write':
        event['writes'].append({'port': 'electric/143', 'version': 1, 'value': '0'})
    elif case == 'bool_id':
        run['events'][0]['id'] = False
    elif case == 'readout_gauss':
        run['checkpoints'][-1]['readout']['gauss_residual'][21] = '1'
    elif case == 'readout_projector':
        run['checkpoints'][-1]['readout']['orbital_norms']['e_c'] = '2'
    elif case == 'zero_variance':
        run['checkpoints'][0]['readout']['charge_variance'][21] = '0'
    elif case.startswith('gauge_'):
        gauge = next(e for e in run['events'] if e['operation'] == 'gauge_transform')
        if case == 'gauge_without_links':
            row = next(w for w in gauge['writes'] if w['port'].startswith('link/') and w['value'] != ['1', '0'])
            row['value'] = ['1', '0']
        else:
            row = next(w for w in gauge['writes'] if w['port'] == 'orbital/e_c/0/0')
            row['value'] = ['3/40', '0']
    elif case == 'corrupt_checkpoint_join':
        run['checkpoints'][-1]['event'] -= 1
    reseal(run)
    with pytest.raises(ValueError):
        v.replay(run, v.geometry_check(receipt['geometry']))


def test_correct_hashes_do_not_certify_mathematics(receipt):
    changed = copy.deepcopy(receipt)
    changed['comparisons']['initial_local_charge_variance'] = '0'
    assert v.validate_custody(changed)['mathematical_replay'] is False
    with pytest.raises(ValueError, match='comparison'):
        v.verify(changed)
