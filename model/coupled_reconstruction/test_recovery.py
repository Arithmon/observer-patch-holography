"""Whole-field local marginals identify derivatives without isolated services."""
import ast
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm

from model.coupled_reconstruction import recovery as r


def operators(d, omega):
    annihilation = np.diag(np.sqrt(np.arange(1, d+8)), 1)
    q = (annihilation+annihilation.T)/np.sqrt(2*omega)
    p = 1j*np.sqrt(omega/2)*(annihilation.T-annihilation)
    return [x[:d, :d] for x in (q, q@q, (p@p).real, q@q@q@q)]


def embed(dimensions, factors):
    matrix = np.array([[1.]])
    for site in reversed(range(len(dimensions))):
        matrix = np.kron(matrix, factors.get(site, np.eye(dimensions[site])))
    return matrix


def pair_matrix(dimensions, first, second, local):
    # Independent matrix-entry construction fixes the joint-probe indexing.
    occupations = list(product(*(range(d) for d in reversed(dimensions))))
    matrix = np.zeros((len(occupations),)*2, complex)
    for a, reverse_left in enumerate(occupations):
        left = tuple(reversed(reverse_left))
        for b, reverse_right in enumerate(occupations):
            right = tuple(reversed(reverse_right))
            if all(left[k] == right[k] for k in range(len(dimensions)) if k not in (first, second)):
                i = left[first]+dimensions[first]*left[second]
                j = right[first]+dimensions[first]*right[second]
                matrix[a, b] = local[i, j]
    return matrix


def projector(dimension, initial, other, sign):
    ket = np.zeros(dimension, complex)
    ket[initial], ket[other] = 1/np.sqrt(2), sign*1j/np.sqrt(2)
    return np.outer(ket, ket.conj())


@pytest.fixture(scope='module')
def field():
    dimensions, omega = (4, 4, 4), 1.7
    quadratic, quartic = [2.9, 3.4, 4.2], [.02, .06, .04]
    edges = {(0, 1): .6, (1, 2): .9}
    q, q2, p2, q4 = operators(4, omega)
    onsite = [.5*p2+.5*a*q2+g*q4 for a, g in zip(quadratic, quartic)]
    a = sum(embed(dimensions, {site: h}) for site, h in enumerate(onsite))
    b = sum(-weight*embed(dimensions, {i: q, j: q}) for (i, j), weight in edges.items())
    # A centered factor cap is sufficient: scalar phases cancel from every
    # projector law.  Take a little extra room above numerical eigennorms.
    width = lambda h: (np.linalg.eigvalsh(h)[-1]-np.linalg.eigvalsh(h)[0])/2
    cap = float(np.ceil(width(a)+width(b)))
    return dict(dimensions=dimensions, omega=omega, onsite=onsite,
                quadratic=quadratic, quartic=quartic, edges=edges, a=a, b=b,
                h=a+b, cap=cap)


def evolution(field, time):
    # One fixed symmetric whole-field step: both matrices act at every probe.
    half = expm(-.5j*time*field['a'])
    return half@expm(-1j*time*field['b'])@half


def blank(field, site=None, occupation=0):
    state = np.zeros(np.prod(field['dimensions']), complex)
    index = 0 if site is None else occupation*np.prod(field['dimensions'][:site], dtype=int)
    state[int(index)] = 1
    return state


def laws(state, effects):
    return np.array([[1-float(np.vdot(state, effect@state).real),
                     float(np.vdot(state, effect@state).real)] for effect in effects])


@pytest.mark.parametrize('site,occupation', list(product(range(3), range(2))))
def test_full_coupled_onsite_marginal_derivative_is_the_local_matrix_element(field, site, occupation):
    effects = [embed(field['dimensions'], {site: projector(4, occupation, occupation+2, sign)})
               for sign in (-1, 1)]
    initial = blank(field, site, occupation)
    observable = (effects[0]-effects[1])/2
    derivative = (1j*np.vdot(initial, (field['h']@observable-observable@field['h'])@initial)).real
    expected = field['onsite'][site][occupation+2, occupation]
    assert derivative == pytest.approx(expected, abs=2e-14)
    times = 1/512, 1/256
    records = np.array([laws(evolution(field, t)@initial, effects) for t in times])
    result = r.coefficient_from_records(records, *times, product_generator_norm_sum_cap=field['cap'])
    assert abs(result['matrix_element']-expected) <= result['radius']
    assert abs(result['matrix_element']-expected) < 1e-8
    assert 'generator_norm_bound' not in result


@pytest.mark.parametrize('pair', list(combinations(range(3), 2)))
def test_full_coupled_joint_marginal_recovers_only_direct_pair_at_first_order(field, pair):
    i, j = pair
    effects = [pair_matrix(field['dimensions'], i, j, projector(16, 0, 5, sign))
               for sign in (-1, 1)]
    initial = blank(field)
    observable = (effects[0]-effects[1])/2
    derivative = (1j*np.vdot(initial, (field['h']@observable-observable@field['h'])@initial)).real
    expected = -field['edges'].get(pair, 0)/(2*field['omega'])
    assert derivative == pytest.approx(expected, abs=2e-14)
    times = 1/512, 1/256
    records = np.array([laws(evolution(field, t)@initial, effects) for t in times])
    result = r.coefficient_from_records(records, *times, product_generator_norm_sum_cap=field['cap'])
    assert abs(result['matrix_element']-expected) <= result['radius']
    assert abs(result['matrix_element']-expected) < 1e-8


def test_marginalization_keeps_excited_spectators_instead_of_postselecting(field):
    site, time = 0, .21
    initial = blank(field)
    state = evolution(field, time)@initial
    local = projector(4, 0, 2, -1)
    marginal = np.vdot(state, embed(field['dimensions'], {site: local})@state).real
    zero = np.diag([1., 0., 0., 0.])
    postselected = np.vdot(state, embed(field['dimensions'], {0: local, 1: zero, 2: zero})@state).real
    assert marginal-postselected > 1e-4
    assert abs(np.vdot(state, state).real-1) < 3e-14


def test_indirect_path_affects_finite_records_without_creating_a_direct_edge(field):
    pair = (0, 2)
    effects = [pair_matrix(field['dimensions'], *pair, projector(16, 0, 5, sign))
               for sign in (-1, 1)]
    initial = blank(field)
    time = .04
    p = laws(evolution(field, time)@initial, effects)
    raw_slope = (p[0, 1]-p[1, 1])/(2*time)
    assert abs(raw_slope) > 1e-5
    small = np.array([laws(evolution(field, t)@initial, effects) for t in (time/2, time)])
    result = r.coefficient_from_records(small, time/2, time, product_generator_norm_sum_cap=field['cap'])
    assert abs(result['matrix_element']) < abs(raw_slope)/100
    assert abs(result['matrix_element']) < result['radius']


def test_whole_product_and_continuous_flow_are_different_calibration_curves(field):
    effects = [embed(field['dimensions'], {0: projector(4, 0, 2, sign)}) for sign in (-1, 1)]
    initial, time = blank(field), .2
    product_law = laws(evolution(field, time)@initial, effects)
    exact_law = laws(expm(-1j*time*field['h'])@initial, effects)
    assert np.max(abs(product_law-exact_law)) > 1e-5
    # Both recover the same generator at zero; their finite records need not
    # agree, and the new certificate never substitutes one for the other.
    for matrix in (evolution(field, 1e-6), expm(-1j*1e-6*field['h'])):
        p = laws(matrix@initial, effects)
        assert (p[0, 1]-p[1, 1])/(2e-6) == pytest.approx(field['onsite'][0][2, 0], abs=1e-8)


def test_real_initial_state_and_real_factors_make_both_probe_differences_odd(field):
    time = .13
    np.testing.assert_allclose(evolution(field, -time), evolution(field, time).conj(), atol=3e-15, rtol=0)
    effects = [pair_matrix(field['dimensions'], 0, 1, projector(16, 0, 5, sign)) for sign in (-1, 1)]
    positive = laws(evolution(field, time)@blank(field), effects)
    negative = laws(evolution(field, -time)@blank(field), effects)
    assert positive[0, 1]-positive[1, 1] == pytest.approx(-(negative[0, 1]-negative[1, 1]), abs=3e-15)


def linear_records(time):
    # Exact linear laws exercise the wrapper's inherited algebra separately
    # from the fully coupled finite-matrix experiments above.
    positions = np.array(list(product((0, 1), repeat=3)))
    permutation = np.array([7, 1, 5, 0, 3, 6, 2, 4])
    positions = positions[permutation]
    adjacency = (np.sum(abs(positions[:, None]-positions[None, :]), axis=-1) == 1).astype(int)
    g, omega = 1/24, 2
    def law(h):
        return [[.5-h*time, .5+h*time], [.5+h*time, .5-h*time]]
    onsite = np.array([[law(np.sqrt(2)*3*g/(2*omega**2)),
                       law(np.sqrt(6)*5*g/(2*omega**2))] for _ in range(8)])
    pairs = list(combinations(range(8), 2))
    bonds = np.array([law(-adjacency[i, j]/(2*omega)) for i, j in pairs])
    return dict(onsite_probabilities=onsite, bond_probabilities=bonds, pairs=pairs,
                fixed_symmetric_steps=1, complete_field_schedule_in_every_trial=True,
                independently_isolated_action_services=False,
                occupation_dimensions=(4,)*8), adjacency


def test_wrapper_recovers_parameters_and_declares_whole_schedule_assumptions():
    u, v = 1/16384, 1/8192
    first, adjacency = linear_records(u)
    second, _ = linear_records(v)
    clock = lambda t: [(1+np.cos(t))/2, (1-np.sin(t))/2]
    result = r.recover(first, second, clock(u), clock(v), omega=2,
                       product_generator_norm_sum_cap=96)
    assert result['quartic'] == pytest.approx(1/24, abs=1e-10)
    assert result['bare_mass_squared'] == pytest.approx(1., abs=1e-9)
    assert np.array_equal(result['adjacency'], adjacency)
    assert result['whole_field_schedule_is_supplied']
    assert not result['single_generator_isolation_used']
    assert result['spectator_outcomes_marginalized']
    assert result['calibration_bounds_compare_product_derivative_at_zero']
    assert not result['separate_calibration_splitting_error_required']
    assert 'isolated_generator_family_is_supplied' not in result
    for row in result['onsite_coefficient_details']:
        for item in row:
            assert item['product_generator_norm_sum_cap'] == 96.
            assert 'generator_norm_bound' not in item
    for item in result['bond_coefficient_details']:
        assert item['product_generator_norm_sum_cap'] == 96.


@pytest.mark.parametrize('steps', [None, True, 0, 2, 1.5])
def test_incompatible_or_missing_product_schedule_is_rejected(steps):
    u, v = 1/16384, 1/8192
    first, _ = linear_records(u)
    second, _ = linear_records(v)
    second['fixed_symmetric_steps'] = steps
    clock = lambda t: [(1+np.cos(t))/2, (1-np.sin(t))/2]
    with pytest.raises(ValueError, match='same positive fixed symmetric step count'):
        r.recover(first, second, clock(u), clock(v), omega=2,
                  product_generator_norm_sum_cap=96)


@pytest.mark.parametrize('key,value', [
    ('complete_field_schedule_in_every_trial', False),
    ('independently_isolated_action_services', True),
    ('occupation_dimensions', (4,)*7),
    ('occupation_dimensions', (8,)+(4,)*7),
    ('occupation_dimensions', None),
])
def test_changed_coupling_access_or_occupation_protocol_is_rejected(key, value):
    u, v = 1/16384, 1/8192
    first, _ = linear_records(u)
    second, _ = linear_records(v)
    second[key] = value
    clock = lambda t: [(1+np.cos(t))/2, (1-np.sin(t))/2]
    with pytest.raises(ValueError, match='both|whole-field'):
        r.recover(first, second, clock(u), clock(v), omega=2,
                  product_generator_norm_sum_cap=96)


def test_changing_step_count_breaks_the_richardson_remainder_despite_same_sum_generator():
    x, z = np.array([[0., 1.], [1., 0.]]), np.diag([1., -1.])
    a, b = .7*x+.2*z, .3*z
    initial = np.array([1., 0.])
    effects = [projector(2, 0, 1, sign) for sign in (-1, 1)]
    def record(time, steps):
        half = expm(-.5j*time*a/steps)
        step = half@expm(-1j*time*b/steps)@half
        return laws(np.linalg.matrix_power(step, steps)@initial, effects)
    # Both schedules have derivative H=a+b and sum of centered factor norms
    # below1.1. Their different cubic terms no longer cancel across times.
    time = .01
    compatible = r.coefficient_from_records(
        np.array([record(time/2, 1), record(time, 1)]), time/2, time,
        product_generator_norm_sum_cap=1.1)
    incompatible = r.coefficient_from_records(
        np.array([record(time/2, 1), record(time, 2)]), time/2, time,
        product_generator_norm_sum_cap=1.1)
    assert abs(compatible['matrix_element']-.7) < compatible['radius']
    assert abs(incompatible['matrix_element']-.7) > 100*incompatible['radius']


def test_new_inverse_never_imports_a_backend_or_source_action():
    tree = ast.parse(Path(r.__file__).read_text())
    imports = [node.module or '' for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    assert not any('source' in name or 'lattice' in name or 'execution' in name for name in imports)
