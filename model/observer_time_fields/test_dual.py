"""Original full binary source oracles for the sparse field dual."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product

import numpy as np
import pytest

from ..observer import native
from . import dual as d

PORTS = (0, 1, 2, 3)
EDGES = ((0, 1), (1, 2), (2, 3))
INCIDENT = ((0, 1),)


@lru_cache(maxsize=None)
def _native_terminal(initial, ticks, cap):
    states = {initial: F(1)}
    for _ in range(ticks):
        waiting, completed = states, {}
        for _ in range(cap):
            following = {}
            for state, mass in waiting.items():
                for seam in EDGES:
                    for coin in (0, 1):
                        output, _ = native.native(state, seam, coin)
                        target = completed if seam in INCIDENT else following
                        target[output] = target.get(output, F(0))+mass/6
            waiting = following
        states = completed
    return states


@pytest.mark.parametrize('degree', (1, 2))
@pytest.mark.parametrize('ticks,cap,bits', ((1, 1, 8), (1, 3, 12), (2, 3, 16)))
def test_every_monomial_coefficient_and_binary_input_against_original_native_law(degree, ticks, cap, bits):
    operator = d.build(PORTS, EDGES, INCIDENT, degree)
    labels, scale = d.states(operator), 1 << bits
    for selected in labels:
        initial = tuple(scale if label == selected else 0 for label in labels)
        lower = d.evolve(operator, initial, ticks, cap, bits)
        exact = []
        for occupied in labels:
            state = tuple(int(p in occupied) for p in PORTS)
            exact.append(sum((mass for final, mass in _native_terminal(state, ticks, cap).items()
                              if all(final[p] for p in selected)), F(0)))
        assert all(F(value, scale) <= target for value, target in zip(lower, exact, strict=True))
        finite_mass = (1-F(2, 3)**cap)**ticks
        assert sum(exact) == finite_mass
        deficit = finite_mass-F(sum(lower), scale)
        assert deficit == sum(target-F(value, scale) for target, value in zip(exact, lower, strict=True))
        assert 0 <= deficit < F(1, 100)
        for state in product((0, 1), repeat=4):
            actual = sum((F(value, scale) for label, value in zip(labels, lower, strict=True)
                          if all(state[p] for p in label)), F(0))
            expected = sum((mass for final, mass in _native_terminal(state, ticks, cap).items()
                            if all(final[p] for p in selected)), F(0))
            assert actual <= expected <= actual+deficit


def test_backward_factor_order_changes_nonstationary_off_root_coefficients():
    operator = d.build((0, 1, 2), ((0, 1), (1, 2)), ((0, 1),))
    correct = d.evolve(operator, (256, 0, 0), ticks=1, cap=2, bits=8)
    assert correct == (96, 80, 16)  # (I+A) B n0
    wrong_forward_order = (96, 96, 0)  # B (I+A) n0
    assert correct != wrong_forward_order
    assert sum(correct) == sum(wrong_forward_order) == 192
    assert correct[0] == wrong_forward_order[0]  # A scalar return check misses this.


def test_exclusion_fixed_pair_partial_mass_and_nonconsecutive_port_labels():
    operator = d.build((4, 9), ((4, 9),), ((4, 9),), degree=2)
    assert d.states(operator) == ((4, 9),)
    assert d.evolve(operator, (129,), 4, cap=5, bits=8) == (129,)
    # Both labels are occupied: swapping them leaves n4*n9 fixed. Two
    # independently hopping labels would introduce spurious collisions.
    outputs = [native.native((1, 0), (0, 1), coin)[0] for coin in (0, 1)]
    assert sum(state[0]*state[1] for state in outputs) == 0
    mean0, mean1 = (F(sum(state[p] for state in outputs), 2) for p in (0, 1))
    assert mean0*mean1 == F(1, 4)  # The product of means is not the pair field.
    singleton = d.build((4, 9), ((4, 9),), ((4, 9),), degree=1)
    assert d.evolve(singleton, (128, 0), 1, cap=1, bits=8) == (64, 64)
    assert d.evolve(singleton, (127, 0), 0, bits=8) == (127, 0)


def test_initial_partial_mass_is_not_normalized_or_reprepared_between_ticks():
    operator = d.build(PORTS, EDGES, INCIDENT, degree=2)
    initial = (129, 0, 0, 0, 0, 0)
    once = d.evolve(operator, initial, 1, cap=4, bits=8)
    twice = d.evolve(operator, once, 1, cap=4, bits=8)
    assert twice == d.evolve(operator, initial, 2, cap=4, bits=8)
    assert sum(twice) < sum(once) < sum(initial)


def test_immutable_public_operator_cannot_poison_cached_maps():
    operator = d.build(PORTS, EDGES, INCIDENT, degree=2)
    initial = (256, 0, 0, 0, 0, 0)
    expected = d.evolve(operator, initial, 1, bits=8)
    assert not hasattr(operator, '__dict__')
    with pytest.raises(AttributeError):
        operator.degree = 1
    with pytest.raises(TypeError):
        d.states(operator)[0][0] = 99
    assert d.evolve(d.build(PORTS, EDGES, INCIDENT, 2), initial, 1, bits=8) == expected


@pytest.mark.parametrize('field,value', (
    ('ports', (False, 1, 2, 3)), ('ports', [0, 1, 2, 3]),
    ('edges', ((1, 0), (1, 2), (2, 3))), ('edges', ((0, 1), (0, 1))),
    ('incident', ((0, 2),)), ('incident', ((0, True),)), ('degree', True), ('degree', F(2))))
def test_original_geometry_types_checked_before_any_cache_lookup(field, value):
    d.evolve(d.build(PORTS, EDGES, INCIDENT, 2), (256, 0, 0, 0, 0, 0), 1, bits=8)
    args = {'ports': PORTS, 'edges': EDGES, 'incident': INCIDENT, 'degree': 2}
    args[field] = value
    with pytest.raises(ValueError):
        d.build(**args)
    with pytest.raises(ValueError):
        d.states(d.Operator(**args))


@pytest.mark.parametrize('initial', ((True, 0, 0, 0), (np.int64(256), 0, 0, 0),
    (256.0, 0, 0, 0), (-1, 0, 0, 0), (257, 0, 0, 0), (1,), [256, 0, 0, 0]))
def test_original_coefficient_vector_not_silently_cast_or_normalized(initial):
    with pytest.raises(ValueError):
        d.evolve(d.build(PORTS, EDGES, INCIDENT), initial, 1, bits=8)


@pytest.mark.parametrize('params', ({'ticks': True}, {'ticks': 17}, {'cap': 0},
    {'cap': 2.0}, {'bits': True}, {'bits': 53}))
def test_native_count_and_integer_arithmetic_limits(params):
    parameters = {'ticks': 1, 'cap': 3, 'bits': 8}
    parameters.update(params)
    with pytest.raises(ValueError):
        d.evolve(d.build(PORTS, EDGES, INCIDENT), (256, 0, 0, 0), **parameters)


def test_overflow_is_refused_before_allocating_the_sparse_dual():
    ports = tuple(range(47))
    edges = tuple((a, b) for a in ports for b in ports if a < b)
    operator = d.build(ports, edges, edges, degree=2)
    with pytest.raises(ValueError, match='int64'):
        d.evolve(operator, (), 1, bits=52)


def test_subclass_dispatch_cannot_hide_geometry_or_coefficient_entries():
    class Trap(tuple):
        def __iter__(self):
            raise AssertionError('untrusted tuple subclass executed')

    class IntegerTrap(int):
        def __hash__(self):
            raise AssertionError('untrusted integer subclass hashed')

    with pytest.raises(ValueError):
        d.build((IntegerTrap(0), 1), ((0, 1),), ((0, 1),))
    with pytest.raises(ValueError):
        d.build(PORTS, Trap(EDGES), INCIDENT)
    with pytest.raises(ValueError):
        d.evolve(d.build(PORTS, EDGES, INCIDENT), Trap((256, 0, 0, 0)), 1, bits=8)


@pytest.fixture(scope='module')
def source_geometry():
    from ..native_action import incident_refinement
    return incident_refinement.geometry()[:3]


def test_actual_carrier_triangle_cover_supplies_full_binary_strict_transfer_bound(source_geometry):
    certificate = d.source_carrier_certificate(*source_geometry)
    assert certificate == {'root_ports': 12, 'internal_seams': 30, 'source_triangles': 20,
        'triangle_edge_multiplicity': 2, 'binary_triangle_states': 8,
        'triangle_count_minimum_eigenvalue': 3, 'triangle_count_maximum_eigenvalue': 6,
        'unnormalized_incident_transfer_lower': 15}
    # Returned metadata cannot poison a subsequent source calculation.
    certificate['unnormalized_incident_transfer_lower'] = 999
    assert d.source_carrier_certificate(*source_geometry)['unnormalized_incident_transfer_lower'] == 15


def test_equal_seam_and_face_censuses_do_not_replace_the_twofold_cover(source_geometry):
    ports, edges, incident = source_geometry
    inside = {edge for edge in edges if max(edge) < 12}
    missing = set(combinations(range(12), 2))-inside
    for removed, added in product(sorted(inside), sorted(missing)):
        changed = inside-{removed} | {added}
        faces = [face for face in combinations(range(12), 3)
                 if set(combinations(face, 2)) <= changed]
        if len(faces) == 20:
            break
    else:
        raise AssertionError('the independent count-preserving wrong graph witness disappeared')
    wrong_edges = tuple(sorted(set(edges)-{removed} | {added}))
    wrong_incident = tuple(sorted(set(incident)-{removed} | {added}))
    assert len(wrong_edges) == len(edges) and len(faces) == 20
    with pytest.raises(ValueError, match='exactly two'):
        d.source_carrier_certificate(ports, wrong_edges, wrong_incident)


def test_root_seams_must_act_on_the_incident_clock(source_geometry):
    ports, edges, incident = source_geometry
    with pytest.raises(ValueError, match='complete incident root carrier'):
        d.source_carrier_certificate(ports, edges, tuple(e for e in incident if e != (0, 1)))
    with pytest.raises(ValueError, match='original nonnegative integer'):
        d.source_carrier_certificate((False, *ports[1:]), edges, incident)


def test_trace_preserving_deterministic_swap_is_not_the_native_triangle(monkeypatch, source_geometry):
    def wrong_native(state, seam, coin):
        final = list(state)
        a, b = seam
        final[a], final[b] = final[b], final[a]
        return tuple(final), (coin, state[a]-final[a])

    # A permutation channel is positive and mass preserving, but dropping
    # the native identity branch invalidates the positive triangle bound.
    monkeypatch.setattr(native, 'native', wrong_native)
    with pytest.raises(ArithmeticError, match='identity plus swap'):
        d.source_carrier_certificate(*source_geometry)
