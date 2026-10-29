"""Independent full-CAR, source-word and complete-reference controls."""
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest
from scipy.linalg import expm

from .field import FieldGeometry, FullFock, evolve, number_pointer, resources, source_word
from .bounds import complete_history_bound, stability_bounds
from ..interacting_field.lattice import TwoFermions
from ..interacting_field import source
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..observer.source_scalar import compiler


def random_columns(space, columns=3, seed=41):
    rng = np.random.default_rng(seed)
    state = rng.normal(size=(space.dimension, columns))+1j*rng.normal(size=(space.dimension, columns))
    return state/np.linalg.norm(state)


def test_full_car_hamiltonian_preserves_all_number_sectors_and_vacuum_filled_limits():
    geometry = FieldGeometry((3, 1, 1))
    space = FullFock(geometry.modes)
    a = [space.annihilator(j).toarray() for j in range(space.modes)]
    for i, j in combinations(range(space.modes), 2):
        assert np.array_equal(a[i]@a[j]+a[j]@a[i], np.zeros((space.dimension,)*2))
        assert np.array_equal(a[i]@a[j].conj().T+a[j].conj().T@a[i], np.zeros((space.dimension,)*2))
    assert np.array_equal(a[0]@a[0].conj().T+a[0].conj().T@a[0], np.eye(space.dimension))
    h = space.hamiltonian(geometry).toarray()
    number = np.diag([i.bit_count() for i in range(space.dimension)])
    assert np.array_equal(h@number, number@h)
    assert np.array_equal(h, h.conj().T)
    assert np.array_equal(h[:, 0], np.zeros(space.dimension))
    assert h[-1, -1] == geometry.volume*geometry.interaction
    assert np.count_nonzero(h[:, -1]) == 1


@pytest.mark.parametrize('p,q,z', [(0, 5, .25+.125j), (1, 4, -.125j), (2, 3, -.25)])
def test_actual_routed_source_gate_matches_independent_all_sector_car_generator(p, q, z):
    space = FullFock(6)
    a, b = space.annihilator(p), space.annihilator(q)
    term = z*(a.conj().T@b)
    expected = expm(-.3j*(term+term.conj().T).toarray())
    raw = gaussian.local_gate('hopping', z, .3).fock_unitary
    actual = space.routed(np.eye(space.dimension), p, q, raw)
    assert np.linalg.norm(actual-expected, 2) < 2e-10
    # An ordinary SWAP loses the occupied-intermediate-mode CAR sign.
    ordinary = np.array([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,1]], complex)
    if q-p > 1:
        wrong = np.eye(space.dimension, dtype=complex)
        route = tuple(range(q-1, p, -1))
        for left in route: wrong = space.adjacent(wrong, left, ordinary)
        wrong = space.adjacent(wrong, p, raw)
        for left in reversed(route): wrong = space.adjacent(wrong, left, ordinary)
        assert np.linalg.norm(wrong-expected, 2) > .01


def test_source_schedule_matches_independent_full_fock_evolution_with_reference():
    geometry = FieldGeometry((2, 2, 1))
    space = FullFock(geometry.modes)
    initial = random_columns(space)
    actual = evolve(geometry, initial, -.125, 4, space=space)
    target = space.exact_evolve(geometry, initial, -.125)
    bound = complete_history_bound(geometry.sides, ((F(-1, 8), 4),))
    assert np.linalg.norm(actual-target) < float(bound)
    assert np.linalg.norm(actual.conj().T@actual-initial.conj().T@initial) < 2e-10
    # Three-, four-, vacuum and filled sectors are all present in this one input.
    norms = {n: np.linalg.norm(actual[[i for i in range(space.dimension) if i.bit_count()==n]])
             for n in (0, 3, 4, space.modes)}
    assert all(norm > 0 for norm in norms.values())


def test_complete_three_dimensional_box_has_all_65536_fock_states_and_source_history():
    geometry = FieldGeometry((2, 2, 2))
    space = FullFock(geometry.modes)
    assert space.dimension == 65536
    initial = random_columns(space, columns=2, seed=91)
    actual = evolve(geometry, initial, 1/32, 2, space=space)
    target = space.exact_evolve(geometry, initial, 1/32)
    bound = complete_history_bound(geometry.sides, ((F(1,32), 2),))
    assert np.linalg.norm(actual-target) < float(bound)
    assert np.linalg.norm(actual.conj().T@actual-initial.conj().T@initial) < 2e-9
    outputs = number_pointer(space, actual, 5)
    assert abs(sum(np.vdot(out,out).real for out in outputs)-1) < 2e-9
    assert np.linalg.norm(sum(out.conj().T@out for out in outputs)-initial.conj().T@initial) < 2e-9


def test_full_fock_source_execution_restricts_to_existing_two_particle_schedule():
    geometry = FieldGeometry((2, 2, 1))
    space = FullFock(geometry.modes)
    old = TwoFermions(geometry.modes)
    old.groups, old.volume, old.interaction = geometry.groups, geometry.volume, geometry.interaction
    old.contact_indices = np.array([old.index(2*x, 2*x+1) for x in range(geometry.volume)])
    rng = np.random.default_rng(59)
    initial = rng.normal(size=(old.dimension, 2))+1j*rng.normal(size=(old.dimension, 2))
    initial /= np.linalg.norm(initial)
    indices = [int((1 << p)+(1 << q)) for p,q in zip(old.first,old.second,strict=True)]
    lifted = np.zeros((space.dimension, 2), complex); lifted[indices] = initial
    actual = evolve(geometry, lifted, .125, 2, space=space)
    target = source.evolve(old, initial, .125, 2)
    assert np.linalg.norm(actual[indices]-target) < 2e-10
    outside = actual.copy(); outside[indices] = 0
    assert np.linalg.norm(outside) < 2e-10
    small_h = old.lift_one_body(geometry.single_h()).toarray()
    small_h[old.contact_indices,old.contact_indices] += geometry.interaction
    h = space.hamiltonian(geometry).toarray()
    assert np.array_equal(h[np.ix_(indices,indices)], small_h)


def test_literal_complete_source_tape_agrees_with_local_contraction_and_resource_count():
    geometry = FieldGeometry((2, 1, 1))
    space = FullFock(geometry.modes)
    tape = source_word(geometry, .125, 1)
    whole = compiler.execute(tape, space.modes)
    actual = evolve(geometry, np.eye(space.dimension), .125, 1, space=space)
    assert np.linalg.norm(whole-actual, 2) < 2e-10
    measured, counted = compiler.resources(tape), resources(geometry, .125, 1)
    assert all(measured[key] == counted[key] for key in source.COUNT_KEYS)
    with pytest.raises(ValueError, match='budget'):
        source_word(geometry, .125, 1, max_pulses=len(tape)-1)


def test_source_pointer_retains_both_outputs_and_reference_across_number_coherences():
    space = FullFock(6)
    initial = random_columns(space)
    outputs = number_pointer(space, initial, 4)
    reference = sum(branch.conj().T@branch for branch in outputs)
    assert len(outputs) == 2
    assert np.linalg.norm(reference-initial.conj().T@initial) < 2e-10
    assert abs(sum(np.vdot(branch,branch).real for branch in outputs)-1) < 2e-10
    expected = sum(np.linalg.norm(initial[i])**2 for i in range(space.dimension) if (i >> 4)&1)
    assert abs(np.vdot(outputs[1],outputs[1]).real-expected) < 2e-10


def test_extensive_stability_and_gibbs_faithfulness_on_full_interacting_carrier():
    geometry = FieldGeometry((2, 2, 1))
    space = FullFock(geometry.modes)
    bounds = stability_bounds(geometry.sides)
    eigenvalues = np.linalg.eigvalsh(space.hamiltonian(geometry).toarray())
    assert eigenvalues[0] >= float(bounds['energy_lower'])-1e-12
    assert eigenvalues[-1] <= float(bounds['energy_upper'])+1e-12
    weights = np.exp(-1.5*eigenvalues); weights /= weights.sum()
    assert weights.min() > 0 and abs(weights.sum()-1) < 1e-14
    assert bounds['generator_norm_sum_upper'] == F(9, 2)


def test_full_fock_bound_charges_volume_and_composes_without_record_multiplier():
    a = complete_history_bound((2, 2, 1), ((F(1, 8), 4),))
    assert complete_history_bound((2, 2, 1), ((F(1, 8), 4),(F(-1, 8),4))) == 2*a
    assert complete_history_bound((3, 3, 3), ((F(1,8),4),)) > a
    assert complete_history_bound((1, 1, 1), ((0,1),)) == 0
    assert complete_history_bound((128,128,128), ((1,1),)) == 1


def test_geometry_and_counts_do_not_require_global_fock_allocation():
    geometry = FieldGeometry((3, 3, 3))
    count = resources(geometry, .1, 2)
    assert geometry.modes == 54 and count['native_pulses'] > 0
    assert count['number_sector_projection_applied'] is False
    with pytest.raises(ValueError, match='explicit full-Fock'):
        FullFock(geometry.modes)


@pytest.mark.parametrize('bad', [True, 0, -1, 1.5, float('nan')])
def test_invalid_split_counts_rejected(bad):
    geometry = FieldGeometry((1,1,1)); space = FullFock(2)
    with pytest.raises(ValueError): evolve(geometry, np.ones(4), .1, bad, space=space)
    with pytest.raises(ValueError): resources(geometry, .1, bad)


def test_invalid_domains_and_inputs_fail_before_execution():
    for sides in [(True,1,1), (0,1,1), (1,2), [1,1,1], (129,1,1)]:
        with pytest.raises(ValueError): FieldGeometry(sides)
    with pytest.raises(ValueError): FieldGeometry((1,1,1), interaction=-1)
    space = FullFock(4)
    for state in [np.full(16, np.nan), np.ones(15), np.ma.array(np.ones(16)), np.zeros((16,0)), np.zeros((16,2**17))]:
        with pytest.raises(ValueError): space.state(state)
    with pytest.raises(ValueError): space.adjacent(np.ones(16), 0, np.ma.array(np.eye(4)))
    with pytest.raises(ValueError): space.routed(np.ones(16), 3, 1, np.eye(4))
    with pytest.raises(ValueError): number_pointer(space, np.ones(16), True)
    with pytest.raises(ValueError): complete_history_bound((1,1,1), ((.1,2),))
    with pytest.raises(ValueError): stability_bounds((1,1,1), speed=True)


def test_contact_derivation_is_genuinely_nonquadratic_on_even_local_observables():
    interacting = FieldGeometry((2,1,1), interaction=1)
    free = FieldGeometry((2,1,1), interaction=0)
    space = FullFock(4)
    delta = space.hamiltonian(interacting)-space.hamiltonian(free)
    hopping = space.annihilator(0).conj().T@space.annihilator(2)
    expected = hopping@(space.number(1)-space.number(3))
    commutator = delta@hopping-hopping@delta
    assert np.array_equal(commutator.toarray(), expected.toarray())
    assert np.linalg.norm(commutator.toarray(), 2) == 1
    # The same occupied transport has different response to a spectator.
    assert commutator[1,4] == 0
    assert abs(commutator[3,6]) == 1


def test_local_boundary_enclosure_is_exact_nontrivial_and_improves_with_padding():
    from .bounds import exponential_upper, local_boundary_bound
    for x in [F(0),F(1,4),F(1),F(3),F(16)]:
        assert float(exponential_upper(x)) >= np.exp(float(x))*(1-1e-14)
    values = [local_boundary_bound(r, 0, 1, F(1,8)) for r in range(2,21)]
    assert all(a >= b for a,b in zip(values,values[1:]))
    assert values[-1] < F(1,1000)
    assert local_boundary_bound(3,0,1,F(-1,8)) == local_boundary_bound(3,0,1,F(1,8))
    assert local_boundary_bound(3,0,1,0) == 0
    assert local_boundary_bound(3,0,1,1,speed=0) == 0
    with pytest.raises(ValueError): local_boundary_bound(2,2,1,F(1,8))
    with pytest.raises(ValueError): local_boundary_bound(2,0,2,F(1,8))
    with pytest.raises(ValueError): local_boundary_bound(4097,0,1,F(1,8))
    with pytest.raises(ValueError): exponential_upper(F(-1))


@pytest.mark.parametrize('interaction', [0., 1., 40.])
def test_locality_and_boundary_comparison_are_uniform_in_onsite_interaction(interaction):
    geometry = FieldGeometry((3,1,1), interaction=interaction)
    space = FullFock(6)
    h = space.hamiltonian(geometry).toarray()
    a = [space.annihilator(j).toarray() for j in range(6)]
    pair = a[0].conj().T@a[1].conj().T
    observable = pair+pair.conj().T  # Even, norm one, and does not commute with onsite U.
    distant = space.number(4).toarray()
    time = .125
    evolution = expm(-1j*time*h)
    evolved = evolution.conj().T@observable@evolution
    lr = .25*np.expm1(48*geometry.speed*time)  # |X|=1, distance=2.
    assert np.linalg.norm(evolved@distant-distant@evolved, 2) <= lr
    boundary = np.zeros_like(h)
    for group in geometry.groups:
        for p,q,z in group:
            if p//2 == 1 and q//2 == 2:
                term = z*a[p].conj().T@a[q]
                boundary += term+term.conj().T
    smaller = expm(-1j*time*(h-boundary))
    restricted = smaller.conj().T@observable@smaller
    difference = np.linalg.norm(evolved-restricted, 2)
    cauchy = geometry.speed*time*.5*np.expm1(48*geometry.speed*time)
    assert 1e-8 < difference <= cauchy
