"""Independent full-CAR controls for local split and history guarantees."""
from fractions import Fraction as F
from itertools import product
from math import factorial

import numpy as np
import pytest
from scipy.linalg import expm

from .bounds import complete_history_bound, local_boundary_bound, stability_bounds
from .field import FieldGeometry, FullFock, evolve, resources
from .local_simulation import (extensive_split_coefficient, factorial_tail_upper,
                              history_split_bound, local_observer_schedule,
                              local_tail_boundary_bound, overlap_split_certificate,
                              source_schedule_counts, split_bound)
from ..common_source.chiral_candidate import gaussian_source
from ..interacting_field import source
from ..observer.source_scalar import compiler


def _commutator(a, b):
    return a@b-b@a


def _terms(geometry, space):
    # Independent complete CAR matrices; no target matrix from the planner.
    ladders = [space.annihilator(j).toarray() for j in range(space.modes)]
    terms = []
    for group in geometry.groups:
        for p, q, z in group:
            value = z*ladders[p].conj().T@ladders[q]
            terms.append(value+value.conj().T)
    terms.extend(geometry.interaction*(space.number(2*x)@space.number(2*x+1)).toarray()
                 for x in range(geometry.volume))
    return terms


@pytest.mark.parametrize('sides', [(2,1,1), (3,1,1), (2,2,1)])
def test_overlap_certificate_bounds_actual_full_fock_nested_commutators(sides):
    geometry = FieldGeometry(sides)
    space = FullFock(geometry.modes)
    terms = _terms(geometry, space)
    actual = 0.
    for i, a in enumerate(terms[:-1]):
        b = sum(terms[i+1:])
        actual += np.linalg.norm(_commutator(a, _commutator(a, b)), 2)/24
        actual += np.linalg.norm(_commutator(b, _commutator(b, a)), 2)/12
    counted = overlap_split_certificate(sides)
    assert 0 < actual <= float(counted['coefficient'])
    assert counted['coefficient'] <= counted['degree_coefficient']
    assert np.allclose(sum(terms), space.hamiltonian(geometry).toarray(), atol=0, rtol=0)


@pytest.mark.parametrize('duration', [F(1,4), F(-1,4)])
def test_tighter_bound_covers_actual_complete_source_unitary_and_number_coherences(duration):
    geometry = FieldGeometry((2,2,1))
    space = FullFock(geometry.modes)
    target = expm(-1j*float(duration)*space.hamiltonian(geometry).toarray())
    rng = np.random.default_rng(391)
    unknown = rng.normal(size=(space.dimension,3))+1j*rng.normal(size=(space.dimension,3))
    unknown /= np.linalg.norm(unknown)
    for steps in (1,2):
        actual = evolve(geometry, np.eye(space.dimension), float(duration), steps)
        bound = split_bound(geometry.sides, duration, steps, refined=True)
        assert np.linalg.norm(actual-target,2) < float(bound)
        assert np.linalg.norm(actual@unknown-target@unknown) < float(bound)
        assert np.linalg.norm(actual.conj().T@actual-np.eye(space.dimension),2) < 1e-9
        # Source execution is on the complete carrier, not merely N=2.
        for n in range(space.modes+1):
            occupied = np.array([i.bit_count()==n for i in range(space.dimension)])
            assert np.linalg.norm(actual[np.ix_(~occupied,occupied)]) < 1e-9


def test_volume_linear_coefficient_removes_the_cubic_volume_penalty():
    for sides in [(2,2,2), (3,3,3), (16,16,16), (10000,10000,10000)]:
        coefficient = extensive_split_coefficient(sides)
        old = stability_bounds(sides)['generator_norm_sum_upper']**3/3
        assert 0 < coefficient < old
        assert coefficient <= F(1463,1024)*np.prod(sides, dtype=object)
    assert overlap_split_certificate((2,2,2))['coefficient'] == F(629,256)
    assert extensive_split_coefficient((2,2,2)) == F(1463,256)
    assert overlap_split_certificate((3,3,3))['coefficient'] == F(49819,4096)
    assert split_bound((2,2,2), F(1,32), 2, refined=True) < (
        complete_history_bound((2,2,2), ((F(1,32),2),))/100)
    assert extensive_split_coefficient((1,1,1)) == 0
    assert extensive_split_coefficient((8,7,6), speed=0, interaction=100) == 0


def test_three_dimensional_complete_fock_execution_meets_the_reduced_budget():
    geometry = FieldGeometry((2,2,2)); space = FullFock(16)
    rng = np.random.default_rng(331)
    initial = rng.normal(size=(65536,2))+1j*rng.normal(size=(65536,2))
    initial /= np.linalg.norm(initial)
    duration = F(1,32)
    actual = evolve(geometry,initial,float(duration),1)
    target = space.exact_evolve(geometry,initial,float(duration))
    bound = split_bound(geometry.sides,duration,1,refined=True)
    assert bound == F(629,8388608)
    assert np.linalg.norm(actual-target) < float(bound)
    assert np.linalg.norm(actual.conj().T@actual-initial.conj().T@initial) < 1e-9


def test_ordered_overlap_bounds_are_exact_and_dominate_many_coefficients():
    for sides in product(range(1,4), repeat=3):
        for c,u in ((F(0),F(9)), (F(1,8),F(0)), (F(3,2),F(7))):
            result = overlap_split_certificate(sides,c,u)
            assert type(result['coefficient']) is F
            assert result['coefficient'] <= extensive_split_coefficient(sides,c,u)
    with pytest.raises(ValueError):
        overlap_split_certificate((17,17,17))


def test_closed_source_counts_match_the_actual_templates_and_routed_counter():
    for kind in ('hopping','fswap'):
        word = gaussian_source.local_word(kind, .25j, -.125) if kind=='hopping' else gaussian_source.local_word(kind)
        count = compiler.resources(word)
        assert (count['single_code_pulses'],count['pair_code_pulses']) == (22,4)
    for interaction in (0,1,7):
        word = source.contact(.125,interaction)[1]['word']
        count = compiler.resources(word)
        assert (count['single_code_pulses'],count['pair_code_pulses']) == (8,2)
    for sides in ((1,1,1),(2,1,1),(1,3,2),(3,2,2)):
        for speed in (F(),F(1,8)):
            geometry = FieldGeometry(sides,speed=float(speed))
            expected = resources(geometry,-.125,3)
            actual = source_schedule_counts(sides,3,speed)
            for key in (*source.COUNT_KEYS,'routing_FSWAP_calls',
                        'completed_hopping_blocks','completed_contact_blocks'):
                assert actual[key] == expected[key]


def test_exact_factorial_tail_enclosure_and_integrated_boundary_improvement():
    for value in (F(),F(1,8),F(3,8),F(1),F(3),F(32)):
        for start in (1,2,5,17,64):
            bound = factorial_tail_upper(value,start)
            lower = sum((value**n/F(factorial(n)) for n in range(start,start+100)),F())
            assert type(bound) is F and lower <= bound
            assert bound <= lower+max(F(1,10**30),lower/F(3))
    for radius in range(1,15):
        bound = local_tail_boundary_bound(radius,0,1,F(1,8))
        assert bound <= local_boundary_bound(radius,0,1,F(1,8))
        assert bound == local_tail_boundary_bound(radius,0,1,F(-1,8))
    assert local_tail_boundary_bound(4,0,1,F(1,8)) < F(1,500)
    assert local_tail_boundary_bound(4,0,1,0) == 0


def test_local_schedule_has_exact_error_budget_and_no_ambient_volume_parameter():
    plan = local_observer_schedule((F(1,8),),F(1,100),refined=True)
    stage, = plan['stages']
    assert (stage['input_radius'],stage['output_radius'],stage['steps']) == (4,0,14)
    assert plan['native_pulses'] == 340803288
    assert plan['half_diamond_bound'] <= F(1,100)
    assert not plan['ambient_volume_parameter_required']
    assert not plan['initial_restricted_state_prepared']
    assert stage['resources']['geometry_shape_admitted']
    assert not stage['resources']['explicit_state_geometry_admitted']
    assert not stage['resources']['global_Fock_array_allocated']
    # Minimal n for this chosen split coefficient and error allocation only.
    k,t,n = stage['split_coefficient'],abs(stage['duration']),stage['steps']
    assert k*t**3/n**2 <= F(1,200) < k*t**3/(n-1)**2


def test_nested_history_uses_larger_earlier_support_and_retains_signed_times():
    times = (F(1,8),F(-1,16))
    plan = local_observer_schedule(times,F(1,100),refined=True)
    first,last = plan['stages']
    assert tuple(s['duration'] for s in plan['stages']) == times
    assert first['input_radius'] > first['output_radius'] == last['input_radius'] > last['output_radius']
    assert first['tested_output_support_sites'] == (2*last['input_radius']+1)**3
    assert not first['resources']['geometry_shape_admitted']
    assert first['coefficient_method'] == 'local degree'
    assert plan['half_diamond_bound'] <= F(1,100)
    assert history_split_bound((2,1,1), ((F(1,4),2),(F(-1,4),2)),refined=True) == (
        2*split_bound((2,1,1),F(1,4),2,refined=True))


def _restrict_purification(value, kept_dimension):
    # Discarded high modes remain explicit environment columns. Retain the
    # final axis as the original untouched reference, not a classical sample.
    old,env,reference = value.shape
    return value.reshape(old//kept_dimension,kept_dimension,env,reference).transpose(
        1,0,2,3).reshape(kept_dimension,-1,reference)


def _joint_density(value):
    field,_,reference = value.shape
    purification = value.transpose(0,2,1).reshape(field*reference,-1)
    return purification@purification.conj().T


def test_complete_local_history_with_intermediate_records_and_entangled_reference():
    # One-sided chain sections of the same open-box action allow a complete
    # small-matrix control. The cubic boundary-count bounds remain upper
    # bounds: deleting transverse bonds only reduces degree and crossings.
    ambient = FieldGeometry((4,1,1)); full = FullFock(8)
    rng = np.random.default_rng(719)
    initial = rng.normal(size=(256,1,2))+1j*rng.normal(size=(256,1,2))
    initial /= np.linalg.norm(initial)
    time,steps = F(1,32),2
    target = full.exact_evolve(ambient,initial.reshape(256,2),float(time)).reshape(256,1,2)
    target_branches = []
    approximate_branches = []
    restricted = _restrict_purification(initial,64)
    three = FieldGeometry((3,1,1)); two = FieldGeometry((2,1,1))
    first = evolve(three,restricted.reshape(64,-1),float(time),steps).reshape(restricted.shape)
    first = _restrict_purification(first,16)
    for flag in (0,1):
        # Same complete, even, number-preserving local pointer on both paths.
        ideal_branch = target.copy()
        ideal_branch[[i for i in range(256) if (i&1)!=flag]] = 0
        ideal_branch = full.exact_evolve(ambient,ideal_branch.reshape(256,-1),float(-time)).reshape(target.shape)
        target_branches.append(_restrict_purification(ideal_branch,4))
        branch = first.copy()
        branch[[i for i in range(16) if (i&1)!=flag]] = 0
        branch = evolve(two,branch.reshape(16,-1),float(-time),steps).reshape(first.shape)
        approximate_branches.append(_restrict_purification(branch,4))
    actual = sum(np.linalg.svd(_joint_density(a)-_joint_density(b),compute_uv=False).sum()/2
                 for a,b in zip(target_branches,approximate_branches,strict=True))
    bound = (local_tail_boundary_bound(2,1,2,time)/2+
             local_tail_boundary_bound(1,0,1,-time)/2+
             split_bound(three.sides,time,steps,refined=True)+
             split_bound(two.sides,-time,steps,refined=True))
    assert 0 < actual < float(bound) < .05
    assert abs(sum(np.vdot(b,b).real for b in approximate_branches)-1) < 1e-9
    target_reference = sum(np.einsum('fer,fes->rs',b.conj(),b) for b in target_branches)
    approximate_reference = sum(np.einsum('fer,fes->rs',b.conj(),b) for b in approximate_branches)
    assert np.linalg.norm(target_reference-approximate_reference) < 1e-9


@pytest.mark.parametrize('bad',[True,0,-1,1.5,float('nan')])
def test_invalid_steps_fail_closed(bad):
    with pytest.raises(ValueError): split_bound((2,1,1),F(1,8),bad)
    with pytest.raises(ValueError): source_schedule_counts((2,1,1),bad)


def test_malformed_domains_and_infeasible_budgets_do_not_produce_a_schedule():
    for epsilon in (True,0,-1,1.1,float('nan')):
        with pytest.raises(ValueError): local_observer_schedule((F(1,8),),epsilon)
    for durations in ([],[.1],[True],None):
        with pytest.raises(ValueError): local_observer_schedule(durations,F(1,10))
    with pytest.raises(ValueError): local_observer_schedule((F(1,8),),F(1,100),max_radius=1)
    with pytest.raises(ValueError): local_observer_schedule((12,),F(1,10))
    with pytest.raises(ValueError): local_observer_schedule((1,),F(1,10),support_sites=2)
    with pytest.raises(ValueError): local_observer_schedule((1,),F(1,10),refined=1)
    with pytest.raises(ValueError): local_tail_boundary_bound(1,1,1,F(1,8))
    with pytest.raises(ValueError): factorial_tail_upper(F(33),2)
    with pytest.raises(ValueError): factorial_tail_upper(F(1),True)
    with pytest.raises(ValueError): history_split_bound((1,1,1),[])
    with pytest.raises(ValueError): split_bound((1,1,1),F(1,8),1,refined=0)
    for speed in (True,-1,.125):
        with pytest.raises(ValueError): extensive_split_coefficient((2,1,1),speed)


def test_zero_duration_and_zero_hopping_do_not_require_a_spatial_buffer():
    for times,c in (((0,0),F(1,8)),((F(1,8),F(-1,4)),F())):
        plan = local_observer_schedule(times,F(1,100),speed=c,max_radius=0)
        assert plan['initial_restriction_radius'] == 0
        assert plan['half_diamond_bound'] == 0
        assert all(s['steps']==1 for s in plan['stages'])
