"""Independent current, sum-rule and gauge controls for both charged candidates."""
import math
import json
import numpy as np
import pytest

from model.full_fock_field.field import FieldGeometry
from model.gauged_field.links import GaugedField
from model.particle_response import dimer
from . import electromagnetism as em


@pytest.mark.parametrize('axis', [0, 1, 2])
def test_full_Fock_lift_matches_independent_existing_two_fermion_law(axis):
    sides = [1, 1, 1]; sides[axis] = 2
    geometry = FieldGeometry(tuple(sides))
    phase = .371
    full = em.operators(geometry, phase)
    two = GaugedField(geometry.sites, angles={(0, 1):phase})
    # Existing two-body basis is lexicographic pairs, not sorted Fock masks.
    masks = [(1 << int(p)) | (1 << int(q)) for p,q in zip(two.first, two.second, strict=True)]
    np.testing.assert_allclose(full['hamiltonian'][np.ix_(masks,masks)],
                               two.hamiltonian.toarray(), atol=1e-15, rtol=0)
    np.testing.assert_allclose(full['current'][np.ix_(masks,masks)],
                               two.current((0,1)).toarray(), atol=1e-15, rtol=0)
    h = 1e-5
    derivative = (em.operators(geometry, phase+h)['hamiltonian']-
                  em.operators(geometry, phase-h)['hamiltonian'])/(2*h)
    np.testing.assert_allclose(full['current'],derivative,atol=2e-11,rtol=0)


@pytest.mark.parametrize('axis', [0, 1, 2])
def test_same_candidates_have_nonzero_current_and_complete_optical_sums(axis):
    sides = [1, 1, 1]; sides[axis] = 2
    report = em.benchmark(FieldGeometry(tuple(sides)))
    assert max(report['operator_residuals'].values()) < 2e-12
    for row in report['candidate_responses']:
        assert row['final_states_retained'] == 16
        assert len(row['transitions']) == 16
        assert {r['final_number'] for r in row['transitions']} == set(range(5))
        # Independent analytic two-level charge-sector result: D offdiagonal=1/2,
        # transition frequency=2t=1/8, J=i[K,D] offdiagonal magnitude=t=1/16.
        assert row['dipole_total_weight'] == pytest.approx(.25,abs=1e-13)
        assert row['current_total_weight'] == pytest.approx(1/256,abs=1e-13)
        assert row['dipole_first_frequency_moment'] == pytest.approx(1/32,abs=1e-13)
        visible = [r for r in row['transitions'] if r['dipole_weight'] > 1e-12]
        assert all(r['final_number'] == row['charge_in_supplied_units_e']+2 for r in visible)
        assert all(abs(r['model_frequency']-.125)<1e-12 for r in visible)
        assert abs(row['current_expectation']) < 1e-13
    assert report['open_graph_cycle_rank'] == 0
    assert report['magnetic_moment_or_gyromagnetic_ratio'] is None
    assert report['QED_identification'].startswith('NOT_ESTABLISHED')


def test_current_must_have_the_correct_sign_and_cannot_be_deleted():
    spectrum = dimer.setup()
    projector = dimer.candidate_projectors(spectrum)['charge_minus_one_hole']
    ops = em.operators()
    # Optical weights detect deletion. The full Ward identity separately
    # fixes current orientation, which squared transition strengths cannot do.
    with pytest.raises(ArithmeticError):
        em.candidate_response(spectrum,projector,ops['dipole'],np.zeros((16,16)),charge=-1)
    wrong = -ops['current']
    with pytest.raises(ArithmeticError):
        em.candidate_response(spectrum,projector,ops['dipole'],wrong,charge=-1)
    residual = np.linalg.norm(1j*(ops['generator']@ops['dipole']-
                                 ops['dipole']@ops['generator'])-wrong,2)
    assert residual > .1


def test_pure_gauge_link_does_not_supply_a_magnetic_moment():
    base = em.operators()
    for phase in (0.,.2,math.pi,3.8):
        changed = em.operators(phase=phase)
        np.testing.assert_allclose(np.linalg.eigvalsh(changed['hamiltonian']),
                                   np.linalg.eigvalsh(base['hamiltonian']),atol=2e-15,rtol=0)
    # Spectral invariance is not absence of electric response.
    assert np.linalg.norm(base['current'],2) > .1


def test_actual_contact_kernel_is_constant_but_lattice_Coulomb_is_not():
    report = em.bare_kernel_comparison()
    np.testing.assert_array_equal(report['contact_fourier_kernel'],[1.,1.,1.])
    assert report['Coulomb_normalized_shape'][0] == 1.
    assert report['Coulomb_normalized_shape'][1] == pytest.approx(1-1/math.sqrt(2))
    assert report['Coulomb_normalized_shape'][2] == pytest.approx(.5)
    assert report['maximum_normalized_shape_difference'] > .7
    # Independently Fourier transform the declared position-space contact.
    potential = np.zeros((8,8,8));potential[0,0,0]=1
    transformed = np.fft.fftn(potential)
    np.testing.assert_array_equal(transformed,np.ones_like(transformed))
    ops = em.operators()
    expected = np.diag([int(mask&3==3)+int(mask&12==12) for mask in range(16)])
    np.testing.assert_array_equal(ops['contact'],expected)
    assert not report['dressed_or_emergent_photon_no_go']
    assert em.bare_kernel_comparison(0)['contact_normalized_shape'] is None


@pytest.mark.parametrize('value',[float('nan'),float('inf'),True,np.ma.array(.2)])
def test_nonfinite_or_ambiguous_external_links_reject(value):
    with pytest.raises(ValueError):em.operators(phase=value)


def test_wrong_candidate_or_projected_carrier_cannot_reuse_the_sum_rule():
    spectrum=dimer.setup();ops=em.operators()
    projector=dimer.candidate_projectors(spectrum)['charge_minus_one_hole']
    with pytest.raises(ValueError):
        em.candidate_response(spectrum,projector,ops['dipole'],ops['current'],charge=1)
    with pytest.raises(ValueError):
        em.candidate_response(spectrum,np.eye(16)/8,ops['dipole'],ops['current'],charge=-1)
    with pytest.raises(ValueError):em.operators(FieldGeometry((1,1,1)))


def test_retained_result_uses_checked_closed_formulas_and_complete_grouping():
    result=em.build()
    assert json.loads(json.dumps(result,allow_nan=False)) == em.build()
    assert result['QED_identification'].startswith('NOT_ESTABLISHED')
    assert result['magnetic_moment_or_gyromagnetic_ratio'] is None
    for candidate in result['candidate_responses']:
        assert candidate['closed_formula_response']['current_total_weight'] == '1/256'
        assert sum(row['final_states'] for row in candidate['complete_grouped_transition_display']) == 16
        assert sum(row['dipole_weight'] for row in candidate['complete_grouped_transition_display']) == .25
    assert 'operator_residuals' not in result
