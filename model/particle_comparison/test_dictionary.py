"""The supplied SI dictionary must act on the complete state, operators and clock."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.particle_response import dimer
from model.particle_comparison import dictionary as d


@pytest.mark.parametrize('axis',[0,1,2])
def test_same_reference_projectors_symmetry_and_assigned_charge_spin(axis):
    x=d.BenchmarkDictionary(axis)
    assert x.dimension==16
    assert np.linalg.matrix_rank(x.reference,tol=1e-10)==1
    q=x.charge_C/float(d.ELEMENTARY_CHARGE_C)
    spin=[s/x.hbar_J_s for s in x.spin_J_s]
    casimir=sum(s@s for s in spin)
    h=x.hamiltonian_J/x.scale_J
    assert np.linalg.norm(h@x.reference,2)<1e-13
    assert np.linalg.norm(q@x.reference,2)<1e-13
    conjugation=x.charge_conjugation
    assert np.linalg.norm(conjugation@q@conjugation.conj().T+q,2)<1e-13
    assert np.linalg.norm(conjugation@h@conjugation.conj().T-h,2)<1e-13
    for name,p in x.projectors.items():
        charge=(-1 if 'minus' in name else 1)*(2 if 'two' in name else 1)
        assert np.linalg.norm(q@p-charge*p,2)<1e-13
        assert np.linalg.norm(casimir@p-(.75 if abs(charge)==1 else 0)*p,2)<1e-13
        if abs(charge)==1:
            assert np.linalg.norm((x.hamiltonian_J/float(d.MEV_TO_J))@p-
                                  float(x.mass_energy_MeV)*p,2)<1e-13
    # The candidates are not used to project away the other ten Fock states.
    assert round(np.trace(sum(x.projectors.values())).real)==6
    assert x.hamiltonian_J.shape==(16,16)


def test_physical_heisenberg_dynamics_and_untouched_reference_match_model_coordinate():
    x=d.BenchmarkDictionary()
    rng=np.random.default_rng(113)
    a=rng.normal(size=(16,16))+1j*rng.normal(size=(16,16))
    a=a+a.conj().T
    duration=.37;seconds=duration*x.seconds_per_model_time
    independent=expm(-1j*x.hamiltonian_J*seconds/x.hbar_J_s)
    model=expm(-1j*x.spectrum.generator*duration)
    assert np.linalg.norm(x.unitary(seconds)-independent,2)<1e-13
    assert np.linalg.norm(x.heisenberg(a,seconds)-model.conj().T@a@model,2)<1e-13
    state=rng.normal(size=(16,3))+1j*rng.normal(size=(16,3));state/=np.linalg.norm(state)
    output=x.unitary(seconds)@state
    assert np.linalg.norm(output.conj().T@output-state.conj().T@state,2)<1e-13
    assert math.isclose(x.spectrum.geometry.speed*x.metres_per_model_length/
                        x.seconds_per_model_time,float(d.LIGHT_SPEED_M_S),rel_tol=2e-16)


@pytest.mark.parametrize('axis',[0,1,2])
def test_local_charge_current_has_continuity_and_independent_peierls_units(axis):
    x=d.BenchmarkDictionary(axis)
    e=float(d.ELEMENTARY_CHARGE_C); scale=e/x.seconds_per_model_time
    q0,q1=x.site_charge_C(0),x.site_charge_C(1)
    assert np.linalg.norm((q0+q1-x.charge_C)/e,2)<1e-14
    current=x.current_A()
    derivative=1j*(x.hamiltonian_J@q0-q0@x.hamiltonian_J)/x.hbar_J_s
    assert np.linalg.norm((derivative+current)/scale,2)<1e-14
    assert np.linalg.norm((x.current_A((1,0))+current)/scale,2)<1e-14
    # Independent directed hopping derivative dH(A)/dA, with A dimensionless.
    from model.interacting_field.lattice import SIGMA
    a=[x.spectrum.space.annihilator(j).toarray() for j in range(4)]
    forward=sum(1j/16*SIGMA[axis][s,t]*(a[s].conj().T@a[2+t])
                for s in range(2) for t in range(2))
    derivative_J=x.scale_J*1j*(forward-forward.conj().T)
    expected=e/x.hbar_J_s*derivative_J
    assert np.linalg.norm((current-expected)/scale,2)<1e-14
    assert np.linalg.norm(current/scale,2)>0
    assert np.linalg.norm((current-current.conj().T)/scale,2)<1e-14


def test_exact_SI_definitions_unit_intervals_and_fully_correlated_calibration():
    x=d.BenchmarkDictionary();r=x.report()
    assert d.PLANCK_J_S==F(662607015,10**42)
    assert d.ELEMENTARY_CHARGE_C==F(1602176634,10**28)
    assert d.LIGHT_SPEED_M_S==299792458
    for key,value in [('seconds_per_model_time',x.seconds_per_model_time),
                      ('metres_per_model_length',x.metres_per_model_length)]:
        interval=r['assigned_units'][key]
        # Float evaluation has a separately declared host error, not an exact
        # interval membership assertion narrower than binary64 spacing.
        assert math.isclose(value,float(F(interval['lo'])),rel_tol=5e-16)
        assert F(interval['lo'])<=F(interval['hi'])
    p=r['paired_charge_one_assignment']
    assert p['calibration_correlation_matrix']==[[1,1],[1,1]]
    assert all(F(v)==x.mass_uncertainty_MeV**2 for row in p['mass_energy_covariance_MeV_squared'] for v in row)
    assert p['calibrated_mass_kg'][0]==p['calibrated_mass_kg'][1]
    assert r['calibration']['independently_predicted_absolute_masses']==0
    assert not r['assigned_units']['native_incident_or_modular_clock_intertwiner_supplied']
    assert F(r['assigned_units']['relative_standard_uncertainty_first_order'])==x.mass_uncertainty_MeV/x.mass_energy_MeV
    resolution=r['dimensionless_resolution_check']
    ratio=resolution['reduced_Compton_length_over_assigned_lattice_spacing']['decimal_enclosure']
    assert F(27,100)<F(ratio['lo'])<F(ratio['hi'])<F(28,100)
    reciprocal=resolution['assigned_lattice_spacing_over_reduced_Compton_length']['decimal_enclosure']
    assert F(36,10)<F(reciprocal['lo'])<F(reciprocal['hi'])<F(37,10)
    assert resolution['universal_no_go_claimed'] is False
    assert r['physical_dictionary_supplied']
    for key in ('physical_identity_derived','physical_clock_independently_calibrated',
                'electromagnetic_dynamics_derived','Lorentz_representation_derived',
                'topological_defect_mapping_derived','OPH_universe_membership_confirmed'):
        assert r[key] is False
    assert r['eligible_CP11_physical_comparisons_added']==0


@pytest.mark.parametrize('axis',[True,-1,3,.0])
def test_wrong_dictionary_axis_rejected(axis):
    with pytest.raises(ValueError):d.BenchmarkDictionary(axis)


def test_malformed_observables_times_and_edges_fail_without_projection():
    x=d.BenchmarkDictionary()
    for a in (np.zeros((4,4)),np.full((16,16),np.nan),np.ma.array(np.eye(16),mask=False)):
        with pytest.raises(ValueError):x.heisenberg(a,0.)
    for t in (True,float('nan'),float('inf'),33*x.seconds_per_model_time):
        with pytest.raises(ValueError):x.unitary(t)
    for edge in ((0,0),(True,0),[0,1]):
        with pytest.raises(ValueError):x.current_A(edge)
