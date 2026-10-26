"""Finite exact geometry identities and independent complete-symbol control."""
import hashlib
import json

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import action,kinetic
from . import spacetime_attachment as attachment


@pytest.fixture(scope="module")
def packet():
    return attachment.report()


def _numeric(matrix):
    return np.array(matrix.evalf(18),dtype=complex)


def test_frozen_chart_uses_actual_native_increment_tensor_and_lapse(packet):
    data=attachment.frozen_geometry()
    assert data["A"][0,0]==1
    assert data["A"][1,1]==s.Rational(5039,181440)
    assert data["A"][1,2]==-s.Rational(5039,362880)
    assert all(value>0 for value in data["A"].eigenvals())
    assert s.simplify(data["N"]**2-s.sqrt(data["A"].det()))==0
    assert packet["frozen_geometry"]["signature"]==[-1,1,1,1]
    assert not packet["frozen_geometry"]["finite_increment_is_continuum_tensor"]


def test_full_actual_chiral_kernel_has_same_frozen_chart_dirac_symbol():
    data=attachment.frozen_geometry()
    momentum=(s.Rational(2,7),s.Rational(-1,3),s.Rational(3,8))
    q=tuple(float(x) for x in attachment.weyl_pullback(momentum))
    kernel=attachment.coordinate_mass(action.mass_kernel())
    actual=kinetic.spin_conjugated_symbol(q,kernel,speed=1.)
    sigma=kinetic.pauli_matrices()
    alpha=[np.block([[np.zeros((2,2)),p],[p,np.zeros((2,2))]]) for p in sigma]
    beta=np.diag([1.,1.,-1.,-1.])
    old=sum((qj*np.kron(a,np.eye(45)) for qj,a in zip(q,alpha,strict=True)),np.zeros((180,180),complex))
    old+=np.kron(beta,_numeric(kernel))
    change=np.kron(_numeric(attachment.clifford_intertwiner()),np.eye(45))
    np.testing.assert_allclose(change@old@change.conj().T,actual,rtol=0,atol=2e-16)


def test_exact_scalar_cone_mass_shell_and_clock_regraduation():
    data=attachment.frozen_geometry()
    p=s.ImmutableMatrix([s.Rational(3,2),-2,s.Rational(1,4)])
    h=attachment.weyl_symbol(p)
    scalar=(p.T*data["A"]*p)[0]
    assert all(s.simplify(x)==0 for x in h*h-scalar*s.eye(2))
    proper_mass=s.Rational(7,3)
    coordinate=attachment.coordinate_mass(proper_mass)
    energy=s.sqrt(scalar+coordinate**2)
    covector=s.ImmutableMatrix([energy,*p])
    assert s.simplify((covector.T*data["g"].inv()*covector)[0]+proper_mass**2)==0
    kappa=s.Rational(5,2)
    omega_new=action.source_coefficients()["Omega"]
    omega_proper=omega_new/(kappa*data["N"])
    mass_new=kappa*data["N"]*proper_mass
    assert s.simplify(mass_new/omega_new-proper_mass/omega_proper)==0


def test_full_incidence_kernel_and_radial_self_adjointness_are_not_erased(packet):
    boundary=packet["radial_boundary"]
    assert boundary["minimal_two_component_deficiency_indices"]==[2,2]
    assert "not self-adjoint" in boundary["scalar_Dirichlet_on_all_spinor_components"]
    assert "Neumann" in boundary["square_of_reflecting_control"]
    domain=packet["finite_domain_obstruction"]
    assert domain["complete_incidence_dimension"]==3570
    assert domain["two_component_Weyl_on_same_scalar_nodes_dimension"]==1440
    assert domain["retained_edge_cycle_zero_modes"]==2130
    assert not domain["complete_unitary_field_identification_possible"]
    assert not packet["scope"]["same_complete_scalar_and_chiral_field_history_proved"]
    assert not packet["same_metric_symbols"]["full_curved_Dirac_square_equals_scalar_B"]


def test_source_custody_is_read_only_and_claim_stays_at_frozen_chart(packet):
    for name,digest in packet["source_bindings"].items():
        assert hashlib.sha256((attachment.ROOT/name).read_bytes()).hexdigest()==digest
    assert packet["Clifford_bridge"]["three_alpha_intertwiners_exact"]
    assert not packet["Clifford_bridge"]["native_word_errors_newly_interval_certified"]
    assert not packet["scope"]["new_full_native_spatial_Dirac_history_executed"]
    assert not packet["scope"]["physical_particles_or_Standard_Model_masses_identified"]
    assert not packet["mass_and_clock"]["actual_native_incident_and_new_clock_histories_identified"]
    json.dumps(packet,allow_nan=False)


def test_exact_input_and_immutable_cached_geometry_boundaries():
    for bad in ((.1,0,0),(s.I,0,0),(True,0,0),(0,0)):
        with pytest.raises(ValueError):
            attachment.weyl_pullback(bad)
    with pytest.raises(TypeError):
        attachment.frozen_geometry()["N"]=1
    with pytest.raises(TypeError):
        attachment.frozen_geometry()["A"][0,0]=1


@pytest.mark.parametrize("law,factor",[("C",s.Rational(7,10)),("E",s.Rational(2,5))])
def test_actual_incident_clock_response_has_exact_strict_norm_and_full_flags(law,factor):
    for deadline in (1,2,3,8):
        data=attachment.native_incident_response(law,deadline)
        response=data["response"]
        assert sum(data["terminal_probabilities"].values())==1
        assert data["incident_probability"]==s.Rational(5,32)
        assert data["remote_probability"]==s.Rational(11,16)
        assert response*s.ones(10,1)==factor*s.ones(10,1)
        assert s.ones(1,10)*response==factor*s.ones(1,10)
        assert all(x>=0 for x in response)
        # Independent singular-value and eigenvalue checks of the exact
        # positive row/column-sum proof, not a floating rank decision.
        array=np.array(response,dtype=float)
        assert abs(np.linalg.norm(array,2)-float(factor))<2e-15
        assert abs(max(abs(np.linalg.eigvals(array)))-float(factor))<3e-15
    short=attachment.native_incident_response(law,1)["response"]
    long=attachment.native_incident_response(law,2)["response"]
    sites=[site for site,_,_ in attachment.collective.fields()]
    assert short[sites.index(7),sites.index(2)]==0
    assert long[sites.index(7),sites.index(2)]==s.Rational(1,540)
    assert short!=long  # one incident tick alone is not the complete record


@pytest.mark.parametrize("law,factor",[("C",s.Rational(7,10)),("E",s.Rational(2,5))])
def test_actual_complete_ramsey_outputs_measure_the_one_tick_contraction(law,factor):
    observed=attachment.record_response_control(law)
    masses={k:s.Rational(v) for k,v in observed["complete_terminal_probabilities"].items()}
    response={k:s.Rational(v) for k,v in observed["raw_signed_Ramsey_contrast_by_terminal_flag"].items()}
    assert masses=={"success":s.Rational(135,512),
                    "selector_timeout":s.Rational(135,512),
                    "deadline_timeout":s.Rational(121,256)}
    assert response["success"]==factor*masses["success"]
    assert s.Rational(observed["successful_conditional_signed_Ramsey_contrast"])==factor
    assert not observed["conditional_normalization_replaces_complete_instrument"]


def test_actual_quadratic_nambu_response_obeys_the_nonzero_intertwiner_gap():
    h=kinetic.spin_conjugated_symbol((.21,-.09,.3))
    assert h.shape==(180,180)
    rng=np.random.default_rng(8107)
    probes=[np.eye(180,10,dtype=complex),
            rng.normal(size=(180,10))+1j*rng.normal(size=(180,10))]
    for time in (0.,.2,1.7):
        unitary=expm(-1j*time*h)
        np.testing.assert_allclose(unitary.conj().T@unitary,np.eye(180),rtol=0,atol=2e-14)
        for law in ("C","E"):
            data=attachment.native_incident_response(law)
            response=np.array(data["response"],dtype=complex)
            bound=float(data["relative_unitary_intertwiner_defect_lower"])
            for probe in probes:
                defect=np.linalg.norm(unitary@probe-probe@response,"fro")/np.linalg.norm(probe,"fro")
                assert defect>=bound-1e-13
            # The universal bound is sharp for a unitary with a fixed vector;
            # a constant matrix-unit response is an independent exact check.
            probe=s.Matrix([[1]*10,[0]*10])
            assert probe-probe*data["response"]==(1-data["operator_norm_and_spectral_radius"])*probe


@pytest.mark.parametrize("law,probability",[("C","3/20"),("E","3/10")])
def test_open_environment_control_is_complete_and_is_not_a_particle_exclusion(law,probability):
    control=attachment.open_response_control(law)
    assert control["phase_flip_probability"]==probability
    assert control["trace_preserving_exact"] and control["charge_populations_preserved_exact"]
    assert not control["full_native_spatial_channel_or_physical_detector_identified"]
    assert not control["source_pulse_dilation_compiled"]
    with pytest.raises(ValueError): attachment.native_incident_response(law,True)
    with pytest.raises(ValueError): attachment.native_incident_response(law,0)
