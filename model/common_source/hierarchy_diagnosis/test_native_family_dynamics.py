"""Independent swap-channel, residual symmetry and projection-memory checks."""
from fractions import Fraction as F
import hashlib
import json

import numpy as np
import pytest
import sympy as s

from . import native_family_dynamics as dynamics
from .. import slow_geometry
from ..chiral_candidate import action
from ...native_modular_reference import local_preparation
from ...observer import native


@pytest.fixture(scope="module")
def packet():
    return dynamics.report()


def _matrix(rows):
    return s.ImmutableMatrix([[s.sympify(x) for x in row] for row in rows])


def _zero(matrix):
    return all(s.simplify(x)==0 for x in matrix)


def test_actual_endpoint_channel_matches_independent_half_swap_including_timeout():
    # At one selector round the nonnegligible terminal probability makes
    # accidental successful-outcome renormalization immediately detectable.
    selection=local_preparation.selectors(1)
    tau=s.Rational(selection["completion_timeout_probability"])
    states=native.states(binary=True);indices={x:i for i,x in enumerate(states)}
    expected=s.eye(66)*tau
    for (a,b),weight in selection["completion_accepted_seam_weights"].items():
        swap=s.zeros(66)
        for column,state in enumerate(states):
            image=list(state);image[a],image[b]=image[b],image[a]
            swap[indices[tuple(image)],column]=1
        expected+=s.Rational(weight)*(s.eye(66)+swap)/2
    actual=dynamics.proposal_operator(1)
    assert actual==expected
    sigma=s.Matrix(list(map(s.Rational,local_preparation.reference())))
    assert actual*sigma==sigma
    assert (actual-tau*s.eye(66))*sigma==(1-tau)*sigma
    branches=dynamics.proposal_branches(1)
    assert sum(flag=="selector_timeout" for _,flag,*rest in branches)==66
    for before in range(66):
        assert sum(p for i,*rest,p in branches if i==before)==1


def test_reference_symmetry_and_full_context_have_different_stabilizers(packet):
    static=packet["static_reference_obstruction"]
    assert static["reference_proper_stabilizer_size"]==2
    assert static["commutant_complex_dimension"]==5
    residual=_matrix(static["half_turn_on_family"])
    assert residual==s.diag(-1,1,-1)
    assert residual*action.family_feature()==action.family_feature()*residual
    response=dynamics.compressed_response()
    assert not _zero(residual*response-response*residual)
    context=packet["full_preparation_context"]
    assert not context["same_caps_after_half_turn"]
    assert context["caps_after_static_half_turn"]==[[0,1,5,7],[0,1,3,11]]
    assert static["weak_CP_Jarlskog_from_only_static_reference_vanishes"]


def test_source_feature_spectral_sandwiches_generate_all_nine_matrix_units(packet):
    feature=action.family_feature();response=dynamics.compressed_response()
    projections=[]
    for i in range(3):
        p=s.eye(3)
        for j in range(3):
            if i!=j:
                p=p*(feature-feature[j,j]*s.eye(3))/(feature[i,i]-feature[j,j])
        projections.append(p.applyfunc(s.simplify))
    units=[]
    for i in range(3):
        for j in range(3):
            if i==j:
                unit=projections[i]
            else:
                assert response[i,j]!=0
                unit=projections[i]*response*projections[j]/response[i,j]
            target=s.zeros(3);target[i,j]=1
            assert _zero(unit-target)
            units.append(s.Matrix(9,1,list(target)))
    assert s.Matrix.hstack(*units).det()!=0
    assert packet["full_preparation_context"]["generated_complex_matrix_algebra_dimension"]==9


def test_nonclosure_and_second_moment_memory_are_visible_in_uncompressed_twelve_port_map(packet):
    c=slow_geometry.source()
    basis=np.array([[float(a)+float(b)*np.sqrt(5) for a,b in row]
                    for row in slow_geometry.raw_frame()])
    left=np.linalg.solve(basis.T@basis,basis.T)
    transition=np.array(_matrix(packet["full_preparation_context"]["actual_T12"]),dtype=float)
    compressed=left@transition@basis
    leakage=(np.eye(12)-basis@left)@transition@basis
    norm=np.trace(leakage.T@leakage)/(basis.T@basis)[0,0]
    assert norm>0
    assert norm==pytest.approx(packet["projection_boundary"]["outside_band_squared_norm_float_diagnostic"],rel=2e-13)
    np.testing.assert_allclose(left@transition@transition@basis-compressed@compressed,
                               leakage.T@leakage/(basis.T@basis)[0,0],rtol=0,atol=5e-16)
    assert not np.allclose(left@transition@transition@basis,compressed@compressed,rtol=0,atol=1e-6)
    assert not packet["projection_boundary"]["autonomous_three_family_repair_dynamics_established"]


def test_flags_source_hashes_and_physical_scope_are_explicit(packet):
    instrument=packet["complete_reduced_instrument"]
    tau=F(instrument["terminal_timeout_probability"])
    assert 0<tau<1 and F(instrument["accepted_probability"])+tau==1
    assert instrument["accepted_seam_coin_current_and_timeout_outcomes_retained"]
    assert instrument["full_diagonal_reference_stationary_exact_for_C_and_E"]
    assert not instrument["fine_auxiliary_selector_tapes_enumerated"]
    assert not packet["scope"]["Yukawa_map_selected"]
    assert not packet["scope"]["response_eigenvalues_identified_as_particle_masses"]
    assert not packet["scope"]["Standard_Model_masses_or_CKM_derived"]
    for name,digest in packet["source_bindings"].items():
        assert hashlib.sha256((dynamics.ROOT/name).read_bytes()).hexdigest()==digest
    json.dumps(packet,allow_nan=False)


@pytest.mark.parametrize("rounds",[True,0,33,1.0])
def test_selector_input_validation_cannot_be_bypassed_by_cached_integer(rounds):
    dynamics.proposal_branches(1)
    with pytest.raises(ValueError):
        dynamics.proposal_branches(rounds)
