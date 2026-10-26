"""Independent source, constrained-energy and stationary-cone controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
import sympy as s
from scipy.linalg import eigh

from . import stationary_selection as q
from ...native_modular_reference import local_preparation
from ...observer import native
from ...repair_scheduler import coins


def numeric(matrix):
    return np.array(matrix.evalf(17), float)


@pytest.fixture(scope="module")
def certificate():
    return q.report()


def test_all_actual_C_E_branches_give_exact_stationary_weight_cone(certificate):
    rows = certificate["primitive_mapping"]["seam_certificates"]
    assert certificate["primitive_mapping"]["input_coin_addresses_checked"] == 2*30*2*66
    assert len(rows) == 30
    assert sum(row["admitted_by_same_reference_stationarity"] for row in rows) == 22
    assert sum(row["crosses_core_boundary"] for row in rows) == 8
    for row in rows:
        a, b = row["seam"]
        expected = F(157, 90000) if (a in q.CORE) != (b in q.CORE) else F()
        assert F(row["reference_Dirichlet_defect"]) == expected
        assert F(row["swap_reference_difference_squared"]) == 4*expected
    # Independent scalar count: nine outer spectators and one other-core
    # spectator, with both swap orientations represented.
    expected = (18*(F(2,225)-F(1,36))**2+2*(F(1,36)-F(2,45))**2)/4
    assert expected == F(157,90000)


@pytest.mark.parametrize("law", ["C", "E"])
def test_full_quantum_reference_channel_not_only_a_mean_statistic(law):
    sigma = np.diag(list(map(float, local_preparation.reference()))).astype(complex)
    weights = q.maximal_stationary_weights()
    output = float(1-sum(weights.values(), F()))*sigma
    for edge, weight in weights.items():
        for coin in (0,1):
            output += float(weight)*coins.apply_coin(sigma, law, coin, seam=edge)
    np.testing.assert_allclose(output, sigma, atol=1e-16, rtol=0)
    cross = q.source()["excluded"][0]
    changed = output+sum((coins.apply_coin(sigma, law, coin, seam=cross) for coin in (0,1)))/64-sigma/64
    assert np.linalg.norm(changed-sigma) > 1e-4
    assert q.stationarity(q.crossing_countercontrol_weights(), law)["stationary"] is False


def test_independent_full_state_KKT_reproduces_all_complement_S_and_Z():
    # Build the source population transition from floor/ceiling, independently
    # of the producer's swap matrix and port graph formula.
    states = native.states(binary=True); indices = {state:i for i,state in enumerate(states)}
    transition = np.eye(66)*10/32
    for a,b in q.source()["admitted"]:
        for before,state in enumerate(states):
            for coin in (0,1):
                total=state[a]+state[b]; following=list(state)
                following[a]=(total+coin)//2; following[b]=(total+1-coin)//2
                transition[indices[tuple(following)],before] += 1/64
    loads=np.array(states,float).T
    gram=loads@loads.T
    # X T = (I-D) X; infer the actual port operator using a full-rank
    # source occupation chart rather than reading the declared graph energy.
    d=np.eye(12)-loads@transition@loads.T@np.linalg.inv(gram)
    geometry=q.family_geometry(); b=numeric(geometry["basis"]); g=float(geometry["gram"])
    constraints=np.vstack((b.T,np.ones((1,12))))
    kkt=np.block([[d,constraints.T],[constraints,np.zeros((4,4))]])
    right=np.vstack((np.zeros((12,3)),g*np.eye(3),np.zeros((1,3))))
    harmonic=np.linalg.solve(kkt,right)[:12]
    result=q.schur_family(q.maximal_stationary_weights())
    np.testing.assert_allclose(harmonic,numeric(result["U"]),atol=3e-13,rtol=0)
    np.testing.assert_allclose(harmonic.T@d@harmonic/g,numeric(result["S0"]),atol=3e-14,rtol=0)
    np.testing.assert_allclose(harmonic.T@harmonic/g,numeric(result["Z"]),atol=3e-13,rtol=0)


def test_omitted_core_internal_edge_lifts_one_zero_but_not_all_three():
    maximal=q.maximal_stationary_weights()
    without_core={edge:weight for edge,weight in maximal.items() if edge!=(0,1)}
    old,new=q.schur_family(without_core),q.schur_family(maximal)
    assert old["rank"] == 1 and new["rank"] == 2
    assert q.stationarity(maximal)["timeout_probability"] == "5/16"
    assert q.stationarity(maximal)["stationary"]
    x=q.family_geometry()["conserved_family_coordinate"]
    assert new["S0"]*x == s.zeros(3,1)
    values=eigh(numeric(new["S0"]),numeric(new["Z"]),eigvals_only=True)
    assert abs(values[0])<1e-14 and np.all(values[1:]>0)
    assert np.linalg.eigvalsh(numeric(new["Z"])-np.eye(3))[0]>-1e-14


def test_entire_compatible_weight_cone_keeps_zero_without_symmetric_weights():
    edges=q.source()["admitted"]
    weights={edge:F(i+1,2*sum(range(1,23))) for i,edge in enumerate(edges)}
    result=q.schur_family(weights)
    assert q.stationarity(weights)["stationary"] and result["rank"] == 2
    assert result["S0"]*q.family_geometry()["conserved_family_coordinate"] == s.zeros(3,1)
    assert not result["S0"].is_diagonal()


@pytest.mark.parametrize("law", ["C","E"])
def test_actual_record_adaptive_histories_preserve_core_charge(law):
    edges=q.source()["admitted"]
    for initial in native.states(binary=True):
        live=[(initial,(),F(1))]; terminal=[]
        for depth in range(3):
            following=[]
            for state,tape,weight in live:
                address=0 if not tape else (7*tape[-1][0]+3*abs(tape[-1][1])+depth)%len(edges)
                edge=edges[address]
                for coin in (0,1):
                    output,(_,current)=native.native(state,edge,coin)
                    label=current if law=="C" else state[edge[0]]-state[edge[1]]
                    record=(*tape,(coin,label))
                    assert sum(output[j] for j in q.CORE)==sum(initial[j] for j in q.CORE)
                    (terminal if label or depth==2 else following).append((output,record,weight/2))
            live=following
        assert sum((weight for _,_,weight in terminal),F())==1


def test_crossing_control_lifts_rank_and_pays_exact_reference_defect():
    weights=q.crossing_countercontrol_weights(); result=q.schur_family(weights)
    assert result["rank"]==3
    assert q.stationarity(weights)["reference_Dirichlet_defect"]==str(F(157,90000*64))
    cross=sum((weight for edge,weight in weights.items() if q.crossing(edge)),F())
    c=q.family_geometry()["conserved_port_vector"]
    assert (c.T*q.port_generator(weights)*c)[0]==18*s.Rational(cross)
    values=eigh(numeric(result["S0"]),numeric(result["Z"]),eigvals_only=True)
    assert 0<values[0] <= float((5-s.sqrt(5))/4)*float(cross)


def test_same_Higgs_route_has_exactly_one_massless_charged_family(certificate):
    row=certificate["same_M15_Higgs_route"]
    assert row["actual_M15_Higgs_current_identities_rechecked"]==72
    assert row["K45_rank"]==28 and row["Delta90_rank"]==56
    assert row["massless_charged_Dirac_fields_counting_color"]==7
    assert row["massive_charged_Dirac_fields_counting_color"]==14
    assert row["massless_left_Weyl_neutrino_fields"]==3
    assert not row["physical_charged_mass_or_vacuum_prediction"]


def test_uniform_reference_is_a_real_countercontrol_outside_fixed_state_domain(certificate):
    row=certificate["different_reference_countercontrol"]
    assert row["actual_total_two_population_graph_connected_states"]==66
    assert row["compressed_family_feature"]=="I3/6"
    assert not row["same_biased_reference_or_its_clock_coefficient_inputs_preserved"]


@pytest.mark.parametrize("invalid", [True,.25,-1,F(-1,4),s.Rational(1,4),"1/4"])
def test_original_weight_domain_rejects_aliases_and_negative_data(invalid):
    with pytest.raises(ValueError):
        q.stationarity({(0,1):invalid})


def test_weight_schema_and_singular_complement_have_explicit_refusal():
    for invalid in ({(False,1):F(1,2)},{(0.,1.):F(1,2)},{(0,2):F(1,2)},{(0,1):2}):
        with pytest.raises(ValueError):q.stationarity(invalid)
    assert q.stationarity({})["stationary"]
    with pytest.raises(ValueError,match="complementary"):
        q.schur_family({})


def test_cached_source_metadata_is_detached_and_receipt_types_are_exact(certificate):
    data=q.source();data["population_maps"].clear();data["rows"][0]["seam"][0]=999
    assert len(q.source()["population_maps"])==60 and q.source()["rows"][0]["seam"][0]!=999
    current=json.loads(json.dumps(certificate))
    q.compare_exact(deepcopy(current),current)
    changed=deepcopy(current);changed["reference"]["binary_total_two_dimension"]=66.0
    with pytest.raises(ValueError,match="type changed"):q.compare_exact(changed,current)
