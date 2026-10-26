"""Full population, complementary-memory and native clock comparisons."""
import json

import numpy as np
import pytest
import sympy as s

from . import memory
from ...observer import native


def _numeric(matrix):
    return np.array(matrix,dtype=float)


def _zero(matrix):
    return all(s.expand(x)==0 or s.simplify(x)==0 for x in matrix)


def _point_preparation():
    states=native.states(binary=True)
    index=next(i for i,x in enumerate(states) if x[0]==1 and x[5]==1)
    result=s.zeros(66,1);result[index]=1
    return result


def test_actual_continuing_accepted_history_retains_each_timeout_without_restart():
    initial=_point_preparation()
    history=memory.proposal_history(3,initial)
    accepted=memory.accepted_state_map()
    tau=1-sum(accepted[:,0]);q=1-tau
    assert history["total_mass"]==1 and history["live_mass"]==q**3
    assert history["live"]==accepted**3*initial
    assert [row["timeout_proposal"] for row in history["terminal_timeouts"]]==[1,2,3]
    for j,row in enumerate(history["terminal_timeouts"]):
        assert row["work"]==tau*accepted**j*initial
        assert sum(row["work"])==tau*q**j
    # An iteration of the one-proposal unread T would hide the absorbing flag.
    assert sum(memory.source_state_map()**3*initial)==1
    assert history["live_mass"]<1


def test_exact_conserved_projectors_are_full_source_components_not_three_guessed_poles():
    fixed=memory.conserved_sectors();d=memory._data()
    assert sorted(map(len,fixed["groups12"]))==[1,1,10]
    assert sorted(map(len,fixed["groups66"]))==[1,10,10,45]
    assert memory.accepted_state_map()*fixed["P66"]==d["q"]*fixed["P66"]
    assert memory.conditional_port_map()*fixed["P12"]==fixed["P12"]
    assert fixed["family_fixed_residue"].rank()==2
    # Normalizing a continuing work state does not remove the retained q^n weight.
    assert d["q"]<1 and d["C66"]*fixed["P66"]==fixed["P66"]


def test_memory_recurrence_reproduces_independent_full_twelve_port_moments():
    d=memory._data();kernel=memory.memory_kernel(4,map_kind="conditional")
    c,b,left=map(_numeric,(d["C"],d["B"],d["left"]))
    r=_numeric(kernel["instantaneous"])
    coefficients=list(map(_numeric,kernel["coefficients"]))
    moments=[np.eye(3)]
    for n in range(5):
        new=r@moments[n]
        for j in range(n):new+=coefficients[n-1-j]@moments[j]
        moments.append(new)
        np.testing.assert_allclose(new,left@np.linalg.matrix_power(c,n+1)@b,rtol=0,atol=1e-14)
    assert np.linalg.norm(moments[2]-r@r)>1e-3


def test_resolvent_matches_an_independent_explicit_complement_elimination():
    d=memory._data();z=2.
    c,b,left,q=map(_numeric,(d["C"],d["B"],d["left"],d["Q"]))
    sigma=left@c@q@np.linalg.solve(z*np.eye(12)-q@c@q,q@c@b)
    expected=np.linalg.inv(z*np.eye(3)-left@c@b-sigma)
    np.testing.assert_allclose(_numeric(memory.family_resolvent(2)),expected,rtol=0,atol=8e-16)
    np.testing.assert_allclose(_numeric(memory.schur_self_energy(2)),sigma,rtol=0,atol=8e-16)


def test_closed_minimal_realization_preserves_exact_source_and_visible_conserved_directions():
    d=memory._data();minimal=memory.minimal_realization()
    assert minimal["dimension"]==11
    v,c,j,o=(minimal[name] for name in ("embedding","conditional_map","input_map","output_map"))
    assert _zero(d["C"]*v-v*c)
    assert _zero(v*j-d["B"])
    assert _zero(o*j-s.eye(3))
    assert _zero(o*c*j-d["left"]*d["C"]*d["B"])
    assert _zero(s.ones(1,12)*v)
    assert minimal["unreachable_total_charge_retained_in_full_source"]


def test_existing_incident_clock_retains_failures_and_unfinished_waits():
    d=memory._data();clock=memory.incident_clock_maps()
    history=memory.finite_tick_history(3,_point_preparation())
    assert history["total_mass"]==1
    assert sum(history["unfinished"])==clock["nonincident_probability"]**3
    for row in history["incident_successes"]:
        j=row["nonincident_waits"]
        assert row["work"]==clock["J66"]*clock["R66"]**j*_point_preparation()
    for row in history["terminal_timeouts"]:
        j=row["nonincident_waits"]
        assert row["work"]==d["tau"]*clock["R66"]**j*_point_preparation()
    assert s.ones(1,12)*(clock["success12"]+clock["failure12"])==s.ones(1,12)
    assert clock["tick_survival"]<1
    # Stopping changes the map nonuniformly: it is not simply another tick size.
    assert not _zero(clock["J12"]*clock["R12"]-clock["R12"]*clock["J12"])


def test_full_population_and_incident_clock_rate_diagnostics_do_not_select_masses():
    d=memory._data();clock=memory.incident_clock_maps()
    accepted=np.linalg.eigvalsh(_numeric(d["C66"]))
    positive=accepted[accepted<1-1e-10]
    assert len(positive)==62 and np.min(positive)>0
    rates=-np.log(positive)
    assert max(rates)/min(rates)<20
    j,r=_numeric(clock["J66"]),_numeric(clock["R66"])
    stopped=np.linalg.solve(np.eye(66)-r.T,j.T).T/float(clock["tick_survival"])
    values=np.linalg.eigvals(stopped)
    assert np.max(abs(values.imag))<2e-14
    assert sum(abs(values.real-1)<1e-10)==4


@pytest.mark.parametrize("bad",[True,-1,33,1.0])
def test_invalid_proposal_horizon_rejected(bad):
    with pytest.raises(ValueError):memory.proposal_history(bad)


def test_resolvent_domain_and_input_normalization_are_explicit():
    with pytest.raises(ValueError):memory.family_resolvent(s.sqrt(2))
    with pytest.raises(ValueError):memory.family_resolvent(1)
    with pytest.raises(ValueError):memory.proposal_history(2,s.zeros(66,1))
    with pytest.raises(ValueError):memory.memory_kernel(True)


def test_stopped_clock_diagnostic_does_not_drop_genuine_complex_or_negative_poles():
    values,residual=memory._clock_poles(np.diag([.3,.8,1.]))
    np.testing.assert_array_equal(values,[.3,.8,1.])
    assert residual==0
    with pytest.raises(ArithmeticError):
        memory._clock_poles(np.array([[0.,-1.],[1.,0.]]))
    with pytest.raises(ArithmeticError):
        memory._clock_poles(-np.eye(2))
