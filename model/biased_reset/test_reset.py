"""Independent kernel/potential, complete-service and hostile controls."""
from collections import Counter
from fractions import Fraction as F
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
import pytest
from . import reset as r
from . import experiment


def independent_core_row(center,weights):
    row=Counter()
    for a in center:
        for seam in r.t.SEAMS:
            if a not in seam:
                continue
            b=next(k for k in seam if k!=a)
            if b not in center:
                target=tuple(sorted((set(center)-{a})|{b}))
                row[target]+=F(weights[b],60*(weights[a]+weights[b]))
    row[center]=1-sum(row.values())
    return row


@pytest.mark.parametrize('fair',[False,True])
def test_actual_native_kernel_and_exact_optimal_potential(fair):
    p,full,_,_,h,maximum,best=r.exact_chain(fair)
    for i,c in enumerate(r.t.LIVE):
        row=independent_core_row(c,p.port_weights)
        assert tuple(row[d] for d in r.t.CORES)==full[i]
        assert h[i]==1+sum(row[d]*h[j] for j,d in enumerate(r.t.LIVE))
        assert 0<h[i]<=maximum
    if fair:
        assert maximum==150 and len(best)==30
        assert all(v==(150 if c in r.t.SEAMS else 120) for c,v in zip(r.t.LIVE,h,strict=True))
    else:
        assert best==((0,1),)
        assert 203<F(maximum)<204


def test_stationary_full_core_reset_rate_and_telescoping_bound():
    cert=r.renewal_certificate()
    pi=list(map(F,cert['stationary_core_probabilities']))
    p,full,q,hazard,h,maximum,_=r.exact_chain()
    assert sum(pi)==1 and all(v>0 for v in pi)
    assert sum(a*b for a,b in zip(pi,hazard,strict=True))==1/maximum
    # A reset may be preemptive/randomized, but every individual increment
    # is bounded by maximum, hence every convex combination is as well.
    for before in (*h,F(0)):
        assert all(after-before<=maximum for after in h)
    values=list(map(F,r.bellman_certificate(64)['terminal_values']))
    assert all(v>=max(F(0),(64-h[i])/maximum) for i,v in enumerate(values))


def test_optimal_reset_does_not_make_transport_clock_invariant():
    cert=r.clock_boundary_certificate()
    assert cert['best_reset_core_is_fixed_clock_invariant'] is True
    assert 0<F(cert['stationary_mass_on_clock_incompatible_codes'])<1
    assert 0<F(cert['stationary_mean_clock_leakage_sine_squared_coefficient'])<=1
    assert F(cert['one_attempt_probability_to_leave_clock_invariant_codes_from_reset_core'])>0
    assert cert['joint_preserved_moving_clock_subspace_obtained'] is False


def test_bellman_horizons_retain_voluntary_actions_and_all_boundary_costs():
    _,_,q,hazard,_,_,_=r.exact_chain()
    one=r.bellman_certificate(1)
    assert list(map(F,one['terminal_values']))==list(hazard)
    two=r.bellman_certificate(2)
    raw=[sum(q[i][j]*hazard[j] for j in range(60))+hazard[i]*(1+min(hazard)) for i in range(60)]
    assert list(map(F,two['terminal_values']))==[min(v,1+min(raw)) for v in raw]
    finite=r.bellman_certificate()
    assert finite['same_boundary_target_proved_optimal_for_all_horizons'] is False
    assert all(not row['voluntary_reset_better_cores'] for row in finite['rows'])
    assert all([0,1] in row['minimizing_initial_cores'] for row in finite['rows'])
    assert all(F(row['maximum_continuation_value_spread'])<1 for row in finite['rows'])


def test_all_raw_slots_reset_only_positive_boundary_flags():
    cert=r.exact_service_certificate()
    counts=cert['census']
    assert counts['all_raw_slots']==43200
    assert counts['boundary_reset_positive']==120
    assert counts['ordinary_continue_positive']==11460
    assert cert['maximum_extra_boundary_reset_pulses']>0
    state=np.arange(512,dtype=complex)/1000
    center,law,seam,_,_=r.t.witness_path()[0]
    for coin in(0,1):
        for raw in range(1<<r.t.circuit.RECORD_BITS[law]):
            images,moved,_=r.t.branch(center,law,seam,coin,raw)
            if not images[0] or moved not in r.t.STOP:
                assert np.array_equal(r.recover_boundary(center,law,seam,coin,raw,state),state)
                assert all(v==0 for v in r.boundary_resources(center,law,seam,coin,raw).values())


def test_source_motion_boundary_reset_and_further_motion_preserve_reference():
    d=r.source_diagnostics()
    assert d['maximum_branch_operator_error']<1e-8
    assert d['maximum_total_mass_error']<1e-8
    assert [row['status'] for row in d['path']]==['continue','reset_continue','continue']
    assert [row['complete_raw_outputs_retained'] for row in d['path']]==[8,16,8]
    assert d['path'][1]['native_core_after'] in [list(c) for c in r.t.STOP]
    assert d['path'][1]['core_after']==[0,1]
    assert all(F(row['unnormalized_prefix_weight_given_selected_seams'])<1 for row in d['path'])


def test_silent_stop_instead_of_reset_fails_source_identity(monkeypatch):
    monkeypatch.setattr(r,'recover_boundary',lambda center,law,seam,coin,raw,state,**kw:state)
    with pytest.raises(ArithmeticError,match='fixed numerical source gate'):
        r.source_diagnostics()


def test_corrupt_emitted_reset_schedule_fails(monkeypatch):
    original=r.t.routing.permutation_edges
    monkeypatch.setattr(r.t.routing,'permutation_edges',lambda p,**kw:original(p,**kw)[1:])
    with pytest.raises(ArithmeticError,match='Gray schedule'):
        r.exact_service_certificate()


@pytest.mark.parametrize('bad',[True,-1,129,1.5])
def test_invalid_horizon_rejected(bad):
    with pytest.raises(ValueError):r.bellman_certificate(bad)


def test_retained_source_custody_and_no_physical_promotion():
    value=json.loads((Path(__file__).parent/'results.json').read_text())
    assert value['source_bindings']==experiment.source_bindings()
    assert value['integrated_gate_accepted'] is False
    assert value['physical_comparison_count_added']==0
    assert value['retained_obstructions']['fixed_reference_clock_invariant_cores']==19
    assert value['retained_obstructions']['reference_like_full_exterior_obtained'] is False
    altered=deepcopy(value);altered['integrated_gate_accepted']=True
    with pytest.raises(ValueError):experiment.compare(value,altered)
