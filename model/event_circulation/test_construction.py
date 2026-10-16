"""Independent jump enumeration, counter laws, preparation and capped records."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import construction as c


def test_event_stationarity_is_not_attempt_stationarity_or_first_event_preparation():
    degree,palm,_,jump = c.embedded()
    assert degree.count(8)==30 and degree.count(10)==36
    assert c.clock._row_times(palm,jump)==palm
    uniform = (F(1,66),)*66
    assert c.clock._row_times(uniform,jump)!=uniform
    assert c.clock._row_times(uniform,jump)!=palm
    assert c.geometry()['single_attempt_Palm_preparation_success']=='5/33'


def test_two_event_joint_law_matches_direct_native_movement_enumeration():
    kernel,_,weights,_ = c.circulation.source()
    index = {state:i for i,state in enumerate(kernel.states)}
    def moves(state):
        result=[]
        for e,seam in enumerate(kernel.seams):
            for coin in (0,1):
                out,record = c.circulation.native.native(state,seam,coin)
                if record[1]:
                    result.append((out,weights[e]*record[1]))
        return result
    movement = {state:moves(state) for state in kernel.states}
    direct = [[F(0)]*66 for _ in range(5)]
    for state in kernel.states:
        first = movement[state]
        for middle,g in first:
            second = movement[middle]
            for final,h in second:
                direct[g+h+2][index[final]] += F(len(first),600)/len(first)/len(second)
    rows,den = c.count_laws((2,))[2]
    assert tuple(tuple(F(x,den) for x in row) for row in rows)==tuple(map(tuple,direct))


def test_complete_counter_mass_variance_and_terminal_state():
    _,palm,_,_ = c.embedded()
    for n,(rows,den) in c.count_laws((0,1,8,32)).items():
        p = [F(sum(row),den) for row in rows]
        assert sum(p)==1
        assert p==list(reversed(p))
        assert sum(F(s-n)**2*value for s,value in enumerate(p))==F(n,10)
        assert tuple(F(sum(row[i] for row in rows),den) for i in range(66))==palm


def test_attempt_and_event_increments_have_different_zero_mean_covariances():
    data=c.geometry()
    assert data['stationary_edge_second_moment_per_event']=='1/30'
    assert F(1,30)*F(5,33)==F(1,198)
    _,palm,marked,_=c.embedded()
    variances=[]
    for i in range(66):
        assert sum(marked[0][i])==sum(marked[2][i])
        variances.append(F(sum(marked[0][i])+sum(marked[2][i]),40))
    assert len(set(variances))>1
    assert sum(p*v for p,v in zip(palm,variances))==F(1,10)


def test_nonuniform_identity_holding_preserves_only_stopped_event_laws():
    data=c.holding_control()
    assert data['all_marked_event_rows_identical']
    assert not data['waiting_costs_and_attempt_records_preserved']
    assert not data['work_evolving_during_holds_covered']


def test_four_outcome_native_instrument_and_all_capped_terminal_states():
    instrument=c.attempt_instrument()
    total=sum(k.conj().T@k for group in instrument.kraus for k in group)
    np.testing.assert_allclose(total,np.eye(66),atol=2e-15,rtol=0)
    data=c.runtime_check()
    assert data['executed_complete_herald_and_three_attempt_records']==128
    assert len(data['first_event_records'])==9
    assert sum(F(row['conditional_probability']) for row in data['first_event_records'])+F(data['conditional_no_event_timeout_probability'])==1
    assert sum(F(row['unconditional_probability']) for row in data['first_event_records'])+F(data['unconditional_no_event_timeout_probability'])+F(data['preparation_failure_probability'])==1
    assert data['complete_stopped_state_residual_diagnostic']<1e-14
    assert data['checkpoint_complete_distance_diagnostic']==0


@pytest.mark.parametrize('horizons',[(),(1,1),(-1,),(True,),(1.5,),(513,)])
def test_invalid_event_horizon_fails(horizons):
    with pytest.raises(ValueError):
        c.count_laws(horizons)
