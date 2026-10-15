"""Native counter laws, conserved gradients and complete stateful execution."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as sp

from . import construction as c


def test_cycle_projection_is_exact_and_gradient_observables_do_not_diffuse():
    data=c.exact_geometry()
    p=sp.Matrix([[sp.Rational(v) for v in row] for row in data['cycle_projection']])
    kernel,_,weights,branches=c.source()
    b=sp.zeros(12,30)
    for e,(a,z) in enumerate(kernel.seams): b[a,e]=-1;b[z,e]=1
    assert p*p==p and p.T==p and p.rank()==19
    assert p*b.T==sp.zeros(30,12)
    # Each vertex gradient integrates to its bounded endpoint load change.
    for i,j,e,coin,current,g in branches:
        assert int(b[0,e])*current==kernel.states[j][0]-kernel.states[i][0]
    assert max(state[0] for state in kernel.states)-min(state[0] for state in kernel.states)==1
    assert sum(weights)==1  # Orientation signs are not all positive.


def test_conditional_zero_drift_but_state_dependent_variance():
    m=c.marked()
    variances=[]
    for i in range(66):
        assert sum(m[0][i])==sum(m[2][i])
        variances.append(F(sum(m[0][i])+sum(m[2][i]),60))
    assert set(variances)=={F(0),F(1,30)}
    assert sum(variances)/66==F(1,66)


def test_complete_two_step_law_matches_independent_native_enumeration():
    kernel,_,weights,_=c.source()
    direct={s:[0]*66 for s in range(-2,3)}
    index={state:i for i,state in enumerate(kernel.states)}
    for state in kernel.states:
        for e,seam in enumerate(kernel.seams):
            for coin in (0,1):
                middle,r=c.native.native(state,seam,coin)
                for f,other in enumerate(kernel.seams):
                    for next_coin in (0,1):
                        final,next_r=c.native.native(middle,other,next_coin)
                        direct[weights[e]*r[1]+weights[f]*next_r[1]][index[final]]+=1
    rows,den=c.count_laws((2,))[2]
    assert den==66*60**2
    assert rows==tuple(tuple(direct[s]) for s in range(-2,3))


def test_complete_three_record_runtime_matches_exact_counter_states():
    output=c.runtime_check()
    assert len(output['complete_records'])==27
    assert sum(F(row['probability']) for row in output['complete_records'])==1
    assert output['checkpoint_complete_distance']==0
    assert output['maximum_complete_output_residual']<1e-14


def test_one_history_variance_scales_without_independent_copies():
    for n,(rows,den) in c.count_laws((1,8,32)).items():
        probabilities=[F(sum(row),den) for row in rows]
        assert sum(probabilities)==1
        assert probabilities==list(reversed(probabilities))
        assert sum(F(s-n)**2*p for s,p in enumerate(probabilities))==F(n,66)
        assert all(sum(row)>=0 for row in rows)


def test_native_instrument_complete_and_preserves_untouched_reference():
    operation=c.instrument()
    summed=sum(k.conj().T@k for group in operation.kraus for k in group)
    np.testing.assert_allclose(summed,np.eye(66),atol=2e-15,rtol=0)
    vector=np.zeros(132,complex);vector[0]=vector[3]=1/np.sqrt(2)
    initial=np.outer(vector,vector.conj())
    p=c.runtime.Program('cycle_reference',66,np.eye(66)/66,
                        (c.runtime.Event('cycle',instruments=(operation,)),))
    result=c.runtime.execute(p,initial=initial)
    total=sum(branch.density for branch in result.branches).reshape(66,2,66,2)
    np.testing.assert_allclose(np.trace(total,axis1=0,axis2=2),np.eye(2)/2,atol=2e-15,rtol=0)


@pytest.mark.parametrize('n',[-1,True,1.5,513,None])
def test_invalid_counter_horizon_is_rejected(n):
    with pytest.raises(ValueError): c.count_laws((n,))


def test_empty_or_duplicate_counter_requests_rejected():
    for counts in ((),(1,1)):
        with pytest.raises(ValueError):c.count_laws(counts)


@pytest.mark.parametrize('value',[float('nan'),float('inf'),-1.,.01])
def test_nonfinite_or_large_diagnostics_fail_closed(value):
    with pytest.raises(ArithmeticError): c.bounded(value,'mutant')
