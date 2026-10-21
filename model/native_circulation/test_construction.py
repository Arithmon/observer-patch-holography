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


def test_axial_readout_uses_the_rank_three_gram_and_every_native_branch():
    row=c.axial_geometry()
    assert row['rank']==3 and row['frame_gram_entries_checked']==144
    assert row['native_branches_checked']==3960
    assert row['all_conditional_drifts_zero'] and row['cycle_orthogonality_exact']
    assert row['covariance_per_attempt_diagonal']=='4/99'
    assert row['off_diagonal_covariances']=='0'
    # The weights' inner products are determined by the recovered Gram alone;
    # no choice of handedness or Cartesian presentation changes this metric.
    weights=c.axial_weights()
    kernel,_,_,_=c.source()
    source_gram=c.transfer.slow_geometry.source().intrinsic_gram_exact()
    gram=sp.Matrix([[sp.Rational(a)+sp.Rational(b)*sp.sqrt(5) for a,b in row]
                    for row in source_gram])
    for i,(a,b) in enumerate(kernel.seams):
        for j,(d,e) in enumerate(kernel.seams):
            expected=gram[a,d]*gram[b,e]-gram[a,e]*gram[b,d]
            assert sp.simplify(weights[:,i].dot(weights[:,j])-expected)==0


def test_axial_counter_accumulates_on_repeated_closed_native_cycles():
    row=c.axial_geometry(); witness=row['closed_cycle']
    initial=tuple(witness['initial_and_final_loads']); state=initial
    kernel,_,_,_=c.source(); weights=c.axial_weights()
    total=sp.zeros(3,1)
    increment=sp.Matrix([sp.sympify(x) for x in witness['axial_counter_increment']])
    assert increment!=sp.zeros(3,1)
    for repeats in range(1,6):
        for event in witness['events']:
            seam=tuple(event['seam'])
            state,record=c.native.native(state,seam,event['coin'])
            assert record[1]==event['current']
            total+=weights[:,kernel.seams.index(seam)]*record[1]
        assert state==initial
        assert (total-repeats*increment).applyfunc(sp.simplify)==sp.zeros(3,1)


def test_axial_geometry_rejects_a_changed_frame_dictionary(monkeypatch):
    original=c.transfer.slow_geometry.raw_frame()
    altered=list(original); altered[1]=original[2]
    monkeypatch.setattr(c.transfer.slow_geometry,'raw_frame',lambda:tuple(altered))
    with pytest.raises(ArithmeticError,match='recovered rank-three Gram'):
        c.axial_weights()


def test_axial_all_initial_state_limit_has_exact_ergodicity_and_radial_controls():
    row = c.axial_radial_limit_control()
    assert row['complete_native_branches_checked'] == 3960
    assert row['configuration_graph_diameter'] == row['Doeblin_attempt_block'] == 4
    assert row['minimum_self_loop_probability'] == '5/6'
    assert row['Doeblin_uniform_reset_probability_lower'] == '11/2160000'
    assert not row['stationary_initial_preparation_required_for_functional_limit']
    assert row['exact_finite_covariance_formula_requires_stationary_preparation']
    assert row['same_zero_counter_one_attempt_radial_second_moment'] == {
        'adjacent': '8/75', 'nonadjacent': '2/15'}
    assert row['stationary_one_attempt_radial_fourth_moment'] == '16/165'
    assert row['Gaussian_one_attempt_radial_fourth_moment_control'] == '80/3267'
    assert not row['physical_spatial_transport_or_independent_product_radius_identified']
    # Independently check the positive four-step transition lower bound on
    # the actual integer kernel; 60^4 fits exactly in the integer datatype.
    kernel, _, _, _ = c.source()
    four = np.linalg.matrix_power(np.asarray(kernel.counts, dtype=np.int64), 4)
    assert four.min() >= 1
    assert np.all(four.sum(axis=1) == 60**4)


def test_axial_radial_moments_use_actual_records_for_every_initial_state():
    kernel, _, _, _ = c.source()
    weights = c.axial_weights()
    norms = [sp.simplify(weights[:, edge].dot(weights[:, edge])) for edge in range(30)]
    squared = []
    fourth = []
    for state in kernel.states:
        moment2 = moment4 = sp.Rational(0)
        for edge, seam in enumerate(kernel.seams):
            for coin in (0, 1):
                _, record = c.native.native(state, seam, coin)
                norm2 = norms[edge]*record[1]**2
                assert norm2 in (0, sp.Rational(4, 5))
                moment2 += norm2/60
                moment4 += norm2**2/60
        squared.append(moment2); fourth.append(moment4)
    assert set(squared) == {sp.Rational(8, 75), sp.Rational(2, 15)}
    assert sum(squared)/66 == sp.Rational(4, 33)
    assert sum(fourth)/66 == sp.Rational(16, 165)


def test_axial_radial_certificate_rejects_nonmartingale_readout(monkeypatch):
    weights = sp.MutableDenseMatrix(c.axial_weights())
    weights[0, 0] += 1
    monkeypatch.setattr(c, 'axial_weights', lambda: weights)
    with pytest.raises(ArithmeticError, match='conditional drift'):
        c.axial_radial_limit_control()


def test_gradient_contamination_destroys_the_axial_martingale(monkeypatch):
    weights=sp.MutableDenseMatrix(c.axial_weights())
    weights[0,0]+=1
    monkeypatch.setattr(c,'axial_weights',lambda:weights)
    with pytest.raises(ArithmeticError,match='martingale difference'):
        c.axial_geometry()


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
