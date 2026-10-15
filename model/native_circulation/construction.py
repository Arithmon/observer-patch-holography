"""Exact current counters and a single-history Gaussian limit.

Native attempts index the history. The uniform scheduler, binary total-two
sector, stationary preparation and chosen cycle readout are declared.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from ..observer import native
from ..primitive_observer import runtime, transfer
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
VALUES = (-1, 0, 1)
SOURCE_PATHS = tuple(dict.fromkeys((
    'model/native_circulation/__init__.py', 'model/native_circulation/__main__.py',
    'model/native_circulation/construction.py', 'model/common_source/protocol.py',
    *transfer.SOURCE_PATHS,
    'Lean/Thermodynamics/GreenKubo.lean',
    'Lean/Thermodynamics/NativeRepairFourCycle.lean',
)))


def bounded(value, name):
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 2e-11:
        raise ArithmeticError(f'{name} failed numerical comparison')
    return value


@lru_cache(maxsize=1)
def source():
    kernel = transfer.native_counts()
    transfer.exact_audit(kernel)
    index = {state: i for i, state in enumerate(kernel.states)}
    edges = tuple(kernel.seams)
    # Apparatus designation uses only the source graph, never measurement data.
    triangle = next(t for t in combinations(range(12), 3)
                    if all(tuple(sorted(e)) in edges for e in combinations(t, 2)))
    oriented = tuple(zip(triangle, triangle[1:]+triangle[:1]))
    weights = tuple(sum(int((a,b)==pair)-int((b,a)==pair) for pair in oriented)
                    for a,b in edges)
    branches = []
    for e,seam in enumerate(edges):
        for coin in (0,1):
            for i,state in enumerate(kernel.states):
                target,record = native.native(state,seam,coin)
                if record != (coin,state[seam[0]]-target[seam[0]]):
                    raise ArithmeticError('native record grammar changed')
                current = record[1]
                branches.append((i,index[target],e,coin,current,weights[e]*current))
    return kernel, triangle, weights, tuple(branches)


def exact_geometry():
    kernel,triangle,weights,branches = source()
    b=sp.zeros(12,30)
    for e,(a,c) in enumerate(kernel.seams):
        b[a,e]=-1; b[c,e]=1
    lap=b*b.T; j=sp.ones(12)/12
    inverse=(lap+j).inv()-j
    projection=sp.eye(30)-b.T*inverse*b
    if b*projection != sp.zeros(12,30) or projection*projection != projection or projection.T != projection:
        raise ArithmeticError('cycle/gradient decomposition failed')
    if sp.trace(projection)!=19 or b*sp.Matrix(weights)!=sp.zeros(12,1):
        raise ArithmeticError('the source cycle sector changed')
    # A current changes loads by B j. Its conditional mean is -B^T x/60.
    conditional=[[F(0) for _ in range(30)] for _ in kernel.states]
    second=[[F(0) for _ in range(30)] for _ in range(30)]
    for i,out,e,coin,current,g in branches:
        conditional[i][e]+=F(current,60)
        second[e][e]+=F(current*current,3960)
        if tuple(kernel.states[out][v]-kernel.states[i][v] for v in range(12)) != tuple(int(b[v,e])*current for v in range(12)):
            raise ArithmeticError('current does not reproduce the actual load change')
    for i,state in enumerate(kernel.states):
        expected=-(b.T*sp.Matrix(state))/60
        if sp.Matrix(conditional[i])!=expected or projection*expected!=sp.zeros(30,1):
            raise ArithmeticError('projected currents are not martingale differences')
    if sp.Matrix(second)!=sp.eye(30)/198:
        raise ArithmeticError('stationary current covariance changed')
    sigma=projection/198
    return {'chosen_triangle':list(triangle), 'oriented_edge_weights':list(weights),
            'incidence_rank':int(b.rank()), 'cycle_rank':19,
            'cycle_projection':[[str(x) for x in row] for row in projection.tolist()],
            'asymptotic_covariance':[[str(x) for x in row] for row in sigma.tolist()],
            'all_native_branch_continuity_checks':len(branches),
            'projected_conditional_mean_exactly_zero_on_all_states':True,
            'stationary_single_edge_second_moment':'1/198',
            'native_count_parameter':'integer uniform native attempts',
            'spatial_dimension_inferred_from_cycle_rank':False}


def marked():
    kernel,_,_,branches=source()
    matrices=[[[0]*66 for _ in range(66)] for _ in VALUES]
    for i,j,_,_,_,g in branches:
        matrices[g+1][i][j]+=1
    if tuple(tuple(sum(matrices[k][i][j] for k in range(3)) for j in range(66)) for i in range(66))!=kernel.counts:
        raise ArithmeticError('marked current maps changed the native transition')
    for i in range(66):
        if sum(g*sum(matrices[g+1][i]) for g in VALUES)!=0:
            raise ArithmeticError('cycle current has a conditional drift')
    return tuple(tuple(tuple(row) for row in m) for m in matrices)


def instrument():
    _,_,weights,branches=source()
    groups={}
    for i,j,e,coin,current,g in branches:
        key=(e,coin,current)
        groups.setdefault(key,np.zeros((66,66),complex))[j,i]=1/math.sqrt(60)
    return runtime.Instrument(kraus=tuple(tuple(groups[key] for key in sorted(groups)
        if weights[key[0]]*key[2]==g) for g in VALUES),visible=True,
        resources=(('native_attempts',1),('circulation_reads',1)),
        binding={'candidate':'C','observable':'oriented triangular circulation from seam and native current',
                 'values':list(VALUES),'all_hidden_labels':'probability summed',
                 'native_channel_semantics':True,'native_pulse_words_emitted':False})


def horizon(value):
    if type(value) is not int or not 0<=value<=512:
        raise ValueError('integer attempt count from zero through512 required')
    return value


def count_laws(horizons=(1,4,16,64,256)):
    """Exact retained counter and terminal-state law; no trajectory cutoff."""
    horizons=tuple(horizon(n) for n in horizons)
    if not horizons or len(set(horizons))!=len(horizons):
        raise ValueError('distinct nonempty attempt counts required')
    m=marked()
    sparse=[[(i,j,count) for i,row in enumerate(matrix) for j,count in enumerate(row) if count]
            for matrix in m]
    rows=[[1]*66] # all66 stationary inputs, integer denominator66
    denominator=66
    saved={}
    for n in range(max(horizons)+1):
        if n in horizons:
            saved[n]=(tuple(tuple(row) for row in rows),denominator)
        if n==max(horizons): break
        updated=[[0]*66 for _ in range(len(rows)+2)]
        for g,transitions in enumerate(sparse):
            for s,row in enumerate(rows):
                out=updated[s+g]
                for i,j,c in transitions:
                    if row[i]: out[j]+=row[i]*c
        rows=updated; denominator*=60
        if sum(map(sum,rows))!=denominator:
            raise ArithmeticError('counter coarsening lost a history or terminal state')
    return saved


def runtime_check():
    operation=instrument()
    program=runtime.Program('one_observer_triangle_current',66,np.eye(66)/66,
        tuple(runtime.Event(f'cycle_read_{n}',instruments=(operation,)) for n in range(3)),
        binding={'stationary_preparation':'I66/66','parameter':'native attempts'})
    output=runtime.execute(program)
    matrices=marked(); residual=0.; records=[]
    aggregate=defaultdict(lambda:np.zeros((66,66),complex))
    for branch in output.branches:
        counts=[1]*66
        word=tuple(VALUES[k] for k in branch.record)
        for k in branch.record:
            counts=[sum(counts[i]*matrices[k][i][j] for i in range(66)) for j in range(66)]
        exact=[F(n,66*60**3) for n in counts]
        expected=np.diag(list(map(float,exact)))
        residual=max(residual,float(np.linalg.norm(branch.density-expected,'fro')))
        aggregate[sum(word)]+=branch.density
        if dict(branch.resources)!={'native_attempts':3,'circulation_reads':3}:
            raise ArithmeticError('recorded circulation resource accounting changed')
        records.append({'word':list(word),'probability':str(sum(exact))})
    law,denominator=count_laws((3,))[3]
    for s,row in enumerate(law):
        expected=np.diag([float(F(value,denominator)) for value in row])
        residual=max(residual,float(np.linalg.norm(aggregate[s-3]-expected,'fro')))
    prefix=runtime.execute(program,stop_after=1)
    restored=runtime.checkpoint_from_json(program,runtime.checkpoint_to_json(prefix))
    continued=runtime.execute(program,checkpoint=restored)
    return {'complete_records':records,'full_terminal_density_dimension':66,
            'maximum_complete_output_residual':bounded(residual,'native current output'),
            'checkpoint_complete_distance':bounded(runtime.complete_distance(output,continued),'checkpoint'),
            'exact_counter_coarsening_checked_against_all_runtime_records':True}


def run():
    geometry=exact_geometry(); laws=count_laws(); summaries=[]
    for n,(rows,denominator) in sorted(laws.items()):
        probabilities=[F(sum(row),denominator) for row in rows]
        moments=[sum((F(s-n)**k*p for s,p in enumerate(probabilities)),F(0)) for k in range(1,5)]
        if moments[0]!=0 or moments[1]!=F(n,66) or moments[2]!=0:
            raise ArithmeticError('native exact counter moments disagree with the martingale identity')
        variance=moments[1]
        empirical=sum(float(p)*math.cos((s-n)/math.sqrt(float(variance))) for s,p in enumerate(probabilities))
        summaries.append({'attempts':n,'complete_counter_law':[str(p) for p in probabilities],
            'score_min':-n,'score_max':n,'mean':str(moments[0]),'variance':str(variance),
            'fourth_moment':str(moments[3]),'excess_kurtosis':str(moments[3]/variance**2-3),
            'standardized_characteristic_at_one_diagnostic':empirical,
            'Gaussian_characteristic_at_one':math.exp(-.5),
            'terminal_counter_joint_integer_counts_sha256':hashlib.sha256(canonical(rows).encode()).hexdigest(),
            'all_scores_and_terminal_states_retained':True})
    m=marked()
    square_correlation=F(0)
    for i in range(66):
        for j in range(66):
            square_correlation+=F((m[0][i][j]+m[2][i][j])*(sum(m[0][j])+sum(m[2][j])),66*60**2)
    if square_correlation!=F(1,1980):
        raise ArithmeticError('the nonindependent squared-current control changed')
    return {'schema':'oph.native_circulation.v1','all_pass':True,
        'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'source_geometry':geometry,'complete_native_execution':runtime_check(),
        'one_history_counter_predictions':summaries,
        'temporal_dependence_control':{'successive_signed_covariance':'0',
            'successive_squared_product':str(square_correlation),'independent_squared_product':str(F(1,66)**2),
            'increments_independent':False},
        'limit':{'triangle_variance_per_native_attempt':'1/66',
            'vector_covariance':'(I-B^T L^+ B)/198','stationary_finite_ergodic_martingale_limit':True,
            'independent_replicas_required':False,'new_Lean_CLT':False,
            'uniform_Berry_Esseen_or_trajectory_error_certificate':False},
        'scope':{'stationary_binary_preparation_uniform_scheduler_and_cycle_readout_supplied':True,
            'physical_clock_spatial_dimension_or_interacting_QFT_derived':False,
            'native_quantum_C_E_selected':False,'new_empirical_agreement':False}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--write',type=Path);mode.add_argument('--check',type=Path)
    args=parser.parse_args(); result=run()
    if args.check:
        compare_receipt(strict_json.load(args.check),result);print('Native circulation receipt reproduced')
    elif args.write:
        with args.write.open('x',encoding='utf-8') as stream:stream.write(canonical(result))
        print(f'Wrote {args.write}')
    else: print(canonical(result),end='')


if __name__=='__main__':main()
