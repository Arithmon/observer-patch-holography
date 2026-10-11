"""Independent exhaustive rational response review; finite declared event menu.

Target 0510578e317de5da0d485d09792dc8bcad2e2f3c. No random sampling.
Reconstructs the icosahedral edge menu from its fixed incidence seed, then
enumerates 66 two-occupancy states x 30 edges x 2 fair-coin outcomes.
This is a review of a named binary process, never full A1-A3 membership.
"""
import hashlib
import json
import platform
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path

SHA = "0510578e317de5da0d485d09792dc8bcad2e2f3c"
FACES = [(0,1,2),(0,3,1),(0,2,4),(0,6,3),(0,4,6),
         (1,5,2),(1,3,7),(1,7,5),(2,8,4),(2,5,8),
         (3,6,9),(3,9,7),(4,10,6),(4,8,10),(5,7,11),
         (5,11,8),(6,10,9),(7,9,11),(8,11,10),(9,10,11)]
EDGES = sorted({tuple(sorted(e)) for f in FACES for e in combinations(f, 2)})
STATES = list(combinations(range(12), 2))
N = 12
assert len(EDGES) == 30 and len(STATES) == 66
lap = [[Q(0) for _ in range(N)] for _ in range(N)]
for a,b in EDGES:
    lap[a][a] += 1; lap[b][b] += 1
    lap[a][b] -= 1; lap[b][a] -= 1


def event_table(weights, unit=Q(1), activity=Q(1)):
    """Global attempt includes retained independent idle branch 1-activity.

    On active attempts, fair local coin chooses either endpoint as the final
    occupied endpoint of an odd edge; one coin preserves the state, one swaps.
    J is signed occupancy transport, in the selected load unit.
    """
    assert sum(weights) == 1 and all(w > 0 for w in weights)
    covariance = [[Q(0) for _ in range(N)] for _ in range(N)]
    triangle_variance = Q(0)
    move_fraction = Q(0)
    triangle = {tuple(sorted(e)) for e in combinations(FACES[0], 2)}
    circulation = {e: (1 if e in [(0,1),(1,2)] else -1) for e in triangle}
    transition = [[Q(0) for _ in STATES] for _ in STATES]
    drift = [[Q(0) for _ in range(N)] for _ in STATES]
    for si,state in enumerate(STATES):
        transition[si][si] += 1-activity
        for (a,b), w in zip(EDGES, weights):
            odd = (a in state) != (b in state)
            for coin in (0,1):
                probability = activity*w/2
                following = state
                current = Q(0)
                if odd and coin:
                    old,new = (a,b) if a in state else (b,a)
                    following = tuple(sorted((set(state)-{old})|{new}))
                    current = unit*(1 if old == a else -1)
                transition[si][STATES.index(following)] += probability
                if current:
                    move_fraction += probability/66
                    covariance[a][a] += probability*current**2/66
                    covariance[b][b] += probability*current**2/66
                    covariance[a][b] -= probability*current**2/66
                    covariance[b][a] -= probability*current**2/66
                    drift[si][a] -= probability*current
                    drift[si][b] += probability*current
                    triangle_variance += probability*(circulation.get((a,b),0)*current)**2/66
    assert all(sum(row) == 1 for row in transition)
    assert all(sum(transition[i][j] for i in range(66)) == 1 for j in range(66))
    assert all(transition[i][j] == transition[j][i] for i in range(66) for j in range(66))
    return covariance,triangle_variance,move_fraction,drift


def multiple(matrix, base, factor):
    return all(matrix[i][j] == factor*base[i][j] for i in range(N) for j in range(N))


def main():
    B,nu,kappa,drift = event_table([Q(1,30)]*30)
    assert multiple(B,lap,Q(1,198)) and nu == Q(1,66) and kappa == Q(5,33)
    for state,row in zip(STATES,drift):
        assert row == [-sum(lap[i][j]*int(j in state) for j in range(N))/60 for i in range(N)]
    out = {"target_revision":SHA,"python":platform.python_version(),
           "event_count":66*30*2,"stationary_state":"uniform binary 2-occupancy shell",
           "seam_schedule":"one uniformly selected edge per total attempt; fair local coin",
           "B_coefficient":"1/198","nu":"1/66","dimensionless_kappa":"5/33",
           "mean_drift_coefficient":"-1/60",
           "ratio_slow_mode_coefficient":"1/6 times (5-sqrt(5))"}
    B2,nu2,k2,_ = event_table([Q(1,30)]*30,unit=Q(2))
    assert multiple(B2,B,4) and nu2 == 4*nu and k2 == kappa
    # gamma is drift rate of the dimensionless occupancy state, unaffected by J units.
    ratio_scale = (k2/nu2)/(kappa/nu)
    assert ratio_scale == Q(1,4)
    squared_activity = 4*k2
    assert (squared_activity/nu2)/(kappa/nu) == 1
    assert 3*Q(4,198)/nu2 == 1
    out["current_scale_2"]={"B_coefficient":"2/99","nu":"2/33",
       "kappa":"5/33","kappa_J":"20/33","L_record_scale":"1",
       "dimensionless_kappa_ratio_scale":"1/4","squared_current_kappa_ratio_scale":"1"}
    Bi,nui,ki,_ = event_table([Q(1,30)]*30,activity=Q(1,2))
    assert multiple(Bi,B,Q(1,2)) and nui == nu/2 and ki == kappa/2
    # Per total attempt gamma also halves, so remeasured kappa produces half the ratio.
    assert ((ki/nui)*Q(1,2))/(kappa/nu) == Q(1,2)
    out["idle_half"]={"B_coefficient":"1/396","nu":"1/132","kappa":"5/66",
        "gamma_scale":"1/2","ratio_remeasured_kappa_scale":"1/2",
        "ratio_frozen_native_kappa_scale":"1"}
    weights = [Q(2 if e == EDGES[0] else 1,31) for e in EDGES]
    Bw,nuw,kw,_ = event_table(weights)
    assert kw == kappa and not multiple(Bw,lap,Bw[0][1]/lap[0][1])
    out["nonuniform_schedule"]={"doubled_edge":EDGES[0],"edge_weights":"2/31 on one edge; 1/31 on other 29",
       "uniform_stationarity_preserved":True,"kappa_preserved":True,
       "B_proportional_to_native_L":False,"nu":str(nuw),
       "consequence":"same stationary shell does not fix geometric response or relaxation under a changed seam schedule"}
    out["scope"]="Exact arithmetic and named-process dependency review; no proof of A1-A3 membership or all-realization exclusion."
    out["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
