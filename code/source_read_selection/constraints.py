"""Producer for exact constrained information-projection controls."""
from fractions import Fraction as F
from math import prod


def group_projection(reference,groups,targets):
    masses = [sum(r for r,g in zip(reference,groups) if g==j) for j in range(len(targets))]
    ratios = [target/mass for target,mass in zip(targets,masses)]
    selected = [reference[i]*ratios[groups[i]] for i in range(len(reference))]
    return {"reference":list(map(str,reference)),"groups":groups,"targets":list(map(str,targets)),
            "ratios":list(map(str,ratios)),"selected":list(map(str,selected))}


def build():
    word = group_projection([F(1,8),F(3,8),F(1,4),F(1,4)],[0,0,1,1],[F(3,4),F(1,4)])
    word.update({"word_order":[[0,0],[0,1],[1,0],[1,1]],
                 "alphabet":[[0,1],[1,2]],"receiver_coefficients":["0","1/4","0","0"],
                 "failure_mass":"7/16","cut_avoidance_mass":"1/8",
                 "product_determinant":"-3/64"})
    full = group_projection([F(1,4),F(1,4),F(1,2)],[0,1,0],[F(1,2),F(1,2)])
    reference = [F(1,4),F(1,4),F(1,2)]
    selected = [F(0),F(1,2),F(1,2)]
    witness = [F(1,4),F(1,2),F(1,4)]
    coarse = lambda p:[p[0]+p[1],p[2]]
    cover = {"history_reference":list(map(str,reference)),"middle_mass":"1/2",
             "coarse_matrix":[[1,1,0],[0,0,1]],"full_matrix":[[1,0,0],[0,1,0],[0,0,1]],
             "coarse_reference":list(map(str,coarse(reference))),
             "coarse_selected":list(map(str,selected)),"full_projection":full,
             "feasible_positive_witness":list(map(str,witness)),
             "constraint_and_coarse_determinant":-1,
             "failure_atom":0,"coarse_failure_reconstruction":{"coefficients":["1","0"],"offset":"-1/2"},
             "coarse_selected_readout":list(map(str,coarse(selected))),
             "coarse_witness_readout":list(map(str,coarse(witness))),
             "full_failure_coefficients":["1","0","0"]}
    r = [F(7,10),F(1,10),F(1,10),F(1,10)]
    p,q = [F(1),F(0),F(0),F(0)],[F(1,4)]*4
    def score(v):
        return prod((x/y)**int(4*x) for x,y in zip(v,r) if x)
    return {"scope":"finite_classical_controls_constraints_reference_and_cover_supplied",
            "correlated_words":word,"same_reference_cover_comparison":cover,
            "nonconvex_control":{"reference":list(map(str,r)),"feasible":[list(map(str,p)),list(map(str,q))],
                                 "selected":list(map(str,p)),"log_score_denominator":4,
                                 "competitor_score_ratio":str(score(q)/score(p)),"convex":False},
            "successful_face":{"reference":word["reference"],"forced_selected":["0","1","0","0"],
                               "excluded_atoms":[0,2,3],"status":"read_success_is_an_added_constraint"}}
