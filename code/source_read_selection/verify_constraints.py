"""Independent rational certificates for full-family KL projections and covers."""
from fractions import Fraction as F
from itertools import permutations

from . import codec


def require(ok,label):
    if not ok:
        raise ValueError(label)


def rationals(values):
    require(type(values) is list,"rational vector")
    result = []
    for value in values:
        require(type(value) is str,"canonical rational")
        try:
            parsed = F(value)
        except (ValueError,ZeroDivisionError) as error:
            raise ValueError("canonical rational") from error
        require(str(parsed)==value,"canonical rational")
        result.append(parsed)
    return result


def distribution(values,faithful=False):
    p = rationals(values)
    require(bool(p) and sum(p)==1 and all(x>=0 for x in p),"probability normalization")
    if faithful:
        require(all(x>0 for x in p),"faithful reference or selected law required")
    return p


def certify_projection(packet):
    r = distribution(packet["reference"],True)
    p = distribution(packet["selected"],True)
    require(len(p)==len(r),"projection dimensions")
    groups = packet["groups"]
    targets = distribution(packet["targets"],True)
    ratios = rationals(packet["ratios"])
    require(type(groups) is list and len(groups)==len(p),"group dimensions")
    require(all(type(g) is int and 0<=g<len(targets) for g in groups),"group index")
    require(set(groups)==set(range(len(targets))) and len(ratios)==len(targets),"complete partition")
    # An independent exact stationarity/Pythagorean certificate: log(p/r)
    # is constant on each group and every group has its prescribed mass.
    # Hence its pairing with every feasible q-p vanishes. Lean proves the
    # resulting KL identity and unique minimum for the entire affine family.
    for group in range(len(targets)):
        indices = [i for i,g in enumerate(groups) if g==group]
        require(sum(p[i] for i in indices)==targets[group],"affine group constraint")
        require(all(p[i]/r[i]==ratios[group] for i in indices),"KL moment certificate")
    return p,r


def determinant(matrix):
    size = len(matrix)
    require(all(len(row)==size for row in matrix),"square matrix")
    answer = F(0)
    for permutation in permutations(range(size)):
        sign = (-1)**sum(permutation[i]>permutation[j] for i in range(size) for j in range(i+1,size))
        term = F(sign)
        for i,j in enumerate(permutation):
            term *= matrix[i][j]
        answer += term
    return answer


def image(matrix,vector):
    return [sum(F(c)*x for c,x in zip(row,vector)) for row in matrix]


def score(p,r,denominator):
    answer = F(1)
    for x,y in zip(p,r):
        exponent = x*denominator
        require(exponent.denominator==1,"integral log-score exponent")
        if x:
            answer *= (x/y)**exponent.numerator
    return answer


def verify(packet):
    require(type(packet) is dict,"projection controls object")
    word = packet["correlated_words"]
    p,r = certify_projection(word)
    codec.equal(word["reference"],["1/8","3/8","1/4","1/4"],"word reference")
    codec.equal(word["groups"],[0,0,1,1],"word constraint grammar")
    codec.equal(word["targets"],["3/4","1/4"],"word constraint values")
    coefficients = []
    for i in range(4):
        a,b = divmod(i,2)
        state = [F(1),F(0),F(0)]
        for left in (a,b):
            mean = sum(state[left:left+2])/2
            state[left:left+2] = [mean,mean]
        coefficients.append(state[2])
    failure = sum(mass for mass,gain in zip(p,coefficients) if gain==0)
    expected_word = {"reference":list(map(str,r)),"selected":list(map(str,p)),"groups":[0,0,1,1],
                     "targets":["3/4","1/4"],"ratios":[str(p[0]/r[0]),str(p[2]/r[2])],
                     "word_order":[[0,0],[0,1],[1,0],[1,1]],"alphabet":[[0,1],[1,2]],
                     "receiver_coefficients":list(map(str,coefficients)),"failure_mass":str(failure),
                     "cut_avoidance_mass":str(p[3]),"product_determinant":str(p[0]*p[3]-p[1]*p[2])}
    codec.equal(word,expected_word,"correlated execution")
    cover = packet["same_reference_cover_comparison"]
    require(set(cover["full_projection"])=={"reference","selected","groups","targets","ratios"},
            "projection fields")
    full,r3 = certify_projection(cover["full_projection"])
    codec.equal(cover["full_projection"]["reference"],["1/4","1/4","1/2"],"history reference")
    codec.equal(cover["full_projection"]["groups"],[0,1,0],"history constraint grammar")
    codec.equal(cover["full_projection"]["targets"],["1/2","1/2"],"history constraint values")
    coarse = [[1,1,0],[0,0,1]]
    selected = [F(0),F(1,2),F(1,2)]
    witness = [F(1,4),F(1,2),F(1,4)]
    cref = image(coarse,r3)
    require(all(x>0 for x in cref) and image(coarse,selected)==cref,"coarse KL zero certificate")
    det = determinant([[0,1,0],*coarse])
    require(det!=0,"cover fails to determine the constrained state")
    require(selected[0]==0<full[0] and witness[0]>0,"support contrast")
    # Identity of linear forms on the affine constraint p1=1/2:
    # p0 = (p0+p1) - 1/2. The negative offset is indispensable here.
    reconstruction = {"coefficients":["1","0"],"offset":"-1/2"}
    require(image(coarse,selected)[0]-F(1,2)==selected[0] and
            image(coarse,witness)[0]-F(1,2)==witness[0],"failure readout")
    expected_cover = {"history_reference":list(map(str,r3)),"middle_mass":"1/2",
                      "coarse_matrix":coarse,"full_matrix":[[1,0,0],[0,1,0],[0,0,1]],
                      "coarse_reference":list(map(str,cref)),"coarse_selected":list(map(str,selected)),
                      "full_projection":cover["full_projection"],"feasible_positive_witness":list(map(str,witness)),
                      "constraint_and_coarse_determinant":int(det),"failure_atom":0,
                      "coarse_failure_reconstruction":reconstruction,
                      "coarse_selected_readout":list(map(str,image(coarse,selected))),
                      "coarse_witness_readout":list(map(str,image(coarse,witness))),
                      "full_failure_coefficients":["1","0","0"]}
    codec.equal(cover,expected_cover,"observer-cover certificate")
    nc = packet["nonconvex_control"]
    nr = [F(7,10),F(1,10),F(1,10),F(1,10)]
    vertex,interior = [F(1),F(0),F(0),F(0)],[F(1,4)]*4
    ratio = score(interior,nr,4)/score(vertex,nr,4)
    require(ratio>1,"nonconvex strict minimum")
    midpoint = [(x+y)/2 for x,y in zip(vertex,interior)]
    require(midpoint not in (vertex,interior),"nonconvex control accidentally convex")
    expected_nc = {"reference":list(map(str,nr)),"feasible":[list(map(str,vertex)),list(map(str,interior))],
                   "selected":list(map(str,vertex)),"log_score_denominator":4,
                   "competitor_score_ratio":str(ratio),"convex":False}
    codec.equal(nc,expected_nc,"nonconvex certificate")
    expected = {"scope":"finite_classical_controls_constraints_reference_and_cover_supplied",
                "correlated_words":expected_word,"same_reference_cover_comparison":expected_cover,
                "nonconvex_control":expected_nc,
                "successful_face":{"reference":list(map(str,r)),"forced_selected":["0","1","0","0"],
                                   "excluded_atoms":[0,2,3],"status":"read_success_is_an_added_constraint"}}
    codec.equal(packet,expected,"complete constrained packet")
    return {"correlated_failure_mass":str(failure),"coarse_selected_failure_mass":"0",
            "full_selected_failure_mass":str(full[0]),"nonconvex_log_score_ratio":str(ratio),
            "support_criterion":"scored_atoms_or_nonnegative_scored_readout_only"}
