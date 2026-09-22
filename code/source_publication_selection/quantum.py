"""Exact matrix controls for the analytic A3 support and process arguments."""
from fractions import Fraction as F

from . import codec


def render(matrix):
    return [[str(x) for x in row] for row in matrix]


def product(left, right):
    return [[sum(a*b for a, b in zip(row, col)) for col in zip(*right)] for row in left]


def case(spec):
    a, b = (F(x, spec[2]) for x in spec[:2])
    r = [[F(1), F(0)], [F(0), F(0)]]
    s = [[a*a, a*b], [a*b, b*b]]
    m = [[(r[i][j]+s[i][j])/2 for j in range(2)] for i in range(2)]
    rs, sr = product(r, s), product(s, r)
    mass, p = b*b, b*b/2
    # Clear rational exponents in the general probability floor, F(sigma)=log 2.
    numerator, denominator = mass.numerator, mass.denominator
    floor_holds = p**numerator >= F(1, 2)**denominator * mass**numerator * (1-mass)**(denominator-numerator)
    return {"unit_vector": [str(a), str(b)], "left": render(r), "right": render(s),
            "selected": render(m), "swap": render([[a, b], [b, -a]]),
            "commutator": render([[rs[i][j]-sr[i][j] for j in range(2)] for i in range(2)]),
            "selected_determinant": str(m[0][0]*m[1][1]-m[0][1]*m[1][0]),
            "selected_eigenvalues": [str((1+a)/2), str((1-a)/2)],
            "witness_failure": str(mass), "selected_failure": str(p),
            "rational_power_floor": floor_holds}


def channels():
    zero = F(0)
    identity = [[F(1), zero], [zero, F(1)]]
    flip = [[zero, F(1)], [F(1), zero]]
    units = [[[F(int(i==a and j==b)) for j in range(2)] for i in range(2)]
             for a in range(2) for b in range(2)]
    families = {"identity": [(F(1), identity)],
                "depolarizing": [(F(1, 2), x) for x in units],
                "dephasing": [(F(1), units[0]), (F(1), units[3])],
                "flip": [(F(1), flip)]}
    output = []
    for name, operators in families.items():
        choi = [[F(0) for _ in range(4)] for _ in range(4)]
        actions = [[[F(0) for _ in range(2)] for _ in range(2)] for _ in units]
        for weight, k in operators:
            vector = [k[j][i] for i in range(2) for j in range(2)]
            for i in range(4):
                for j in range(4):
                    choi[i][j] += weight*vector[i]*vector[j]/2
            for index, x in enumerate(units):
                image = product(product(k, x), list(zip(*k)))
                for i in range(2):
                    for j in range(2):
                        actions[index][i][j] += weight*image[i][j]
        output.append({"name": name, "choi": render(choi),
                       "matrix_unit_images": [render(x) for x in actions],
                       "uniform_bit_joint": [[str(choi[2*i+j][2*i+j]) for j in range(2)] for i in range(2)],
                       "mismatch": str(choi[1][1]+choi[2][2])})
    return output


def packet():
    return {"segments": [case(spec) for spec in codec.QUANTUM_CASES],
            "face": {"selected": ["1/2", "1/2", "0"], "forbidden_mass": "0", "rank": 2},
            "singular_reference": {"reference": ["1", "0"], "finite_face_rank": 1,
                                   "discarded_state_failure": "1"},
            "channels": channels(), "unconstrained_process_optimizer": "depolarizing",
            "zero_mismatch_process_optimizer": "dephasing",
            "flip_witness_failure_floor": "1/4"}
