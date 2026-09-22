"""Independent scalar matrix replay; no producer, optimizer or floating point."""
from fractions import Fraction as F

from . import codec


def matrix(rows):
    return tuple(tuple(F(x) for x in row) for row in rows)


def multiply(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def stringify(a):
    return [[str(x) for x in row] for row in a]


def assert_equal(a, b, name):
    if a != b:
        raise ValueError(name)


def segment(spec):
    a, b = F(spec[0], spec[2]), F(spec[1], spec[2])
    if not (0 < a < 1 and 0 < b < 1 and a*a+b*b == 1):
        raise ValueError("rational unit vector")
    basis = matrix([[1, 0], [0, 0]])
    swap = matrix([[a, b], [b, -a]])
    ident = matrix([[1, 0], [0, 1]])
    assert_equal(multiply(swap, swap), ident, "involution")
    # Construct the second state by conjugating the first, not its outer formula.
    second = multiply(multiply(swap, basis), swap)
    assert_equal(multiply(second, second), second, "pure state")
    assert_equal(sum(second[i][i] for i in range(2)), 1, "normalization")
    midpoint = tuple(tuple((basis[i][j]+second[i][j])/2 for j in range(2)) for i in range(2))
    assert_equal(multiply(multiply(swap, midpoint), swap), midpoint, "optimizer symmetry")
    det = midpoint[0][0]*midpoint[1][1]-midpoint[0][1]*midpoint[1][0]
    if not (midpoint[0][0] > 0 and det > 0):
        raise ValueError("positive definite midpoint")
    # Verify characteristic roots from trace and product, without eigensolvers.
    eigenvalues = ((1+a)/2, (1-a)/2)
    assert_equal(sum(eigenvalues), sum(midpoint[i][i] for i in range(2)), "spectral trace")
    assert_equal(eigenvalues[0]*eigenvalues[1], det, "spectral product")
    commutator = tuple(tuple(sum(basis[i][k]*second[k][j]-second[i][k]*basis[k][j]
                                  for k in range(2)) for j in range(2)) for i in range(2))
    if not any(x for row in commutator for x in row):
        raise ValueError("commuting control")
    bad = second[1][1]
    n, d = bad.numerator, bad.denominator
    # Equivalent integer inequality after all denominators are cleared.
    lhs = midpoint[1][1]**n
    rhs = bad**n * (1-bad)**(d-n) / 2**d
    floor = lhs.numerator*rhs.denominator >= rhs.numerator*lhs.denominator
    if not floor:
        raise ValueError("quantitative failure floor")
    return {"unit_vector": [str(a), str(b)], "left": stringify(basis), "right": stringify(second),
            "selected": stringify(midpoint), "swap": stringify(swap),
            "commutator": stringify(commutator), "selected_determinant": str(det),
            "selected_eigenvalues": list(map(str, eigenvalues)), "witness_failure": str(bad),
            "selected_failure": str(midpoint[1][1]), "rational_power_floor": floor}


def process(name):
    images = []
    for a in range(2):
        for b in range(2):
            if name == "identity":
                out = [[F(int(i==a and j==b)) for j in range(2)] for i in range(2)]
            elif name == "flip":
                out = [[F(int(i==1-a and j==1-b)) for j in range(2)] for i in range(2)]
            elif name == "depolarizing":
                out = [[F(int(a==b and i==j), 2) for j in range(2)] for i in range(2)]
            elif name == "dephasing":
                out = [[F(int(a==b and i==j==a)) for j in range(2)] for i in range(2)]
            else:
                raise ValueError("channel")
            images.append(matrix(out))
    choi = tuple(tuple(images[2*(i//2)+(j//2)][i%2][j%2]/2 for j in range(4)) for i in range(4))
    for a in range(2):
        for b in range(2):
            assert_equal(sum(choi[2*a+j][2*b+j] for j in range(2)), F(int(a==b), 2), "trace preservation")
            assert_equal(sum(choi[2*j+a][2*j+b] for j in range(2)), F(int(a==b), 2), "unitality")
    # Exact Choi positivity: orthogonal support projectors with positive eigenvalues.
    eigenvalue = F(1, {"identity": 1, "flip": 1, "depolarizing": 4, "dephasing": 2}[name])
    assert_equal(multiply(choi, choi), tuple(tuple(eigenvalue*x for x in row) for row in choi), "Choi projector spectrum")
    assert_equal(choi, tuple(zip(*choi)), "Choi self-adjointness")
    assert_equal(sum(choi[i][i] for i in range(4)), 1, "Choi normalization")
    return {"name": name, "choi": stringify(choi), "matrix_unit_images": [stringify(x) for x in images],
            "uniform_bit_joint": [[str(choi[2*i+j][2*i+j]) for j in range(2)] for i in range(2)],
            "mismatch": str(choi[1][1]+choi[2][2])}


def expected():
    channels = [process(name) for name in ("identity", "depolarizing", "dephasing", "flip")]
    # The full channel reference is feasible; restricting it to the zero-mismatch
    # support gives diag(1/2,0,0,1/2), also trace preserving. Positivity forces
    # every zero-mismatch state's support into this same rank-two subspace.
    selected = matrix(channels[1]["choi"])
    good = (0, 3)
    mass = sum(selected[i][i] for i in good)
    restricted = tuple(tuple(selected[i][j]/mass if i in good and j in good else F(0)
                             for j in range(4)) for i in range(4))
    assert_equal(restricted, matrix(channels[2]["choi"]), "record-face optimizer")
    assert_equal(F(channels[3]["mismatch"]), 1, "deterministic failure witness")
    # exp(-D(pure || I/4))=1/4, versus actual selected mismatch 1/2.
    if F(channels[1]["mismatch"]) < F(1, 4):
        raise ValueError("process witness floor")
    reference3 = (F(1, 3),)*3
    allowed = (True, True, False)
    face_mass = sum(x for x, keep in zip(reference3, allowed) if keep)
    face = [x/face_mass if keep else F(0) for x, keep in zip(reference3, allowed)]
    singular = (F(1), F(0))
    discarded = (F(0), F(1))
    return {"segments": [segment(spec) for spec in codec.QUANTUM_CASES],
            "face": {"selected": list(map(str, face)), "forbidden_mass": str(face[2]),
                     "rank": sum(x>0 for x in face)},
            "singular_reference": {"reference": list(map(str, singular)),
                                   "finite_face_rank": sum(x>0 for x in singular),
                                   "discarded_state_failure": str(discarded[1])},
            "channels": channels, "unconstrained_process_optimizer": "depolarizing",
            "zero_mismatch_process_optimizer": "dephasing",
            "flip_witness_failure_floor": str(F(1, 4))}
