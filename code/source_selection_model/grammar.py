"""Complete affine state coordinates and every omitted-coordinate witness.

Sparse entries use (central block, row, column, real/imaginary part).
The positivity argument is analytic: every witness below is a Hermitian
matrix unit or a difference of two diagonal units, with operator norm one.
"""
from fractions import Fraction as Q


def coordinates():
    out = [[p, a, a, 0] for p in range(12) for a in range(6)]
    out.pop()  # normalization supplies the anchor diagonal
    out.extend([p, a, b, part] for p in range(12) for a in range(6)
               for b in range(a+1, 6) for part in range(2))
    return out


def verify_coordinates(rows):
    if type(rows) is not list or len(rows) != 431:
        raise ValueError("complete affine grammar count")
    expected = {(p, a, b, part) for p in range(12) for a in range(6)
                for b in range(a, 6) for part in range(2)
                if (a != b or part == 0) and (p, a, b, part) != (11, 5, 5, 0)}
    actual = set()
    for row in rows:
        if type(row) is not list or len(row) != 4 or any(type(x) is not int for x in row):
            raise ValueError("coordinate schema")
        key = tuple(row)
        if key not in expected or key in actual:
            raise ValueError("omitted or duplicated coordinate")
        actual.add(key)
    if actual != expected:
        raise ValueError("coordinate coverage")
    # Each coordinate gets a distinct admissible invisible intervention after
    # it is deleted. The checker evaluates every other coordinate, including
    # normalization, rather than trusting a stored rank or True flag.
    for p, a, b, part in rows:
        key = p, a, b, part
        perturbation = {key: Q(1)}
        if a == b:
            perturbation[11, 5, 5, 0] = Q(-1)
        else:
            perturbation[p, b, a, part] = Q(-1 if part == 1 else 1)
        if sum(v for (_, i, j, k), v in perturbation.items() if i == j and k == 0) != 0:
            raise ValueError("intervention trace")
        if any(perturbation.get(tuple(other), 0) != int(tuple(other) == key) for other in rows):
            raise ValueError("omission intervention visibility")
        if any(perturbation.get((q, j, i, k), 0) != v*(-1 if k else 1)
               for (q, i, j, k), v in perturbation.items()):
            raise ValueError("intervention Hermiticity")
    return {"normalized_state_dimension": len(rows), "omission_witnesses": len(rows),
            "reference_eigenvalue": "1/72", "perturbation_size": "1/288",
            "perturbed_eigenvalue_lower_bound": str(Q(1, 72)-Q(1, 288))}


def controls():
    return verify_coordinates(coordinates())
