"""Exact finite controls. Analytic quantifiers are proved in DERIVATION.md."""

from itertools import product

import sympy as sp


def require(condition):
    if not condition:
        raise ValueError("exact RG principle certificate failed")


def clock_certificate():
    phi = (1 + sp.sqrt(5)) / 2
    vertices = set()
    for a, b in product((-1, 1), repeat=2):
        v = (0, a, b * phi)
        for k in range(3):
            vertices.add(v[k:] + v[:k])
    rotations = [sp.diag(1, -1, -1), sp.diag(-1, 1, -1),
                 sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])]
    for rotation in rotations:
        require(rotation.det() == 1 and rotation.T * rotation == sp.eye(3))
        require({tuple(rotation * sp.Matrix(v)) for v in vertices} == vertices)
    a, b, c, d, e, f = sp.symbols("a b c d e f")
    hessian = sp.Matrix([[a, d, e], [d, b, f], [e, f, c]])
    equations = [x for rotation in rotations for x in rotation.T * hessian * rotation - hessian]
    solution = sp.linsolve(equations, (a, b, c, d, e, f))
    require(solution == sp.FiniteSet((c, c, c, 0, 0, 0)))
    # Homogeneous polyhedral norm: unequal normalized values imply nonroundness.
    # A quadratic expansion of a homogeneous squared norm would make it exactly quadratic.
    axis_squared = phi**2
    diagonal_squared = (1 + phi)**2 / 3
    anisotropy = sp.simplify(axis_squared - diagonal_squared)
    require(anisotropy == sp.Rational(1, 3))
    tau, kappa, radius, count = sp.symbols("tau kappa radius count", positive=True)
    refined = sp.expand(count * (tau * radius / count + kappa * (radius / count)**2))
    require(sp.simplify(refined - (tau * radius + kappa * radius**2 / count)) == 0)
    require(sp.limit(refined, count, sp.oo) == tau * radius)
    # Entropy susceptibility in a two-level slice, with natural logarithms.
    x = sp.symbols("x", real=True)
    relative_entropy = (sp.Rational(1, 2) + x) * sp.log(1 + 2*x) + (sp.Rational(1, 2) - x) * sp.log(1 - 2*x)
    require(sp.diff(relative_entropy, x).subs(x, 0) == 0)
    require(sp.diff(relative_entropy, x, 2).subs(x, 0) / 2 == 2)
    return {"ports": len(vertices), "invariant_quadratic_dimension": 1,
            "polyhedral_equal_radius_squared_gap": str(anisotropy),
            "refined_overhead": "kappa*radius^2/count", "entropy_quadratic_coefficient": 2}


def algebra_certificate(d=12):
    if type(d) is not int or not 2 <= d <= 12:
        raise ValueError("reference dimension must be an integer in [2,12]")
    # Sparse Kraus matrix units on C^d tensor M_d: row=(p,p), column=(p,b).
    kraus = [((p, p), (p, b)) for p in range(d) for b in range(d)]
    identity = {(p, b): 0 for p in range(d) for b in range(d)}
    for _, column in kraus:
        identity[column] += 1
    require(all(value == 1 for value in identity.values()))
    checked = 0
    for p, a, b in product(range(d), repeat=3):
        # Evaluate sum K (e_p tensor E_ab) K* from sparse row/column matching.
        output = {}
        for row, column in kraus:
            if column == (p, a) and column == (p, b):
                key = (row, row)
                output[key] = output.get(key, 0) + 1
        expected = {((p, p), (p, p)): 1} if a == b else {}
        require(output == expected)
        checked += 1
    # Off-diagonal buffer effect distinguishes phase; diagonal port copying erases it.
    plus = sp.Matrix([[1, 1], [1, 1]]) / 2
    minus = sp.Matrix([[1, -1], [-1, 1]]) / 2
    require(sp.trace(plus * plus) == 1 and sp.trace(plus * minus) == 0)
    return {"central_ports": d, "kraus_operators": d*d,
            "matrix_units_checked": checked, "trace_preserving": True,
            "phase_effect_probabilities": [1, 0]}


def resource_certificate():
    rows = []
    for n in (2, 3, 4, 5, 8, 13, 21, 32, 64):
        q, m = n*n, n//2
        receivers, reads = (q - 2*m)**3, (2*m + 1)**3
        require(3*m*m < n*n)
        require(8 * receivers * reads >= n**9)
        require(8 * n * receivers * reads >= n**10)
        require(reads >= n**3)
        rows.append({"n": n, "q": q, "interior_receivers": receivers,
                     "reads_per_receiver": reads, "read_incidences_per_horizon": n*receivers*reads})
    return rows
