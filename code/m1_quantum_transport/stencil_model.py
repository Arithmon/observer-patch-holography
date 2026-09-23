"""Candidate general-stencil velocity and exact support-polytope evidence."""

import itertools
import sympy as s

from .model import I, IDENTITY, PAULI, digest, matrix


def geometries():
    cube = list(itertools.product((-1, 1), repeat=3))
    return {
        "tetrahedron": [p for p in cube if p[0]*p[1]*p[2] == 1],
        "cube": cube,
        "octahedron": [p for p in itertools.product((-1, 0, 1), repeat=3) if sum(x*x for x in p) == 1],
        "integer_ball_2": [p for p in itertools.product(range(-2, 3), repeat=3) if 0 < sum(x*x for x in p) <= 4],
    }


def polytopes():
    out = {}
    for name, points in geometries().items():
        facets = set()
        for triangle in itertools.combinations(points, 3):
            system = s.Matrix(triangle)
            if system.det() == 0:
                continue
            normal = tuple(system.inv()*s.ones(3, 1))
            if all(sum(a*b for a, b in zip(normal, point)) <= 1 for point in points):
                facets.add(normal)
        radius2 = max(sum(x*x for x in p) for p in points)
        inradius2 = min(1/sum(x*x for x in n) for n in facets)
        out[name] = dict(points=[list(p) for p in sorted(points)],
                         facets=[[str(x) for x in n] for n in sorted(facets)],
                         inradius_squared=str(inradius2), circumradius_squared=str(radius2),
                         geometric_ratio_squared=str(inradius2/radius2),
                         nonzero_displacements=len(points), cap_bound=str(1-s.Rational(2, len(points))))
    return out


def velocities():
    rows, derivatives = [], []
    directions = [p for p in itertools.product((-1, 0, 1), repeat=3) if any(p)]
    for x, y in itertools.product(range(8), repeat=2):
        steps = [s.cos(s.pi*k/4)*IDENTITY-I*s.sin(s.pi*k/4)*PAULI[i]
                 for i, k in enumerate((x, y, 0))]
        u = steps[2]*steps[1]*steps[0]
        partials = [steps[2]*steps[1]*(-I*PAULI[0])*steps[0],
                    steps[2]*(-I*PAULI[1])*steps[1]*steps[0],
                    (-I*PAULI[2])*steps[2]*steps[1]*steps[0]]
        generators = [(I*u.conjugate().T*d).applyfunc(s.expand) for d in partials]
        derivatives.append([[x, y, 0], [matrix(a) for a in generators]])
        for n in directions:
            a = sum((n[i]*generators[i] for i in range(3)), s.zeros(2))
            h = sum(abs(v) for v in n)
            lo, hi = h*IDENTITY+a, h*IDENTITY-a
            rows.append([[x, y, 0], list(n), -h, h, matrix(a),
                         [str(s.expand(z)) for z in (lo[0, 0], lo[1, 1], lo.det(), hi[0, 0], hi[1, 1], hi.det())]])
    return dict(momentum_grid=8, momenta=64, directions=26, comparisons=len(rows),
                derivatives_sha256=digest(derivatives), loewner_rows_sha256=digest(rows))


def candidate():
    return dict(polytopes=polytopes(), velocities=velocities(), tilt_controls=tilt_controls())


def tilt_controls():
    """Differentiate complete unitaries on both sides of the positive-radius scope."""
    variables = s.symbols('kx ky kz', real=True)
    origin = dict.fromkeys(variables, 0)
    result = {}
    for name, drift, scale in (("stationary", 0, 0), ("translation", 1, 0),
                               ("tilted_isotropic", 1, 3)):
        u = s.exp(-I*drift*variables[0])*IDENTITY
        for i, k in enumerate(variables):
            u = (s.cos(scale*k)*IDENTITY-I*s.sin(scale*k)*PAULI[i])*u
        zero = u.subs(origin)
        generators = [I*zero.conjugate().T*u.diff(k).subs(origin) for k in variables]
        tilt = s.Matrix([s.trace(a)/2 for a in generators])
        velocity = s.Matrix([[s.trace(a*p)/2 for p in PAULI] for a in generators])
        least = min(velocity.singular_values())
        radius = least-s.sqrt(tilt.dot(tilt))
        if scale:
            points = [(drift+scale*x, scale*y, scale*z) for x, y, z in
                      itertools.product((-1, 1), repeat=3)]
        else:
            points = [(drift, 0, 0)]
        result[name] = dict(period="1", internal_dimension=2,
                            generators=[matrix(a) for a in generators],
                            tilt=[str(x) for x in tilt], velocity=[[str(x) for x in row] for row in velocity.tolist()],
                            nonzero_displacements=sum(any(p) for p in points),
                            least_singular_value=str(least),
                            positive_centered_radius=str(radius) if radius > 0 else None,
                            cap_applies=bool(radius > 0))
    return result
