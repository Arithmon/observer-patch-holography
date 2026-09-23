"""Independent Fourier-kernel differentiation and exact supporting-face reconstruction."""

from fractions import Fraction as F
import itertools

from .check import need, digest, Q, Z, J, eye, mat, add, scale, mul, dagger, encoded, native_word


def support_polytope(points):
    need(type(points) is list and points, "missing displacement support")
    need(all(type(p) is list and len(p) == 3 and all(type(x) is int for x in p) for p in points),
         "support must have exact integer coordinates")
    points = sorted(set(tuple(p) for p in points))
    facets = set()
    for a, b, c in itertools.combinations(points, 3):
        u, v = [b[i]-a[i] for i in range(3)], [c[i]-a[i] for i in range(3)]
        normal = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
        if not any(normal):
            continue
        for sign in (-1, 1):
            n = [sign*x for x in normal]
            h = sum(n[i]*a[i] for i in range(3))
            if all(sum(n[i]*p[i] for i in range(3)) <= h for p in points):
                need(h > 0, "support must surround the origin in three dimensions")
                facets.add(tuple(F(x, h) for x in n))
    need(bool(facets), "support has no full-dimensional interior")
    radius2 = max(sum(x*x for x in p) for p in points)
    inradius2 = min(1/sum(x*x for x in n) for n in facets)
    count = sum(any(p) for p in points)
    need(count >= 4, "three-dimensional cone needs at least four displacements")
    return dict(points=[list(p) for p in points], facets=[[str(x) for x in n] for n in sorted(facets)],
                inradius_squared=str(inradius2), circumradius_squared=str(radius2),
                geometric_ratio_squared=str(inradius2/radius2), nonzero_displacements=count,
                cap_bound=str(1-F(2, count)))


def validate_isotropic_support(points, period_squared, field_speed_squared, native_speed_squared):
    values = (period_squared, field_speed_squared, native_speed_squared)
    need(all(type(v) is str for v in values), "clock and squared speeds require rational strings")
    t2, v2, c2 = map(F, values)
    need(t2 > 0 and c2 > 0 and v2 > 0, "positive clock and speeds required")
    result = support_polytope(points)
    need(F(result["circumradius_squared"]) <= c2*t2, "actual displacement exceeds native clock")
    need(v2*t2 <= F(result["inradius_squared"]), "isotropic field ball escapes finite support")
    need(v2 <= c2*F(result["cap_bound"])**2, "angular resource bound violated")
    return result


def validate_design(points, weights):
    """Check a weighted equal-radius integer design before normalizing its directions."""
    need(type(points) is list and type(weights) is list and len(points) == len(weights), "design shape")
    need(all(type(p) is list and len(p) == 3 and all(type(x) is int for x in p) for p in points), "design coordinates")
    need(len(points) >= 4 and len(set(map(tuple, points))) == len(points), "distinct minimal design outcomes")
    need(all(type(w) is str for w in weights), "exact design weights required")
    w = list(map(F, weights))
    need(all(x > 0 for x in w) and sum(w) == 1, "positive normalized design weights")
    radius2 = sum(x*x for x in points[0])
    need(radius2 > 0 and all(sum(x*x for x in p) == radius2 for p in points), "equal nonzero flight speed")
    for i in range(3):
        need(sum(a*p[i] for a, p in zip(w, points)) == 0, "design has a tilt")
        for j in range(3):
            need(sum(a*p[i]*p[j] for a, p in zip(w, points)) == (F(radius2, 3) if i == j else 0),
                 "design does not saturate isotropic velocity")
    return True


def polytopes():
    cube = [list(p) for p in itertools.product((-1, 1), repeat=3)]
    tetra = [p for p in cube if p[0]*p[1]*p[2] == 1]
    octa = [[0, 0, sign] for sign in (-1, 1)]+[[0, sign, 0] for sign in (-1, 1)]+[[sign, 0, 0] for sign in (-1, 1)]
    ball = [[x, y, z] for x in range(-2, 3) for y in range(-2, 3) for z in range(-2, 3)
            if 0 < x*x+y*y+z*z <= 4]
    return {name: support_polytope(points) for name, points in
            (("tetrahedron", tetra), ("cube", cube), ("octahedron", octa), ("integer_ball_2", ball))}


def velocities():
    raw = native_word((0, 1, 2))["kernel"]
    kernel = [(b, [[Z(Q(F(z[0])), Q(F(z[1]))) for z in row] for row in a]) for b, a in raw]
    cosine = [Q(1), Q(0, F(1, 2)), Q(), Q(0, F(-1, 2)), Q(-1), Q(0, F(-1, 2)), Q(), Q(0, F(1, 2))]
    phases = [Z(cosine[k], -cosine[(k-2) % 8]) for k in range(8)]
    rows, derivatives = [], []
    directions = [p for p in itertools.product((-1, 0, 1), repeat=3) if any(p)]
    for x, y in itertools.product(range(8), repeat=2):
        u = mat([[0, 0], [0, 0]])
        moments = [mat([[0, 0], [0, 0]]) for _ in range(3)]
        for b, coefficient in kernel:
            contribution = scale(coefficient, phases[(b[0]*x+b[1]*y) % 8])
            u = add(u, contribution)
            for i in range(3):
                moments[i] = add(moments[i], scale(contribution, b[i]))
        need(mul(dagger(u), u) == eye(2), "finite-stencil unitarity")
        generators = [mul(dagger(u), moment) for moment in moments]
        derivatives.append([[x, y, 0], [encoded(a) for a in generators]])
        for n in directions:
            a = mat([[0, 0], [0, 0]])
            for i in range(3):
                a = add(a, scale(generators[i], n[i]))
            need(dagger(a) == a, "native velocity must be Hermitian")
            h = max(sum(n[i]*b[i] for i in range(3)) for b, _ in kernel)
            l = min(sum(n[i]*b[i] for i in range(3)) for b, _ in kernel)
            lower, upper = add(a, scale(eye(2), -l)), add(scale(eye(2), h), scale(a, -1))
            minors = []
            for bound in (lower, upper):
                for value in (bound[0][0], bound[1][1], bound[0][0]*bound[1][1]-bound[0][1]*bound[1][0]):
                    need(value.im == Q() and value.re.b == 0 and value.re.a >= 0, "Loewner support-plane bound")
                    minors.append(str(value.re.a))
            rows.append([[x, y, 0], list(n), l, h, encoded(a), minors])
    return dict(momentum_grid=8, momenta=64, directions=26, comparisons=len(rows),
                derivatives_sha256=digest(derivatives), loewner_rows_sha256=digest(rows))


def reconstruct():
    return dict(polytopes=polytopes(), velocities=velocities())
