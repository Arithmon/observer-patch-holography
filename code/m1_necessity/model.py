"""Construct candidates. The independent checker does not import this module."""

from itertools import product
from math import isqrt
from fractions import Fraction


def parameters(n):
    if type(n) is not int or not 1 <= n <= 24:
        raise ValueError("reference n must be an integer in [1,24]")
    r = 2**n
    return r, r*r


def normalized_flight(w, radius):
    square = sum(x*x for x in w)
    return tuple((1 if x >= 0 else -1)*isqrt(radius*radius*x*x//square) for x in w)


def sparse(n):
    radius, _ = parameters(n)
    vectors = {(0, 0, 0)}
    # Enumerate six faces; a set removes shared edges and corners.
    for axis in range(3):
        for sign in (-1, 1):
            for a, b in product(range(-n, n+1), repeat=2):
                w = [a, b]
                w.insert(axis, sign*n)
                vectors.add(normalized_flight(w, radius))
    for axis in range(3):
        for bit in range(n+1):
            for sign in (-1, 1):
                v = [0, 0, 0]
                v[axis] = sign*2**bit
                vectors.add(tuple(v))
    return sorted(vectors)


def stencil(n, family):
    radius, _ = parameters(n)
    if family == "sparse":
        return sparse(n)
    if n > 4:
        raise ValueError("expanded reference stencils are bounded to n<=4")
    if family == "dense":
        return [v for v in product(range(-radius, radius+1), repeat=3)
                if sum(x*x for x in v) <= radius*radius]
    if family == "axial":
        result = set(sparse(n))
        for axis, length in product(range(3), range(1, radius+1)):
            for sign in (-1, 1):
                v = [0, 0, 0]
                v[axis] = sign*length
                result.add(tuple(v))
        return sorted(result)
    raise ValueError("unknown stencil family")


def route(n, displacement):
    radius, _ = parameters(n)
    d = tuple(displacement)
    if len(d) != 3 or any(type(x) is not int for x in d):
        raise ValueError("integer displacement required")
    maximum = max(map(abs, d))
    if not maximum:
        return []
    w = tuple((1 if x >= 0 else -1)*((2*n*abs(x)+maximum)//(2*maximum)) for x in d)
    flight = normalized_flight(w, radius)
    repeats = isqrt(sum(x*x for x in d))//radius
    result = [flight]*repeats
    for axis in range(3):
        residual = d[axis]-repeats*flight[axis]
        sign = 1 if residual >= 0 else -1
        whole, rest = divmod(abs(residual), radius)
        v = [0, 0, 0]
        v[axis] = sign*radius
        result.extend([tuple(v)]*whole)
        for bit in range(n):
            if rest & (1 << bit):
                v = [0, 0, 0]
                v[axis] = sign*(1 << bit)
                result.append(tuple(v))
    return result


def triangular(t):
    return (t+1)*(t+2)//2 if t >= 0 else 0


def rectangle_below(t, nx, ny):
    return triangular(t)-triangular(t-nx)-triangular(t-ny)+triangular(t-nx-ny)


def cut_count(q, vectors, diagonal=False):
    total = 0
    for x, y, z in vectors:
        projection = x+y if diagonal else x
        if projection <= 0:
            continue
        if not diagonal:
            total += x*(q-abs(y))*(q-abs(z))
        else:
            lx, ly = max(0, -x), max(0, -y)
            nx, ny = q-abs(x), q-abs(y)
            upper, lower = q-1-lx-ly, q-projection-1-lx-ly
            total += (rectangle_below(upper, nx, ny)-rectangle_below(lower, nx, ny))*(q-abs(z))
    return total


def cuts(n):
    radius, q = parameters(n)
    vectors = sparse(n)
    extra = [v for v in vectors if sum(x != 0 for x in v) > 1]
    axis_x = q*q*radius*(radius+1)//2
    axis_d = q*q*radius*(radius+1)-q*radius*(radius+1)*(2*radius+1)//3
    result = {"sparse_degree": len(vectors), "axial_degree": 6*radius+len(extra)+1,
              "sparse_coordinate": cut_count(q, vectors),
              "sparse_diagonal": cut_count(q, vectors, True),
              "axial_coordinate": axis_x+cut_count(q, extra),
              "axial_diagonal": axis_d+cut_count(q, extra, True)}
    if n <= 4:
        dense = stencil(n, "dense")
        result.update(dense_degree=len(dense), dense_coordinate=cut_count(q, dense),
                      dense_diagonal=cut_count(q, dense, True))
    return result


def action(n, family):
    radius,_ = parameters(n)
    if family == "axial":
        vectors = [v for v in sparse(n) if sum(x != 0 for x in v) > 1]
        square_axis = radius*(radius+1)*(2*radius+1)//3
        fourth_axis = radius*(radius+1)*(2*radius+1)*(3*radius*radius+3*radius-1)//15
        matrix = [[square_axis if i == j else 0 for j in range(3)] for i in range(3)]
        fourth = 3*fourth_axis
        degree = 6*radius+len(vectors)
        parity = [2*radius*t for t in (1,2,3)]
    else:
        vectors = [v for v in stencil(n,family) if any(v)]
        matrix, fourth, degree = [[0]*3 for _ in range(3)],0,len(vectors)
        parity = [0,0,0]
    for v in vectors:
        for i,j in product(range(3),repeat=2):
            matrix[i][j] += v[i]*v[j]
        fourth += sum(x*x for x in v)**2
        for t in (1,2,3):
            parity[t-1] += 1-(-1)**(sum(v[:t]) % 2)
    moment = sum(matrix[i][i] for i in range(3))
    return {"degree":degree, "second_moments":matrix, "fourth_norm_moment":fourth,
            "alpha_times_h_squared":str(Fraction(6,moment)),
            "tick_squared_eigenvalues":[str(Fraction(6*radius*radius*p,moment)) for p in parity],
            "tick_squared_mean_eigenvalue":str(Fraction(6*radius*radius*degree,moment))}


def graph(n, family, layers=3, intervention=None):
    _, q = parameters(n)
    if n not in (1, 2) or layers != 3:
        raise ValueError("full graph reference uses n=1,2 and three transitions")
    sites = list(product(range(q), repeat=3))
    vectors = stencil(n, family)
    parents = {x: [tuple(x[i]-v[i] for i in range(3)) for v in vectors
                   if all(0 <= x[i]-v[i] < q for i in range(3))] for x in sites}
    states = [{x: 1+(x[0]+2*x[1]+3*x[2]) % 7 for x in sites}]
    for layer in range(layers+1):
        if layer:
            previous = states[-1]
            states.append({x: 1+sum(previous[y] for y in parents[x]) for x in sites})
        if intervention is not None and layer == intervention[0]:
            states[-1][tuple(intervention[1])] += 1
    return sites, parents, states
