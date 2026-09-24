"""Independent reference definitions, integer predicates and outgoing graph replay."""

from fractions import Fraction as F
from itertools import product, permutations
import hashlib
import json


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key")
        result[key] = value
    return result


def strict_load(path):
    from pathlib import Path
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=no_duplicates,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def level(n):
    require(type(n) is int and 1 <= n <= 24, "invalid reference level")
    return 1 << n, 1 << (2*n)


def vector(v):
    require(type(v) in (tuple, list) and len(v) == 3 and
            all(type(x) is int for x in v), "nonintegral vector")
    return tuple(v)


def floor_ratio_root(top, bottom, cap):
    # Squared inequalities, with no float or producer square-root routine.
    lo, hi = 0, cap+1
    while hi-lo > 1:
        middle = (lo+hi)//2
        if middle*middle*bottom <= top:
            lo = middle
        else:
            hi = middle
    require(lo*lo*bottom <= top < (lo+1)**2*bottom, "floor-root bracket")
    return lo


def expected_sparse(n):
    radius, _ = level(n)
    answer = {(0, 0, 0)}
    # Partition the cube boundary by its first maximal coordinate.
    for first in range(3):
        ranges = [range(-n+1, n) if i < first else
                  (-n, n) if i == first else range(-n, n+1) for i in range(3)]
        for w in product(*ranges):
            square = sum(x*x for x in w)
            answer.add(tuple((-1 if x < 0 else 1)*floor_ratio_root(
                radius*radius*x*x, square, radius) for x in w))
    for axis in range(3):
        length = 1
        while length <= radius:
            for sign in (-1, 1):
                v = tuple(sign*length if i == axis else 0 for i in range(3))
                answer.add(v)
            length *= 2
    return answer


def check_sparse(n, submitted):
    radius, _ = level(n)
    require(type(submitted) is list, "stencil list required")
    values = [vector(v) for v in submitted]
    require(len(set(values)) == len(values), "duplicate read vector")
    require(set(values) == expected_sparse(n), "incomplete or substituted sparse stencil")
    require(len(values) <= 24*n*n+6*n+9, "degree bound")
    require(all(sum(x*x for x in v) <= radius*radius for v in values), "superluminal edge")
    # Check generators of the full cubic symmetry, not just vector counts.
    chosen = set(values)
    for x, y, z in values:
        require((-x,y,z) in chosen and (y,z,x) in chosen and (y,x,z) in chosen,
                "lost cubic stencil symmetry")
    return chosen


def check_route(n, displacement, steps, allowed=None):
    radius, _ = level(n)
    d = vector(displacement)
    require(type(steps) is list, "route list required")
    accepted = expected_sparse(n) if allowed is None else allowed
    place = [0, 0, 0]
    points = [place.copy()]
    for raw in steps:
        v = vector(raw)
        require(v in accepted and v != (0,0,0), "illegal or uncharged route step")
        place = [a+b for a,b in zip(place,v)]
        points.append(place)
    require(tuple(place) == d, "route endpoint mismatch")
    square = sum(x*x for x in d)
    count = len(steps)
    require(square <= (count*radius)**2, "outer cone bound")
    left = (count-3*n-2)*radius
    factor = F(1)+F(3,n)+F(3,radius)
    require(left <= 0 or left*left <= factor*factor*square, "constructive time bound")
    # Check the physical boundary consequence too: every submitted point must
    # lie in the proved Euclidean tube around the closed endpoint segment.
    # Compare against ((3/n+3/R)*sqrt(square)+2*R)^2 by rational squaring.
    slope = F(3,n)+F(3,radius)
    for point in points:
        dot = sum(x*y for x,y in zip(point,d))
        if square and 0 < dot < square:
            distance_squared = F(sum(x*x for x in point))-F(dot*dot,square)
        else:
            end = d if square and dot >= square else (0,0,0)
            distance_squared = F(sum((x-y)**2 for x,y in zip(point,end)))
        excess = distance_squared-slope*slope*square-4*radius*radius
        require(excess <= 0 or excess*excess <= 16*radius*radius*slope*slope*square,
                "route leaves the proved boundary tube")
    return {"steps": count, "points_sha256": digest(points)}


def ramp_rows(lx, hx, ly, hy, lower, upper):
    # Count the clipped affine row widths by their breakpoints. This is not
    # the producer's inclusion-exclusion formula for triangular numbers.
    cuts = {lx, hx+1}
    for threshold in (lower-ly, lower-hy, upper-ly, upper-hy):
        cuts.update(k for k in range(threshold-1, threshold+3) if lx <= k <= hx+1)
    def width(x):
        return max(0, min(hy,upper-x)-max(ly,lower-x)+1)
    total = 0
    edges = sorted(cuts)
    for a, stop in zip(edges,edges[1:]):
        b = stop-1
        numerator = (b-a+1)*(width(a)+width(b))
        require(numerator % 2 == 0, "nonintegral row sum")
        if b > a:
            middle = (a+b)//2
            require((width(middle)-width(a))*(b-a) == (width(b)-width(a))*(middle-a),
                    "missed clipped-row breakpoint")
        total += numerator//2
    return total


def cut(q, vectors, diagonal=False):
    total = 0
    for raw in vectors:
        x,y,z = vector(raw)
        projection = x+y if diagonal else x
        if projection <= 0:
            continue
        lx,hx = max(0,-x),min(q-1,q-1-x)
        ly,hy = max(0,-y),min(q-1,q-1-y)
        if diagonal:
            rows = ramp_rows(lx,hx,ly,hy,q-projection,q-1)
        else:
            left = max(lx,q//2-x)
            right = min(hx,q//2-1)
            rows = max(0,right-left+1)*(hy-ly+1)
        total += rows*(q-abs(z))
    return total


def expected_cuts(n, vectors):
    radius,q = level(n)
    chosen = check_sparse(n,vectors)
    extra = {v for v in chosen if sum(x != 0 for x in v) > 1}
    # Closed axial sums, independently factored from the published formulas.
    x0 = (radius*(radius+1)//2)*q*q
    y0 = 2*q*(q*(radius*(radius+1)//2)-radius*(radius+1)*(2*radius+1)//6)
    result = {"sparse_degree":len(chosen), "axial_degree":len(extra)+6*radius+1,
              "sparse_coordinate":cut(q,chosen), "sparse_diagonal":cut(q,chosen,True),
              "axial_coordinate":x0+cut(q,extra), "axial_diagonal":y0+cut(q,extra,True)}
    degree = 24*n*n+6*n+9
    require(result["sparse_coordinate"] <= degree*radius*q*q, "sparse cut upper bound")
    ratio = F(result["axial_diagonal"],2*result["axial_coordinate"])
    b = 1-F(2*radius+1,3*q)
    require(b/(1+F(2*degree,radius)) <= ratio <= b+F(2*degree,radius),
            "anisotropy ratio enclosure")
    if n <= 4:
        dense = set()
        for x,y in product(range(-radius,radius+1),repeat=2):
            remaining = radius*radius-x*x-y*y
            if remaining < 0:
                continue
            for z in range(-radius,radius+1):
                if z*z <= remaining:
                    dense.add((x,y,z))
        result.update(dense_degree=len(dense),dense_coordinate=cut(q,dense),
                      dense_diagonal=cut(q,dense,True))
        require(32*result["dense_coordinate"] >= radius**4*q*q, "dense cut lower bound")
        require(chosen <= dense, "sparse family not a dense subgraph")
    return result


def check_action(n, family, submitted):
    radius,_ = level(n)
    chosen = expected_sparse(n)
    axis_part = family == "axial"
    if axis_part:
        chosen = {v for v in chosen if sum(x != 0 for x in v) > 1}
    elif family == "dense":
        require(n <= 4, "dense action reference scope")
        chosen = {v for v in product(range(-radius,radius+1),repeat=3)
                  if sum(x*x for x in v) <= radius*radius}
    else:
        require(family == "sparse", "unknown action")
    chosen.discard((0,0,0))
    matrix = [[sum(v[i]*v[j] for v in chosen) for j in range(3)] for i in range(3)]
    fourth = sum(sum(x*x for x in v)**2 for v in chosen)
    parity = [sum(2 for v in chosen if sum(v[:t]) % 2) for t in (1,2,3)]
    count = len(chosen)
    if axis_part:
        # Polynomial sums are independently checked by finite differences in tests.
        s2 = F(radius**3,3)+F(radius**2,2)+F(radius,6)
        s4 = F(radius**5,5)+F(radius**4,2)+F(radius**3,3)-F(radius,30)
        for i in range(3):
            matrix[i][i] += int(2*s2)
        fourth += int(6*s4)
        count += 6*radius
        parity = [x+2*radius*t for t,x in enumerate(parity,1)]
    moment = sum(matrix[i][i] for i in range(3))
    require(moment > 0, "zero moment")
    require(all(3*matrix[i][j] == (moment if i == j else 0)
                for i,j in product(range(3),repeat=2)), "nonisotropic action")
    require(fourth <= radius*radius*moment, "Taylor moment bound")
    expected = {"degree":count, "second_moments":matrix, "fourth_norm_moment":fourth,
                "alpha_times_h_squared":str(F(6,moment)),
                "tick_squared_eigenvalues":[str(F(6*radius*radius*p,moment)) for p in parity],
                "tick_squared_mean_eigenvalue":str(F(6*radius*radius*count,moment))}
    require(F(expected["tick_squared_mean_eigenvalue"]) >= 6, "trace stability obstruction")
    require(max(map(F,expected["tick_squared_eigenvalues"])) > 4,
            "reference parity modes do not witness instability")
    require(canonical(submitted) == canonical(expected), "action certificate disagrees with actual stencil")
    return expected


def graph_replay(n, family, layers=3, intervention=None):
    radius,q = level(n)
    require(n in (1,2) and layers == 3, "full graph reference scope")
    if family == "dense":
        offsets = {v for v in product(range(-radius,radius+1),repeat=3)
                   if sum(x*x for x in v) <= radius*radius}
    else:
        offsets = expected_sparse(n)
        if family == "axial":
            offsets.update(tuple(j if i == axis else 0 for i in range(3))
                           for axis in range(3) for j in range(-radius,radius+1))
        else:
            require(family == "sparse", "unknown graph family")
    sites = list(product(range(q),repeat=3))
    index = {v:i for i,v in enumerate(sites)}
    outgoing = [[] for _ in sites]
    for i,x in enumerate(sites):
        for v in offsets:
            target = tuple(a+b for a,b in zip(x,v))
            if target in index:
                outgoing[i].append(index[target])
    states = [[1+(x[0]+2*x[1]+3*x[2]) % 7 for x in sites]]
    for layer in range(layers+1):
        if layer:
            new = [1]*len(sites)
            for i,value in enumerate(states[-1]):
                for j in outgoing[i]:
                    new[j] += value
            states.append(new)
        if intervention is not None and layer == intervention[0]:
            states[-1][index[tuple(intervention[1])]] += 1
    return sites,outgoing,states


def direct_cut_pairs(n, family, diagonal=False):
    """Full small-box pair enumeration; no cut formula is used."""
    sites,outgoing,_ = graph_replay(n,family)
    _,q = level(n)
    threshold = 2*q-1 if diagonal else q-1
    projection = lambda x: 2*(x[0]+x[1] if diagonal else x[0])
    count = 0
    for i,targets in enumerate(outgoing):
        if projection(sites[i]) >= threshold:
            continue
        for j in targets:
            if projection(sites[j]) > threshold:
                count += 1
    return count
