"""Projective cross-ratio data and oriented spherical-cap reconstruction.

Floating calculations reject unresolved configurations. A supplied side witness
chooses cap orientation; circle data alone determine only an unoriented plane.
"""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Complex, Integral, Real

import numpy as np


ROUND = 256 * np.finfo(float).eps


def _real_array(value, shape_tail):
    a = np.asarray(value)
    if a.dtype.kind not in 'iuf' or a.ndim != len(shape_tail):
        raise ValueError('finite real coordinates are required')
    if any(size is not None and a.shape[i] != size for i, size in enumerate(shape_tail)):
        raise ValueError('coordinate shape mismatch')
    a = np.asarray(a, dtype=float)
    if not np.isfinite(a).all():
        raise ValueError('coordinates must be finite')
    return a


def sphere_points(points):
    points = _real_array(points, (None, 3))
    if len(points) == 0 or np.max(np.abs(np.linalg.norm(points, axis=1) - 1)) > 1e-10:
        raise ValueError('points must lie on the unit sphere; no implicit normalization')
    return points


def _projective(z):
    if isinstance(z, (bool, np.bool_)) or not isinstance(z, Complex) or np.ndim(z) != 0:
        raise ValueError('a complex projective coordinate is required')
    try:
        z = complex(z)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('invalid projective coordinate') from exc
    if np.isnan(z.real) or np.isnan(z.imag) or np.isinf(z.imag):
        raise ValueError('invalid projective coordinate')
    if np.isinf(z.real):
        if z.imag != 0 or z.real < 0:
            raise ValueError('infinity uses the canonical coordinate +inf+0j')
        return np.array([1, 0], dtype=complex)
    scale = max(1, abs(z.real), abs(z.imag))
    return np.array([z / scale, 1 / scale], dtype=complex)


def _det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def _ratio(num, den):
    if abs(num) + abs(den) <= ROUND:
        raise ValueError('cross ratio is degenerate or numerically unresolved')
    if den == 0:
        return complex(np.inf, 0)
    if abs(den) <= ROUND * max(abs(num), abs(den)):
        raise ValueError('finite cross ratio is too close to a pole to resolve')
    z = num / den
    if not np.isfinite(z):
        raise ValueError('cross ratio overflow')
    return complex(z)


def cross_ratio(z1, z2, z3, z4):
    """(z1-z3)(z2-z4)/((z1-z4)(z2-z3)), including exact infinity."""
    values = (z1, z2, z3, z4)
    vectors = list(map(_projective, values))
    for i in range(4):
        for j in range(i):
            bracket = _det(vectors[i], vectors[j])
            if ((bracket == 0 and values[i] != values[j])
                    or 0 < abs(bracket) <= ROUND):
                raise ValueError('cross-ratio points are too close to resolve')
    a, b, c, d = vectors
    ac, bd, ad, bc = _det(a, c), _det(b, d), _det(a, d), _det(b, c)
    return _ratio(ac * bd, ad * bc)


def mobius_normalize(z, g1, g2, g3):
    return cross_ratio(z, g2, g1, g3)


def stereographic(p):
    """North-pole stereographic coordinate with a stable second chart."""
    x, y, z = sphere_points(np.asarray(p)[None, :])[0]
    if x == 0 and y == 0 and z > 0:
        return complex(np.inf, 0)
    if z > 0:
        coordinate = complex((1 + z) / complex(x, -y))
        if not np.isfinite(coordinate):
            raise ValueError('finite stereographic coordinate overflow')
        return coordinate
    return complex(x, y) / (1 - z)


def _gauge(gauge, size):
    if (not isinstance(gauge, (list, tuple)) or len(gauge) != 3
            or any(not isinstance(i, Integral) or isinstance(i, (bool, np.bool_)) for i in gauge)
            or len(set(gauge)) != 3 or any(i < 0 or i >= size for i in gauge)):
        raise ValueError('gauge must name three distinct valid indices')
    return tuple(gauge)


def cross_ratio_receipts(points, gauge):
    """Coordinate-based fixture producer; returns only invariant scalar data."""
    points = sphere_points(points)
    gauge = _gauge(gauge, len(points))
    zs = [stereographic(p) for p in points]
    triple = [zs[i] for i in gauge]
    if any(abs(_det(_projective(triple[i]), _projective(triple[j]))) <= ROUND
           for i in range(3) for j in range(i)):
        raise ValueError('gauge points coincide or are unresolved')
    values = [mobius_normalize(z, *triple) for z in zs]
    # Structural pole at the gauge point is exact, including in the finite chart.
    values[gauge[0]], values[gauge[1]], values[gauge[2]] = 0j, 1 + 0j, complex(np.inf, 0)
    return reconstruct_from_cross_ratios(values, gauge)


def reconstruct_from_cross_ratios(values, gauge):
    """Read supplied normalized cross ratios; never access source coordinates.

    Entry i is CR(z_i,g2;g1,g3). A labelled distinct configuration is determined
    by these scalars up to a common Mobius transformation.
    """
    data = np.asarray(values)
    if data.ndim != 1 or data.dtype.kind not in 'iufc' or len(data) < 3:
        raise ValueError('a numeric vector of cross-ratio receipts is required')
    gauge = _gauge(gauge, len(data))
    projective = [_projective(z) for z in data]
    result = np.asarray(data, dtype=complex).copy()
    if result[gauge[0]] != 0 or result[gauge[1]] != 1 or result[gauge[2]] != complex(np.inf, 0):
        raise ValueError('gauge receipts must be exactly 0, 1, infinity')
    if any(abs(_det(projective[i], projective[j])) <= ROUND
           for i in range(len(data)) for j in range(i)):
        raise ValueError('distinct receipt points coincide or are unresolved')
    return result


def sphere_from_cross_ratios(values, gauge):
    """Unit-sphere representative in the supplied gauge, computed from ratios."""
    coordinates = reconstruct_from_cross_ratios(values, gauge)
    points = []
    for coordinate in coordinates:
        u, v = _projective(coordinate)
        uv = u * np.conj(v)
        scale = abs(u)**2 + abs(v)**2
        points.append([2 * uv.real / scale, 2 * uv.imag / scale,
                       (abs(u)**2 - abs(v)**2) / scale])
    return np.asarray(points)


@dataclass(frozen=True)
class CapFit:
    normal: np.ndarray
    max_incidence_residual: float
    plane_residual: float
    third_singular_value: float
    side_margin: float

    def direction_error_bound(self, input_operator_error: float) -> float:
        """Sine-angle bound to a true rank-three plane, given ||A-A_true||_2.

        The data-error bound is an input, not estimated from fitting residuals.
        This bounds the Euclidean unit plane direction, before Lorentz scaling.
        """
        e = _budget(input_operator_error, upper=None)
        if e >= self.third_singular_value:
            raise ValueError('input error does not resolve the plane rank')
        return min(1.0, (self.plane_residual + e) / (self.third_singular_value - e))


def _budget(value, upper=0.1):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not np.isfinite(value) or value < 0 or (upper is not None and value > upper)):
        raise ValueError('a finite nonnegative residual budget is required')
    return float(value)


def fit_cap(boundary_points, *, interior_point, max_residual=1e-10):
    """Fit a spacelike plane and select its side using an explicit witness.

    All boundary samples must meet the stated normalized incidence budget.
    No spectral flooring, radius clipping, outlier removal or sign default.
    """
    points = sphere_points(boundary_points)
    if len(points) < 3:
        raise ValueError('at least three distinct boundary points are required')
    budget = _budget(max_residual)
    witness = sphere_points(_real_array(interior_point, (3,))[None, :])[0]
    a = np.column_stack((-np.ones(len(points)), points))
    _, singular, vh = np.linalg.svd(a, full_matrices=len(points) < 4)
    if singular[2] <= ROUND * singular[0]:
        raise ValueError('circle plane is degenerate or unresolved')
    u = vh[-1]
    q = float(np.dot(u[1:], u[1:]) - u[0] ** 2)
    if q <= ROUND:
        raise ValueError('cap plane is not resolvably spacelike')
    n = u / np.sqrt(q)
    residual = float(np.max(np.abs(a @ n)))
    if not np.isfinite(n).all() or not np.isfinite(residual) or residual > budget:
        raise ValueError('boundary samples exceed the declared circle residual budget')
    side = float(np.dot(n[1:], witness) - n[0])
    if abs(side) <= budget + ROUND * np.linalg.norm(n):
        raise ValueError('side witness is on or unresolved from the boundary')
    if side < 0:
        n = -n
    n.setflags(write=False)
    return CapFit(n, residual, float(np.linalg.norm(a @ u)), float(singular[2]), abs(side))


def produced_cap_normal(boundary_points, *, interior_point, max_residual=1e-10):
    return fit_cap(boundary_points, interior_point=interior_point,
                   max_residual=max_residual).normal


def cap_from_cross_ratios(values, gauge, *, boundary_indices, interior_index,
                          max_residual=1e-10):
    """Consume scalar receipts and labelled inside/boundary data in one frame."""
    points = sphere_from_cross_ratios(values, gauge)
    if (not isinstance(boundary_indices, (list, tuple)) or len(boundary_indices) < 3
            or any(not isinstance(i, Integral) or isinstance(i, (bool, np.bool_))
                   or not 0 <= i < len(points) for i in boundary_indices)
            or len(set(boundary_indices)) != len(boundary_indices)
            or not isinstance(interior_index, Integral)
            or isinstance(interior_index, (bool, np.bool_))
            or not 0 <= interior_index < len(points) or interior_index in boundary_indices):
        raise ValueError('distinct valid boundary indices and a separate inside index are required')
    return fit_cap(points[list(boundary_indices)], interior_point=points[interior_index],
                   max_residual=max_residual)
