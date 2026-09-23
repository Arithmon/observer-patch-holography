"""Produce compact cut and moment evidence for explicitly specified stencils."""

from fractions import Fraction as F
from itertools import product
from math import isqrt

import numpy as np

CONFIGS = ((12, 4, 1, 2), (16, 4, 1, 2), (24, 4, 2, 2))
LABELS = ('plane', 'residue', 'diagonal', 'checkerboard', 'ball',
          'localized', 'connected', 'permuted_half')


def scale(t):
    if type(t) is not int or not 1 <= t <= 12:
        raise ValueError('reference scale must be an integer in [1,12]')
    q, m, k, r = (2**(e*t) for e in (16, 8, 6, 7))
    return dict(t=t, q=q, m=m, k=k, r=r, radius=m*k,
                tick_times_c_over_L=str(F(m*k, q)),
                compactness_error_factor=str(F(m*q, r**4)),
                added_area_factor=str(F(r**4, q*q)),
                alias_frequency_squared_scale=str(F(r**5, k**5)),
                spectral_saturation=str(F(q*q, m*m*k*k)))


def ball_columns(radius, q, spacing=1):
    """Exact complete vertical columns, without a radius-cubed archive."""
    if type(radius) is not int or not 0 <= radius <= 128:
        raise ValueError('column reference radius must be in [0,128]')
    degree = flux = moment = fourth = axis2 = axis4 = cut = 0
    for x, y in product(range(-radius, radius+1), repeat=2):
        residual = radius*radius-x*x-y*y
        if residual < 0:
            continue
        z = isqrt(residual)
        count = 2*z+1
        s2 = z*(z+1)*(2*z+1)//3
        s4 = z*(z+1)*(2*z+1)*(3*z*z+3*z-1)//15
        degree += count
        axis2 += x*x*count
        axis4 += x**4*count
        moment += (x*x+y*y)*count+s2
        fourth += (x*x+y*y)**2*count+2*(x*x+y*y)*s2+s4
        if x > 0:
            flux += x*count
            cut += spacing*x*(q-spacing*abs(y))*((2*z+1)*q-spacing*z*(z+1))
    return dict(degree=degree, flux=flux, moment=moment, fourth=fourth,
                axis2=axis2, axis4=axis4, coordinate_cut=cut)


def moment_case(kind):
    p = scale(1)
    q, m, k = p['q'], p['m'], p['k']
    radius = {'raw': 0, 'critical': k, 'repaired': p['r']}[kind]
    coarse = ball_columns(k, q, m)
    fine = ball_columns(radius, q)
    extra = [2**b for b in range(8) if 2**b > radius]
    degree = coarse['degree']+fine['degree']-1+6*len(extra)
    moment = m*m*coarse['moment']+fine['moment']+6*sum(s*s for s in extra)
    fourth = m**4*coarse['fourth']+fine['fourth']+6*sum(s**4 for s in extra)
    alpha = F(6*q*q, moment)
    theta_lo, theta_hi = F(6, m), F(44, 7*m)  # 3 < pi < 22/7
    lower = max(F(0), alpha*(theta_lo**2*fine['axis2']/2-theta_hi**4*fine['axis4']/24))
    upper = alpha*(theta_hi**2*fine['axis2']/2+4*len(extra))
    coordinate = coarse['coordinate_cut']+fine['coordinate_cut']+q*q*sum(extra)
    alias = (2*q**3//m)*(fine['flux']+sum(extra))
    return dict(kind=kind, q=q, m=m, k=k, r=radius, degree=degree,
                second_norm_moment=moment, fourth_norm_moment=fourth,
                coarse=coarse, fine=fine, remaining_axis_lengths=extra,
                coordinate_cut=coordinate, periodic_residue_cut=alias,
                pi_times_residue_entropy=str(F(4*alias, q**4)),
                uniform_alpha_L_squared=str(alpha),
                one_alias_eigenvalue_L_squared_bounds=[str(lower), str(upper)])


def stencil(config, family):
    q, m, k, radius = config
    if config not in CONFIGS or family not in ('T', 'U'):
        raise ValueError('unknown finite reference graph')
    chosen = {tuple(m*x for x in v) for v in product(range(-k,k+1), repeat=3)
              if sum(x*x for x in v) <= k*k}
    for b in range(m.bit_length()-1):
        for axis, sign in product(range(3), (-1, 1)):
            chosen.add(tuple(sign*2**b if i == axis else 0 for i in range(3)))
    if family == 'U':
        chosen.update(v for v in product(range(-radius,radius+1), repeat=3)
                      if sum(x*x for x in v) <= radius*radius)
    return sorted(chosen)


def labels(config):
    q, m, _, _ = config
    x, y, z = np.indices((q,q,q))
    ball = 16*((x-q//2)**2+(y-q//2)**2+(z-q//2)**2) < q*q
    residue = x % m < m//2
    localized = ball & residue
    return dict(plane=x<q//2, residue=residue, diagonal=(x+y+z)%q<q//2,
                checkerboard=(x+y+z)%2==0, ball=ball, localized=localized,
                connected=localized | (ball & (y==q//2) & (z==q//2)),
                permuted_half=(137*(x*q*q+y*q+z)+29)%(q**3)<q**3//2)


def graph_cuts(config, family, topology):
    if topology not in ('periodic', 'clipped'):
        raise ValueError('unknown topology')
    q = config[0]
    vectors = stencil(config, family)
    masks = labels(config)
    arrays = np.stack([masks[name] for name in LABELS])
    counts = np.zeros(len(LABELS), dtype=np.int64)
    for v in vectors:
        if topology == 'periodic':
            other = np.roll(arrays, tuple(-x for x in v), axis=(1,2,3))
            counts += np.count_nonzero(arrays != other, axis=(1,2,3))
        else:
            left = (slice(None),)+tuple(slice(max(0,-x),min(q,q-x)) for x in v)
            right = (slice(None),)+tuple(slice(max(0,x),min(q,q+x)) for x in v)
            counts += np.count_nonzero(arrays[left] != arrays[right], axis=(1,2,3))
    if any(int(x)%2 for x in counts):
        raise ValueError('directed crossing parity')
    reads = q**3*len(vectors) if topology == 'periodic' else sum(
        (q-abs(x))*(q-abs(y))*(q-abs(z)) for x,y,z in vectors)
    return masks, dict(zip(LABELS, (int(x)//2 for x in counts))), reads
