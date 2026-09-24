"""Produce the entropy-preserving sparse family with a different radius."""

from itertools import product
from math import isqrt
from fractions import Fraction

from .model import cut_count


def parameters(t):
    if type(t) is not int or not 1 <= t <= 8:
        raise ValueError('balanced reference t must be in [1,8]')
    return 2**(8*t), 2**(4*t), 2**(3*t)


def nearest(top, bottom):
    return (1 if top >= 0 else -1)*((2*abs(top)+bottom)//(2*bottom))


def route(t, displacement):
    q,m,k = parameters(t)
    d = tuple(displacement)
    if len(d) != 3 or any(type(x) is not int for x in d):
        raise ValueError('integer displacement required')
    z = tuple(nearest(x,m) for x in d)
    square = sum(x*x for x in z)
    length = isqrt(square)
    if length*length < square:
        length += 1
    count = (length+k-3)//(k-2)
    steps = []
    previous = (0,0,0)
    for j in range(1,count+1):
        point = tuple(m*nearest(j*x,count) for x in z)
        steps.append(tuple(x-y for x,y in zip(point,previous)))
        previous = point
    residual = tuple(x-y for x,y in zip(d,previous))
    for axis in range(3):
        for bit in range(4*t):
            if abs(residual[axis]) & (1 << bit):
                steps.append(tuple((1 if residual[axis]>=0 else -1)*2**bit if i==axis else 0
                                   for i in range(3)))
    return steps


def certificate(t):
    if t not in (1,2,3):
        raise ValueError('balanced cut catalog is t=1,2,3')
    q,m,k = parameters(t)
    coordinate=diagonal=moment=fourth=degree=flux=0
    # Sum entire vertical columns, not a radius-cubed vector archive.
    for x,y in product(range(-k,k+1),repeat=2):
        residual = k*k-x*x-y*y
        if residual < 0:
            continue
        z = isqrt(residual)
        count = 2*z+1
        s2 = z*(z+1)*(2*z+1)//3
        s4 = z*(z+1)*(2*z+1)*(3*z*z+3*z-1)//15
        degree += count
        moment += m*m*(count*(x*x+y*y)+s2)
        fourth += m**4*(count*(x*x+y*y)**2+2*(x*x+y*y)*s2+s4)
        z_width = count*q-m*z*(z+1)
        if x>0:
            flux += x*count
            coordinate += m*x*(q-m*abs(y))*z_width
        if t<=2 and x+y>0:
            diagonal += cut_count(q,[(m*x,m*y,0)],True)//q*z_width
    # The zero offset counts as a read, but contributes no action moment.
    degree += 24*t
    moment += 6*sum(2**(2*b) for b in range(4*t))
    fourth += 6*sum(2**(4*b) for b in range(4*t))
    coordinate += q*q*(m-1)
    result = dict(q=q,spacing=m,coarse_radius=k,radius=m*k,degree=degree,
                  coordinate=coordinate,coarse_flux=flux,second_norm_moment=moment,
                  fourth_norm_moment=fourth)
    result['stabilized']={
        'coarse_edge_weight_times_h_squared':str(Fraction(3,moment)),
        'nearest_extra_weight_times_h_squared':'1/2',
        'fourth_error_coefficient_times_inverse_h_squared':str((Fraction(fourth,moment)+1)/8),
        'naive_alias_modes':m**3,
        'naive_alias_eigenvalue_times_L_squared_upper':str(Fraction(288*t*q*q,moment)),
        'stabilized_nonzero_alias_eigenvalue_times_L_squared_lower':8*q*q//(m*m)}
    if t<=2:
        diagonal += 2*q*sum(2**b*(q-2**b) for b in range(4*t))
        result['diagonal'] = diagonal
    return result
