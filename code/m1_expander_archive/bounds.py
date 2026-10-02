"""Exact contraction, tail coefficients and a constructive refinement schedule."""

from fractions import Fraction
import math


def pair(x):
    return [x.numerator, x.denominator]


def schedule(q):
    if type(q) is not int or q < 2:
        raise ValueError('integer refinement q>=2')
    ell = (q-1).bit_length()+1
    target = 2048*ell
    m = math.isqrt(target)
    m += int(m*m < target)
    levels = (16*ell-1).bit_length()
    return dict(q=q, logarithm_ceiling=ell, m=m, bits=m*m,
                levels=levels, fault_power=2**levels)


def candidate():
    d = 8**12
    width, depth = d+2, 2+d+d*(d-1)//2
    tails = []
    for v, m in ((5, 2), (9, 3), (16, 4), (24, 6)):
        tails.append(dict(locations=v, required=m,
                          coefficients=[(-1)**(j-m)*math.comb(j-1, m-1) for j in range(m, v+1)]))
    return dict(power=12, spectral_square=pair(Fraction(25, 32)**12),
                old_error_fraction=[1, 16], fault_fraction=[1, 64],
                contraction=[4, 9], output_fraction=pair(Fraction(4, 9)*Fraction(5, 64)+Fraction(1, 64)),
                degree=d, width_per_bit=width, depth=depth, service_locations_per_bit=width*depth,
                archive_threshold_denominator_factor=3072, export_threshold_denominator=768,
                archive_exponential_denominator=32, export_exponential_denominator=8,
                schedule=[schedule(2**r) for r in (1, 3, 6, 10, 20, 40)],
                tails=tails, quantum_decay=32, archive_decay=64, synthesis_decay=24,
                inventory=[4, 8, 1], depth_exponents=[1, 8, 1], active_volume=[5, 16, 2],
                diagnostic_bits=[5, 17, 2], storage_time=[6, 25, 3],
                joint_decay=20, accounting_decay=16)
