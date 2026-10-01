"""Integer/fraction proof controls; no producer formulas are imported."""

from fractions import Fraction
import math
from .format import keys, need, integer


def fraction(row):
    need(type(row) is list and len(row) == 2 and all(type(x) is int for x in row)
         and row[1] > 0, 'rational scalar')
    x = Fraction(*row)
    need(row == [x.numerator, x.denominator], 'reduced rational')
    return x


def verify(row):
    keys(row, 'power spectral_square old_error_fraction fault_fraction contraction output_fraction '
              'degree width_per_bit depth service_locations_per_bit archive_threshold_denominator_factor '
              'export_threshold_denominator archive_exponential_denominator export_exponential_denominator '
              'schedule tails quantum_decay archive_decay synthesis_decay inventory depth_exponents '
              'active_volume diagnostic_bits storage_time joint_decay accounting_decay')
    integer(row['power'], 12, 12)
    gap = Fraction(1)
    d = 1
    for _ in range(12):
        gap *= Fraction(50, 64)
        d *= 8
    need(fraction(row['spectral_square']) == gap and gap < Fraction(1, 16), 'powered gap')
    alpha, tau = fraction(row['old_error_fraction']), fraction(row['fault_fraction'])
    need(alpha == Fraction(1, 16) and tau == Fraction(1, 64), 'fixed live error contract')
    contraction = Fraction(1, 16)/Fraction(3, 8)**2
    need(fraction(row['contraction']) == contraction, 'spectral contraction')
    output = contraction*(alpha+tau)+tau
    need(fraction(row['output_fraction']) == output and output < alpha and alpha+tau < Fraction(1, 2),
         'refresh closes and majority stays live')
    depth = 2+d+(d*d-d)//2
    for key, value in [('degree', d), ('width_per_bit', d+2), ('depth', depth),
                       ('service_locations_per_bit', (d+2)*depth),
                       ('archive_threshold_denominator_factor', 3072),
                       ('export_threshold_denominator', 768),
                       ('archive_exponential_denominator', 32), ('export_exponential_denominator', 8),
                       ('quantum_decay', 32), ('archive_decay', 64), ('synthesis_decay', 24)]:
        integer(row[key], value, value)
    # e<3: at the rational thresholds, each tail base is <=1/16.
    # Convexity of 1/x on [1,2] gives log(2)>2/3, hence log(16)>8/3.
    need(Fraction(8, 3)/64-Fraction(1, 3072) > Fraction(1, 32)
         and Fraction(8, 3)/16-Fraction(1, 768) > Fraction(1, 8),
         'fixed-threshold exponential constants')
    need(type(row['schedule']) is list and len(row['schedule']) == 6, 'complete schedule catalog')
    for item, r in zip(row['schedule'], (1, 3, 6, 10, 20, 40)):
        keys(item, 'q logarithm_ceiling m bits levels fault_power')
        for x in item.values():
            integer(x, 1, 2**50)
        need(item['q'] == 2**r and item['logarithm_ceiling'] == r+1, 'dyadic schedule input')
        m, k = item['m'], item['levels']
        need((m-1)**2 < 2048*(r+1) <= m*m == item['bits'], 'integer archive ceiling')
        need(2**(k-1) < 16*(r+1) <= 2**k == item['fault_power'], 'integer protection ceiling')
        need(2*item['fault_power'] >= 32*(r+1), 'fixed theta=1/4 quantum error')
        need(item['bits'] >= 2048*(r+1), 'archive refinement decay using log(2)<1')
    need(type(row['tails']) is list and len(row['tails']) == 4, 'complete marked-tail controls')
    for item, (v, m) in zip(row['tails'], ((5, 2), (9, 3), (16, 4), (24, 6))):
        keys(item, 'locations required coefficients')
        integer(item['locations'], v, v)
        integer(item['required'], m, m)
        need(type(item['coefficients']) is list and len(item['coefficients']) == v-m+1
             and all(type(x) is int for x in item['coefficients']), 'integer tail coefficients')
        # Solve the triangular path-multiplicity equations rather than using
        # the producer's closed formula.
        solved = []
        for k in range(m, v+1):
            solved.append(1-sum(math.comb(k, m+j)*a for j, a in enumerate(solved)))
        need(item['coefficients'] == solved, 'independent inclusion-exclusion coefficients')
    vectors = {'inventory': [4, 8, 1], 'depth_exponents': [1, 8, 1]}
    vectors['active_volume'] = [a+b for a, b in zip(vectors['inventory'], vectors['depth_exponents'])]
    vectors['diagnostic_bits'] = [vectors['active_volume'][0], vectors['active_volume'][1]+1,
                                  vectors['active_volume'][2]]
    vectors['storage_time'] = [a+b for a, b in zip(vectors['diagnostic_bits'], vectors['depth_exponents'])]
    for key, value in vectors.items():
        need(type(row[key]) is list and all(type(x) is int for x in row[key]) and row[key] == value,
             'charged lifetime '+key)
    error = min(row['quantum_decay']-vectors['active_volume'][0],
                row['archive_decay']-vectors['active_volume'][0],
                row['synthesis_decay']-vectors['inventory'][0])
    integer(row['joint_decay'], error, error)
    integer(row['accounting_decay'], error-4, error-4)
