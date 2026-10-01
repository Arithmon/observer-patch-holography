"""Independent Pascal/fraction replay of coherent tails and the combined budget."""

from fractions import Fraction
import math
from .format import close, integer, keys, need, number


def pascal(n):
    rows = [[1]]
    for _ in range(n):
        old = [0]+rows[-1]+[0]
        rows.append([old[i]+old[i+1] for i in range(len(old)-1)])
    return rows


def fraction(row):
    need(type(row) is list and len(row) == 2 and all(type(x) is int for x in row)
         and row[0] >= 0 and row[1] > 0 and math.gcd(*row) == 1, 'canonical positive rational')
    return Fraction(*row)


def verify(row):
    keys(row, 'tails archive_tail quantum_level quantum_fault_power amplitude_decay_denominator '
         'active_locations quantum_q_degree synthesis_request_locations synthesis_per_gate_decay '
         'joint_decay accounting_decay resources model archive_normalization')
    catalog = ((1, 1, 2), (5, 2, 8), (9, 3, 16), (13, 4, 32), (24, 6, 64))
    need(type(row['tails']) is list and len(row['tails']) == len(catalog), 'complete coherent-tail catalog')
    for item, (v, m, den) in zip(row['tails'], catalog):
        keys(item, 'locations required_faults amplitude_denominator coefficient_by_fault_count '
             'coherent_bad_majorant binomial_majorant')
        for key, value in dict(locations=v, required_faults=m, amplitude_denominator=den).items():
            integer(item[key], value, value)
        combinations = pascal(v)
        coefficients = []
        for k in range(v+1):
            value = 0
            for j in range(m, k+1):
                value += (-1 if (j-m) % 2 else 1)*combinations[k][j]*combinations[j-1][m-1]
            need(value == int(k >= m), 'marked-set identity')
            coefficients.append(value)
        need(type(item['coefficient_by_fault_count']) is list
             and all(type(x) is int for x in item['coefficient_by_fault_count'])
             and item['coefficient_by_fault_count'] == coefficients, 'exact emitted path multiplicities')
        numerator, upper = Fraction(0), Fraction(0)
        # Expand both polynomials, independently of the producer's factored
        # binomial upper expression. Each coefficient dominates separately.
        for j in range(m, v+1):
            a = combinations[v][j]*combinations[j-1][m-1]
            b = combinations[v][m]*combinations[v-m][j-m]
            need(a <= b, 'termwise norm majorant')
            numerator += Fraction(a, den**j)
            upper += Fraction(b, den**j)
        need(fraction(item['coherent_bad_majorant']) == numerator
             and fraction(item['binomial_majorant']) == upper
             and numerator <= upper, 'coherent norm tail including multiplicities')
    need(type(row['archive_tail']) is list and len(row['archive_tail']) == 5, 'archive asymptotic controls')
    for item, exponent in zip(row['archive_tail'], (8, 16, 64, 256, 1024)):
        keys(item, 'log2_q r word_length log_tail_bound')
        integer(item['log2_q'], exponent, exponent)
        r = exponent+1
        integer(item['r'], r, r)
        integer(item['word_length'], 4*r+1, 4*r+1)
        q = 1 << exponent
        log_u = math.log(math.log(2*q)**4)-math.log(q)/2
        expected = math.log(q)*5+math.log(math.log(2*q))*8+math.exp(log_u)+log_u*(r+1)
        close(number(item['log_tail_bound']), expected, 'amplitude archive tail')
    for name, expected in dict(quantum_level=6, quantum_fault_power=64, amplitude_decay_denominator=2,
                               quantum_q_degree=5-64//2, synthesis_per_gate_decay=24,
                               joint_decay=24-4, accounting_decay=24-4-4).items():
        integer(row[name], expected, expected)
    for name, expected in dict(active_locations=[5, 16], synthesis_request_locations=[4, 8]).items():
        need(type(row[name]) is list and all(type(x) is int for x in row[name]) and row[name] == expected,
             'combined program bound '+name)
    resource = dict(active_inventory=[4, 8], depth=[1, 8], active_volume=[5, 16],
                    diagnostic_bits=[5, 17], physical_storage_volume=[6, 25])
    keys(row['resources'], 'active_inventory depth active_volume diagnostic_bits physical_storage_volume')
    for name, expected in resource.items():
        need(type(row['resources'][name]) is list and all(type(x) is int for x in row['resources'][name])
             and row['resources'][name] == expected, 'complete noisy-record lifetime '+name)
    need([a+b for a, b in zip(resource['active_inventory'], resource['depth'])] == resource['active_volume']
         and [a+b for a, b in zip(resource['diagnostic_bits'], resource['depth'])] == resource['physical_storage_volume'],
         'lifetime product identities')
    need(row['model'] == 'full-interface Markov quantum generators and local classical Poisson faults',
         'composed noise-model boundary')
    need(row['archive_normalization'] == 'C=Cprime=1 illustrative only; not a native threshold or onset',
         'normalized examples are not hardware thresholds')
