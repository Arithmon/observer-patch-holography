"""Norm-tail and lifetime evidence for composition with noisy public records."""

from fractions import Fraction
import math


def candidate():
    tails = []
    for locations, required, denominator in ((1, 1, 2), (5, 2, 8), (9, 3, 16), (13, 4, 32), (24, 6, 64)):
        eta = Fraction(1, denominator)
        coefficients = [sum((-1)**(j-required)*math.comb(k, j)*math.comb(j-1, required-1)
                            for j in range(required, k+1)) for k in range(locations+1)]
        tail = sum(math.comb(locations, j)*math.comb(j-1, required-1)*eta**j
                   for j in range(required, locations+1))
        upper = math.comb(locations, required)*eta**required*(1+eta)**(locations-required)
        tails.append(dict(locations=locations, required_faults=required, amplitude_denominator=denominator,
                          coefficient_by_fault_count=coefficients,
                          coherent_bad_majorant=[tail.numerator, tail.denominator],
                          binomial_majorant=[upper.numerator, upper.denominator]))
    budgets = []
    for power in (8, 16, 64, 256, 1024):
        # q=2^power, C=C'=1: a normalized example, not a measured threshold.
        r = power+1
        log_q, log_2q = power*math.log(2), (power+1)*math.log(2)
        log_u = 4*math.log(log_2q)-log_q/2
        log_bound = 5*log_q+8*math.log(log_2q)+math.exp(log_u)+(r+1)*log_u
        budgets.append(dict(log2_q=power, r=r, word_length=4*r+1, log_tail_bound=log_bound))
    return dict(tails=tails, archive_tail=budgets,
                quantum_level=6, quantum_fault_power=64, amplitude_decay_denominator=2,
                active_locations=[5, 16], quantum_q_degree=-27,
                synthesis_request_locations=[4, 8], synthesis_per_gate_decay=24,
                joint_decay=20, accounting_decay=16,
                resources=dict(active_inventory=[4, 8], depth=[1, 8], active_volume=[5, 16],
                               diagnostic_bits=[5, 17], physical_storage_volume=[6, 25]),
                model='full-interface Markov quantum generators and local classical Poisson faults',
                archive_normalization='C=Cprime=1 illustrative only; not a native threshold or onset')
