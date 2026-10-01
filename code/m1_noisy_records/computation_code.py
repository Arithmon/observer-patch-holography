"""Explicit punctured RM(5,11) CSS code, with compact generator custody."""

import hashlib
import math


def rows(variables, degree):
    return [sum(1 << x for x in range(1 << variables) if x & monomial == monomial)
            for monomial in range(1 << variables) if monomial.bit_count() <= degree]


def candidate():
    m, r = 11, 5
    generators = rows(m, r)
    raw = b''.join(x.to_bytes((1 << m)//8, 'little') for x in generators)
    return dict(variables=m, degree=r, extended_length=1 << m,
                extended_dimension=sum(math.comb(m, i) for i in range(r+1)),
                extended_distance=1 << (m-r), quantum_length=(1 << m)-1,
                quantum_dimension=1, quantum_distance=(1 << (m-r))-1,
                corrects=((1 << (m-r))-2)//2,
                procedure_spread_bound=4, correction_then_gate_spread_bound=16,
                tolerated_faults_per_rectangle=1,
                generator_sha256=hashlib.sha256(raw).hexdigest())
