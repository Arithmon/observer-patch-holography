"""Majority-logic decoding of the actual punctured RM(5,11) record code."""

import numpy as np

M, R = 11, 5
LENGTH = (1 << M)-1
MONOMIALS = tuple(mask for degree in range(R, -1, -1)
                  for mask in range(1 << M) if mask.bit_count() == degree)
POINTS = np.arange(1 << M)


def polynomial(monomials):
    word = np.zeros(1 << M, dtype=np.uint8)
    for mask in monomials:
        if type(mask) is not int or not 0 <= mask < 1 << M or mask.bit_count() > R:
            raise ValueError('RM(5,11) monomial')
        word ^= (POINTS & mask) == mask
    return word


def decode(received):
    if (type(received) is not list or len(received) != LENGTH
            or any(type(x) is not int or x not in (0, 1) for x in received)):
        raise ValueError('2047 ordinary bits required')
    candidates = {}
    for fill in (0, 1):
        residual = np.array(received+[fill], dtype=np.uint8)
        decoded = np.zeros(1 << M, dtype=np.uint8)
        coefficients = []
        for mask in MONOMIALS:
            axes = tuple(M-1-j for j in range(M) if mask >> j & 1)
            derivatives = np.bitwise_xor.reduce(residual.reshape((2,)*M), axis=axes)
            if 2*int(derivatives.sum()) > derivatives.size:
                value = (POINTS & mask) == mask
                decoded ^= value
                residual ^= value
                coefficients.append(mask)
        distance = int(np.count_nonzero(decoded[:-1] != received))
        if distance <= 31:
            candidates[tuple(coefficients)] = dict(logical=int(decoded[-1]), distance=distance,
                                                  coefficients=coefficients, failure=False)
    if len(candidates) > 1:
        raise ValueError('violated punctured code distance')
    return next(iter(candidates.values()), dict(logical=0, distance=None, coefficients=[], failure=True))


def specifications():
    rows = []
    polynomials = ((), (0,), (31,), (0, 3, 37, 1024, 1539), (1, 2, 4, 8, 16, 1023 & 31))
    for index, monomials in enumerate(polynomials):
        for errors in (0, 1, 15, 30, 31):
            positions = sorted({(73*j+101*index) % LENGTH for j in range(errors)})
            rows.append((list(monomials), positions, 'within_radius'))
    support = [x for x in range(LENGTH) if x & 31 == 31]
    rows.append(([], support[:32], 'wrong_label_control'))
    rows.append(([], sorted(73*j % LENGTH for j in range(32)), 'retained_failure_control'))
    return rows


def candidate():
    result = []
    for monomials, errors, kind in specifications():
        received = polynomial(monomials)[:-1].astype(int).tolist()
        for i in errors:
            received[i] ^= 1
        result.append(dict(monomials=monomials, errors=errors, kind=kind, decoded=decode(received)))
    return result
