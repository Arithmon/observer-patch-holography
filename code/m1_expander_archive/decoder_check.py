"""Independent packed-bit finite differences and complete radius rejection."""

from functools import lru_cache
from .format import keys, need, integer


@lru_cache(maxsize=1)
def masks():
    # Each integer is a truth table; no NumPy reductions or producer routines.
    rows = []
    for degree in reversed(range(6)):
        for monomial in range(2048):
            if monomial.bit_count() == degree:
                cube = sum(1 << x for x in range(2048) if x & monomial == 0)
                truth = sum(1 << x for x in range(2048) if x & monomial == monomial)
                rows.append((monomial, cube, truth))
    return rows


def reference(received):
    need(type(received) is int and 0 <= received < 1 << 2047, 'punctured word integer')
    candidates = {}
    for fill in (0, 1):
        residual, answer, coefficients = received | (fill << 2047), 0, []
        for monomial, cube, truth in masks():
            derivative = residual
            for j in range(11):
                if monomial & (1 << j):
                    derivative ^= derivative >> (1 << j)
            votes = (derivative & cube).bit_count()
            if 2*votes > cube.bit_count():
                answer ^= truth
                residual ^= truth
                coefficients.append(monomial)
        distance = ((answer ^ received) & ((1 << 2047)-1)).bit_count()
        if distance <= 31:
            candidates[tuple(coefficients)] = dict(logical=answer >> 2047, distance=distance,
                                                  coefficients=coefficients, failure=False)
    need(len(candidates) <= 1, 'unique radius-31 punctured decoding')
    return next(iter(candidates.values()), dict(logical=0, distance=None, coefficients=[], failure=True))


def verify(rows):
    need(type(rows) is list and len(rows) == 27, 'complete actual-code decoder catalog')
    bases = ([], [0], [31], [0, 3, 37, 1024, 1539], [1, 2, 4, 8, 16, 31])
    catalog = [(list(base), sorted((73*j+101*i) % 2047 for j in range(count)), 'within_radius')
               for i, base in enumerate(bases) for count in (0, 1, 15, 30, 31)]
    catalog += [([], [x for x in range(2047) if x & 31 == 31][:32], 'wrong_label_control'),
                ([], sorted(73*j % 2047 for j in range(32)), 'retained_failure_control')]
    for row, (monomials, errors, kind) in zip(rows, catalog):
        keys(row, 'monomials errors kind decoded')
        need(type(row['monomials']) is list and row['monomials'] == monomials
             and all(type(x) is int for x in row['monomials']), 'frozen polynomial')
        need(type(row['errors']) is list and row['errors'] == errors
             and all(type(x) is int for x in row['errors']), 'frozen error set')
        need(row['kind'] == kind, 'decoder control classification')
        encoded = sum((sum((x & mask) == mask for mask in monomials) % 2) << x for x in range(2048))
        received = (encoded & ((1 << 2047)-1)) ^ sum(1 << x for x in errors)
        expected = reference(received)
        result = row['decoded']
        keys(result, 'logical distance coefficients failure')
        integer(result['logical'], 0, 1)
        need(type(result['failure']) is bool and type(result['coefficients']) is list
             and all(type(x) is int for x in result['coefficients']), 'typed decoder result')
        need(result['distance'] is None or type(result['distance']) is int, 'distance type')
        need(result == expected, 'independent full 2047-bit decoding')
        if kind == 'within_radius':
            need(not result['failure'] and result['logical'] == encoded >> 2047
                 and result['distance'] == len(errors), 'correct radius-31 logical record')
        elif kind == 'wrong_label_control':
            need(not result['failure'] and result['logical'] == 1 and result['distance'] == 31,
                 '32-error wrong-label witness must survive')
        else:
            need(result['failure'] and result['logical'] == 0, 'complete fallback, not postselection')
