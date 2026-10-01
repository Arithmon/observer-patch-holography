"""Actual parent decoder composed with arbitrary code-to-M6 basis faults."""

import numpy as np
from scipy.linalg import expm
from m1_fixed_noise.algebra import decoder
from .interfaces import reduction
from .format import pack


def integer_maps():
    rows = decoder()
    maps = np.zeros((64, 2, 128), dtype=np.int64)
    for s, row in enumerate(rows):
        for logical, entries in enumerate(row['kraus_rows']):
            for index, sign in entries:
                maps[s, logical, index] = sign
    return maps[0].T.copy(), maps


def on_site(matrix, site, code):
    # Physical codeword bits are little endian, as in the parent decoder.
    out = np.zeros_like(code, dtype=np.result_type(matrix.dtype, code.dtype))
    for before in range(128):
        source = (before >> site) & 1
        for target in range(2):
            after = (before & ~(1 << site)) | (target << site)
            out[after] += matrix[target, source]*code[before]
    return out


def basis_census():
    code, maps = integer_maps()
    reductions = [k.real.astype(np.int64) for group in reduction(2) for k in group]
    hist = {}
    total = 0
    for site in range(7):
        for a in range(6):
            for b in range(2):
                error = np.zeros((6, 2), dtype=np.int64)
                error[a, b] = 1
                for reduction_map in reductions:
                    logical = maps@on_site(reduction_map@error, site, code)
                    for row in logical:
                        if not np.array_equal(row, np.eye(2, dtype=np.int64)*row[0, 0]):
                            raise ValueError('not a corrected scalar branch')
                        n = int(row[0, 0])
                        hist[str(n)] = hist.get(str(n), 0)+1
                        total += 1
    return dict(sites=7, full_interface_basis_operators_per_site=12,
                reduction_kraus_branches=5, syndrome_branches=64,
                branch_identities=total, scalar_numerator_histogram=hist,
                scalar_denominator=8)


def noise_operators(kind):
    embedding = np.eye(6, 2, dtype=complex)
    if kind == 'erasure':
        a, b = np.zeros((6, 2), complex), np.zeros((6, 2), complex)
        a[4, 0], b[4, 1] = 1, 1
        return [a, b]
    if kind == 'leak_damping':
        a = embedding.copy()
        a[1, 1] = np.sqrt(.75)
        b = np.zeros((6, 2), complex)
        b[4, 1] = .5
        return [a, b]
    if kind == 'coherent_leak':
        h = np.zeros((6, 6))
        h[1, 4] = h[4, 1] = .43
        return [expm(-1j*h)@embedding]
    raise ValueError('correction noise catalog')


def corrected(maps, code, operators, site, reductions):
    ks = []
    for a in reductions:
        for error in operators:
            ks.extend(maps@on_site(a@error, site, code))
    return ks


def choi(ks):
    return sum(np.outer(k.T.reshape(-1), k.T.reshape(-1).conj()) for k in ks)/2


def candidate():
    code, maps = integer_maps()
    code, maps = code/np.sqrt(8), maps/np.sqrt(8)
    groups = reduction(2)
    reductions = sum(groups, [])
    psi = np.array([np.sqrt(.6), 1j*np.sqrt(.4)])
    effect = np.array([[.35, .12+.05j], [.12-.05j, .75]])
    cases = []
    for kind in ('erasure', 'leak_damping', 'coherent_leak'):
        ops = noise_operators(kind)
        for site in range(7):
            ks = corrected(maps, code, ops, site, reductions)
            valid = corrected(maps, code, ops, site, groups[0])
            output = sum(np.outer(k@psi, (k@psi).conj()) for k in ks)
            valid_output = sum(np.outer(k@psi, (k@psi).conj()) for k in valid)
            acceptance = np.trace(valid_output).real
            cases.append(dict(kind=kind, site=site, noise_operators=[pack(k) for k in ops],
                              recovered_choi=pack(choi(ks)),
                              raw_carrier_acceptance=float(acceptance),
                              corrected_accounting=float(np.trace(effect@output).real),
                              product_abort_accounting=float(np.trace(effect@valid_output).real+1-acceptance)))
    two = np.array([(-1)**((j & 3).bit_count()) for j in range(128)])[:, None]*code
    return dict(basis_census=basis_census(), cases=cases,
                two_site_phase_choi=pack(choi(maps@two)))
