"""Independent exact parity construction and complete logical-channel oracle."""

import math
import numpy as np
from .format import close, integer, keys, need, number, unpack


def integer_code():
    checks = [[j for j in range(7) if (j+1) >> bit & 1] for bit in range(3)]
    code = np.zeros((128, 2), dtype=np.int64)
    for word in range(128):
        if all(sum((word >> j) & 1 for j in positions) % 2 == 0 for positions in checks):
            code[word, word.bit_count() % 2] = 1
    maps = np.zeros((64, 2, 128), dtype=np.int64)
    for s in range(64):
        flip = 1 << (s//8-1) if s//8 else 0
        phase = 1 << (s % 8-1) if s % 8 else 0
        for word in range(128):
            original = word ^ flip
            maps[s, :, word] = code[original]*(-1)**((original & phase).bit_count() % 2)
    return code, maps


def local_action(op, site, code):
    # Tensor-axis contraction is independent of the producer's integer-index
    # rewriting. Reverse the axis because the canonical code is little endian.
    tensor = code.reshape((2,)*7+(2,))
    moved = np.moveaxis(tensor, 6-site, 0)
    changed = np.tensordot(op, moved, axes=(1, 0))
    return np.moveaxis(changed, 0, 6-site).reshape(128, 2)


def exact_basis_check(code, maps):
    histogram = {}
    for site in range(7):
        for a in range(6):
            for b in range(2):
                for branch in range(5):
                    op = np.zeros((2, 2), dtype=np.int64)
                    if branch == 0 and a < 2:
                        op[a, b] = 1
                    elif branch and a == branch+1:
                        op[0, b] = 1
                    images = np.einsum('sij,jk->sik', maps, local_action(op, site, code))
                    for image in images:
                        scalar = int(image[0, 0])
                        need(image[0, 1] == image[1, 0] == 0 and image[1, 1] == scalar,
                             'every full-interface basis branch preserves the logical qubit')
                        histogram[str(scalar)] = histogram.get(str(scalar), 0)+1
    return histogram


def reduced_operators(full):
    # Each environment branch is retained separately. Coherently adding the
    # four leakage rows would be a different, generally non-TP channel.
    result = [(0, op[:2]) for op in full]
    for row in range(2, 6):
        for op in full:
            k = np.zeros((2, 2), complex)
            k[0] = op[row]
            result.append((1, k))
    return result


def verify(row):
    keys(row, 'basis_census cases two_site_phase_choi')
    code, maps = integer_code()
    close(code.T@code, 8*np.eye(2), 'code normalization')
    close(sum(k.T@k for k in maps), 8*np.eye(128), 'all syndrome spaces')
    histogram = exact_basis_check(code, maps)
    expected = dict(sites=7, full_interface_basis_operators_per_site=12,
                    reduction_kraus_branches=5, syndrome_branches=64,
                    branch_identities=7*12*5*64, scalar_numerator_histogram=histogram,
                    scalar_denominator=8)
    keys(row['basis_census'], 'sites full_interface_basis_operators_per_site reduction_kraus_branches '
         'syndrome_branches branch_identities scalar_numerator_histogram scalar_denominator')
    need(row['basis_census'] == expected and all(type(v) is int for k, v in row['basis_census'].items()
                                              if k != 'scalar_numerator_histogram')
         and all(type(v) is int for v in row['basis_census']['scalar_numerator_histogram'].values()),
         'exact complete basis census')
    code, maps = code/np.sqrt(8), maps/np.sqrt(8)
    need(type(row['cases']) is list and len(row['cases']) == 21, 'all site/noise cases')
    bell = np.array([1, 0, 0, 1])/np.sqrt(2)
    psi = np.array([math.sqrt(.6), 1j*math.sqrt(.4)])
    effect = np.array([[.35, .12+.05j], [.12-.05j, .75]])
    ideal_accounting = float(np.vdot(psi, effect@psi).real)
    for index, item in enumerate(row['cases']):
        keys(item, 'kind site noise_operators recovered_choi raw_carrier_acceptance '
             'corrected_accounting product_abort_accounting')
        kind, site = ('erasure', 'leak_damping', 'coherent_leak')[index//7], index % 7
        need(item['kind'] == kind, 'correction case order')
        integer(item['site'], site, site)
        need(type(item['noise_operators']) is list and len(item['noise_operators']) == (1 if kind == 'coherent_leak' else 2),
             'complete physical noise operators')
        operators = [unpack(k, (6, 2)) for k in item['noise_operators']]
        frozen = []
        if kind == 'erasure':
            for b in range(2):
                op = np.zeros((6, 2), complex)
                op[4, b] = 1
                frozen.append(op)
            acceptance = 0.
        elif kind == 'leak_damping':
            first, second = np.eye(6, 2, dtype=complex), np.zeros((6, 2), complex)
            first[1, 1], second[4, 1] = math.sqrt(.75), .5
            frozen = [first, second]
            acceptance = .875
        else:
            first = np.eye(6, 2, dtype=complex)
            first[1, 1], first[4, 1] = math.cos(.43), -1j*math.sin(.43)
            frozen = [first]
            acceptance = 1-math.sin(.43)**2/2
        for a, b in zip(operators, frozen):
            close(a, b, 'actual full-interface fault')
        close(sum(k.conj().T@k for k in operators), np.eye(2), 'noise TP before correction')
        actual = np.zeros((4, 4), complex)
        for _, k in reduced_operators(operators):
            for branch in maps@local_action(k, site, code):
                vector = branch.T.reshape(-1)/np.sqrt(2)
                actual += np.outer(vector, vector.conj())
        close(actual, np.outer(bell, bell), 'complete reference-preserving leakage correction')
        close(unpack(item['recovered_choi'], (4, 4)), actual, 'emitted recovered Choi state')
        close(number(item['raw_carrier_acceptance']), acceptance, 'raw carrier rejection is not logical failure')
        close(number(item['corrected_accounting']), ideal_accounting, 'corrected positive accounting')
        close(number(item['product_abort_accounting']), acceptance*ideal_accounting+1-acceptance,
              'obsolete whole-carrier abort penalty')
    # This decoder corrects one arbitrary error, not every possible fault.
    z_bell = np.array([1, 0, 0, -1])/np.sqrt(2)
    two = code.copy()
    for j in range(128):
        two[j] *= (-1)**((j & 3).bit_count() % 2)
    output = np.zeros((4, 4), complex)
    for k in maps@two:
        vector = k.T.reshape(-1)/np.sqrt(2)
        output += np.outer(vector, vector.conj())
    close(output, np.outer(z_bell, z_bell), 'two-site uncorrectable witness')
    close(unpack(row['two_site_phase_choi'], (4, 4)), output, 'actual two-site channel')
