"""Finite native conservation from faithful invariance and a Dirichlet bound.

The channel forgets a new attempt's records only to form its Heisenberg
pullback. Its conservation conclusion also applies to source-only observables
in the complete classical-feedback instrument. It supplies no stationary
density on a growing tape, physical time or active quantum recovery.
"""
from __future__ import annotations

from fractions import Fraction as F
import json

import numpy as np

from ..local_repair_family import local_bias as local, mixing
from ..observer import quantum_mixing


DIMENSION = 66


def _parameters(parameters):
    """Use the owning finite-family admission and detach exact payloads."""
    ports, seams, cn, cd, rn, rd = mixing._parameters(parameters)
    return local.Parameters(ports, seams, F(cn, cd), F(rn, rd))


def _weights(parameters, c_weights):
    parameters = _parameters(parameters)
    if c_weights is None:
        return parameters, (parameters.c_weight,) * len(local.SEAMS)
    if type(c_weights) is not tuple or len(c_weights) != len(local.SEAMS):
        raise ValueError('one exact C-mixture weight per native seam required')
    # The frozen global mixture is a special case. The theorem additionally
    # permits different source-independent law choices at each seam.
    return parameters, tuple(mixing._fraction(weight) for weight in c_weights)


def branches(parameters, c_weights=None):
    """Complete weighted partial injections, including zero laws and hold."""
    parameters, weights = _weights(parameters, c_weights)
    result = []
    for edge, seam_weight, c_weight in zip(local.SEAMS, parameters.seam_weights,
                                         weights, strict=True):
        chance = parameters.repair_probability * F(seam_weight, sum(parameters.seam_weights))
        for law, selection in (('C', c_weight), ('E', 1-c_weight)):
            for coin, raw, pairs in local.raw_maps(law, edge):
                probability = chance * selection * local.coin_probability(parameters, edge, coin)
                result.append((probability, pairs))
    result.append((1-parameters.repair_probability, tuple((i, i) for i in range(DIMENSION))))
    return tuple(result)


def faithful_channel_certificate(parameters, c_weights=None):
    """Exact completeness and faithful invariance of the actual H66 maps."""
    parameters, weights = _weights(parameters, c_weights)
    pi = local.reference(parameters)
    effects, output = [F()] * DIMENSION, [F()] * DIMENSION
    raw = branches(parameters, weights)
    for probability, pairs in raw:
        if len({i for i, _ in pairs}) != len(pairs) or len({j for _, j in pairs}) != len(pairs):
            raise ArithmeticError('native branch is not a partial injection')
        for i, j in pairs:
            effects[i] += probability
            output[j] += probability*pi[i]
    if effects != [F(1)]*DIMENSION or tuple(output) != pi or min(pi) <= 0:
        raise ArithmeticError('complete source channel lost its faithful invariant density')
    return {'dimension': DIMENSION, 'raw_branches_including_zero_and_hold': len(raw),
            'positive_branches': sum(bool(w and pairs) for w, pairs in raw),
            'complete': True, 'faithful_invariant_reference': True,
            'smallest_reference_weight': str(min(pi)),
            'largest_reference_weight': str(max(pi))}


def extra_dephasing_certificate():
    """C-E is diagonal and nonnegative on H66 matrix units per seam/coin.

    This is a Hilbert-Schmidt quadratic-form comparison of forgetting
    channels, not an assertion that C-E is a positive or CP map.
    """
    checked, differences = 0, 0
    for edge in local.SEAMS:
        maps = {law: local.raw_maps(law, edge) for law in ('C', 'E')}
        for coin in (0, 1):
            actions = {}
            for law in ('C', 'E'):
                values = {}
                for branch_coin, _, pairs in maps[law]:
                    if branch_coin != coin:
                        continue
                    for i, out_i in pairs:
                        for j, out_j in pairs:
                            key = (i, j)
                            if key in values:
                                raise ArithmeticError('one input unit appears in multiple raw slots')
                            values[key] = (out_i, out_j)
                actions[law] = values
            if any(actions['C'].get(key) != value for key, value in actions['E'].items()):
                raise ArithmeticError('E does not refine C on this complete quantum branch')
            for key in actions['C'].keys()-actions['E'].keys():
                if actions['C'][key] != key:
                    raise ArithmeticError('extra C coherence is not an unchanged matrix unit')
                differences += 1
            checked += DIMENSION**2
    return {'matrix_unit_coin_seam_checks': checked,
            'additional_C_matrix_unit_terms': differences,
            'C_minus_E_diagonal_nonnegative_in_Hilbert_Schmidt_basis': True,
            'C_minus_E_positive_map_claimed': False}


def uniform_gap_certificate():
    """A rational gap from existing exact population and coherence controls."""
    fair = local.Parameters(c_weight=F(1))
    minor = mixing.population_minorization(fair, block=20)
    population = mixing.probability_kernel(fair)
    if any(population[i][j] != population[j][i]
           for i in range(DIMENSION) for j in range(DIMENSION)):
        raise ArithmeticError('the fair population comparison is not symmetric')
    if not 0 < minor['beta'] <= 1:
        raise ArithmeticError('no positive population minorization')
    # Actual off-diagonal channel blocks are symmetric, entrywise
    # nonnegative and have constant column sums at most5/6. Their spectral
    # maximum is consequently at most5/6 by the absolute row-sum bound.
    blocks, units = 0, 0
    for sector in quantum_mixing.sectors(binary=True):
        if sector.difference == quantum_mixing.ZERO:
            continue
        matrix = quantum_mixing.native_block('C', sector.difference, binary=True)
        if matrix != tuple(zip(*matrix, strict=True)):
            raise ArithmeticError('fair native coherence block is not symmetric')
        if any(x < 0 for row in matrix for x in row) or max(map(sum, matrix)) > F(5, 6):
            raise ArithmeticError('fair native coherence block exceeds the proved ceiling')
        blocks += 1
        units += len(sector.pairs)
    gap = min(minor['beta']/minor['block'], F(1, 6))
    return {'population_block': minor['block'],
            'exact_population_minorization': str(minor['beta']),
            'coherence_blocks': blocks, 'offdiagonal_matrix_units': units,
            'coherence_gap_lower': '1/6', 'uniform_gap_lower': str(gap),
            'gap_is_rational_lower_bound_not_fitted_eigenvalue': True}


def gap_bound(parameters, c_weights=None):
    """Weighted GNS Poincare bound, uniform over the per-seam law mixture."""
    parameters, _ = _weights(parameters, c_weights)
    pi = local.reference(parameters)
    smallest = min(F(weight, sum(parameters.seam_weights)) * local.coin_probability(parameters, edge, coin)
                   for edge, weight in zip(local.SEAMS, parameters.seam_weights, strict=True)
                   for coin in (0, 1))
    multiplier = min(pi)/max(pi) * 60 * parameters.repair_probability * smallest
    base = F(uniform_gap_certificate()['uniform_gap_lower'])
    return multiplier * base


def seam_GNS_certificate(parameters):
    """Check each actual seam map is a positive GNS self-adjoint contraction.

    ``positive`` concerns its Hilbert-space quadratic form, separately from
    the complete positivity of the Heisenberg channel. Exact native matrix
    units split into one- or two-dimensional blocks. Their weighted Gram
    matrices and identity-minus-map matrices have nonnegative principal
    minors. No floating spectrum or fitted similarity transform is used.
    """
    parameters = _parameters(parameters)
    pi = local.reference(parameters)
    diagonal_blocks, two_blocks = 0, 0
    for edge in local.SEAMS:
        for law in ('C', 'E'):
            entries = {}
            for coin, _, pairs in local.raw_maps(law, edge):
                weight = local.coin_probability(parameters, edge, coin)
                for i, out_i in pairs:
                    for j, out_j in pairs:
                        key = (i*DIMENSION+j, out_i*DIMENSION+out_j)
                        entries[key] = entries.get(key, F())+weight
            neighbors = {}
            for i, j in entries:
                if i != j:
                    neighbors.setdefault(i, set()).add(j)
                    neighbors.setdefault(j, set()).add(i)
            seen = set()
            for i in range(DIMENSION**2):
                if i in seen:
                    continue
                component = {i} | neighbors.get(i, set())
                if len(component) > 2 or any(not neighbors.get(j, set()) <= component for j in component):
                    raise ArithmeticError('a native seam has a larger than two-unit GNS block')
                seen.update(component)
                indices = sorted(component)
                gram = [[pi[a % DIMENSION]*entries.get((a, b), F())
                         for b in indices] for a in indices]
                complement = [[(pi[a % DIMENSION] if a == b else F())-gram[u][v]
                               for v, b in enumerate(indices)] for u, a in enumerate(indices)]
                for matrix in (gram, complement):
                    if any(matrix[a][b] != matrix[b][a] for a in range(len(indices)) for b in range(len(indices))):
                        raise ArithmeticError('native map violates weighted self-adjointness')
                    if any(matrix[a][a] < 0 for a in range(len(indices))):
                        raise ArithmeticError('native GNS quadratic form has a negative diagonal')
                    if len(indices) == 2 and matrix[0][0]*matrix[1][1] < matrix[0][1]**2:
                        raise ArithmeticError('native GNS block has a negative principal minor')
                diagonal_blocks += len(indices) == 1
                two_blocks += len(indices) == 2
    return {'seam_law_maps': 2*len(local.SEAMS),
            'one_dimensional_blocks': diagonal_blocks,
            'two_dimensional_blocks': two_blocks,
            'matrix_units_checked': diagonal_blocks+2*two_blocks,
            'weighted_self_adjoint_positive_contractions': True,
            'physical_time_or_stationary_archive_claimed': False}


def positive_diagonal_Kraus_span(parameters, c_weights=None):
    """An invertible positive diagonal lies in the supported Kraus span.

    The integer sum rescales each nonzero scalar-weighted native branch.
    It is an algebraic span witness, not an available coherent operation.
    """
    counts = [0]*DIMENSION
    number = 0
    for probability, pairs in branches(parameters, c_weights):
        if probability and pairs and all(i == j for i, j in pairs):
            number += 1
            for i, _ in pairs:
                counts[i] += 1
    if min(counts) <= 0:
        raise ArithmeticError('supported diagonal Kraus sum is not strictly positive')
    return {'diagonal_branches': number,
            'minimum_positive_diagonal': min(counts),
            'maximum_positive_diagonal': max(counts),
            'coherent_sum_implemented': False}


def apply_heisenberg(parameters, observable, c_weights=None):
    """Numerical check on the exact admitted channel; no fitted spectrum."""
    matrix = np.asarray(observable, dtype=complex)
    if matrix.shape != (DIMENSION, DIMENSION) or not np.isfinite(matrix).all():
        raise ValueError('one finite complex H66 observable required')
    result = np.zeros_like(matrix)
    for probability, pairs in branches(parameters, c_weights):
        if not probability or not pairs:
            continue
        inputs, outputs = zip(*pairs, strict=True)
        result[np.ix_(inputs, inputs)] += float(probability)*matrix[np.ix_(outputs, outputs)]
    return result


def dirichlet_control(parameters, observable, c_weights=None):
    """Independent commutator versus channel trace pairing diagnostics."""
    parameters, weights = _weights(parameters, c_weights)
    matrix = np.asarray(observable, dtype=complex)
    evolved = apply_heisenberg(parameters, matrix, weights)
    pi = np.asarray([float(x) for x in local.reference(parameters)])
    left = float(np.sum(pi*np.sum(matrix.conj()*(matrix-evolved), axis=0)).real)
    right = 0.0
    for probability, pairs in branches(parameters, weights):
        if not probability or not pairs:
            continue
        inputs, outputs = zip(*pairs, strict=True)
        # [X,M]=XM-MX, with M_out,in=1. Avoid dense matrix products.
        commutator = np.zeros_like(matrix)
        commutator[:, inputs] += matrix[:, outputs]
        commutator[outputs, :] -= matrix[inputs, :]
        right += float(probability)*float(np.sum(pi*np.sum(abs(commutator)**2, axis=0)))/2
    mean = np.dot(pi, np.diag(matrix))
    centered = matrix-mean*np.eye(DIMENSION)
    variance = float(np.sum(pi*np.sum(abs(centered)**2, axis=0)))
    residual = float(np.sum(pi*np.sum(abs(matrix-evolved)**2, axis=0)))
    return {'dirichlet_trace_pairing': left,
            'half_weighted_commutator_sum': right,
            'weighted_variance': variance,
            'weighted_fixed_point_residual_squared': residual}


def certificate():
    canonical = mixing.canonical_parameters()
    base = uniform_gap_certificate()
    gamma = F(base['uniform_gap_lower'])
    return {'schema': 'oph.native_fixed_algebra.v1',
            'channel': faithful_channel_certificate(canonical),
            'seam_GNS_contractions': seam_GNS_certificate(canonical),
            'positive_diagonal_Kraus_span': positive_diagonal_Kraus_span(canonical),
            'extra_E_dephasing': extra_dephasing_certificate(),
            'fair_reference_gap': base,
            'canonical_GNS_gap_lower': str(gap_bound(canonical)),
            'whole_bounded_family_GNS_gap_lower': str(gamma/10240),
            'fixed_source_observable_algebra': 'complex scalars on H66',
            'source_only_classical_feedback_conservation_requires_every_conditional_map': True,
            'exogenous_non_IID_product_contraction_uses_common_endpoint_weights': True,
            'joint_record_feedback_mixing_claimed': False,
            'complete_joint_history_stationarity_assumed': False,
            'active_quantum_feedback_or_history_charges_excluded': False,
            'time_parameter': 'native attempt count',
            'physical_particle_mass_or_clock_identified': False}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
