"""Faithful conservation and contraction for both native exchange lifts.

The negative lift includes its declared record-selected phase service.
One step here is a COMPLETE chosen-character channel application, not a
claim of equal pulse cost, physical duration, or arbitrary coherent recovery.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import json

import numpy as np

from . import fixed
from ..exchange_sector import statistics
from ..local_repair_family import local_bias as local, mixing


DIMENSION = 66


def _character(character):
    if type(character) is not int or character not in (-1, 1):
        raise ValueError('integer exchange character +1 or -1 required')
    return character


def _rows(law, seam, character):
    """Complete physical raw slots, retaining zero and abort addresses."""
    _character(character)
    rows = []
    for coin in (0, 1):
        for raw in range(4 if law == 'C' else 8):
            row = statistics.branch(law, seam, coin, raw, character)
            if ((row.law, row.seam, row.coin, row.raw) != (law, seam, coin, raw)
                    or len(row.inputs) != len(row.outputs) or len(row.inputs) != len(row.signs)
                    or len(set(row.inputs)) != len(row.inputs)
                    or len(set(row.outputs)) != len(row.outputs)):
                raise ArithmeticError('signed native branch lost its complete partial-injection contract')
            for i, j, sign in zip(row.inputs, row.outputs, row.signs, strict=True):
                if (type(i) is not int or type(j) is not int or not 0 <= i < 66 or not 0 <= j < 66
                        or type(sign) is not int or sign not in (-1, 1)):
                    raise ArithmeticError('invalid signed native column')
                source, target = statistics.PAIRS[i], statistics.PAIRS[j]
                # Independent ordered-Fock parity, not a scalar-per-record gauge.
                between = sum(seam[0] < p < seam[1] for p in source)
                expected = character**between if source != target else 1
                if sign != expected:
                    raise ArithmeticError('column phase differs from the selected exchange character')
            rows.append(row)
    return tuple(rows)


@dataclass(frozen=True)
class WeightedBranch:
    probability: F
    columns: tuple[tuple[int, int, int], ...]
    label: tuple


def branches(parameters, character=-1, c_weights=None):
    parameters, weights = fixed._weights(parameters, c_weights)
    _character(character)
    result = []
    for seam, seam_weight, c_weight in zip(local.SEAMS, parameters.seam_weights, weights, strict=True):
        chance = parameters.repair_probability*F(seam_weight, sum(parameters.seam_weights))
        for law, choice in (('C', c_weight), ('E', 1-c_weight)):
            for row in _rows(law, seam, character):
                probability = chance*choice*local.coin_probability(parameters, seam, row.coin)
                result.append(WeightedBranch(probability, tuple(zip(row.inputs, row.outputs, row.signs)),
                                             (law, seam, row.coin, row.raw)))
    result.append(WeightedBranch(1-parameters.repair_probability,
                                 tuple((i, i, 1) for i in range(66)), ('hold',)))
    return tuple(result)


def faithful_channel_certificate(parameters, character=-1, c_weights=None):
    parameters, weights = fixed._weights(parameters, c_weights)
    pi = local.reference(parameters)
    effects, output = [F()]*66, [F()]*66
    raw = branches(parameters, character, weights)
    for row in raw:
        for i, j, sign in row.columns:
            effects[i] += row.probability*sign**2
            output[j] += row.probability*sign**2*pi[i]
    if effects != [F(1)]*66 or tuple(output) != pi or min(pi) <= 0:
        raise ArithmeticError('selected exchange lift lost completeness or the common faithful reference')
    return {'character': character, 'raw_slots_including_zero_and_hold': len(raw),
            'complete': True, 'same_faithful_diagonal_reference': True,
            'smallest_reference_weight': str(min(pi)), 'largest_reference_weight': str(max(pi))}


def algebra_certificate(rows):
    """Exact effect separation and signed transition forest, without ranks.

    Rows are actual supported weighted branches from ``branches``. This
    helper also permits removing seams as a controlled hypothesis failure.
    No absolute-value substitution of the quantum channel is used.
    """
    if type(rows) not in (tuple, list) or len(rows) > 721:
        raise ValueError('a bounded native branch inventory required')
    for row in rows:
        if (type(row) is not WeightedBranch or type(row.probability) is not F
                or type(row.probability.numerator) is not int or type(row.probability.denominator) is not int
                or row.probability.denominator <= 0
                or not 0 <= row.probability.numerator <= row.probability.denominator
                or type(row.columns) is not tuple or len(row.columns) > 66):
            raise ValueError('a valid exact weighted partial injection required')
        if (any(type(column) is not tuple or len(column) != 3
                or any(type(x) is not int for x in column)
                or not 0 <= column[0] < 66 or not 0 <= column[1] < 66
                or column[2] not in (-1, 1) for column in row.columns)
                or len({i for i, _, _ in row.columns}) != len(row.columns)
                or len({j for _, j, _ in row.columns}) != len(row.columns)):
            raise ValueError('unit-modulus injective signed columns required')
    rows = tuple(row for row in rows if row.probability and row.columns)
    if not rows:
        raise ValueError('supported branches required')
    supports = [frozenset(i for i, _, _ in row.columns) for row in rows]
    unresolved = [(i, j) for i in range(66) for j in range(i+1, 66)]
    selected = []
    while unresolved:
        scores = [sum((i in support) != (j in support) for i, j in unresolved) for support in supports]
        choice = max(range(len(rows)), key=lambda i: scores[i])
        if scores[choice] == 0:
            raise ValueError('supported effects do not separate all native source addresses')
        selected.append(choice)
        unresolved = [(i, j) for i, j in unresolved if (i in supports[choice]) == (j in supports[choice])]
    signatures = tuple(sum((i in supports[j]) << k for k, j in enumerate(selected)) for i in range(66))
    if len(set(signatures)) != 66:
        raise ArithmeticError('effect signatures fail to construct every primitive projector')
    neighbors = {i: {} for i in range(66)}
    for branch, row in enumerate(rows):
        for i, j, phase in row.columns:
            if i != j:
                # Multiplication by this NONZERO sign recovers a matrix unit;
                # it does not assert equivalence of the complete instruments.
                neighbors[i].setdefault(j, (branch, i, j, phase))
                neighbors[j].setdefault(i, (branch, i, j, phase))
    remaining, components, forest = set(range(66)), [], []
    while remaining:
        root = min(remaining)
        reached, queue = {root}, [root]
        for parent in queue:
            for child, witness in sorted(neighbors[parent].items()):
                if child not in reached:
                    reached.add(child)
                    queue.append(child)
                    branch, source, target, phase = witness
                    forest.append({'parent': parent, 'child': child, 'branch': rows[branch].label,
                                   'source': source, 'target': target, 'nonzero_phase': phase,
                                   'use_adjoint': parent != source})
        components.append(sorted(reached))
        remaining -= reached
    return {'dimension': 66, 'separating_effects': [rows[i].label for i in selected],
            'effect_signatures': signatures, 'components': components, 'transition_forest': forest,
            'generated_star_algebra_dimension': sum(len(component)**2 for component in components),
            'commutant_dimension': len(components), 'scalar_record_gauge_equivalence_claimed': False}


def _entries(rows, weight):
    result = {}
    for row in rows:
        probability = weight(row.coin)
        for i, out_i, si in zip(row.inputs, row.outputs, row.signs, strict=True):
            for j, out_j, sj in zip(row.inputs, row.outputs, row.signs, strict=True):
                key = (66*i+j, 66*out_i+out_j)
                result[key] = result.get(key, F())+probability*si*sj
    return {key: value for key, value in result.items() if value}


def matrix_unit_signs():
    """Signs of a GNS isometry; this is NOT a positive quantum operation."""
    return tuple(tuple((-1)**sum(d < q for q in set(left) & set(right)
                                 for d in set(left) ^ set(right))
                       for right in statistics.PAIRS) for left in statistics.PAIRS)


def matrix_unit_transform(observable):
    """Apply the mathematical intertwiner S, without treating it as a channel."""
    matrix = np.asarray(observable, dtype=complex)
    if matrix.shape != (66, 66) or not np.isfinite(matrix).all():
        raise ValueError('one finite complex H66 observable required')
    return np.asarray(matrix_unit_signs())*matrix


def similarity_certificate():
    """Exact branchwise similarity, with an explicit failure of positivity."""
    phases = matrix_unit_signs()
    if any(phases[i][i] != 1 or phases[i][j] != phases[j][i]
           or phases[i][j] not in (-1, 1) for i in range(66) for j in range(66)):
        raise ArithmeticError('matrix-unit transform lost its unital adjoint-preserving isometry')
    checks, raw_slots = 0, 0
    for seam in local.SEAMS:
        for law in ('C', 'E'):
            for plus, minus in zip(_rows(law, seam, 1), _rows(law, seam, -1), strict=True):
                if (plus.inputs, plus.outputs, plus.abort) != (minus.inputs, minus.outputs, minus.abort):
                    raise ArithmeticError('the two character branches have different native supports')
                for i, out_i, si in zip(minus.inputs, minus.outputs, minus.signs, strict=True):
                    for j, out_j, sj in zip(minus.inputs, minus.outputs, minus.signs, strict=True):
                        if si*sj != phases[i][j]*phases[out_i][out_j]:
                            raise ArithmeticError('signed branch does not obey the matrix-unit similarity')
                        checks += 1
                raw_slots += 1
    # The sign triangle has a negative quadratic form. Thus S cannot be
    # silently used as a positive preparation, physical gauge, or readout.
    addresses = [statistics.INDEX[pair] for pair in ((0, 1), (0, 2), (1, 2))]
    witness = (1, -1, 1)
    expectation = sum(F(witness[i]*witness[j]*phases[a][b], 9)
                      for i, a in enumerate(addresses) for j, b in enumerate(addresses))
    if expectation != F(-1, 3):
        raise ArithmeticError('nonpositive-intertwiner counterexample changed')
    return {'complete_raw_slots': raw_slots, 'branch_matrix_units_checked': checks,
            'identity': 'Phi_minus = S Phi_plus S, S^2=I, for every labeled branch',
            'same_spectrum_for_every_fixed_selector': True,
            'GNS_isometry_for_every_faithful_diagonal_reference': True,
            'all_diagonal_source_observables_fixed_by_S': True,
            'S_is_positive_or_a_physical_channel': False,
            'positive_input_control_pairs': [[0, 1], [0, 2], [1, 2]],
            'positive_input': 'equal coherent superposition of the three addresses',
            'output_witness_vector': list(witness),
            'negative_normalized_output_expectation': str(expectation)}


def seam_GNS_certificate(parameters, character=-1):
    """Check the actual signed blocks and both weighted positivity minors."""
    parameters = fixed._parameters(parameters)
    pi = local.reference(parameters)
    one, two, negative_entries = 0, 0, 0
    for seam in local.SEAMS:
        for law in ('C', 'E'):
            entries = _entries(_rows(law, seam, character),
                               lambda coin: local.coin_probability(parameters, seam, coin))
            negative_entries += sum(value < 0 for value in entries.values())
            neighbors = {}
            for i, j in entries:
                if i != j:
                    neighbors.setdefault(i, set()).add(j)
                    neighbors.setdefault(j, set()).add(i)
            seen = set()
            for i in range(66**2):
                if i in seen:
                    continue
                component = {i} | neighbors.get(i, set())
                if len(component) > 2 or any(not neighbors.get(j, set()) <= component for j in component):
                    raise ArithmeticError('signed seam map has a larger than two-unit block')
                seen.update(component)
                indices = sorted(component)
                gram = [[pi[a % 66]*entries.get((a, b), F()) for b in indices] for a in indices]
                complement = [[(pi[a % 66] if a == b else F())-gram[u][v]
                               for v, b in enumerate(indices)] for u, a in enumerate(indices)]
                for matrix in (gram, complement):
                    if (any(matrix[u][v] != matrix[v][u] for u in range(len(indices)) for v in range(len(indices)))
                            or any(matrix[u][u] < 0 for u in range(len(indices)))
                            or (len(indices) == 2 and matrix[0][0]*matrix[1][1] < matrix[0][1]**2)):
                        raise ArithmeticError('signed channel failed a weighted self-adjoint contraction minor')
                one += len(indices) == 1
                two += len(indices) == 2
    return {'character': character, 'matrix_units_checked': one+2*two,
            'seam_law_maps': 60, 'negative_matrix_entries_retained': negative_entries,
            'weighted_self_adjoint_positive_contractions': True}


def extra_dephasing_certificate(character=-1):
    checks, differences = 0, 0
    for seam in local.SEAMS:
        rows = {law: _rows(law, seam, character) for law in ('C', 'E')}
        for coin in (0, 1):
            entries = {law: _entries(tuple(r for r in rows[law] if r.coin == coin), lambda _: F(1))
                       for law in ('C', 'E')}
            if any(entries['C'].get(key) != value for key, value in entries['E'].items()):
                raise ArithmeticError('signed E does not refine signed C')
            for key in entries['C'].keys()-entries['E'].keys():
                if key[0] != key[1] or entries['C'][key] != 1:
                    raise ArithmeticError('signed C-E is not diagonal nonnegative')
                differences += 1
            checks += 66**2
    return {'character': character, 'matrix_unit_coin_seam_checks': checks,
            'additional_C_matrix_unit_terms': differences,
            'C_minus_E_diagonal_nonnegative_quadratic_form': True,
            'C_minus_E_positive_map_claimed': False}


def uniform_gap_certificate(character=-1):
    """Prove the signed coherence bound by verified entrywise domination.

    The signed matrix itself remains in the certificate. Absolute values
    give a norm UPPER bound only after symmetry, population equality and
    the positive comparison row sums have been checked independently.
    """
    _character(character)
    actual, positive = {}, {}
    for seam in local.SEAMS:
        for target, chi in ((actual, character), (positive, 1)):
            for key, value in _entries(_rows('C', seam, chi), lambda _: F(1, 60)).items():
                target[key] = target.get(key, F())+value
    if any(actual.get((j, i), F()) != value for (i, j), value in actual.items()):
        raise ArithmeticError('signed fair channel is not symmetric')
    row_bounds = [F() for _ in range(66**2)]
    for (i, j), value in positive.items():
        if value < 0 or abs(actual.get((i, j), F())) > value:
            raise ArithmeticError('signed fair channel exceeds its positive comparison entry')
        row_bounds[i] += value
    if any(abs(value) > positive.get(key, F()) for key, value in actual.items()):
        raise ArithmeticError('signed channel created an unbounded matrix-unit transition')
    population = mixing.probability_kernel(local.Parameters(c_weight=F(1)))
    for i in range(66):
        for j in range(66):
            if actual.get((67*i, 67*j), F()) != population[i][j]:
                raise ArithmeticError('signed channel changed the fair population kernel')
    for (i, j) in actual:
        if (i//66 == i % 66) != (j//66 == j % 66):
            raise ArithmeticError('population and coherence spaces ceased to be invariant')
    ceiling = max(row_bounds[i] for i in range(66**2) if i//66 != i % 66)
    if ceiling > F(5, 6):
        raise ArithmeticError('positive coherence row bound exceeds5/6')
    base = fixed.uniform_gap_certificate()
    return {**base, 'character': character, 'absolute_coherence_row_sum_upper': str(ceiling),
            'actual_signed_matrix_is_symmetric': True, 'population_kernel_unchanged': True,
            'entrywise_absolute_domination_verified': True,
            'signed_channel_replaced_by_absolute_values': False}


def gap_bound(parameters, character=-1, c_weights=None):
    _character(character)
    # The checked signed fair gap and signed C-E comparison yield the SAME
    # scalar coefficient comparison as the existing positive-lift proof.
    return fixed.gap_bound(parameters, c_weights)


def apply_heisenberg(parameters, observable, character=-1, c_weights=None):
    matrix = np.asarray(observable, dtype=complex)
    if matrix.shape != (66, 66) or not np.isfinite(matrix).all():
        raise ValueError('one finite complex H66 observable required')
    result = np.zeros_like(matrix)
    for row in branches(parameters, character, c_weights):
        if not row.probability or not row.columns:
            continue
        inputs, outputs, signs = zip(*row.columns, strict=True)
        phases = np.outer(signs, signs)
        result[np.ix_(inputs, inputs)] += float(row.probability)*phases*matrix[np.ix_(outputs, outputs)]
    return result


def dirichlet_control(parameters, observable, character=-1, c_weights=None):
    parameters, weights = fixed._weights(parameters, c_weights)
    matrix = np.asarray(observable, dtype=complex)
    evolved = apply_heisenberg(parameters, matrix, character, weights)
    pi = np.array(list(map(float, local.reference(parameters))))
    left = float(np.sum(pi*np.sum(matrix.conj()*(matrix-evolved), axis=0)).real)
    right = 0.0
    for row in branches(parameters, character, weights):
        if not row.probability or not row.columns:
            continue
        inputs, outputs, signs = zip(*row.columns, strict=True)
        commutator = np.zeros_like(matrix)
        commutator[:, inputs] += matrix[:, outputs]*np.array(signs)[None, :]
        commutator[outputs, :] -= np.array(signs)[:, None]*matrix[inputs, :]
        right += float(row.probability)*float(np.sum(pi*np.sum(abs(commutator)**2, axis=0)))/2
    centered = matrix-np.dot(pi, np.diag(matrix))*np.eye(66)
    return {'dirichlet_trace_pairing': left, 'half_weighted_commutator_sum': right,
            'weighted_variance': float(np.sum(pi*np.sum(abs(centered)**2, axis=0))),
            'weighted_fixed_point_residual_squared': float(np.sum(pi*np.sum(abs(matrix-evolved)**2, axis=0)))}


def certificate():
    parameters = mixing.canonical_parameters()
    characters = {}
    for character in (1, -1):
        algebra = algebra_certificate(branches(parameters, character))
        if algebra['commutant_dimension'] != 1 or algebra['generated_star_algebra_dimension'] != 66**2:
            raise ArithmeticError('the selected lift did not generate the full source algebra')
        characters[str(character)] = {
            'faithful_channel': faithful_channel_certificate(parameters, character),
            'source_algebra': algebra, 'seam_GNS': seam_GNS_certificate(parameters, character),
            'extra_E_dephasing': extra_dephasing_certificate(character),
            'fair_gap': uniform_gap_certificate(character),
            'canonical_gap_lower': str(gap_bound(parameters, character)),
            'whole_bounded_family_gap_lower': str(F(fixed.uniform_gap_certificate()['uniform_gap_lower'])/10240)}
    return {'schema': 'oph.signed_native_fixed_algebra.v1', 'characters': characters,
            'exact_nonpositive_matrix_unit_similarity': similarity_certificate(),
            'fixed_source_observables': 'complex scalars on H66 for either lift',
            'arbitrary_untouched_reference_fixed_algebra': 'I_source tensor B(reference)',
            'classical_feedback_conservation_requires_all_conditional_inputs': True,
            'exogenous_character_and_selector_switching_uses_same_reference': True,
            'joint_feedback_or_growing_archive_mixing_claimed': False,
            'arbitrary_coherent_recovery_excluded': False,
            'time_parameter': 'complete selected-character channel applications',
            'negative_character_includes_declared_record_selected_phase_service': True,
            'equal_pulse_cost_or_physical_duration_claimed': False,
            'physical_spin_species_mass_or_CP03_acceptance_claimed': False}


if __name__ == '__main__':
    print(json.dumps(certificate(), indent=2))
