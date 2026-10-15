"""Recover relative agreement energy from the actual native mean response.

This attachment concerns the existing uniform-seam calibration, with fair
placement coins. It derives a normalized diagonal symmetrizer and seam
coefficients from that calibration's response. It does not identify the
native load with a scalar quantum field or its event count with physical
time. Concentrated preparation and access to load means are audit inputs;
the production observer's public interface contains committed currents.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F

from ..observer import native
from ..common_source import slow_geometry


SOURCE_PATHS = (
    'model/observer/native.py', 'model/common_source/slow_geometry.py',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    'Lean/Geometry/SourceActionMeasure.lean',
    'Lean/SourceLaw/WeightedSeams.lean',
    'Lean/SourceLaw/RepairAdmissibility.lean',
    'Lean/SourceLaw/RepairConvergence.lean',
    'Lean/SourceLaw/VacuumIdentity.lean',
)


def _zero(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def _mv(matrix, vector):
    return tuple(sum((a*b for a, b in zip(row, vector, strict=True)), F(0))
                 for row in matrix)


def _mm(left, right):
    return tuple(tuple(sum((left[i][k]*right[k][j] for k in range(len(right))), F(0))
                       for j in range(len(right[0]))) for i in range(len(left)))


def _mean(state, seams):
    result = [F(0)]*len(state)
    weight = F(1, 2*len(seams))
    for seam in seams:
        for coin in (0, 1):
            output, _ = native.native(state, seam, coin)
            for i, value in enumerate(output):
                result[i] += weight*value
    return tuple(result)


def recover_response():
    """Infer all columns from allowed total-two inputs, then check all inputs."""
    seams = tuple(slow_geometry.source().seams())
    if len(seams) != 30 or len(set(seams)) != 30:
        raise ArithmeticError('the declared carrier seam dictionary changed')
    states = native.states()
    if len(states) != 78 or any(len(state) != 12 for state in states):
        raise ArithmeticError('the declared native total-two carrier changed')
    columns = tuple(tuple(value/2 for value in _mean(
        tuple(2*int(i == j) for i in range(12)), seams)) for j in range(12))
    transition = tuple(tuple(columns[j][i] for j in range(12)) for i in range(12))
    for state in states:
        if _mean(state, seams) != _mv(transition, state):
            raise ArithmeticError('the actual native mean response is not linear on this carrier')
    response = tuple(tuple(F(i == j)-transition[i][j] for j in range(12))
                     for i in range(12))
    return transition, response, seams


def relative_action(response):
    """Solve positive diagonal reciprocity, independently of native scheduling.

    Connected reciprocal support determines every capacity ratio. The final
    normalization is a convention; common rescaling preserves the response.
    This is the concrete hypothesis map of normalized_symmetrizer_unique.
    """
    if (not isinstance(response, (tuple, list)) or not response
            or any(not isinstance(row, (tuple, list)) or len(row) != len(response)
                   for row in response)
            or any(type(value) not in (int, F) for row in response for value in row)):
        raise ValueError('a nonempty square exact rational response is required')
    matrix = tuple(tuple(map(F, row)) for row in response)
    n = len(matrix)
    if (any(sum(row) != 0 for row in matrix)
            or any(matrix[i][j] > 0 for i in range(n) for j in range(n) if i != j)
            or any((matrix[i][j] == 0) != (matrix[j][i] == 0)
                   for i in range(n) for j in range(n))):
        raise ValueError('a conservative response with reciprocal negative support is required')
    capacity = {0: F(1)}
    pending, tree = [0], []
    while pending:
        i = pending.pop(0)
        for j in range(n):
            if j != i and matrix[i][j] != 0:
                value = capacity[i]*matrix[i][j]/matrix[j][i]
                if j in capacity:
                    if capacity[j] != value:
                        raise ValueError('response fails positive diagonal reciprocity')
                else:
                    capacity[j] = value
                    pending.append(j)
                    tree.append((i, j))
    if len(capacity) != n:
        raise ValueError('response support must be connected for normalized uniqueness')
    total = sum(capacity.values())
    masses = tuple(capacity[i]/total for i in range(n))
    stiffness = tuple(tuple(masses[i]*matrix[i][j] for j in range(n)) for i in range(n))
    if (any(value <= 0 for value in masses)
            or any(stiffness[i][j] != stiffness[j][i] for i in range(n) for j in range(n))):
        raise ArithmeticError('the recovered diagonal symmetrizer is invalid')
    seams = tuple((i, j, -stiffness[i][j]) for i in range(n) for j in range(i+1, n)
                  if stiffness[i][j] != 0)
    independent = _zero(n)
    for i, j, conductance in seams:
        independent[i][i] += conductance
        independent[j][j] += conductance
        independent[i][j] -= conductance
        independent[j][i] -= conductance
    if tuple(map(tuple, independent)) != stiffness:
        raise ArithmeticError('seam energy does not reproduce the recovered response')
    return {'capacities': masses, 'stiffness': stiffness, 'seams': seams,
            'spanning_tree': tuple(tree)}


def energy(readings, attachment):
    """Half the recovered seam quadratic form; exact diagnostic readings only."""
    values = tuple(readings)
    if (len(values) != len(attachment['capacities'])
            or any(type(value) not in (int, F) for value in values)):
        raise ValueError('one exact rational reading per recovered node is required')
    return sum((c*(values[i]-values[j])**2/2 for i, j, c in attachment['seams']), F(0))


def _weighted_repair_check(response, attachment):
    masses = attachment['capacities']
    out = _zero(len(masses))
    rates = []
    for i, j, c in attachment['seams']:
        rate = c*(1/masses[i]+1/masses[j])
        rates.append(rate)
        # Derive the mean loss directly from capacity-weighted agreement.
        out[i][i] += rate*masses[j]/(masses[i]+masses[j])
        out[i][j] -= rate*masses[j]/(masses[i]+masses[j])
        out[j][j] += rate*masses[i]/(masses[i]+masses[j])
        out[j][i] -= rate*masses[i]/(masses[i]+masses[j])
    if tuple(map(tuple, out)) != response:
        raise ArithmeticError('weighted-seam theorem does not attach to the actual response')
    return rates


def _diagonal_quantum_extension(seams):
    """Prove diagonal output laws ignore input coherences for both instruments.

    Within each record branch the map from input configuration to output is
    injective. Therefore L_r^* D L_r is diagonal for every diagonal D. The
    verified classical-basis response then extends by linearity to every
    quantum input, including reductions of arbitrary reference correlations.
    This does not identify the complete C/E instruments, which differ.
    """
    states = native.states()
    checked = 0
    for law in ('C', 'E'):
        for seam in seams:
            seen = set()
            for coin in (0, 1):
                for state in states:
                    output, record = native.native(state, seam, coin)
                    label = record if law == 'C' else (coin, state[seam[0]]-state[seam[1]])
                    address = (label, output)
                    if address in seen:
                        raise ArithmeticError('native recorded branch is not injective on input configurations')
                    seen.add(address)
                    checked += 1
    return {'recorded_branch_addresses_checked': checked,
            'laws': ['C', 'E'], 'configuration_carrier_dimension': len(states),
            'every_record_branch_injective': True,
            'diagonal_output_observable_pullbacks_are_diagonal': True,
            'mean_response_extends_to_all_quantum_and_reference_inputs': True,
            'complete_quantum_instruments_identified_with_each_other': False}


def _noise_control():
    state = (1, 0, 1)+(0,)*9
    currents = tuple(native.native(state, (0, 1), coin)[1][1] for coin in (0, 1))
    mean = F(sum(currents), 2)
    variance = sum((F(value)-mean)**2 for value in currents)/2
    fourth = sum((F(value)-mean)**4 for value in currents)/2
    gaussian = 3*variance**2
    if (sorted(currents) != [0, 1] or variance != F(1, 4)
            or fourth != F(1, 16) or fourth == gaussian):
        raise ArithmeticError('the actual atomic-noise counterexample changed')
    return {'input_loads': list(state), 'selected_seam': [0, 1],
            'coin_order_currents': list(currents), 'mean_current': str(mean),
            'variance': str(variance), 'centered_fourth_moment': str(fourth),
            'Gaussian_fourth_moment_at_same_variance': str(gaussian),
            'native_conditional_noise_is_Gaussian': False}


def _coherent_control(seams):
    """Actual C/E native sums, averaged over the same fair uniform scheduler."""
    states = native.states(binary=True)
    a = tuple(int(i in (2, 3)) for i in range(12))
    u = tuple(int(i in (0, 2)) for i in range(12))
    if a not in states or u not in states:
        raise ArithmeticError('coherent witness left the binary carrier')
    probabilities = {}
    for law in ('C', 'E'):
        output = defaultdict(F)
        for seam in seams:
            for coin in (0, 1):
                for x in (a, u):
                    for y in (a, u):
                        ox, rx = native.native(x, seam, coin)
                        oy, ry = native.native(y, seam, coin)
                        same = (rx == ry if law == 'C' else
                                x[seam[0]]-x[seam[1]] == y[seam[0]]-y[seam[1]])
                        if same:
                            output[ox, oy] += F(1, 4*len(seams))
        if sum(value for (x, y), value in output.items() if x == y) != 1:
            raise ArithmeticError('the complete C/E output is not normalized')
        plus = sum(output[x, y]/2 for x in (a, u) for y in (a, u))
        probabilities[law] = [str(plus), str(1-plus)]
    if probabilities['C'] == probabilities['E']:
        raise ArithmeticError('same mean response no longer separates complete quantum laws')
    return {'preparation': '(|a>+|u>)/sqrt(2)', 'a_loads': list(a), 'u_loads': list(u),
            'seam_scheduler': 'uniform over the same thirty seams',
            'read': 'projector onto the prepared ket and its full complement',
            'complete_binary_probabilities': probabilities,
            'same_classical_population_transition': True,
            'same_relative_mean_action_selects_quantum_instrument': False,
            'coherent_preparation_and_read_supplied': True}


def _stationary_vacuum_control(attachment, seams):
    """Compare actual stationary load covariance with a proposed vacuum reading.

    The tested interpretation identifies centered loads with canonical field
    coordinates and the reconstructed capacities with inertia. It is an
    additional candidate identification, not an assumption of the runtime.
    Even after removing the conserved constant direction, no common hbar
    makes this covariance the ground-state covariance of that quadratic law.
    """
    states = native.states(binary=True)
    count = len(states)
    outgoing = defaultdict(int)
    for state in states:
        for seam in seams:
            for coin in (0, 1):
                outgoing[native.native(state, seam, coin)[0]] += 1
    if set(outgoing) != set(states) or any(outgoing[state] != 2*len(seams) for state in states):
        raise ArithmeticError('the actual binary uniform reference is not stationary')
    means = tuple(F(sum(state[i] for state in states), count) for i in range(12))
    covariance = tuple(tuple(F(sum(state[i]*state[j] for state in states), count)-means[i]*means[j]
                             for j in range(12)) for i in range(12))
    projector = tuple(tuple(F(i == j)-F(1, 12) for j in range(12)) for i in range(12))
    if covariance != tuple(tuple(F(5, 33)*value for value in row) for row in projector):
        raise ArithmeticError('the exact stationary native covariance changed')
    stiffness = attachment['stiffness']
    product = _mm(_mm(covariance, stiffness), covariance)
    # For the hypothetical same-capacity quadratic Hamiltonian, the Gaussian
    # ground relation is C V C = hbar**2 T/4 on the mean-free sector, with
    # inverse inertia T = 12*Pi. This algebraic criterion avoids square roots.
    inverse_inertia = tuple(tuple(12*value for value in row) for row in projector)
    adjacent = next((i, j) for i in range(12) for j in range(i+1, 12) if stiffness[i][j] < 0)
    nonadjacent = next((i, j) for i in range(12) for j in range(i+1, 12) if stiffness[i][j] == 0)
    required = tuple(4*product[i][j]/inverse_inertia[i][j] for i, j in (adjacent, nonadjacent))
    if required[0] <= 0 or required[1] != 0:
        raise ArithmeticError('the exact incompatible-vacuum entries changed')
    return {'stationary_reference': 'uniform law on all 66 binary total-two states',
            'stationarity_branches_checked': count*2*len(seams),
            'mean_load_per_port': '1/6',
            'covariance': '(5/33)*(I-J_all_ones/12)',
            'covariance_diagonal': str(covariance[0][0]),
            'covariance_off_diagonal': str(covariance[0][1]),
            'conserved_constant_mode_removed': True,
            'hypothetical_field_identification': 'centered native loads are canonical q; inertia is recovered M',
            'ground_relation_tested': 'C*K_action*C = (hbar^2/4)*M_inverse on the mean-free sector',
            'adjacent_entry': list(adjacent), 'nonadjacent_entry': list(nonadjacent),
            'hbar_squared_required_by_adjacent_entry': str(required[0]),
            'hbar_squared_required_by_nonadjacent_entry': str(required[1]),
            'one_common_hbar_matches_ground_relation': False,
            'excludes_other_observer_fields_or_states': False}


def _spectrum_check(response):
    """Factor the exact recovered response, without input spectral projectors."""
    from sympy import Matrix, Poly, factor_list
    integer_matrix = [[60*value for value in row] for row in response]
    if any(value.denominator != 1 for row in integer_matrix for value in row):
        raise ArithmeticError('the recovered response is not the expected rational operator')
    polynomial = Matrix([[int(value) for value in row] for row in integer_matrix]).charpoly()
    factor, terms = factor_list(polynomial.as_expr(), polynomial.gen)
    recovered = sorted((tuple(int(value) for value in Poly(term, polynomial.gen).all_coeffs()), int(power))
                       for term, power in terms)
    if factor != 1 or recovered != [((1, -10, 20), 3), ((1, -6), 5), ((1, 0), 1)]:
        raise ArithmeticError('the exact spectrum of the source-recovered action changed')
    return {'matrix': '60*R, recovered from actual native outputs',
            'characteristic_polynomial_coefficients': [int(value) for value in polynomial.all_coeffs()],
            'factorization': 'lambda*(lambda-6)^5*(lambda^2-10*lambda+20)^3',
            'smallest_nonzero_response_eigenvalue': '(5-sqrt(5))/60',
            'slow_response_dimension': 3,
            'physical_spatial_dimension_selected': False}


def run():
    transition, response, seams = recover_response()
    attachment = relative_action(response)
    masses = attachment['capacities']
    rates = _weighted_repair_check(response, attachment)
    if (masses != (F(1, 12),)*12 or set(rates) != {F(1, 30)}
            or {c for _, _, c in attachment['seams']} != {F(1, 720)}
            or {tuple(sorted(seam)) for seam in seams}
                != {(i, j) for i, j, _ in attachment['seams']}):
        raise ArithmeticError('the recovered native relative action changed')
    # The independently recovered stiffness has the same exact null direction
    # as conservation. It is positive semidefinite, not positively grounded.
    if _mv(attachment['stiffness'], (1,)*12) != (F(0),)*12:
        raise ArithmeticError('the constant native zero mode changed')
    witness = (1,)+(0,)*11
    sample_energy = energy(witness, attachment)
    quadratic = sum(F(witness[i])*value for i, value in
                    enumerate(_mv(attachment['stiffness'], witness)))/2
    if sample_energy != quadratic or sample_energy != F(1, 288):
        raise ArithmeticError('the recovered seam quadratic form changed')
    return {'schema': 'oph.primitive_observer.native_relative_action.v1', 'all_pass': True,
            'response_recovery': {'independent_total_two_inputs': 12,
                                  'full_native_inputs_checked': 78,
                                  'binary_native_inputs_included': 66,
                                  'native_branches_checked': 4680,
                                  'response_rows': [[str(x) for x in row] for row in response]},
            'relative_action': {'normalized_capacities': list(map(str, masses)),
                                'connected_spanning_tree': [list(edge) for edge in attachment['spanning_tree']],
                                'seam_conductances': [[i, j, str(c)] for i, j, c in attachment['seams']],
                                'recovered_seam_rates': list(map(str, rates)),
                                'stiffness': 'L/720', 'energy': 'sum_seams (f_a-f_b)^2/1440',
                                'unit_port_reading_energy': str(sample_energy),
                                'normalized_capacity_unique': True,
                                'absolute_action_scale_selected': False,
                                'constant_zero_mode': True},
            'theorem_attachment': {
                'relative_capacity': 'OPH.SourceActionMeasure.normalized_symmetrizer_unique',
                'mean_drift': 'OPH.SourceLaw.repairGenerator_apply',
                'energy': 'OPH.SourceLaw.stiffness_quadratic',
                'exact_hypotheses_checked': ['connected reciprocal support', 'positive capacities',
                                           'unit capacity sum', 'all reciprocity equations',
                                           'distinct seam endpoints', 'nonnegative conductances',
                                           'actual all-input native mean equals weighted-seam generator'],
                'new_Lean_proof_generated': False},
            'noise_control': _noise_control(), 'quantum_law_control': _coherent_control(seams),
            'diagonal_quantum_extension': _diagonal_quantum_extension(seams),
            'stationary_vacuum_control': _stationary_vacuum_control(attachment, seams),
            'recovered_response_spectrum': _spectrum_check(response),
            'history_register_boundary': {
                'same_spatial_mean_drift_attached': True,
                'actual_native_precision_on_time_slabs_supplied': False,
                'all_small_slab_admissibility_verified': False,
                'uniform_time_seam_rate_2_over_dt_supplied_by_native': False,
                'positive_grounding_supplied_by_native': False,
                'native_noise_is_the_declared_Gaussian_Langevin_noise': False,
                'vacuum_covariance_theorem_instantiated': False,
                'inertia_countercontrol': 'J=2M is diagonal and repair admissible; its time-seam rate is 4/dt',
                'clock_premise_needed_to_exclude_countercontrol': True},
            'scope': {'carrier_architecture_and_uniform_scheduler_supplied': True,
                      'placement_coin_fairness_supplied': True,
                      'load_means_are_working_state_audits': True,
                      'production_feedback_observer_has_this_uniform_schedule': False,
                      'parameter': 'native attempt count in the declared uniform calibration',
                      'relative_mean_action_derived_from_actual_native_response': True,
                      'physical_field_action_or_interaction_coefficients_selected': False,
                      'quantum_inertia_vacuum_or_physical_clock_selected': False,
                      'empirical_outcomes_received_or_fitted': False}}
