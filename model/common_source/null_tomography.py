"""Exact native-port decoder for a declared Lorentzian tensor-charge map.

The directions come from the existing slow-response carrier; Lorentzian
lifting and physical stress/geometry measurements remain supplied interfaces.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import Counter

from ..observer import native
from .slow_geometry import (_add, _dot, _encode_vector, raw_frame,
                            registered_binding, source)


# Symmetric four-tensor coordinates; each off-diagonal entry occurs twice
# when contracted. The declared raw-coordinate metric is diag(-r^2,1,1,1),
# r^2=phi+2. Its normalized chart is the usual diag(-1,1,1,1).
TENSOR_COMPONENTS = ((0, 0), (0, 1), (0, 2), (0, 3), (1, 1),
                     (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))


def null_tomography_design():
    """Exact charges B(k_p,k_p) for the CONDITIONAL lift k_p=(1,raw_p).

    The carrier supplies raw_p, not a Lorentzian physical interpretation or
    a stress measurement. This map tests the finite directional requirement
    of EinsteinBranch.Tensor on that declared lift.
    """
    c = source()
    return [[c.q5_mul(c.q5(1 if i == j else 2), c.q5_mul(k[i], k[j]))
             for i, j in TENSOR_COMPONENTS]
            for axis in raw_frame() for k in [(c.q5(1),)+axis]]


def _null_tomography_decoder():
    """Closed inverse on the B00=0 slice, from exact icosahedral moments."""
    c = source()
    inverse_r2 = c.q5_inv(c.q5_add(c.GOLDEN, c.q5(2)))
    inverse_r4 = c.q5_mul(inverse_r2, inverse_r2)
    rows = []
    for i, j in TENSOR_COMPONENTS[1:]:
        if i == 0:
            rows.append([c.q5_mul(c.q5(F(1, 8)), c.q5_mul(inverse_r2, u[j-1]))
                         for u in raw_frame()])
        else:
            rows.append([c.q5_sub(
                c.q5_mul(c.q5(F(5, 8)), c.q5_mul(inverse_r4, c.q5_mul(u[i-1], u[j-1]))),
                c.q5_mul(c.q5(F(int(i == j), 8)), inverse_r2)) for u in raw_frame()])
    return rows


def decode_null_charges(charges):
    """Return the unique B00=0 representative, rejecting incompatible data.

    Twelve arbitrary exact readings need not be charges of one symmetric
    tensor. No projection or least-squares fit replaces this compatibility
    test. Entries are exact (a,b) pairs representing a+b*sqrt(5).
    """
    if (not isinstance(charges, (tuple, list)) or len(charges) != 12
            or any(not isinstance(q, (tuple, list)) or len(q) != 2
                   or any(type(v) not in (int, F) for v in q) for q in charges)):
        raise ValueError('twelve exact Q(sqrt5) charge pairs are required')
    values = tuple((F(a), F(b)) for a, b in charges)
    candidate = (source().q5(0),)+tuple(_dot(row, values) for row in _null_tomography_decoder())
    if tuple(_dot(row, candidate) for row in null_tomography_design()) != values:
        raise ValueError('the twelve readings are incompatible with a common symmetric tensor')
    return candidate


def _max_row_sum(matrix):
    """The exact induced infinity norm, also valid on real input errors."""
    c = source()
    largest = c.q5(0)
    for row in matrix:
        total = c.q5(0)
        for value in row:
            total = c.q5_add(total, value if c.q5_sign(value) >= 0 else c.q5_sub(c.q5(0), value))
        if c.q5_sign(c.q5_sub(total, largest)) > 0:
            largest = total
    return largest


def null_tomography_certificate():
    """Finite exact inverse/kernel certificate on the native twelve ports.

    All matrix identities are checked over Q(sqrt5); hence the displayed
    inverse also proves the statement for arbitrary real tensor components.
    This is an executable algebraic certificate, not an additional Lean proof.
    """
    c = source()
    zero, one = c.q5(0), c.q5(1)
    frame = raw_frame()
    r2 = c.q5_add(c.GOLDEN, c.q5(2))
    if c.q5_sign(r2) <= 0:
        raise AssertionError('the declared metric requires a positive spatial radius')
    inverse_r2 = c.q5_inv(r2)
    gram = c.intrinsic_gram_exact()
    if any(c.q5_mul(_dot(u, v), inverse_r2) != gram[i][j]
           for i, u in enumerate(frame) for j, v in enumerate(frame)):
        raise AssertionError('null frame differs from the native response Gram')
    design, decoder = null_tomography_design(), _null_tomography_decoder()
    sliced = [row[1:] for row in design]
    identity9 = [[c.q5(int(i == j)) for j in range(9)] for i in range(9)]
    if c.q5_matmul(decoder, sliced) != identity9:
        raise AssertionError('native null-charge decoder is not an exact left inverse')
    metric = [c.q5_sub(zero, r2), zero, zero, zero, one, zero, zero, one, zero, one]
    if any(_dot(row, metric) != zero for row in design):
        raise AssertionError('declared carrier lifts are not null')
    # Rank >=9 by the left inverse, <=9 by this nonzero metric kernel.
    # The full inverse verifies that there are no other invisible components.
    quotient = [[(inverse_r2 if i in (3, 6, 8) else zero)]+identity9[i]
                for i in range(9)]
    if c.q5_matmul(decoder, design) != quotient:
        raise AssertionError('full tensor decoder does not quotient by precisely the metric')
    projector = c.q5_matmul(sliced, decoder)
    partner = [next(j for j, v in enumerate(frame) if _add(u, v) == (zero,)*3)
               for u in frame]
    even = [[c.q5(F(int(i == j)+int(partner[i] == j), 2)) for j in range(12)]
            for i in range(12)]
    slow = c.slow_band_projector_exact()
    if projector != [[c.q5_add(even[i][j], slow[i][j]) for j in range(12)] for i in range(12)]:
        raise AssertionError('tensor charges do not match the native constant/slow/quadrupole bands')
    residual = [[c.q5_sub(c.q5(int(i == j)), projector[i][j]) for j in range(12)]
                for i in range(12)]
    lap = [[c.q5(int(x)) for x in row] for row in c.laplacian()]
    if c.q5_matmul(lap, residual) != [[c.q5_mul(c.q5(5, 1), x) for x in row] for row in residual]:
        raise AssertionError('incompatible charges are not precisely the native fast band')
    merged = c.q5_matmul(even, design)
    # Six independent spatial components survive antipodal averaging; the
    # three mixed time/space entries vanish, in addition to the metric kernel.
    spatial_decoder = decoder[3:]
    if any(row[j] != zero for row in merged for j in (1, 2, 3)):
        raise AssertionError('antipodal averaging should erase time/space mixed charges')
    identity6 = [[c.q5(int(i == j)) for j in range(6)] for i in range(6)]
    if c.q5_matmul(spatial_decoder, [row[4:] for row in merged]) != identity6:
        raise AssertionError('antipodal merged charges lost a spatial tensor component')
    return {'schema': 'oph.common_source.native_null_tomography.v1',
            'registered_binding': registered_binding(),
            'tensor_component_order': [list(pair) for pair in TENSOR_COMPONENTS],
            'declared_raw_metric_diagonal_Qsqrt5': _encode_vector((metric[0], one, one, one)),
            'design_Qsqrt5': [_encode_vector(row) for row in design],
            'B00_zero_decoder_Qsqrt5': [_encode_vector(row) for row in decoder],
            'exact_design_rank': 9, 'exact_kernel_dimension': 1,
            'exact_charge_compatibility_constraints': 3,
            'charge_image': 'constant + native slow rank3 + antipodal-even tracefree rank5',
            'incompatible_reading_band': 'native fast rank3, Laplacian eigenvalue 5+sqrt(5)',
            'antipodal_merged_rank': 6, 'antipodal_merged_kernel_dimension': 4,
            'canonical_tensor_max_error_per_charge_max_error_Qsqrt5': c.q5_str(_max_row_sum(decoder)),
            'compatibility_residual_max_error_per_charge_max_error_Qsqrt5': c.q5_str(_max_row_sum(residual)),
            'checks': {'decoder_left_inverse_exact': True, 'full_kernel_is_metric_exact': True,
                       'charge_image_matches_native_bands_exact': True,
                       'antipodal_loss_exact': True},
            'scope': {'conditional_lorentz_lift': True,
                      'physical_null_directions_derived': False,
                      'directional_stress_or_geometry_measurements_constructed': False,
                      'einstein_equation_derived': False}, 'all_pass': True}


def local_occupancy(initial_bit, port, incident_records):
    """Read one bit using its disclosed initial value and own signed currents.

    Each record is (oriented native seam, signed current). No later source
    state is an input. Access to the initial bit and a complete local tape
    are the protocol contract; these are not inferred from a hidden load.
    """
    if type(initial_bit) is not int or initial_bit not in (0, 1):
        raise ValueError('a disclosed binary initial occupancy is required')
    if type(port) is not int or port not in range(12):
        raise ValueError('a native port is required')
    value = initial_bit
    for seam, current in incident_records:
        if (not isinstance(seam, (list, tuple)) or len(seam) != 2
                or any(type(i) is not int for i in seam)
                or tuple(sorted(seam)) not in source().seams() or port not in seam
                or type(current) is not int or current not in (-1, 0, 1)):
            raise ValueError('an incident signed native current is required')
        value += current*(int(port == seam[1])-int(port == seam[0]))
        if value not in (0, 1):
            raise ValueError('the local record is incompatible with binary occupancy')
    return value


def _binary_initial(initial):
    if (not isinstance(initial, (list, tuple)) or len(initial) != 12
            or any(type(x) is not int or x not in (0, 1) for x in initial)
            or sum(initial) != 2):
        raise ValueError('a binary total-two preparation is required')
    return tuple(initial)


def execute_record_word(initial, word):
    """Execute a declared seam/coin word and return only disclosed records."""
    initial = _binary_initial(initial)
    state, local, events = initial, [[] for _ in range(12)], []
    for seam, coin in word:
        if (not isinstance(seam, (tuple, list)) or len(seam) != 2
                or any(type(i) is not int for i in seam)
                or tuple(sorted(seam)) not in source().seams()
                or type(coin) is not int or coin not in (0, 1)):
            raise ValueError('a native seam and fair-coin outcome are required')
        state, (_, current) = native.native(state, seam, coin)
        for port in seam:
            local[port].append((tuple(seam), current))
        events.append((tuple(seam), coin, current))
        # Hidden states are used only to audit the implemented readout.
        if tuple(local_occupancy(initial[i], i, local[i]) for i in range(12)) != state:
            raise AssertionError('local records disagree with the executed native state')
    return {'initial_bits': initial, 'events': tuple(events),
            'local_records': tuple(tuple(tape) for tape in local),
            'readings': tuple(local_occupancy(initial[i], i, local[i]) for i in range(12))}


@lru_cache(maxsize=1)
def _binary_record_kernel():
    states = native.states(binary=True)
    index = {state: i for i, state in enumerate(states)}
    rows = []
    for state in states:
        counts = Counter()
        for seam in source().seams():
            for coin in (0, 1):
                output = execute_record_word(state, [(seam, coin)])
                counts[index[output['readings']]] += 1
        rows.append(tuple((i, F(n, 60)) for i, n in sorted(counts.items())))
    return tuple(rows)


def native_record_law(initial, attempts):
    """Exact joint law of the twelve local readouts at one common horizon.

    The horizon counts IID uniform native seam attempts, including zero
    currents. Each branch retains its committed records; summing them here
    computes the endpoint readout law, not a replacement history law.
    """
    initial = _binary_initial(initial)
    if type(attempts) is not int or attempts < 0:
        raise ValueError('a nonnegative native attempt count is required')
    states = native.states(binary=True)
    probabilities = [F(int(state == initial)) for state in states]
    for _ in range(attempts):
        following = [F(0)]*len(states)
        for p, row in zip(probabilities, _binary_record_kernel(), strict=True):
            for index, weight in row:
                following[index] += p*weight
        probabilities = following
    mean = tuple(sum(p*state[i] for p, state in zip(probabilities, states, strict=True))
                 for i in range(12))
    return tuple(probabilities), mean


def _charge_residual():
    c = source()
    projector = c.q5_matmul([row[1:] for row in null_tomography_design()], _null_tomography_decoder())
    return [[c.q5_sub(c.q5(int(i == j)), projector[i][j]) for j in range(12)] for i in range(12)]


def _power(value, n):
    result = source().q5(1)
    while n:
        if n % 2:
            result = source().q5_mul(result, value)
        value = source().q5_mul(value, value)
        n //= 2
    return result


WITHHELD_PORTS = (3, 7, 11)
FIT_PORTS = tuple(i for i in range(12) if i not in WITHHELD_PORTS)


def _nine_port_inverse():
    """Exact inverse of the fixed nine-port design; no withheld inputs."""
    c = source()
    rows = [null_tomography_design()[i][1:]+[c.q5(int(i == j)) for j in FIT_PORTS]
            for i in FIT_PORTS]
    for j in range(9):
        pivot = next(i for i in range(j, 9) if rows[i][j] != c.q5(0))
        rows[j], rows[pivot] = rows[pivot], rows[j]
        inverse = c.q5_inv(rows[j][j])
        rows[j] = [c.q5_mul(inverse, value) for value in rows[j]]
        for i in range(9):
            if i != j:
                factor = rows[i][j]
                rows[i] = [c.q5_sub(a, c.q5_mul(factor, b))
                           for a, b in zip(rows[i], rows[j], strict=True)]
    return [row[9:] for row in rows]


def withheld_charge_predictions(nine_charges):
    """Predict ports3,7,11 from ports0,1,2,4,5,6,8,9,10 in that order."""
    if (not isinstance(nine_charges, (tuple, list)) or len(nine_charges) != 9
            or any(not isinstance(q, (tuple, list)) or len(q) != 2
                   or any(type(v) not in (int, F) for v in q) for q in nine_charges)):
        raise ValueError('nine exact Q(sqrt5) charge pairs are required')
    values = tuple((F(a), F(b)) for a, b in nine_charges)
    tensor = (source().q5(0),)+tuple(_dot(row, values) for row in _nine_port_inverse())
    return tuple(_dot(null_tomography_design()[i], tensor) for i in WITHHELD_PORTS)


def independent_history_sampling_bound(trials, tolerance):
    """Joint twelve-mean failure bound, without assuming independent ports.

    Independent repetitions of the complete preparation/history are needed.
    Binary readouts have variance <=1/4; Chebyshev plus a twelve-port union
    bound gives P(max_i |sample_mean_i-mean_i| >= tolerance) <=3/(N*tolerance^2).
    No serial samples from one native trajectory are counted as independent.
    """
    if (type(trials) is not int or trials <= 0
            or type(tolerance) not in (int, F) or tolerance <= 0):
        raise ValueError('positive independent-trial count and exact error tolerance are required')
    return min(F(1), F(3, trials)/F(tolerance)**2)


def interval_tensor_estimate(intervals, native_attempts=0):
    """Propagate twelve unconditional occupation intervals without postselection.

    Complete resolved/unknown record laws may supply these intervals. Their
    origin is a caller contract, not inferred here. Coefficient uncertainty
    bounds R*m, the canonical compatible component of the unknown mean m;
    m's fast component has the separately declared native bound. Midpoint
    compatibility never certifies exact compatibility of the unknown mean.
    """
    if (not isinstance(intervals, (tuple, list)) or len(intervals) != 12
            or any(not isinstance(row, (tuple, list)) or len(row) != 2
                   or any(type(x) not in (int, F) for x in row)
                   or not 0 <= row[0] <= row[1] <= 1 for row in intervals)):
        raise ValueError('twelve exact unconditional occupation intervals in [0,1] are required')
    if not sum(lo for lo, _ in intervals) <= 2 <= sum(hi for _, hi in intervals):
        raise ValueError('occupation intervals exclude the conserved total two')
    if type(native_attempts) is not int or native_attempts < 0:
        raise ValueError('a nonnegative uniform native attempt count is required')
    c = source()
    midpoint = tuple(c.q5((F(lo)+F(hi))/2) for lo, hi in intervals)
    radius = max((F(hi)-F(lo))/2 for lo, hi in intervals)
    tensor = (c.q5(0),)+tuple(_dot(row, midpoint) for row in _null_tomography_decoder())
    residual = tuple(_dot(row, midpoint) for row in _charge_residual())
    norm = _max_row_sum([[x] for x in residual])
    coefficient_error = c.q5_mul(c.q5(1, F(-1, 5)), c.q5(radius))
    residual_error = c.q5_mul(c.GOLDEN, c.q5(radius))
    lower = c.q5_sub(norm, residual_error)
    if c.q5_sign(lower) < 0:
        lower = c.q5(0)
    upper = c.q5_add(norm, residual_error)
    fast_bound = c.q5_mul(c.q5(F(1, 4), F(1, 20)),
                          _power(c.q5(F(11, 12), F(-1, 60)), native_attempts))
    compatible = all(x == c.q5(0) for x in residual)
    if compatible and decode_null_charges(midpoint) != tensor:
        raise AssertionError('an exactly compatible midpoint has inconsistent decoding')
    return {'schema': 'oph.common_source.interval_tensor.v1',
            'declared_uniform_native_attempts': native_attempts,
            'midpoint_readings_Qsqrt5': _encode_vector(midpoint),
            'canonical_tensor_estimate_Qsqrt5': _encode_vector(tensor),
            'maximum_interval_halfwidth': str(radius),
            'canonical_component_error_upper_Qsqrt5': c.q5_str(coefficient_error),
            'midpoint_compatibility_residual_Qsqrt5': _encode_vector(residual),
            'midpoint_exactly_tensor_compatible': compatible,
            'compatibility_residual_error_upper_Qsqrt5': c.q5_str(residual_error),
            'unknown_mean_residual_norm_interval_Qsqrt5': _encode_vector((lower, upper)),
            'native_fast_residual_upper_Qsqrt5': c.q5_str(fast_bound),
            'native_bound_excluded_by_intervals': c.q5_sign(c.q5_sub(lower, fast_bound)) > 0,
            'withheld_error_upper_under_native_bound_Qsqrt5': c.q5_str(c.q5_mul(
                c.q5(2, 2), c.q5_add(fast_bound, c.q5(radius)))),
            'scope': {'input_intervals_are_unconditional_caller_contract': True,
                      'unknown_mean_exact_tensor_compatibility_asserted': False,
                      'native_bound_pass_is_sufficient_box_feasibility_test': False,
                      'finite_sample_uncertainty_added_automatically': False}}


def native_record_tensor_certificate():
    """Record-to-tensor attachment with an explicit fast-band error budget."""
    c = source()
    null_tomography_certificate()
    residual = _charge_residual()
    states = native.states(binary=True)
    fast = c.q5(F(55, 60), F(-1, 60))
    mean_operator = [[c.q5(F(int(i == j))-F(int(x), 60)) for j, x in enumerate(row)]
                     for i, row in enumerate(c.laplacian())]
    for initial, row in zip(states, _binary_record_kernel(), strict=True):
        observed = [c.q5(sum(p*states[j][i] for j, p in row)) for i in range(12)]
        if observed != [_dot(r, [c.q5(x) for x in initial]) for r in mean_operator]:
            raise AssertionError('actual record means do not implement I-L/60')
    if c.q5_matmul(residual, mean_operator) != [[c.q5_mul(fast, x) for x in row] for row in residual]:
        raise AssertionError('native transport does not give the declared fast residual decay')
    bound = c.q5(F(1, 4), F(1, 20))
    attained, compatible = False, []
    for state in states:
        r = tuple(_dot(row, [c.q5(x) for x in state]) for row in residual)
        compatible.append(all(x == c.q5(0) for x in r))
        for x in r:
            absolute = x if c.q5_sign(x) >= 0 else c.q5_sub(c.q5(0), x)
            if c.q5_sign(c.q5_sub(absolute, bound)) > 0:
                raise AssertionError('the uniform binary-sector fast residual bound failed')
            attained = attained or absolute == bound
    if sum(compatible) != 6 or not attained:
        raise AssertionError('exact binary compatible preparations or sharp residual bound changed')
    # Every rational odd vector with zero fast part also has zero slow part:
    # field conjugation exchanges the two projectors, whose sum is odd parity.
    slow = c.slow_band_projector_exact()
    if residual != [[(a, -b) for a, b in row] for row in slow]:
        raise AssertionError('the rational-preparation conjugation control failed')
    inverse = _nine_port_inverse()
    sliced = [null_tomography_design()[i][1:] for i in FIT_PORTS]
    if c.q5_matmul(inverse, sliced) != [[c.q5(int(i == j)) for j in range(9)] for i in range(9)]:
        raise AssertionError('the fixed nine-port design cannot predict the withheld directions')
    prediction = c.q5_matmul([null_tomography_design()[i][1:] for i in WITHHELD_PORTS], inverse)
    withheld_operator = [[prediction[k][FIT_PORTS.index(i)] if i in FIT_PORTS else c.q5(-int(i == port))
                          for i in range(12)] for k, port in enumerate(WITHHELD_PORTS)]
    if any(x != c.q5(0) for row in c.q5_matmul(withheld_operator, null_tomography_design()) for x in row):
        raise AssertionError('withheld errors do not annihilate compatible tensor charges')
    withheld_amplification = _max_row_sum(withheld_operator)
    examples = []
    for label, pair in (('antipodal', (0, 3)), ('adjacent', (0, 1))):
        initial = tuple(int(i in pair) for i in range(12))
        initial_residual = tuple(_dot(row, [c.q5(x) for x in initial]) for row in residual)
        for n in (0, 1, 2, 16, 64):
            law, mean = native_record_law(initial, n)
            charges = tuple(c.q5(x) for x in mean)
            actual = tuple(_dot(row, charges) for row in residual)
            expected = tuple(c.q5_mul(_power(fast, n), x) for x in initial_residual)
            if actual != expected or sum(law) != 1:
                raise AssertionError('actual record law does not match native residual evolution')
            predicted = withheld_charge_predictions([charges[i] for i in FIT_PORTS])
            differences = tuple(c.q5_sub(predicted[j], charges[i]) for j, i in enumerate(WITHHELD_PORTS))
            if label == 'antipodal' and any(x != c.q5(0) for x in actual):
                raise AssertionError('the antipodal preparation lost exact compatibility')
            if label == 'antipodal':
                decode_null_charges(charges)
                if any(x != c.q5(0) for x in differences):
                    raise AssertionError('withheld native means do not match the common tensor')
            else:
                try:
                    decode_null_charges(charges)
                except ValueError:
                    pass
                else:
                    raise AssertionError('an incompatible finite native mean was silently projected')
                if all(x == c.q5(0) for x in differences):
                    raise AssertionError('withheld directions failed to detect the native fast band')
            examples.append({'preparation': label, 'native_attempts': n,
                             'endpoint_readout_law': list(map(str, law)),
                             'mean_readings': list(map(str, mean)),
                             'canonical_tensor_estimate_Qsqrt5': _encode_vector(
                                 (c.q5(0),)+tuple(_dot(row, charges) for row in _null_tomography_decoder())),
                             'compatibility_residual_Qsqrt5': _encode_vector(actual),
                             'withheld_errors_Qsqrt5': _encode_vector(differences),
                             'uniform_residual_bound_Qsqrt5': c.q5_str(c.q5_mul(bound, _power(fast, n)))})
    return {'schema': 'oph.common_source.native_record_tensor.v1',
            'native_local_record_branches_checked': len(states)*60,
            'binary_preparations_exactly_compatible': sum(compatible),
            'fast_decay_per_uniform_attempt_Qsqrt5': c.q5_str(fast),
            'sharp_binary_initial_fast_bound_Qsqrt5': c.q5_str(bound),
            'rational_exact_compatibility_requires_antipodal_even_mean': True,
            'fit_ports': list(FIT_PORTS), 'withheld_ports': list(WITHHELD_PORTS),
            'withheld_error_amplification_Qsqrt5': c.q5_str(withheld_amplification),
            'examples': examples,
            'scope': {'readout_uses_initial_local_bit_and_own_currents': True,
                      'uniform_independent_seam_schedule_supplied': True,
                      'endpoint_means_are_not_complete_history_laws': True,
                      'arbitrary_individual_readings_exactly_compatible': False,
                      'physical_stress_tensor_or_lorentz_lift_derived': False}, 'all_pass': True}


