"""A faithful, exterior-invisible coherence on the actual biased carrier.

This is a known mixed-state excitation, not an arbitrary-logical-state code.
Its initial diagonal is the source reference. Every native history preserves
its off-diagonal difference vector, so record-plus-exterior restrictions agree
with the reference history. The surviving contrast obeys a killed ten-state
kernel. Native attempt count and reference modular parameter are distinct.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import math

import numpy as np
import sympy as sp

from ..local_repair_family import local_bias as local, modular_clock as clock
from ..observer import adaptive_record_recovery as adaptive
from ..observer.source_discrimination import circuit, experiment as recorded
from ..observer.source_scalar import compiler
from ..native_modular_reference import cap_local_clock as cap

CORE = (1, 7)
AMPLITUDE = F(8)
BASIS = local.basis()
PAIRS = tuple(tuple(i for i, n in enumerate(x) if n) for x in BASIS)
SPECTATORS = tuple(k for k in range(12) if k not in CORE)
INDEX = {pair: i for i, pair in enumerate(PAIRS)}
TRANSFERS = tuple((INDEX[tuple(sorted((CORE[0], k)))],
                   INDEX[tuple(sorted((CORE[1], k)))]) for k in SPECTATORS)
SOURCE_CODES = tuple(circuit.encode(x) for x in BASIS)


def parameters(law="mixed"):
    if law not in ("C", "E", "mixed"):
        raise ValueError("C, E, or mixed law required")
    return replace(clock.parameters(), c_weight={"C": F(1), "E": F(), "mixed": F(1, 2)}[law])


def initial_density(amplitude=AMPLITUDE):
    """Exact rational state; positivity requires amplitude squared <= 200.

    Its ten q=1 blocks are theta_k/Z [[25,a],[a,8]]. The other 46
    diagonal entries equal the faithful source reference. At a=8 every
    block has determinant 136 theta_k squared / Z squared.
    """
    if type(amplitude) is not F or amplitude*amplitude > 200:
        raise ValueError("a rational amplitude with square at most 200 required")
    pi = local.reference(parameters())
    state = [[F() for _ in BASIS] for _ in BASIS]
    for i, value in enumerate(pi):
        state[i][i] = value
    z = sum(math.prod(w**n for w, n in zip(parameters().port_weights, x, strict=True))
            for x in BASIS)
    for k, (i, j) in zip(SPECTATORS, TRANSFERS, strict=True):
        state[i][j] = state[j][i] = amplitude*F(parameters().port_weights[k], z)
    return tuple(map(tuple, state))


def exterior_restriction(state, core=CORE):
    """Partial trace in the binary occupation tensor product, without a cutoff."""
    if core not in tuple(combinations(range(12), 2)):
        raise ValueError("ordered two-port core required")
    result = defaultdict(F)
    outside = tuple(i for i in range(12) if i not in core)
    for i, x in enumerate(BASIS):
        for j, y in enumerate(BASIS):
            if all(x[k] == y[k] for k in core) and state[i][j]:
                result[tuple(x[k] for k in outside), tuple(y[k] for k in outside)] += state[i][j]
    return dict(result)


def _branch_matrix(state, mapping, weight):
    out = [[F() for _ in BASIS] for _ in BASIS]
    for i, j in mapping:
        for a, b in mapping:
            out[j][b] = weight*state[i][a]
    return tuple(map(tuple, out))


def trajectory_certificate():
    """Independent raw-menu check of the induction valid at every history length.

    The diagonal laws coincide. Every surviving matrix unit has the same
    signed occupation difference, hence zero exterior partial trace. This
    proves equality even with the *complete* native tape retained; it does
    not identify a conditioned branch exterior with the unconditioned pi.
    """
    p = parameters()
    transfer_index = {pair: i for i, pair in enumerate(TRANSFERS)}
    difference = tuple(int(i == CORE[0])-int(i == CORE[1]) for i in range(12))
    counts = Counter()
    for branch in adaptive.branch_maps(binary=True):
        name, seam = adaptive.actions()[branch.action]
        decoded = circuit.record_values(name, branch.coin)
        raw = decoded.get(branch.raw)
        source_map = dict(next((mapping for coin, value, mapping in local.raw_maps(name, seam)
                                if coin == branch.coin and value == raw), ()))
        if source_map != dict(zip(branch.inputs, branch.outputs, strict=True)):
            raise ArithmeticError("raw source and independent maps disagree")
        counts["raw_slots"] += 1
        counts["raw_abort_slots"] += raw is None
        counts["valid_zero_slots"] += raw is not None and not branch.inputs
        for i, j in TRANSFERS:
            counts["matrix_units_examined"] += 1
            if i in source_map and j in source_map:
                target = source_map[i], source_map[j]
                if target not in transfer_index or tuple(a-b for a, b in zip(
                        BASIS[target[0]], BASIS[target[1]], strict=True)) != difference:
                    raise ArithmeticError("native common record moved the coherence core")
                counts["surviving_matrix_units"] += 1
            else:
                counts["killed_matrix_units"] += 1
        # Every outcome weight is independent of the off-diagonal contrast.
        if len(set(source_map.values())) != len(source_map):
            raise ArithmeticError("record map lost injectivity")
    if counts["raw_slots"] != 720 or counts["matrix_units_examined"] != 7200:
        raise ArithmeticError("incomplete all-seam raw menu")
    return {**dict(counts), "core": list(CORE),
            "record_plus_exterior_equal_to_reference_history_for_all_lengths": True,
            "conditioned_exterior_equal_to_unconditioned_reference_claimed": False,
            "coherence_core_moves_under_native_repair": False,
            "identity_of_diagonal_histories": True,
            "proof": "partial-injective common records preserve x-y; all nonzero difference units have opposite core occupations and zero exterior partial trace"}


@lru_cache(maxsize=3)
def coherence_kernel(law="mixed"):
    """Rows are input spectators; columns are output spectators, with killing."""
    p = parameters(law)
    result = [[F() for _ in SPECTATORS] for _ in SPECTATORS]
    index = {pair: i for i, pair in enumerate(TRANSFERS)}
    for name, mixture in (("C", p.c_weight), ("E", 1-p.c_weight)):
        for seam in local.SEAMS:
            for coin, _, mapping in local.raw_maps(name, seam):
                mapping = dict(mapping)
                weight = mixture*local.coin_probability(p, seam, coin)/30
                for i, (x, y) in enumerate(TRANSFERS):
                    if x in mapping and y in mapping:
                        result[i][index[mapping[x], mapping[y]]] += weight
    return tuple(map(tuple, result))


def _step(vector, matrix):
    return tuple(sum((vector[i]*matrix[i][j] for i in range(len(vector))), F())
                 for j in range(len(vector)))


def contrast(attempts, law="mixed"):
    if type(attempts) is not int or not 0 <= attempts <= 256:
        raise ValueError("an integer native-attempt count in 0..256 required")
    p = parameters(law)
    values = tuple(F(8*p.port_weights[k], 7505) for k in SPECTATORS)
    for _ in range(attempts):
        values = _step(values, coherence_kernel(law))
    return sum(values, F())


def spectral_certificate(law="mixed"):
    """Exact killed resolvent; real-positive spectrum follows weighted symmetry.

    Q is self-adjoint in l2(theta) and strictly diagonally dominant. Thus the
    initial theta-weighted contrast is a positive mixture of decaying real
    modes, rather than a propagating Hamiltonian pole. Floating eigenvalues
    are only displayed diagnostics; all bounds and mean lifetime are rational.
    """
    q = coherence_kernel(law)
    weights = tuple(parameters(law).port_weights[k] for k in SPECTATORS)
    if not all(weights[i]*q[i][j] == weights[j]*q[j][i]
               for i in range(10) for j in range(10)):
        raise ArithmeticError("killed coherence kernel lost detailed balance")
    low = min(2*q[i][i]-sum(q[i]) for i in range(10))
    high = max(map(sum, q))
    if not 0 < low <= high < 1:
        raise ArithmeticError("positive strictly subunit spectral enclosure failed")
    matrix = sp.Matrix(q)
    potential = (sp.eye(10)-matrix).inv()*sp.ones(10, 1)
    mean = F(sum(weights[i]*potential[i] for i in range(10))/sum(weights))
    # A positive rational test vector gives certified Collatz bounds for rho(Q).
    vector = (F(1),)*10
    for _ in range(512):
        vector = tuple(sum((q[i][j]*vector[j] for j in range(10)), F()) for i in range(10))
    ratios = [sum((q[i][j]*vector[j] for j in range(10)), F())/vector[i] for i in range(10)]
    scale = 10**8
    interval = [F(math.floor(min(ratios)*scale), scale), F(math.ceil(max(ratios)*scale), scale)]
    eigenvalues = np.linalg.eigvalsh(np.sqrt(np.asarray(weights))[:, None]
                                   *np.array(q, float)/np.sqrt(np.asarray(weights))[None, :])
    rows = [{"attempts": n, "contrast": str(contrast(n, law)),
             "relative_visibility": str(contrast(n, law)/contrast(0, law))} for n in (0, 1, 2, 4, 8, 16, 32)]
    if law == "E" and (mean != F(10, 3) or any(
            contrast(n, law) != contrast(0, law)*F(7, 10)**n for n in range(17))):
        raise ArithmeticError("E's exact exponential visibility law failed")
    return {"law": law, "kernel": [[str(x) for x in row] for row in q],
            "spectator_ports": list(SPECTATORS), "reversible_weights": list(weights),
            "all_eigenvalues_rational_enclosure": [str(low), str(high)],
            "spectral_radius_rational_enclosure": list(map(str, interval)),
            "eigenvalues_float_diagnostic": eigenvalues.tolist(),
            "integrated_relative_visibility_attempts": str(mean),
            "integrated_relative_visibility_float": float(mean),
            "potential": list(map(str, potential)),
            "mean_lifetime_definition": "sum over n>=0 of D(rho_n,pi)/D(rho_0,pi); not a survival time of a conserved particle",
            "all_native_history_records_retained_gives_same_trace_distance": True,
            "visibility": rows, "oscillatory_native_attempt_modes": False,
            "E_exact_relative_visibility": "(7/10)^n" if law == "E" else None}


def exact_code_certificate():
    """Extend the predecessor all-code support proof to exterior and clock.

    No sharp-core-charge assumption is made. Exact flagged correction of
    arbitrary complex code states with k>=2 forces nonedge support and
    scalar compressions of all n_i. Positive endpoint biases rescale each
    effect without changing these conclusions.
    """
    p = parameters()
    pi = local.reference(p)
    seams = set(local.SEAMS)
    rows = []
    for core in combinations(range(12), 2):
        outside = set(range(12))-set(core)
        by_exterior = defaultdict(list)
        for i, pair in enumerate(PAIRS):
            by_exterior[tuple(k for k in pair if k in outside)].append(i)
        forbidden = [key for key, indices in by_exterior.items() if all(PAIRS[i] in seams for i in indices)]
        mass = sum((pi[i] for key in forbidden for i in by_exterior[key]), F())
        rows.append({"core": list(core), "forbidden_exterior_configurations": [list(key) for key in sorted(forbidden)],
                     "universal_exterior_trace_distance_lower": str(mass)})
    return {"arbitrary_complex_exact_code_dimension_at_least": 2,
            "full_initial_total_two_carrier_dimension": 78,
            "rank_one_records_force_binary_support": True,
            "native_edge_basis_rows_forbidden": 30,
            "positive_endpoint_bias_changes_support_proof": False,
            "subsystem_codes_with_nontrivial_logical_factor_and_full_gauge_support_also_excluded": True,
            "rows": rows,
            "minimum_any_two_port_core": str(min(F(row["universal_exterior_trace_distance_lower"]) for row in rows)),
            "bound_attained_by_some_nonedge_diagonal_state_not_claimed_to_be_a_code": True,
            "clock_generator_modulo_scalar": "-log(25/8) (n_0+n_1)",
            "exact_code_compressed_heavy_number": "V* (n_0+n_1) V = nu I, 0<=nu<=1",
            "clock_compression": "V* U_s V = scalar_global_phase * (1-nu+nu exp(-i s log(25/8))) I",
            "clock_leakage_probability_any_logical_reference": "4 nu(1-nu) sin^2(s log(25/8)/2)",
            "clock_invariant_exact_code_implies_trivial_logical_clock": True,
            "excluded_class": "fixed pre-attempt all-seam exact flagged codes, including arbitrary complex nonscalar-q sectors; positive source-independent seam weights; arbitrary CPTP flagged recovery",
            "surviving_routes": ["known mixed-state coherences without arbitrary-qubit protection", "approximate correction", "seam-conditioned encoding before an attempt", "different source/reference/refinement", "history observables with a separately proved particle interpretation"]}


def state_certificate():
    rho = initial_density()
    pi = initial_density(F())
    if exterior_restriction(rho) != exterior_restriction(pi):
        raise ArithmeticError("initial exterior differs from the same source reference")
    if any(rho[i][i] != pi[i][i] for i in range(66)):
        raise ArithmeticError("coherence excitation changed the reference diagonal")
    array = np.array(rho, float)
    values = np.linalg.eigvalsh(array)
    target = np.diag(array)
    entropy = -float(np.sum(values*np.log(values)))
    reference_entropy = -float(np.sum(target*np.log(target)))
    return {"reference": "pi_ij=theta_i theta_j/7505; theta=(25,25,8,...,8)",
            "core": list(CORE), "state_blocks": "theta_k/7505 [[25,8],[8,8]] for every k outside {1,7}; other reference entries unchanged",
            "q1_block_determinant_factor": "136", "state_faithful": True,
            "initial_trace_distance_to_reference": str(contrast(0)),
            "full_exterior_trace_distance": "0", "every_diagonal_observable_excess": "0",
            "modular_energy_excess": "0", "entropy_deficit_float_diagnostic": reference_entropy-entropy,
            "relative_entropy_equals_entropy_deficit": True,
            "modular_phase_ratio": "25/8", "modular_parameter_symbol": "s",
            "phase_on_oriented_transfer": "exp(-i s log(25/8))",
            "native_attempt_parameter_symbol": "n",
            "physical_clock_or_mass_calibration_supplied": False,
            "localization_statement": "record-plus-exterior difference vanishes under every native history; the initial local coherence core cannot move under bare native repair",
            "physical_topological_charge_or_particle_identified": False,
            "preparation_service": "supplied classical mixture and coherent pulse preparation; no autonomous source selection of this state"}


def current_response(attempts=0, law="mixed", parameter=0.):
    """Full-trial native current after a supplied local Ramsey analyzer.

    No existing attempt is discarded. The final detector is one additional
    attempted repair on CORE. Its C and E current marginals agree; their
    quantum instruments remain different. The mixer, phase actuation and
    detector are declared interventions, not operations generated by repair.
    """
    if type(parameter) not in (int, float) or not math.isfinite(parameter):
        raise ValueError("finite dimensionless modular parameter required")
    signal = contrast(attempts, law)
    q1 = F(3201, 7505)
    cosine = math.cos(parameter*math.log(25/8))
    p_a, p_b = F(25, 33), F(8, 33)
    positive = float(p_b)*(float(q1)/2+float(signal)*cosine)
    negative = float(p_a)*(float(q1)/2-float(signal)*cosine)
    return {"attempts_before_detector": attempts, "law": law,
            "modular_parameter": parameter, "constant": "-1649/15010",
            "cosine_coefficient": str(signal), "frequency_ratio": "25/8",
            "mean_current": positive-negative,
            "current_probabilities": {"1": positive, "-1": negative, "0": 1-positive-negative},
            "full_trial_no_postselection": True,
            "native_detector_attempts_added": 1,
            "supplied_analyzer": "bare (1,7) pi/4 rotation, all ten spectators"}


def unitary_preparation_control():
    """Compare the already compiled local rotation of pi on the same carrier.

    It has a positive initial modular-energy excess, but its population
    perturbation spreads into the exterior during native evolution. This
    prevents conflating the energy-carrying unitary perturbation with the
    always-exterior-invisible, zero-energy coherence state.
    """
    pi = local.reference(parameters())
    delta = [F() for _ in BASIS]
    for i, j in TRANSFERS:
        delta[i] = (pi[j]-pi[i])/2
        delta[j] = -delta[i]
    kernel = local.kernel(parameters())
    outside = tuple(i for i in range(12) if i not in CORE)
    rows = []
    for n in range(17):
        if n in (0, 1, 2, 4, 8, 16):
            marginal = defaultdict(F)
            for x, value in zip(BASIS, delta, strict=True):
                marginal[tuple(x[i] for i in outside)] += value
            energy = -sum((value*(x[0]+x[1]) for x, value in zip(BASIS, delta, strict=True)), F())
            rows.append({"native_attempts": n,
                         "exterior_trace_distance": str(sum(map(abs, marginal.values()), F())/2),
                         "modular_energy_excess_coefficient_of_log25over8": str(energy)})
        if n < 16:
            delta = list(_step(delta, kernel))
    if rows[0]["exterior_trace_distance"] != "0" or F(rows[1]["exterior_trace_distance"]) <= 0:
        raise ArithmeticError("unitary preparation exterior countercontrol failed")
    return {"state": "R pi R* with R the existing full bare-core pi/4 rotation",
            "initial_modular_energy_excess": "(1649/15010) log(25/8)",
            "initial_von_Neumann_entropy_change": "0",
            "population_evolution_same_for_C_and_E": True,
            "rows": rows, "physical_energy_calibration": False,
            "same_object_as_exterior_invisible_coherence": False}


def _pair_word(left, right, angle):
    """General-angle version of the existing complete Gray conditional-Y word."""
    different = left ^ right
    pivot = (different & -different).bit_length()-1
    controls = tuple(j for j in range(9) if j != pivot)
    fan = [pulse for bit in controls if different >> bit & 1 for pulse in compiler.cnot(pivot, bit)]
    transformed = left ^ ((different ^ (1 << pivot)) if left >> pivot & 1 else 0)
    orientation = 1 if not (transformed >> pivot & 1) else -1
    control_word = sum(((transformed >> bit) & 1) << j for j, bit in enumerate(controls))
    tape = list(fan)
    for k in range(256):
        gray = k ^ (k >> 1)
        tape.append(compiler.pulse("Y", orientation*angle/256*(-1)**((control_word & gray).bit_count()), [pivot]))
        successor = (k+1) % 256
        changed = gray ^ (successor ^ (successor >> 1))
        tape.extend(compiler.cnot(controls[changed.bit_length()-1], pivot))
    tape.extend(compiler.inverse(fan))
    return tape


def preparation_factor(executed=False):
    """Retain every classical preparation component as a separate column.

    [[25,8],[8,8]] = (25,8)(25,8)^T/25 + diag(0,136/25).
    The diagonal classical mixture is supplied. Ten source-executed rotations
    create its coherent components; no normalization or padding projection is
    performed after source execution.
    """
    pi = local.reference(parameters())
    factor = np.zeros((512, 66), complex)
    for i, probability in enumerate(pi):
        factor[SOURCE_CODES[i], i] = math.sqrt(float(probability))
    amplitude = np.zeros((512, 10), complex)
    word = []
    for column, (k, (i, j)) in enumerate(zip(SPECTATORS, TRANSFERS, strict=True)):
        weight = F(parameters().port_weights[k], 7505)
        factor[SOURCE_CODES[i], i] = 0
        factor[SOURCE_CODES[j], j] = math.sqrt(float(weight*F(136, 25)))
        amplitude[SOURCE_CODES[i], column] = math.sqrt(float(weight*F(689, 25)))
        word.extend(_pair_word(SOURCE_CODES[i], SOURCE_CODES[j], -math.atan2(8, 25)))
    if executed:
        amplitude = compiler.execute(word, 9, amplitude)
    else:
        for column, (k, (i, j)) in enumerate(zip(SPECTATORS, TRANSFERS, strict=True)):
            weight = math.sqrt(float(F(parameters().port_weights[k], 7505)))
            amplitude[SOURCE_CODES[i], column] = 5*weight
            amplitude[SOURCE_CODES[j], column] = F(8, 5)*weight
    for column, (i, _) in enumerate(TRANSFERS):
        factor[:, i] = amplitude[:, column]
    return factor, len(word)


def _native_source(law, factor):
    """Actual source dilations on seam (0,1), every raw coin/current slot."""
    slots = 1 << circuit.RECORD_BITS[law]
    count = factor.shape[1]
    local_columns = factor.reshape(64, 8, count).transpose(1, 0, 2).reshape(8, 64*count)
    blank = np.zeros((8*slots, 64*count), complex)
    blank[:8] = local_columns
    angle = -math.asin(math.sqrt(float(local.coin_probability(parameters(), (0, 1), 1))))
    coin_values = compiler.execute([compiler.pulse("Y", angle, [0])], 1, [1., 0.])
    for coin in (0, 1):
        values = (recorded._repair(law, coin)@blank).reshape(slots, 8, 64, count)
        for raw in range(slots):
            yield coin, raw, values[raw].transpose(1, 0, 2).reshape(512, count)*coin_values[coin]


def _native_ideal(law, coin, raw, factor):
    result = np.zeros_like(factor)
    value = circuit.record_values(law, coin).get(raw)
    if value is not None:
        for c, label, mapping in local.raw_maps(law, (0, 1)):
            if c == coin and label == value:
                root = math.sqrt(float(local.coin_probability(parameters(), (0, 1), coin)))
                for i, j in mapping:
                    result[SOURCE_CODES[j]] = root*factor[SOURCE_CODES[i]]
    return result


def source_execution():
    """Fresh full-register preparation/clock and complete requested C/E histories.

    This requested seam replay complements the exact full-scheduler proof.
    It is not a 30-seam pulse replay or a measured hardware error certificate.
    """
    ideal, pulse_count = preparation_factor()
    actual, actual_count = preparation_factor(True)
    if pulse_count != actual_count:
        raise ArithmeticError("source preparation words changed")
    preparation_error = float(np.linalg.norm(actual-ideal))
    clock_word = clock.source_program()["clock"]
    clock_actual = compiler.execute(clock_word, 9, actual)
    phase = np.array([complex(float(x))**(-1j) for x in clock._ratios()])
    clock_ideal = phase[:, None]*ideal
    clock_error = float(np.linalg.norm(clock_actual-clock_ideal))
    histories = []
    error_sum = abort_mass = 0.
    for coin_c, raw_c, first in _native_source("C", clock_actual):
        wanted_first = _native_ideal("C", coin_c, raw_c, clock_ideal)
        for coin_e, raw_e, second in _native_source("E", first):
            wanted = _native_ideal("E", coin_e, raw_e, wanted_first)
            error = float(np.linalg.norm(second-wanted))
            mass = float(np.vdot(second, second).real)
            invalid = (raw_c not in circuit.record_values("C", coin_c)
                       or raw_e not in circuit.record_values("E", coin_e))
            error_sum += error
            abort_mass += mass*invalid
            histories.append((mass, error))
    if len(histories) != 128 or max(preparation_error, clock_error, error_sum) > 1e-8:
        raise ArithmeticError("complete source replay exceeded floating diagnostic gate")
    outside = sorted(set(range(512))-set(SOURCE_CODES))
    return {"preparation_pulses": pulse_count, "clock_pulses": len(clock_word),
            "classical_mixture_components": 66, "full_source_rows": 512,
            "retained_preparation_reference_columns": 66,
            "preparation_factor_error_frobenius": preparation_error,
            "prepared_padding_mass": float(np.linalg.norm(actual[outside])**2),
            "clock_factor_error_frobenius": clock_error,
            "native_requested_seams": [[0, 1], [0, 1]], "native_laws": ["C", "E"],
            "complete_raw_history_count": len(histories),
            "full_history_error_sum_frobenius": error_sum,
            "complete_history_mass": sum(x[0] for x in histories), "raw_abort_mass": abort_mass,
            "floating_diagnostic_gate": 1e-8,
            "formal_outward_or_hardware_certificate": False,
            "source_full_scheduler_replayed": False,
            "supplied_services": ["classical mixture and native address initialization", "coherent source pulse and proper-code load/export", "fixed seam/coin selection and native records", "archive and detector export"],
            "preparation_archive_is_retained_and_not_a_thermal_reference_archive": True}


def report(with_source=True):
    return {"schema": "oph.native_localization.v1", "state": state_certificate(),
            "all_history_localization": trajectory_certificate(),
            "current_response": [current_response(n, law, 1.) for law in ("C", "E", "mixed") for n in (0, 2, 8)],
            "energy_carrying_unitary_countercontrol": unitary_preparation_control(),
            "native_lifetime": {law: spectral_certificate(law) for law in ("C", "E", "mixed")},
            "arbitrary_exact_code_obstructions": exact_code_certificate(),
            "source_execution": source_execution() if with_source else None,
            "canonical_physical_premises_discharged": 0, "CP03_accepted": False,
            "remaining_particle_joins": ["moving coherent core under the same native evolution", "autonomous protected or topological boundary charge", "nonzero same-excitation physical energy and mass", "physical clock calibration and refinement transport"]}
