"""Information in the actual native tape and its retained quantum output.

For a complete finite instrument, D(input||reference) = D(tape||reference
tape) + mean conditional D + information loss. The loss is nonnegative by
data processing; it is not a heat, area, or physical entropy-production law.
Native branch translations let us evaluate the quantum blocks in their
actual image supports without allocating empty endpoint/padding rows.

Numerical entropies include EVERY strictly positive computed eigenvalue,
however small. A tolerance only guards negative eigensolver roundoff; it does
not infer support or rank. These floating diagnostics are not outward entropy
error certificates. Analytic identities, exact native support counts and
Fraction completeness checks do not depend on the numerical eigensolver.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from itertools import product
import math

import numpy as np

from ..native_equilibration import history_entropy as population
from ..native_modular_reference import local_preparation
from ..observer import native
from .protocol import _original_scalars


ZERO_TOLERANCE = 1e-12
CHECK_TOLERANCE = 2e-10


def _density(value, dimension):
    _original_scalars(value)
    if np.ma.isMaskedArray(value):
        raise ValueError("masked density operators are inadmissible")
    raw = np.asarray(value)
    if (raw.shape != (dimension, dimension) or raw.dtype.kind not in "iufc"
            or not np.isfinite(raw).all()):
        raise ValueError("a finite density matrix with the declared shape required")
    result = np.array(raw, dtype=complex, copy=True)
    # PSD requires each original diagonal population to be nonnegative. This
    # check has no eigensolver roundoff and must not erase a small negative
    # supplied population under the spectral diagnostic's tolerance.
    if np.any(np.diag(raw).real < 0):
        raise ValueError("nonnegative original diagonal populations required")
    if np.linalg.norm(result-result.conj().T, 'fro') > ZERO_TOLERANCE:
        raise ValueError("a Hermitian density matrix required")
    if abs(np.trace(result).real-1) > ZERO_TOLERANCE:
        raise ValueError("a normalized density matrix required")
    _entropy(result)
    return result


def _entropy(rho):
    """Floating entropy diagnostic; every positive eigenvalue contributes."""
    values = np.linalg.eigvalsh((rho+rho.conj().T)/2)
    if values[0] < -ZERO_TOLERANCE:
        raise ValueError("a positive semidefinite density matrix required")
    return -math.fsum(float(x)*math.log(float(x)) for x in values if x > 0)


def _shannon(weights):
    return -math.fsum(float(p)*math.log(float(p)) for p in weights if p)


def _reference_marginal(rho, dimension, reference_dimension):
    return np.trace(rho.reshape(dimension, reference_dimension,
                                dimension, reference_dimension), axis1=0, axis2=2)


def _native_blocks(paths, rho, reference_dimension):
    """Exact partial injections, with all quantum matrix elements retained.

    Each returned block uses the ordered physical endpoint basis of that
    tape's image. Permuting it to input order is a unitary change of basis,
    not a trace, a diagonalization of the state, or a conditional reset.
    """
    states = native.states(binary=True)
    sigma = local_preparation.reference()
    rows = population.fibres(paths)
    beta = population.flag_probability(paths)
    for (tape, _), paths in rows.items():
        ordered = sorted(paths, key=lambda path: states.index(path.endpoint))
        indices = [path.initial*reference_dimension+r
                   for path in ordered for r in range(reference_dimension)]
        chance = ordered[0].chance/beta
        block = float(chance)*rho[np.ix_(indices, indices)]
        reference_diagonal = np.array([float(chance*sigma[path.initial])/reference_dimension
                                       for path in ordered for _ in range(reference_dimension)])
        yield tape, ordered, block, reference_diagonal


def native_history_ledger(rho, *, law="E", prefix=3, deadline=3,
                          stop_on_current=True, success_only=False,
                          reference_dimension=1):
    """Actual adaptive native history relative to sigma_native tensor I_R/r.

    Arbitrary mixed/coherent inputs on the binary carrier and 1..6 untouched
    reference levels are admitted. Spectral values are floating diagnostics.
    Success-only C is the existing explicitly flagged state-independent channel;
    the default retains all timeout paths. No branch is normalized away.
    """
    if type(reference_dimension) is not int or not 1 <= reference_dimension <= 6:
        raise ValueError("one through six finite reference levels required")
    if type(success_only) is not bool:
        raise ValueError("success_only must be an explicit Boolean")
    # Validate before the shared population enumerator's cache: its keys can
    # identify 3, 3.0 and Fraction(3), and E does not otherwise use deadline.
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("C or E required")
    if type(prefix) is not int or not 0 <= prefix <= 5:
        raise ValueError("integer native prefix from zero through five required")
    if type(deadline) is not int or not 2 <= deadline <= 16:
        raise ValueError("integer census deadline from two through sixteen required")
    if type(stop_on_current) is not bool:
        raise ValueError("stop_on_current must be an explicit Boolean")
    dimension = 66
    rho = _density(rho, dimension*reference_dimension)
    paths = population.histories(law, prefix, deadline, stop_on_current)
    if success_only:
        paths = tuple(path for path in paths if path.terminal == "success")
    beta = population.flag_probability(paths)
    sigma = np.array(list(map(float, local_preparation.reference())))
    reference_diagonal = np.repeat(sigma/reference_dimension, reference_dimension)
    input_entropy = _entropy(rho)
    input_cross_entropy = -float(np.dot(np.diag(rho).real, np.log(reference_diagonal)))
    input_relative = input_cross_entropy-input_entropy
    probabilities, reference_probabilities = [], []
    conditional, conditional_entropy, counting, attempts = [], [], [], []
    reference_output = np.zeros((reference_dimension, reference_dimension), complex)
    charge_resolved = True
    for tape, row, block, target in _native_blocks(paths, rho, reference_dimension):
        p, q = float(np.trace(block).real), float(target.sum())
        probabilities.append(p)
        reference_probabilities.append(q)
        reference_output += _reference_marginal(block, len(row), reference_dimension)
        charge_resolved &= len({path.endpoint for path in row}) == len(row)
        charge_resolved &= len({local_preparation.reference()[path.initial] for path in row}) == 1
        if p > 0:
            state = block/p
            entropy = _entropy(state)
            cross = -float(np.dot(np.diag(state).real, np.log(target/q)))
            conditional.append(p*(cross-entropy))
            conditional_entropy.append(p*entropy)
            counting.append(p*math.log(len(row)))
            attempts.append(p*len(tape))
    if abs(math.fsum(probabilities)-1) > CHECK_TOLERANCE:
        raise ArithmeticError("complete native quantum history lost probability")
    marginal = _reference_marginal(rho, dimension, reference_dimension)
    reference_error = float(np.linalg.norm(reference_output-marginal, 'fro'))
    record_relative = math.fsum(p*math.log(p/q) for p, q in
                                zip(probabilities, reference_probabilities, strict=True) if p)
    conditional_relative = math.fsum(conditional)
    output_relative = record_relative+conditional_relative
    loss = input_relative-output_relative
    if loss < -CHECK_TOLERANCE or reference_error > CHECK_TOLERANCE:
        raise ArithmeticError("native information ledger lost data processing or its reference")
    tape_entropy = _shannon(probabilities)
    output_entropy = tape_entropy+math.fsum(conditional_entropy)
    mean_attempts = math.fsum(attempts)
    coin_cost = mean_attempts*math.log(2)+math.log(float(beta))
    quantum_functional = output_entropy-coin_cost-math.log(reference_dimension)
    counted_functional = tape_entropy+math.fsum(counting)-coin_cost
    # This comparison uses the same faithful 66-state reference as the
    # existing population theorem; the protected code below has another one.
    native_cross = input_cross_entropy-math.log(reference_dimension)
    envelope_residual = quantum_functional-(native_cross-output_relative)
    if abs(envelope_residual) > CHECK_TOLERANCE:
        raise ArithmeticError("native partial-injection quantum entropy identity failed")
    return {"input_relative_entropy": input_relative,
            "record_relative_entropy": record_relative,
            "conditional_quantum_relative_entropy": conditional_relative,
            "output_relative_entropy": output_relative,
            "information_loss": loss,
            "input_entropy": input_entropy, "output_entropy": output_entropy,
            "record_entropy": tape_entropy, "mean_native_attempts": mean_attempts,
            "coin_and_flag_cost": coin_cost,
            "quantum_counted_functional": quantum_functional,
            "population_counted_functional": counted_functional,
            "charge_resolved": charge_resolved,
            "quantum_envelope_residual": envelope_residual,
            "reference_marginal_error": reference_error,
            "reference_dimension": reference_dimension,
            "success_probability_exact": str(beta),
            "terminal_record_count": len(probabilities),
            "entropy_values_are_outward_error_certificates": False,
            "record_probabilities": probabilities}


def code_stopped_weights(law="C", maximum=3):
    """Exact six-code native tree: stop on first nonzero current or at cap.

    Every event is corrected before continuation. Seam selection depends on
    the committed previous coin and value. The full tape includes the chosen
    seam, coin, raw label and terminal stopping decision; no timeout is dropped.
    """
    from . import record_code as code
    if (type(law) is not str or law not in ("C", "E")
            or type(maximum) is not int or not 1 <= maximum <= 4):
        raise ValueError("C/E and a one-through-four native-attempt cap required")
    live = [((), F(1))]
    terminal = []
    for depth in range(maximum):
        following = []
        for tape, weight in live:
            edge = code.SEAMS[0 if not tape else (7*tape[-1][2]+tape[-1][3]) % len(code.SEAMS)]
            for coin in (0, 1):
                for raw in range(4 if law == "C" else 8):
                    p = code.branch_weight(law, edge, coin, raw)
                    if not p:
                        continue  # Exactly zero on the declared code, retained by the source API.
                    next_tape = (*tape, (*edge, coin, raw))
                    item = (next_tape, weight*p)
                    current = {0: 0, 1: 1, 2: -1}[raw % 4]
                    (terminal if current or depth+1 == maximum else following).append(item)
        live = following
    if live or sum((p for _, p in terminal), F()) != 1:
        raise ArithmeticError("stopped code tree is incomplete")
    return tuple(terminal)


def protected_code_ledger(rho, *, law="C", maximum=3, reference_dimension=1):
    """Entropy of the exact actively corrected six-code history.

    The reference is I_6/6 tensor I_R/r ON THE CODE, not the faithful native
    population Gibbs state. Native branch isometries are checked exactly by
    record_code; this analytic transport evaluates their finite stopped tree.
    It does not claim to execute the compiled recovery pulse words here.
    """
    from . import record_code as code
    if type(reference_dimension) is not int or not 1 <= reference_dimension <= 6:
        raise ValueError("one through six finite reference levels required")
    if (type(law) is not str or law not in ("C", "E")
            or type(maximum) is not int or not 1 <= maximum <= 4):
        raise ValueError("C/E and a one-through-four native-attempt cap required")
    rho = _density(rho, 6*reference_dimension)
    # Derive and check the six-column Gram scalar on every actual source seam.
    for edge in code.SEAMS:
        for coin in (0, 1):
            for raw in range(4 if law == "C" else 8):
                images = code.branch_images(law, edge, coin, raw)
                counts = list(map(len, images))
                flat = [address for row in images for address in row]
                if (len(images) != 6 or len(set(counts)) != 1
                        or len(set(flat)) != len(flat)
                        or code.branch_weight(law, edge, coin, raw) != F(counts[0], 12)):
                    raise ArithmeticError("actual native code branch lost its scalar Gram")
    terminal = code_stopped_weights(law, maximum)
    weights = [p for _, p in terminal]
    mean = sum((p*len(tape) for tape, p in terminal), F())
    one_event = (math.log(12)-F(5, 6)*math.log(5) if law == "C"
                 else math.log(3)+F(2, 3)*math.log(2))
    tape_entropy = _shannon(weights)
    if abs(tape_entropy-float(mean)*one_event) > CHECK_TOLERANCE:
        raise ArithmeticError("bounded stopping lost the native tape entropy chain rule")
    input_entropy = _entropy(rho)
    conditional_relative = math.log(6*reference_dimension)-input_entropy
    reference_entropy = _entropy(_reference_marginal(rho, 6, reference_dimension))
    return {"logical_dimension": 6, "reference_dimension": reference_dimension,
            "reference": "maximally mixed on the code times maximally mixed reference",
            "input_entropy": input_entropy,
            "record_entropy": tape_entropy,
            "output_entropy": tape_entropy+input_entropy,
            "input_conditional_entropy": input_entropy-reference_entropy,
            "output_conditional_entropy": tape_entropy+input_entropy-reference_entropy,
            "input_relative_entropy": conditional_relative,
            "record_relative_entropy": 0.,
            "conditional_quantum_relative_entropy": conditional_relative,
            "information_loss": 0.,
            "one_attempt_record_entropy": one_event,
            "mean_attempts_exact": str(mean),
            "terminal_records": len(terminal),
            "coin_only_subtraction_excess": tape_entropy-float(mean)*math.log(2),
            "transport": "exact checked branch isometries and record-selected correction",
            "compiled_source_tree_executed_by_this_function": False}


def passive_code_information(law="C", word=((0, 1), (2, 3))):
    """Exact Bell/reference loss when the same six-code is left uncorrected.

    Compare the encoded maximally entangled input with I_6/6 tensor I_6/6.
    Both have the same local input and hence exactly the same complete tape.
    Each native branch remains a partial injection, so its logical Gram is
    diagonal with independently counted surviving supports. Its conditional
    entanglement entropy is the Shannon entropy of these six support counts.
    This evaluates the full joint quantum output, including its reference;
    the tape alone and BOTH local input marginals cannot detect the loss.
    """
    from . import record_code as code
    from ..observer import record_recovery as recovery
    if (type(law) is not str or law not in ("C", "E")
            or type(word) is not tuple or not 1 <= len(word) <= 4
            or any(type(edge) is not tuple or len(edge) != 2
                   or any(type(port) is not int for port in edge)
                   or edge not in code.SEAMS for edge in word)):
        raise ValueError("C/E and a native seam word of one through four events required")
    groups = defaultdict(lambda: [set() for _ in range(6)])
    for logical, support in enumerate(code.SUPPORTS):
        for initial in support:
            for coins in product((0, 1), repeat=len(word)):
                state, tape = initial, ()
                for edge, coin in zip(word, coins, strict=True):
                    state, record = recovery.recorded_step(law, state, edge, coin)
                    tape = (*tape, record)
                if state in groups[tape][logical]:
                    raise ArithmeticError("a native history ceased to be injective")
                groups[tape][logical].add(state)
    total, entropy = F(), 0.
    loss_form = defaultdict(F)
    population.v._log_add(loss_form, F(1), 6)
    unequal = 0
    rows = []
    for tape, supports in sorted(groups.items()):
        counts = tuple(map(len, supports))
        flat = [state for support in supports for state in support]
        if len(set(flat)) != len(flat):
            raise ArithmeticError("different logical columns share one native branch image")
        p = F(sum(counts), 36*2**len(word))
        total += p
        entropy += float(p)*_shannon(F(n, sum(counts)) for n in counts)
        for n in counts:
            if n:
                population.v._log_add(loss_form, p*F(n, sum(counts)), F(n, sum(counts)))
        unequal += len(set(counts)) != 1
        rows.append({"record": tape, "logical_support_counts": counts,
                     "probability_exact": str(p)})
    if total != 1:
        raise ArithmeticError("passive native code history is incomplete")
    loss = math.log(6)-entropy
    loss_form = population.v._clean(loss_form)
    denominator = math.lcm(*(coefficient.denominator for coefficient in loss_form.values()))
    powers = {int(prime): int(coefficient*denominator) for prime, coefficient in loss_form.items()}
    left = math.prod(prime**power for prime, power in powers.items() if power > 0)
    right = math.prod(prime**(-power) for prime, power in powers.items() if power < 0)
    return {"input_relative_entropy": 2*math.log(6),
            "record_relative_entropy": 0.,
            "conditional_quantum_relative_entropy": math.log(6)+entropy,
            "information_loss": loss,
            "information_loss_prime_logs_exact": {str(prime): str(value) for prime, value in loss_form.items()},
            "positive_information_loss_exact": left > right,
            "mean_conditional_entanglement_entropy": entropy,
            "non_scalar_branch_grams": unequal,
            "complete_branch_probability_exact": str(total),
            "branches": rows}


def field_reference_information(law="C"):
    """Same-protocol native/field information ledger, including both timeouts.

    Uses the source-executed field amplitudes. The comparison input is the
    product of its work and reference marginals. The two inputs therefore
    have identical probabilities for every work-only record, while their
    conditional joint outputs need not agree. Native code records and selector
    addresses are retained by their complete exact factor laws, as in the
    field runtime; they are not coherently averaged or success conditioned.
    """
    from . import record_field as field
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("C or E required")
    experiment, arrays = field.detector_experiment()
    initial = arrays["initial"]
    if initial.shape != (field.DIMENSION, 2):
        raise ArithmeticError("the same field/reference input is required")
    gram = initial.conj().T@initial
    if np.linalg.norm(gram-np.eye(2)/2, 'fro') > CHECK_TOLERANCE:
        raise ArithmeticError("the declared field input is not its two-level Bell state")
    accepted = F(15, 16)
    paths = [("first_selector_timeout", F(1, 16), initial)]
    for first in (0, 1):
        paths.append((f"second_selector_timeout_{first}", accepted/16,
                      arrays[f"after_first_{first}"]))
        for second in (0, 1):
            paths.append((f"terminal_{first}_{second}", accepted**2,
                          arrays[f"terminal_{first}_{second}"]))
    conditional, mass = 0., 0.
    marginal = np.zeros((2, 2), complex)
    rows = []
    for label, weight, state in paths:
        branch_gram = state.conj().T@state
        norm = float(np.trace(branch_gram).real)
        p = float(weight)*norm
        mass += p
        marginal += float(weight)*branch_gram
        entanglement = _entropy(branch_gram/norm) if norm > 0 else 0.
        relative = math.log(2)+entanglement if norm > 0 else 0.
        conditional += p*relative
        rows.append({"terminal": label, "selector_weight_exact": str(weight),
                     "probability": p,
                     "conditional_reference_entanglement_entropy": entanglement,
                     "conditional_joint_relative_entropy": relative})
    if abs(mass-1) > CHECK_TOLERANCE or np.linalg.norm(marginal-gram, 'fro') > CHECK_TOLERANCE:
        raise ArithmeticError("the complete field entropy audit lost a timeout or reference")
    loss = 2*math.log(2)-conditional
    expected = float(accepted)*math.log(2)
    if abs(loss-expected) > CHECK_TOLERANCE:
        raise ArithmeticError("destructive field read differs from its complete information loss")
    # The exact native branch laws supply the full extra record entropy.
    factors = field.record_factors(law)
    entropy_by_seam = [_shannon(F(row["weight_exact"]) for row in edge["local_record_factor"])
                       for edge in factors["accepted_seams"]]
    if max(entropy_by_seam)-min(entropy_by_seam) > CHECK_TOLERANCE:
        raise ArithmeticError("field native record entropy depends on the selected seam")
    mean_selectors = 1+accepted
    mean_repaired_cycles = accepted+accepted**2
    external_record_entropy = float(mean_selectors)*math.log(32)
    external_record_entropy += field.MODES*float(mean_repaired_cycles)*entropy_by_seam[0]
    return {"input_relative_entropy": 2*math.log(2),
            "record_relative_entropy": 0.,
            "conditional_quantum_relative_entropy": conditional,
            "information_loss": loss,
            "information_loss_exact": "(15/16)*log(2)",
            "first_read_probability_exact": str(accepted),
            "mean_selector_attempts_exact": str(mean_selectors),
            "mean_corrected_native_cycles_exact": str(mean_repaired_cycles),
            "selector_and_native_record_entropy": external_record_entropy,
            "complete_record_entropy": external_record_entropy+expected,
            "complete_joint_output_entropy": external_record_entropy+expected,
            "native_randomness_subtracted_quantum_entropy": expected,
            "native_recovery_information_loss": 0.,
            "field_read_information_loss": loss,
            "second_read_additional_information_loss": 0.,
            "complete_reference_marginal_error": float(np.linalg.norm(marginal-gram, 'fro')),
            "complete_quantum_terminal_blocks": rows,
            "input_reference": "product of actual initial work and untouched-reference marginals",
            "physical_heat_or_gravitational_entropy_identification": False}


def report():
    """Compact diagnostic receipt for the three explicitly distinct references."""
    sigma = np.diag(list(map(float, local_preparation.reference())))
    states = native.states(binary=True)

    def pair(left, right):
        vector = np.zeros(66, complex)
        for sites, amplitude in ((left, 1.), (right, 1j)):
            state = tuple(int(i in sites) for i in range(12))
            vector[states.index(state)] = amplitude/math.sqrt(2)
        return np.outer(vector, vector.conj())

    def compact(row):
        return {key: value for key, value in row.items() if key not in ('record_probabilities', 'branches')}

    bell = np.eye(6, dtype=complex).ravel()/math.sqrt(6)
    return {"schema": "oph.common_source.record_information.v1",
            "same_native_population_reference": {
                "complete_C_with_timeouts": compact(native_history_ledger(sigma, law='C')),
                "complete_E": compact(native_history_ledger(sigma)),
                "surviving_spectator_coherence": compact(native_history_ledger(pair((9, 10), (10, 11)))),
                "destroyed_core_coherence": compact(native_history_ledger(pair((0, 2), (1, 2))))},
            "protected_code_reference": {law: protected_code_ledger(
                np.outer(bell, bell.conj()), law=law, reference_dimension=6) for law in ('C', 'E')},
            "passive_two_event_code_control": {law: compact(passive_code_information(law))
                                               for law in ('C', 'E')},
            "same_field_experiment_reference": {law: field_reference_information(law)
                                                for law in ('C', 'E')},
            "scope": {"relative_entropy_chain_includes_conditional_quantum_outputs": True,
                      "native_coin_cost_varies_with_the_actual_stopping_law": True,
                      "native_timeout_quantum_outputs_retained": True,
                      "code_reference_equals_full_native_Gibbs_reference": False,
                      "field_reads_are_information_preserving": False,
                      "numerical_entropies_are_outward_error_certificates": False,
                      "heat_area_or_physical_stress_identified": False}}
