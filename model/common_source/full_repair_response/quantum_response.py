"""Complete finite C/E coherence diagnostics, without identifying damping as mass.

The actual finite completion selector supplies accepted seam weights and a
TERMINAL timeout. Native Kraus supports are reused unchanged. Signed-difference
blocks avoid a dense4356-dimensional eigensolve. The emitted path iterator
retains every native work record, timeout and detector outcome; selector-private
fine records have already been traced out by the existing selector contract.
No new native pulse compiler, charged CAR attachment or physical clock is added.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math

import numpy as np
from scipy.sparse import coo_matrix

from ...observer import native, quantum_mixing
from ...repair_scheduler import coins
from ...native_modular_reference import local_preparation as preparation
from ...native_clock import regional_preparation as incident


STATES = native.states(binary=True)
INDEX = {state: i for i, state in enumerate(STATES)}
DIMENSION = len(STATES)
TOLERANCE = 2e-10


def _law(law):
    if type(law) is not str or law not in ("C", "E"):
        raise ValueError("actual native C or E law required")
    return law


def _rounds(value):
    if type(value) is not int or not 1 <= value <= 32:
        raise ValueError("a finite selector round count1..32 required")
    return value


def _freeze(value):
    raw = np.asarray(value)
    return np.frombuffer(raw.tobytes(), dtype=raw.dtype).reshape(raw.shape)


def selector(selector_rounds=32):
    """Actual reduced selector; timeout remains a terminal outcome."""
    packet = preparation.selectors(_rounds(selector_rounds))
    weights = packet["completion_accepted_seam_weights"]
    timeout = F(packet["completion_timeout_probability"])
    if sum(weights.values(), F())+timeout != 1 or not 0 < timeout < 1:
        raise ArithmeticError("finite selector omitted a terminal outcome")
    if any(set(seam).intersection(preparation.CORE) for seam in weights):
        raise ArithmeticError("completion selector touched its protected core")
    return dict(weights), timeout, 1-timeout


@dataclass(frozen=True)
class Event:
    record: tuple[int, int, int, int]
    inputs: tuple[int, ...]
    outputs: tuple[int, ...]
    weight: F


@lru_cache(maxsize=8)
def _events(law, selector_rounds):
    _law(law)
    weights, _, _ = selector(selector_rounds)
    result = []
    for seam, weight in sorted(weights.items()):
        for coin, label, mapping in coins.coin_maps(law, seam):
            inputs, outputs = zip(*mapping, strict=True)
            result.append(Event((*seam, coin, label), tuple(inputs), tuple(outputs), weight/2))
    return tuple(result)


@lru_cache(maxsize=1)
def _sectors():
    result = {}
    for indices in preparation.difference_blocks():
        left, right = divmod(indices[0], DIMENSION)
        difference = tuple(a-b for a, b in zip(STATES[left], STATES[right], strict=True))
        result[difference] = tuple(indices)
    if sum(map(len, result.values())) != DIMENSION**2:
        raise ArithmeticError("the native matrix-unit partition is incomplete")
    return result


@lru_cache(maxsize=8)
def _kernel(law, selector_rounds):
    events = _events(_law(law), _rounds(selector_rounds))
    _, timeout, survival = selector(selector_rounds)
    columns = [defaultdict(F) for _ in range(DIMENSION**2)]
    effects = [F() for _ in STATES]
    block_of = {index: difference for difference, indices in _sectors().items() for index in indices}
    for event in events:
        mapping = tuple(zip(event.inputs, event.outputs, strict=True))
        for i, _ in mapping:
            effects[i] += event.weight
        for i, a in mapping:
            for j, b in mapping:
                row, column = DIMENSION*a+b, DIMENSION*i+j
                if block_of[row] != block_of[column]:
                    raise ArithmeticError("a committed Kraus branch moved signed charge difference")
                columns[column][row] += event.weight
    if effects != [survival]*DIMENSION:
        raise ArithmeticError("accepted selector survival is state dependent")
    for column, entries in enumerate(columns):
        for row, value in entries.items():
            if columns[row].get(column, F()) != value:
                raise ArithmeticError("the actual accepted channel is not Hilbert-Schmidt self-adjoint")
    rows, cols, values = [], [], []
    for column, entries in enumerate(columns):
        for row, value in entries.items():
            rows.append(row); cols.append(column); values.append(float(value))
    sparse = coo_matrix((values, (rows, cols)), shape=(DIMENSION**2, DIMENSION**2)).tocsr()
    sigma = preparation.reference()
    stationary = [F() for _ in STATES]
    for j, probability in enumerate(sigma):
        for row, value in columns[(DIMENSION+1)*j].items():
            a, b = divmod(row, DIMENSION)
            if a != b:
                raise ArithmeticError("native completion created coherence from the reference")
            stationary[a] += value*probability
    if stationary != [survival*value for value in sigma]:
        raise ArithmeticError("actual accepted channel does not preserve the specified reference")
    return tuple(tuple(sorted(entries.items())) for entries in columns), sparse, timeout, survival


def apply_channel(operator, law, *, selector_rounds=32, conditional=False):
    """Accepted A, or C=A/q, tensored with an untouched finite reference.

    C is explicitly the state-independent survivor-normalized diagnostic.
    Neither function recycles terminal timeouts. Use the path iterator to
    retain the complete stopped history rather than iterating A+tau*Id.
    """
    if type(conditional) is not bool or np.ma.isMaskedArray(operator):
        raise ValueError("plain finite operator and Boolean conditional flag required")
    raw = np.asarray(operator)
    if (raw.ndim != 2 or raw.shape[0] != raw.shape[1] or raw.shape[0] % DIMENSION
            or not DIMENSION <= len(raw) <= 4*DIMENSION or raw.dtype.kind not in "iufc"
            or not np.isfinite(raw).all()):
        raise ValueError("finite native operator with1..4 untouched reference levels required")
    _, sparse, _, survival = _kernel(_law(law), _rounds(selector_rounds))
    reference = len(raw)//DIMENSION
    columns = raw.reshape(DIMENSION, reference, DIMENSION, reference).transpose(0, 2, 1, 3)
    output = sparse@columns.reshape(DIMENSION**2, reference**2)
    if conditional:
        output /= float(survival)
    return output.reshape(DIMENSION, DIMENSION, reference, reference).transpose(0, 2, 1, 3).reshape(raw.shape)


def exact_block(law, difference, *, selector_rounds=32):
    """Exact Fraction block of the survivor-normalized channel C=A/q."""
    if (type(difference) is not tuple or len(difference) != 12
            or any(type(value) is not int for value in difference) or difference not in _sectors()):
        raise ValueError("an attained native signed-difference tuple required")
    columns, _, _, survival = _kernel(_law(law), _rounds(selector_rounds))
    indices = _sectors()[difference]
    dictionaries = {j: dict(columns[j]) for j in indices}
    return tuple(tuple(dictionaries[j].get(i, F())/survival for j in indices) for i in indices)


@dataclass(frozen=True)
class QuantumBlock:
    difference: tuple[int, ...]
    indices: tuple[int, ...]
    matrix: np.ndarray
    values: np.ndarray
    vectors: np.ndarray
    gershgorin_lower: F


@lru_cache(maxsize=8)
def reconstruction(law, selector_rounds=32):
    """All4356 matrix units, with eigensolves of size at most66."""
    _law(law); _rounds(selector_rounds)
    blocks, lower, fixed = [], F(1), 0
    digest = hashlib.sha256()
    for difference, indices in _sectors().items():
        exact = exact_block(law, difference, selector_rounds=selector_rounds)
        bound = min(row[i]-sum(abs(x) for j, x in enumerate(row) if i != j) for i, row in enumerate(exact))
        if bound <= 0:
            raise ArithmeticError("the actual quantum block lacks the claimed strict positive lower bound")
        lower = min(lower, bound)
        matrix = np.asarray(exact, float)
        values, vectors = np.linalg.eigh(matrix)
        if values[0] <= 0 or values[-1] > 1+TOLERANCE:
            raise ArithmeticError("native positive transfer reconstruction failed")
        fixed += int(np.count_nonzero(abs(values-1) < TOLERANCE))
        digest.update(json.dumps([difference, [[str(x) for x in row] for row in exact]],
                                 separators=(",", ":")).encode())
        blocks.append(QuantumBlock(difference, indices, _freeze(matrix), _freeze(values), _freeze(vectors), bound))
    if fixed != 6:
        raise ArithmeticError("protected-core fixed algebra no longer has dimension6")
    return tuple(blocks), lower, digest.hexdigest()


def sector_response(law, difference, *, selector_rounds=32):
    exact = exact_block(law, difference, selector_rounds=selector_rounds)
    values = np.linalg.eigvalsh(np.asarray(exact, float))
    return {"law": law, "signed_difference": list(difference), "dimension": len(exact),
            "conditional_transfer_eigenvalues": values.tolist(),
            "Euclidean_log_rates_per_surviving_work_proposal": [-math.log(x) for x in values if x < 1-TOLERANCE],
            "native_coherent_oscillation_phases": [0.]*len(values),
            "physical_mass_interpretation_established": False}


def modular_structure(selector_rounds=32):
    """Exact modular frequencies and commutation for every matrix unit."""
    sigma = preparation.reference()
    ratios = Counter(sigma[i]/sigma[j] for i in range(DIMENSION) for j in range(DIMENSION))
    checked = {}
    for law in ("C", "E"):
        columns, _, _, _ = _kernel(law, selector_rounds)
        count = 0
        for column, entries in enumerate(columns):
            i, j = divmod(column, DIMENSION)
            for row, _ in entries:
                a, b = divmod(row, DIMENSION)
                if sigma[a]/sigma[b] != sigma[i]/sigma[j]:
                    raise ArithmeticError("native completion does not commute with the same modular flow")
                count += 1
        checked[law] = count
    if set(ratios) != {F(1), F(25, 8), F(8, 25), F(8, 5), F(5, 8), F(5), F(1, 5)}:
        raise ArithmeticError("the full native modular frequency menu changed")
    return {"reference": "the same ideal protected-core local_preparation.reference()",
            "all_matrix_units": DIMENSION**2,
            "exact_weight_ratio_multiplicities": {str(x): n for x, n in sorted(ratios.items())},
            "modular_frequencies": ["0", "+/-log(25/8)", "+/-log(8/5)", "+/-log5"],
            "positive_frequency_max_over_min": math.log(5)/math.log(8/5),
            "every_accepted_channel_and_terminal_identity_commutes_with_modular_flow": True,
            "nonzero_operator_entries_checked": checked,
            "protected_core_algebra_modular_flow": "identity",
            "commutation_selects_modular_flow_as_particle_dynamics": False,
            "a_calibrated_clock_supplies_charged_field_transport": False}


def transport_obstruction():
    """Per-branch invariant, hence valid with classical feedback/memory too."""
    groups = 0
    for law in ("C", "E"):
        for seam in quantum_mixing.seams():
            for _, _, mapping in coins.coin_maps(law, seam):
                displacements = {tuple(b-a for a, b in zip(STATES[i], STATES[o], strict=True)) for i, o in mapping}
                if len(displacements) != 1:
                    raise ArithmeticError("a committed native branch changes the signed-difference support")
                groups += 1
    annihilation_addresses = 0
    for state in STATES:
        for port, occupied in enumerate(state):
            if occupied:
                target = list(state); target[port] -= 1
                if tuple(target) in INDEX:
                    raise ArithmeticError("bare annihilation unexpectedly acts inside total load2")
                annihilation_addresses += 1
    return {"all30_seams_and_both_laws_checked": True, "partial_injection_branches_checked": groups,
            "every_committed_branch_preserves_signed_difference": True,
            "positive_and_negative_difference_sites_propagate": False,
            "only_common_occupation_on_zero_difference_sites_can_move": True,
            "classical_record_memory_or_feedback_evades_this_invariant": False,
            "scope": "compositions of the actual committed C/E repair branches; extra coherent controls are different operations",
            "carrier_total_native_load": 2,
            "bare_odd_port_annihilation_restriction": "P_N=2*c_port*P_N=2=0",
            "annihilation_addresses_leaving_carrier": annihilation_addresses,
            "full12_mode_CAR_irreducible_dimension": 2**12,
            "carrier_dimension": DIMENSION,
            "full12_mode_unital_CAR_representation_on_this_carrier_possible": DIMENSION % (2**12) == 0,
            "relative_core_charge_coherence_is_already_an_odd_or_SM_charged_field": False,
            "additional_encoded90_mode_CAR_field_identified_with_native66_carrier": False,
            "larger_network_or_composite_field_constructions_excluded": False}


def _probe(outside_port=5):
    if type(outside_port) is not int or not 0 <= outside_port < 12 or outside_port in preparation.CORE:
        raise ValueError("one port outside the protected core required")
    first, other = preparation.CORE
    address = lambda ports: INDEX[tuple(int(p in ports) for p in range(12))]
    return address((first, other)), address((other, outside_port)), address((first, outside_port))


def coherence_probe(law, *, selector_rounds=32, outside_port=5):
    x, y, _ = _probe(outside_port)
    columns, _, timeout, survival = _kernel(_law(law), _rounds(selector_rounds))
    weights, _, _ = selector(selector_rounds)
    degree = sum((w for seam, w in weights.items() if outside_port in seam), F())/survival
    eigenvalue = 1-degree/(2 if law == "C" else 1)
    column = DIMENSION*x+y
    if dict(columns[column]) != {column: survival*eigenvalue}:
        raise ArithmeticError("the native charged-relative coherence is not the predicted scalar eigenmode")
    return {"law": law, "native_matrix_unit": [x, y], "outside_port": outside_port,
            "relative_core_charge_change": 1, "total_native_charge_change": 0,
            "operator_is_even_global_transfer_not_odd_CAR": True,
            "conditional_native_eigenvalue_exact": str(eigenvalue),
            "accepted_eigenvalue_exact": str(survival*eigenvalue),
            "terminal_timeout_probability_exact": str(timeout),
            "decay_log_per_surviving_work_proposal": -math.log(float(eigenvalue)),
            "native_coherent_frequency_per_work_proposal": 0.,
            "modular_phase_weight_ratio_exact": str(preparation.reference()[x]/preparation.reference()[y]),
            "if_separately_composed_with_modular_flow": "response(n,s)=lambda^n*exp(+i*log(8/5)*s) for U_s=sigma^(i*s)",
            "modular_flow_composition_is_selected_native_continuation": False}


@lru_cache(maxsize=1)
def seam_positivity_certificate():
    """Replay the existing exact local-channel polynomial/source checks."""
    certificate = incident.exact_audit(preparation.CORE)
    if any(not row["exact_unital_trace_preserving_positive_self_adjoint"]
           for row in certificate["laws"].values()):
        raise ArithmeticError("a native local channel lost its positive transfer identity")
    return certificate


def incident_clock_response(law, *, selector_rounds=32, outside_port=5, deadline=3):
    """Exact coherence of the actual finite stopped incident clock.

    Accepted J is incident on CAPS[0]; R is accepted but remote. A selector
    timeout ends immediately; deadline R**h also ends. Both remain present.
    This retains the source remote-count memory rather than pretending that
    one observer tick equals one global proposal.
    """
    if type(deadline) is not int or not 1 <= deadline <= 64:
        raise ValueError("a finite incident-clock deadline1..64 required")
    x, y, _ = _probe(outside_port)
    weights, timeout, _ = selector(selector_rounds)
    cap = set(preparation.CAPS[0])
    p_j = sum((w for seam, w in weights.items() if cap.intersection(seam)), F())
    p_r = sum(weights.values(), F())-p_j
    j = r = F()
    for event in _events(_law(law), _rounds(selector_rounds)):
        mapping = dict(zip(event.inputs, event.outputs, strict=True))
        if x in mapping and y in mapping:
            if (mapping[x], mapping[y]) != (x, y):
                raise ArithmeticError("selected clock coherence unexpectedly moved")
            if cap.intersection(event.record[:2]):
                j += event.weight
            else:
                r += event.weight
    if not 0 <= j <= p_j or not 0 <= r <= p_r or not p_j+p_r+timeout == 1:
        raise ArithmeticError("clocked coherence or probability decomposition failed")
    successes, selector_timeouts = [], []
    for remote in range(deadline):
        successes.append({"remote_proposals": remote, "probability_exact": str(p_j*p_r**remote),
                          "coherence_multiplier_exact": str(j*r**remote)})
        selector_timeouts.append({"remote_proposals": remote, "probability_exact": str(timeout*p_r**remote),
                                 "coherence_multiplier_exact": str(timeout*r**remote)})
    waiting = sum((p_r**k for k in range(deadline)), F())
    quantum_waiting = sum((r**k for k in range(deadline)), F())
    success_probability = p_j*waiting
    timeout_probability = timeout*waiting+p_r**deadline
    if success_probability+timeout_probability != 1:
        raise ArithmeticError("stopped observer clock omitted timeout mass")
    return {"law": law, "clock_cap": list(preparation.CAPS[0]), "outside_port": outside_port,
            "deadline_work_proposals": deadline,
            "incident_probability_per_proposal_exact": str(p_j),
            "remote_probability_per_proposal_exact": str(p_r),
            "incident_coherence_multiplier_exact": str(j), "remote_coherence_multiplier_exact": str(r),
            "accepted_tick_by_remote_count": successes,
            "selector_timeout_by_remote_count": selector_timeouts,
            "deadline_timeout": {"remote_proposals": deadline, "probability_exact": str(p_r**deadline),
                                 "coherence_multiplier_exact": str(r**deadline)},
            "accepted_clock_tick_probability_exact": str(success_probability),
            "total_terminal_timeout_probability_exact": str(timeout_probability),
            "accepted_tick_coherence_multiplier_exact": str(j*quantum_waiting),
            "conditional_tick_coherence_eigenvalue_exact": str(j*quantum_waiting/success_probability),
            "new_autonomous_oscillation_frequency": 0.,
            "finite_clock_map": "S_h=J*f_h(R), f_h(R)=sum_(k=0)^(h-1)R^k; D_h=tau*f_h(R)+R^h",
            "positive_spectrum_proof": "Each actual seam channel is HS-positive. Thus J,R>=0, norm(R)<=p_R<1 and f_h(R)>0. S_h is similar to f_h(R)^(1/2)*J*f_h(R)^(1/2)>=0, so it has only real nonnegative eigenvalues.",
            "source_clock_reconstruction_has_Lorentzian_particle_poles": False,
            "terminal_flags_or_remote_count_memory_discarded_as_success": False}


def complex_pole_decomposition(law, *, selector_rounds=32, outside_port=5):
    """Separate attenuation from the phase of the SAME reference's flow.

    This is an explicit conditional composition: apply one unit s=1 of
    U_s=sigma^(i*s) after each chosen proposal/tick map. It is not a proof
    that committed native repair chooses that coherent continuation.
    """
    probe = coherence_probe(law, selector_rounds=selector_rounds, outside_port=outside_port)
    clock = incident_clock_response(law, selector_rounds=selector_rounds, outside_port=outside_port)
    ratio = F(probe["modular_phase_weight_ratio_exact"])
    omega = math.log(float(ratio))
    rows = []
    for name, value in (("conditional_work_proposal", probe["conditional_native_eigenvalue_exact"]),
                        ("conditional_stopped_incident_tick", clock["conditional_tick_coherence_eigenvalue_exact"])):
        attenuation = F(value)
        pole = float(attenuation)*np.exp(1j*omega)
        rows.append({"response_parameter": name, "attenuation_exact": str(attenuation),
                     "complex_multiplier": [float(pole.real), float(pole.imag)],
                     "damping_minus_log_modulus": -math.log(float(attenuation)),
                     "phase_per_unit_modular_step": float(np.angle(pole)),
                     "modular_frequency": omega,
                     "damping_is_identified_with_frequency_or_rest_mass": False})
    return {"law": law, "outside_port": outside_port, "modular_step": 1,
            "same_reference_ratio_exact": str(ratio), "complex_response_formula": "z=lambda*exp(i*omega*s)",
            "pole_decomposition": "-log|z|=-log(lambda); arg(z)=omega*s modulo2*pi",
            "rows": rows, "native_uncomposed_oscillation_frequency": 0.,
            "phase_does_not_come_from_the_real_native_decay_rate": True,
            "one_modular_advance_per_chosen_tick_is_a_composition_prescription": True,
            "physical_particle_clock_or_selected_continuation_established": False}


def joint_modular_spectrum(law, *, selector_rounds=32):
    """All4356 modes split into exact sigma-frequency classes before eigensolve."""
    sigma = preparation.reference()
    blocks, _, _ = reconstruction(law, selector_rounds)
    groups = defaultdict(list)
    for block in blocks:
        positions = defaultdict(list)
        for local, index in enumerate(block.indices):
            i, j = divmod(index, DIMENSION)
            positions[sigma[i]/sigma[j]].append(local)
        for ratio, local in positions.items():
            submatrix = block.matrix[np.ix_(local, local)]
            groups[ratio].extend(np.linalg.eigvalsh(submatrix).tolist())
    rows = []
    for ratio, values in sorted(groups.items()):
        rows.append({"weight_ratio_exact": str(ratio), "operator_modes": len(values),
                     "modular_frequency_log_ratio": math.log(float(ratio)),
                     "native_attenuation_minimum": min(values), "native_attenuation_maximum": max(values),
                     "joint_multiplier_if_modular_flow_composed": "lambda*exp(i*s*log(weight_ratio))"})
    if sum(row["operator_modes"] for row in rows) != DIMENSION**2:
        raise ArithmeticError("joint modular/native response discarded operator modes")
    return {"rows": rows, "source_operator_modes_retained": DIMENSION**2,
            "native_map_and_modular_flow_simultaneously_resolved": True,
            "C_E_or_clock_normalization_can_change_attenuation_without_new_sigma_frequencies": True,
            "same_sigma_frequency_menu_is_a_Standard_Model_mass_spectrum": False}


def _positive_ldl(matrix):
    """Exact rational positive-pivot LDL certificate, without symbolic rank."""
    size = len(matrix)
    lower = [[F(i == j) for j in range(size)] for i in range(size)]
    pivots = []
    for i in range(size):
        pivot = matrix[i][i]-sum((lower[i][k]**2*pivots[k] for k in range(i)), F())
        if pivot <= 0:
            raise ArithmeticError("the finite rational gap certificate has a nonpositive pivot")
        pivots.append(pivot)
        for j in range(i+1, size):
            lower[j][i] = (matrix[j][i]-sum((lower[j][k]*lower[i][k]*pivots[k]
                                            for k in range(i)), F()))/pivot
    if any(matrix[i][j] != sum((lower[i][k]*pivots[k]*lower[j][k]
                                for k in range(min(i, j)+1)), F())
           for i in range(size) for j in range(size)):
        raise ArithmeticError("exact rational LDL reconstruction failed")
    digest = hashlib.sha256(json.dumps([[str(x) for x in row] for row in lower],
                                       separators=(",", ":")).encode()).hexdigest()
    return {"dimension": size, "positive_pivots_exact": list(map(str, pivots)),
            "lower_triangular_factor_sha256": digest, "exact_LDL_reconstruction_checked": True}


@lru_cache(maxsize=2)
def operator_log_spread_certificate(selector_rounds=32):
    """Finite exact instance bound for conditional proposal maps, not masses.

    Six fixed-core matrix sectors have dimensions45,10,10,10,10,1. For each
    nonconstant sector use Q=(e_i-e_last) to certify Q*(I-C-gI)Q>0, g=1/100.
    All other4270 matrix units have strict row-sum loss>g. This proves the
    nonfixed gap without importing an exclusion-gap theorem or diagonalizing
    the full operator space.
    """
    _rounds(selector_rounds)
    core, outside = preparation.CORE, tuple(p for p in range(12) if p not in preparation.CORE)
    core_sectors = defaultdict(list)
    for i, x in enumerate(STATES):
        for j, y in enumerate(STATES):
            if all(x[p] == y[p] for p in outside):
                key = tuple(x[p] for p in core), tuple(y[p] for p in core)
                core_sectors[key].append(DIMENSION*i+j)
    if sorted(map(len, core_sectors.values())) != [1, 10, 10, 10, 10, 45]:
        raise ArithmeticError("the exact protected-core sector decomposition changed")
    columns_c, _, _, survival = _kernel("C", selector_rounds)
    columns_e, _, _, _ = _kernel("E", selector_rounds)
    dictionaries = {law: [dict(column) for column in columns]
                    for law, columns in (("C", columns_c), ("E", columns_e))}
    gap, covered, certificates = F(1, 100), set(), []
    for key, indices in sorted(core_sectors.items()):
        selected = set(indices); covered.update(indices)
        matrices = {}
        for law in ("C", "E"):
            if any(set(dictionaries[law][index])-selected for index in indices):
                raise ArithmeticError("native repair leaks out of a fixed-core matrix sector")
            matrices[law] = tuple(tuple(dictionaries[law][column].get(row, F())/survival
                                        for column in indices) for row in indices)
        if matrices["C"] != matrices["E"] or any(sum(row, F()) != 1 for row in matrices["C"]):
            raise ArithmeticError("the protected matrix sectors lost their common stochastic action")
        size, transfer = len(indices), matrices["C"]
        generator = [[F(i == j)-transfer[i][j] for j in range(size)] for i in range(size)]
        # Q spans the zero-sum subspace; Q*Q has2 on diagonal and1 elsewhere.
        last = size-1
        compressed = [[generator[i][j]-generator[i][last]-generator[last][j]+generator[last][last]
                       -gap*(1+int(i == j)) for j in range(last)] for i in range(last)]
        certificate = _positive_ldl(compressed)
        certificate.update({"core_ket_bra": [list(key[0]), list(key[1])], "source_sector_dimension": size,
                            "certified_form": "Q^T[(I-C)-(1/100)I]Q, Q_i=e_i-e_last"})
        certificates.append(certificate)
    if len(covered) != 86:
        raise ArithmeticError("the fixed-core certificate lost matrix units")
    others = set(range(DIMENSION**2))-covered
    rows = {}
    for law in ("C", "E"):
        # Self-adjointness makes column sums equal the corresponding row sums.
        maximum = max(sum(dictionaries[law][index].values(), F())/survival for index in others)
        if maximum >= 1-gap:
            raise ArithmeticError("a remaining quantum block violates its exact killing bound")
        _, lower, _ = reconstruction(law, selector_rounds)
        if lower <= F(7, 50):
            raise ArithmeticError("the source positivity floor no longer proves the logarithmic upper bound")
        rows[law] = {"other_matrix_units": len(others), "maximum_other_row_sum_exact": str(maximum),
                     "positive_transfer_floor_exact": str(lower), "floor_exceeds_7_over50_exact": True}
    exp_two_lower = sum((F(2)**k/math.factorial(k) for k in range(6)), F())
    if not exp_two_lower > F(50, 7):
        raise ArithmeticError("the exact exponential comparison failed")
    return {"selector_rounds": selector_rounds, "six_fixed_sectors": certificates,
            "fixed_operator_dimension_exact": 6, "other_sectors": rows,
            "nonfixed_conditional_transfer_upper_bound_exact": "99/100",
            "nonzero_log_rate_lower_strict": "1/100", "log_rate_upper_strict": "2",
            "log_rate_spread_upper_strict": "200",
            "upper_log_proof": "lambda>7/50; exp(2)>sum_(k=0)^5 2^k/k!=109/15>50/7; hence-log(lambda)<2",
            "lower_log_proof": "nonfixed lambda<=99/100 implies-log(lambda)>=-log(99/100)>1/100",
            "exact_exp_two_Taylor_lower": str(exp_two_lower),
            "imported_exclusion_gap_theorem_required": False,
            "existing_1_over600_source_bound_modified": False,
            "scope": "all nonzero HS operator-log rates of C=A/q per surviving completion work proposal, both native C/E laws",
            "bound_claimed_for_all_incident_clock_or_coherent_arc_spectra": False,
            "these_operator_logs_are_physical_masses": False}


@lru_cache(maxsize=8)
def _incident_sparse(law, selector_rounds):
    cap = set(preparation.CAPS[0])
    rows, cols, values = ([], []), ([], []), ([], [])
    for event in _events(law, selector_rounds):
        kind = int(not cap.intersection(event.record[:2]))  #0=incident,1=remote
        for i, a in zip(event.inputs, event.outputs, strict=True):
            for j, b in zip(event.inputs, event.outputs, strict=True):
                rows[kind].append(DIMENSION*a+b); cols[kind].append(DIMENSION*i+j)
                values[kind].append(float(event.weight))
    return tuple(coo_matrix((values[k], (rows[k], cols[k])), shape=(DIMENSION**2, DIMENSION**2)).tocsr()
                 for k in (0, 1))


def full_incident_clock_spectrum(law, *, selector_rounds=32):
    """All quantum modes of the geometric remote-wait stopped-clock map.

    This is S=J(I-R)^-1 and terminalD=tau(I-R)^-1, WITHOUT a finite remote
    deadline. The selector itself is still finite; its timeout is terminal.
    It is distinct from S_h in incident_clock_response and the complete short
    histories. No terminal timeout is restarted as an ordinary repair event.
    """
    _law(law); _rounds(selector_rounds)
    weights, timeout, _ = selector(selector_rounds)
    p_j = sum((weight for seam, weight in weights.items() if set(preparation.CAPS[0]).intersection(seam)), F())
    p_r = 1-timeout-p_j
    success = p_j/(p_j+timeout)
    incident, remote = _incident_sparse(law, selector_rounds)
    eigenvalues, guaranteed_zero, maximum_residual = [], 0, 0.
    for indices in _sectors().values():
        selected = list(indices)
        j = incident[selected, :][:, selected].toarray()
        r = remote[selected, :][:, selected].toarray()
        size = len(indices)
        m = np.eye(size)-r
        values, vectors = np.linalg.eigh(m)
        if values[0] < float(p_j+timeout)-TOLERANCE:
            raise ArithmeticError("remote-wait resolvent lost its source positive lower bound")
        half = (vectors*np.sqrt(values))@vectors.T
        invhalf = (vectors/np.sqrt(values))@vectors.T
        b = invhalf@j@invhalf/float(success)
        stopped = np.linalg.solve(m.T, j.T).T/float(success)
        maximum_residual = max(maximum_residual, float(np.max(abs(invhalf@stopped@half-b))))
        result = np.linalg.eigvalsh(b)
        if result[0] < -TOLERANCE or result[-1] > 1+TOLERANCE:
            raise ArithmeticError("positive stopped-clock similarity has invalid spectrum")
        eigenvalues.extend(result.tolist())
        if not np.any(j):
            guaranteed_zero += size  # Exact support zero: positive terms cannot cancel.
    values = np.array(eigenvalues)
    positive_nonfixed = values[(values > TOLERANCE)&(values < 1-TOLERANCE)]
    rates = -np.log(positive_nonfixed)
    return {"law": law, "operator_modes_retained": len(values), "largest_block_inverse": 66,
            "geometric_remote_wait_resolvent": "S=J(I-R)^-1; D=tau(I-R)^-1",
            "remote_wait_has_finite_deadline": False,
            "finite_selector_timeout_is_terminal": True,
            "conditional_tick_success_probability_exact": str(success),
            "remote_probability_per_proposal_exact": str(p_r),
            "finite_deadline_complete_history_differs_with_probability_at_most": "p_R^h",
            "similar_positive_operator": "(I-R)^(-1/2) J (I-R)^(-1/2) divided by tick-success probability",
            "positive_similarity_residual": maximum_residual,
            "fixed_modes_at_numerical_tolerance": int(np.count_nonzero(abs(values-1) < TOLERANCE)),
            "zero_transfer_modes_at_numerical_tolerance": int(np.count_nonzero(abs(values) < TOLERANCE)),
            "exact_zero_modes_from_identically_zero_incident_blocks": guaranteed_zero,
            "zero_transfer_log": "not finite; these erased modes are retained, not omitted as particles",
            "minimum_resolved_positive_transfer_eigenvalue": float(positive_nonfixed.min()),
            "maximum_resolved_nonfixed_transfer_eigenvalue": float(positive_nonfixed.max()),
            "resolved_finite_log_rate_spread_diagnostic": float(rates.max()/rates.min()),
            "native_coherent_oscillation_phases": "zero for every positive eigenvalue; zero eigenvalues erase the mode",
            "floating_spectrum_outward_certified": False,
            "Lorentzian_particle_mass_identification_established": False}


@dataclass(frozen=True)
class QuantumPath:
    records: tuple[tuple, ...]
    terminal_flag: str
    accepted_proposals: int
    pointer: int
    state: np.ndarray


def _initial_probe(outside_port):
    x, y, z = _probe(outside_port)
    state = np.zeros((DIMENSION, 2), complex)
    state[x, 0] = state[y, 0] = .5
    state[z, 1] = 1/math.sqrt(2)
    return state


def _endpoint(state, outside_port):
    x, y, _ = _probe(outside_port)
    swapped = np.zeros_like(state); swapped[x], swapped[y] = state[y], state[x]
    support = np.zeros_like(state); support[x], support[y] = state[x], state[y]
    common = state/math.sqrt(2)+(.5-1/math.sqrt(2))*support
    return common+swapped/2, common-swapped/2


def iter_history_outputs(law, *, steps=2, selector_rounds=32, outside_port=5):
    """Every output of a short stopped CQ history, including zeros/timeouts.

    State arrays are native-work rows by two UNTOUCHED quantum reference
    columns. The supplied positive binary endpoint effect is(I+S_xy)/2;
    its reference factor is identity. All actual native work records survive.
    This executes native Kraus semantics, not new M6 pulse words.
    """
    if type(steps) is not int or not 0 <= steps <= 2:
        raise ValueError("complete tree diagnostic supports0..2 work proposals")
    events = _events(_law(law), _rounds(selector_rounds))
    _, timeout, _ = selector(selector_rounds)
    initial = _initial_probe(outside_port)
    def terminal(state, records, flag, accepted):
        for pointer, output in enumerate(_endpoint(state, outside_port)):
            yield QuantumPath(records, flag, accepted, pointer, output)
    def visit(state, records, depth):
        if depth == steps:
            yield from terminal(state, records, "accepted", depth)
            return
        # A timeout terminates. It is never reinserted into a later update.
        yield from terminal(math.sqrt(float(timeout))*state,
                            (*records, ("selector_timeout",)), "selector_timeout", depth)
        for event in events:
            output = np.zeros_like(state)
            output[list(event.outputs)] = math.sqrt(float(event.weight))*state[list(event.inputs)]
            yield from visit(output, (*records, event.record), depth+1)
    yield from visit(initial, (), 0)


def history_diagnostic(law, *, steps=2, selector_rounds=32, outside_port=5):
    probe = coherence_probe(law, selector_rounds=selector_rounds, outside_port=outside_port)
    _, timeout, survival = selector(selector_rounds)
    eigenvalue = F(probe["conditional_native_eigenvalue_exact"])
    counts, probabilities = Counter(), defaultdict(list)
    reference = np.zeros((2, 2), complex)
    digest = hashlib.sha256()
    zero = 0
    for path in iter_history_outputs(law, steps=steps, selector_rounds=selector_rounds, outside_port=outside_port):
        probability = float(np.vdot(path.state, path.state).real)
        counts[path.terminal_flag] += 1
        probabilities[path.terminal_flag, path.pointer].append(probability)
        reference += path.state.conj().T@path.state
        zero += probability == 0
        # Only exact record grammar belongs in an exact replay digest. Floating
        # states are recomputed and compared through diagnostics below; their
        # raw bytes depend on roundoff/platform and are not state-byte custody.
        digest.update(json.dumps([path.records, path.terminal_flag, path.accepted_proposals, path.pointer],
                                 separators=(",", ":")).encode("utf-8"))
    masses = {flag: math.fsum([*probabilities[flag, 0], *probabilities[flag, 1]]) for flag in counts}
    plus = math.fsum(values for (flag, pointer), rows in probabilities.items() if pointer == 0 for values in rows)
    minus = math.fsum(values for (flag, pointer), rows in probabilities.items() if pointer == 1 for values in rows)
    response = (survival*eigenvalue)**steps+timeout*sum(((survival*eigenvalue)**k for k in range(steps)), F())
    expected_mean = float(response)/2  # Initial native S_xy mean is exactly1/2.
    expected_masses = {"accepted": float(survival**steps), "selector_timeout": float(1-survival**steps)}
    residual = max(abs(plus+minus-1), abs((plus-minus)-expected_mean),
                   float(np.max(abs(reference-np.eye(2)/2))),
                   *(abs(masses.get(flag, 0)-value) for flag, value in expected_masses.items()))
    if residual > TOLERANCE:
        raise ArithmeticError("complete native coherence history lost weight, reference or response")
    return {"work_proposals_upper": steps, "selector_rounds": selector_rounds,
            "complete_endpoint_paths": sum(counts.values()), "paths_by_terminal_flag": dict(counts),
            "zero_probability_paths_retained": zero, "probability_by_terminal_flag": masses,
            "accepted_probability_exact": str(survival**steps),
            "timeout_probability_exact": str(1-survival**steps),
            "pointer_probabilities": [plus, minus], "initial_relative_charge_coherence_mean": .5,
            "expected_complete_pointer_mean": expected_mean,
            "complete_response_and_reference_residual": residual,
            "endpoint_state_shape": [DIMENSION, 2], "endpoint_record_stream_sha256": digest.hexdigest(),
            "endpoint_record_digest_fields": ["records", "terminal_flag", "accepted_proposals", "pointer"],
            "floating_endpoint_states_recomputed_not_byte_hashed": True,
            "floating_endpoint_checks": "aggregate branch probability, pointer response and untouched-reference diagnostics; no endpoint state-byte custody asserted",
            "work_resets": 0, "postselection": False, "terminal_timeouts_recycled": False,
            "all_native_work_records_and_both_detector_outcomes_retained": True,
            "selector_private_fine_tape_in_comparison": False,
            "source_execution_level": "actual native partial-injection CP semantics, not newly compiled M6 words",
            "endpoint_coherent_preparation_and_readout_are_supplied": True,
            "physical_clock_or_mass_identification": False}


def report(*, selector_rounds=32):
    weights, timeout, survival = selector(selector_rounds)
    laws = {}
    for law in ("C", "E"):
        blocks, lower, digest = reconstruction(law, selector_rounds)
        values = np.concatenate([block.values for block in blocks])
        nonfixed = values[values < 1-TOLERANCE]
        rates = -np.log(nonfixed)
        laws[law] = {"all_matrix_units": DIMENSION**2,
                     "block_dimensions": dict(sorted(Counter(len(block.indices) for block in blocks).items())),
                     "exact_conditional_block_sha256": digest,
                     "exact_positive_operator_lower_bound": str(lower), "operator_upper_bound": "1",
                     "protected_fixed_algebra_dimension": 6,
                     "minimum_conditional_transfer_eigenvalue": float(values.min()),
                     "maximum_nonfixed_conditional_transfer_eigenvalue": float(nonfixed.max()),
                     "minimum_nonzero_Euclidean_log_rate": float(rates.min()),
                     "maximum_Euclidean_log_rate": float(rates.max()),
                     "Euclidean_log_rate_spread_diagnostic": float(rates.max()/rates.min()),
                     "all_native_transfer_eigenvalues_real_and_positive": True,
                     "native_unitary_carrier_dynamics_established": False,
                     "operator_Hilbert_space_log": "H_HS=-log(C); C^n=exp(-n*H_HS)",
                     "accepted_map_log": "-log(A)=-log(q)*I_HS-log(C); the common shift is selector survival, not a species mass",
                     "common_terminal_survival_log_shift": -math.log1p(-float(timeout)),
                     "log_parameter": "surviving discrete work-proposal count, not a calibrated Lorentzian clock",
                     "reflection_positive_two_point_bond_form": "sum_ij <O_i,C^(n_i+n_j+1)O_j>=<V,CV>>=0",
                     "full_OS_or_local_Lorentzian_field_reconstruction_established": False,
                     "fractional_transfer_powers_assumed_CP": False,
                     "relative_charge_probe": coherence_probe(law, selector_rounds=selector_rounds),
                     "finite_incident_clock_coherence": [incident_clock_response(law, selector_rounds=selector_rounds,
                                                                                    outside_port=port) for port in (5, 3)],
                     "complex_native_plus_modular_response": [complex_pole_decomposition(law, selector_rounds=selector_rounds,
                                                                                         outside_port=port) for port in (5, 3)],
                     "full_joint_modular_spectrum": joint_modular_spectrum(law, selector_rounds=selector_rounds),
                     "full_geometric_remote_wait_quantum_clock_spectrum": full_incident_clock_spectrum(law, selector_rounds=selector_rounds),
                     "complete_quantum_history": history_diagnostic(law, selector_rounds=selector_rounds)}
    return {"schema": "oph.full_repair_response.quantum_response.v1",
            "native_dimension": DIMENSION, "selector_rounds": selector_rounds,
            "accepted_selector_probability_exact": str(survival),
            "terminal_selector_timeout_exact": str(timeout),
            "accepted_seam_weights_exact": [[list(seam), str(value)] for seam, value in sorted(weights.items())],
            "conditional_channel_is_explicit_survivor_diagnostic": True,
            "actual_stopped_history_never_conditions_away_timeouts": True,
            "same_ideal_reference_target_used": True, "finite_reference_preparation_replayed": False,
            "laws": laws, "modular_structure": modular_structure(selector_rounds),
            "exact_all_seam_positive_transfer_certificate": seam_positivity_certificate(),
            "charged_field_identification_boundary": transport_obstruction(),
            "exact_conditional_operator_log_spread_certificate": operator_log_spread_certificate(selector_rounds),
            "coherent_alternatives": {
                "existing_weighted_arc_walk": "separate coherent_arc law; coherent event reflection/reversal is not the committed C/E channel",
                "existing_modular_clock": "native modular source flow supplies phases of the same reference, but does not select a particle kinetic operator",
                "damping_alone_is_a_particle_mass": False,
            },
            "empirical_mass_inputs": [], "physical_particles_or_mass_ratios_derived": False,
            "floating_spectral_and_history_checks_outward_certified": False}
