"""One continuing charged interferometer, propagation probe and quantum clock.

The native-source backend contracts actual source words. The continuous
backend is an independent target-matrix comparison. All work/reference outputs
remain unnormalized. Neither backend prepares work again between its reads.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import math

import numpy as np

from . import history, model
from ...observer.source_scalar import compiler


QUBITS = 10  # eight fermion work bits, two clock bits
HADAMARD = np.array([[1., 1.], [1., -1.]]) / math.sqrt(2)
WEAK_STRENGTH = F(1, 9)  # declared apparatus contrast; no target-mass dependence


@dataclass
class Result:
    branches: dict[tuple[int, ...], np.ndarray]
    aborts: dict[tuple, np.ndarray]
    prefix_probabilities: dict[tuple[int, ...], float]
    resources: dict
    specification: dict

    def probabilities(self):
        return {key: norm(value) for key, value in self.branches.items()}


def norm(state):
    return float(np.vdot(state, state).real)


def _resource_add(left, right):
    for key in ("native_pulses", "proper_code_service_events", "fresh_pointer_blanks", "committed_classical_bits"):
        left[key] = left.get(key, 0) + int(right.get(key, 0))


def _source_preparation(interaction):
    vector, metadata = model.charge_transfer_preparation(interaction=interaction)
    p0, p1 = map(F, model.ideal_clock_calibration()["populations"])
    word = history.preparation_word(vector)
    for clock, purifier in ((8, 10), (9, 11)):
        word += compiler.remap(history.preparation_word(np.sqrt([float(p0), float(p1)])), [clock])
        word += compiler.cnot(clock, purifier)
        word += [compiler.pulse("H", math.pi/2, [clock])]
    state = np.zeros((1 << 12, 1), complex)
    state[0, 0] = 1
    prepared, preparation_resources = history.execute_word(word, 12, state)
    # High purifier bits become untouched reference columns, never detector data.
    prepared = prepared[:, 0].reshape(4, 1 << QUBITS).T
    # Only expose counters actually aggregated over the whole history.
    # Detailed compiler counters remain explicitly preparation-only.
    resources = {
        "preparation_source_word": preparation_resources,
        "native_pulses": preparation_resources["native_pulses"],
        "proper_code_service_events": preparation_resources["proper_code_service_events"],
        "fresh_pointer_blanks": 0, "committed_classical_bits": 0,
        "initial_blank_logical_qubits": 12, "work_preparations": 1,
    }
    resources["clock_purification_qubits"] = 2
    resources["work_and_clock_qubits"] = QUBITS
    return prepared, vector, resources, metadata


def _ideal_preparation(interaction):
    vector, metadata = model.charge_transfer_preparation(interaction=interaction)
    p0, p1 = map(F, model.ideal_clock_calibration()["populations"])
    clock = HADAMARD @ np.diag(np.sqrt([float(p0), float(p1)]))
    state = np.zeros((1024, 4), complex)
    for r0 in range(2):
        for r1 in range(2):
            # Clock bit8 is low within the two-clock register.
            state[:, r0+2*r1] = np.kron(np.kron(clock[:, r1], clock[:, r0]), vector)
    return state, vector, {"work_preparations": 1}, metadata


def _ideal_diagonal(state, effect, targets):
    roots = (np.sqrt(effect), np.sqrt(1-effect))
    return tuple(model.apply_subsystem(state, np.diag(root), targets, QUBITS) for root in roots)


def _ideal_projector(state, vector):
    projector = np.outer(vector, vector.conj())
    plus = model.apply_subsystem(state, projector, range(8), QUBITS)
    return plus, state-plus


def _ideal_pauli(state, bit, axis):
    pauli = {"X": np.array([[0., 1.], [1., 0.]]),
             "Y": np.array([[0., -1j], [1j, 0.]])}[axis]
    changed = model.apply_subsystem(state, pauli, [bit], QUBITS)
    return (state+changed)/2, (state-changed)/2


def execute(*, interaction=True, interval=F(1, 4), steps=64, kick=True,
            propagation=True, second_interval_interaction=None, backend="source"):
    """Execute one prepared history with five logical records and all aborts.

    At s=t/2: weak charge read. At s=t: neutral charge-transfer return read.
    Apply a site phase, continue the SAME work for t, then read site occupation
    and both clock quadratures. Clock phase at 2t measures the interval used by
    the earlier return record. Earlier records are never overwritten.
    """
    if type(interaction) is not bool or type(kick) is not bool or type(propagation) is not bool:
        raise ValueError("Boolean law/intervention controls required")
    if second_interval_interaction is not None and type(second_interval_interaction) is not bool:
        raise ValueError("optional Boolean later-interaction control required")
    if type(interval) not in (int, F) or not 0 < interval <= F(1, 2):
        raise ValueError("exact interval in (0,1/2] required by the phase-window protocol")
    if type(steps) is not int or not 1 <= steps <= 4096:
        raise ValueError("positive split count at most4096 required")
    if backend == "source":
        state, vector, resources, metadata = _source_preparation(interaction)
    elif backend == "continuous":
        state, vector, resources, metadata = _ideal_preparation(interaction)
    else:
        raise ValueError("source or independent continuous backend required")
    branches, aborts = {(): state}, {}
    prefix_probabilities = {}
    max_scratch_return = 0.

    def advance(duration, *, work_propagates=True, interaction_control=None):
        nonlocal branches
        outputs = {}
        cost = None
        for record, value in branches.items():
            if work_propagates:
                enabled = interaction if interaction_control is None else interaction_control
                output, cost = model.evolve(value, duration, steps, interaction=enabled, backend=backend)
            else:
                # Pin the work only after the committed mass read. Both clocks
                # still advance by the same interval; this is an intervention.
                output = value.copy()
                for bit in (8, 9):
                    phase = -model.parameters().clock*float(duration)/2
                    if backend == "source":
                        output = compiler.execute([compiler.pulse("Z", phase, [bit])], QUBITS, output)
                    else:
                        output = model.apply_subsystem(output, np.diag(np.exp(1j*phase*np.array([1., -1.]))), [bit], QUBITS)
                cost = {"native_pulses": 2, "proper_code_service_events": 4}
            outputs[record] = output
        branches = outputs
        _resource_add(resources, cost or {})  # cost of one realized path

    def read(kind, *, effect=None, targets=None, bit=None, axis=None):
        nonlocal branches, max_scratch_return
        outputs = {}
        cost = None
        for record, value in branches.items():
            if backend == "source":
                if kind == "diagonal":
                    result = history.diagonal_read(value, effect, targets=targets)
                elif kind == "projector":
                    result = history.projector_read(value, vector, targets=range(8))
                else:
                    result = history.pauli_read(value, bit, axis=axis)
                quantum = result.outputs
                cost = result.resources
                for raw, failed in zip(result.failure_records, result.failure_outputs, strict=True):
                    aborts[(*record, "scratch_abort", *raw)] = failed
                max_scratch_return = max(max_scratch_return, result.scratch_return_error)
            elif kind == "diagonal":
                quantum = _ideal_diagonal(value, effect, targets)
            elif kind == "projector":
                quantum = _ideal_projector(value, vector)
            else:
                quantum = _ideal_pauli(value, bit, axis)
            for outcome, output in enumerate(quantum):
                outputs[(*record, outcome)] = output
        branches = outputs
        _resource_add(resources, cost or {})

    advance(interval/2)
    charge_z = np.array([1. if mask.bit_count() == 3 else -1. for mask in range(16)])
    read("diagonal", effect=(1+float(WEAK_STRENGTH)*charge_z)/2, targets=range(4))
    advance(interval/2)
    read("projector")
    prefix_probabilities = {record: norm(value) for record, value in branches.items()}

    # A gauge-neutral local kick creates a nonstationary density response.
    # Its application is an apparatus setting; it is not an additional force
    # silently added to the autonomous interval Hamiltonian.
    angle = math.pi/2 if kick else 0.
    kick_word = [compiler.pulse("I", -angle/2, [0]), compiler.pulse("Z", angle/2, [0])]
    for record in branches:
        if backend == "source":
            branches[record] = compiler.execute(kick_word, QUBITS, branches[record])
        else:
            branches[record] = model.apply_subsystem(branches[record], np.diag([1., np.exp(-1j*angle)]), [0], QUBITS)
    _resource_add(resources, compiler.resources(kick_word) if backend == "source" else {})
    advance(interval, work_propagates=propagation, interaction_control=second_interval_interaction)
    read("diagonal", effect=np.array([0., 1.]), targets=[1])
    read("pauli", bit=8, axis="X")
    read("pauli", bit=9, axis="Y")
    resources["maximum_scratch_return_residual"] = max_scratch_return
    resources["maximum_source_qubits_with_scratch"] = 13
    resources["maximum_logical_qubits_including_untouched_reference"] = 15
    if backend == "source":
        resources["proper_code_load_events"] = resources["native_pulses"]
        resources["proper_code_export_events"] = resources["native_pulses"]
        if resources["proper_code_service_events"] != 2*resources["native_pulses"]:
            raise ArithmeticError("inconsistent whole-history source cost ledger")
    resources["physical_service_error_certified"] = False
    specification = {
        "schema": "oph.common_source.mass_response.history.v1", "backend": backend,
        "law": "source-dependent four-mode spinless SSH ring with quartic density interaction",
        "interaction_enabled": interaction, "source_parameters": model.parameters().record(),
        "interval_exact": str(F(interval)), "split_count_per_interval": steps,
        "kick_enabled": kick, "second_interval_work_propagation": propagation,
        "second_interval_interaction": interaction if second_interval_interaction is None else second_interval_interaction,
        "kick_site": 0, "receiver_site": 1, "receiver_graph_distance": 1,
        "preparation": metadata, "weak_strength_exact": str(WEAK_STRENGTH),
        "weak_visibility": math.sqrt(1-float(WEAK_STRENGTH)**2),
        "clock": model.ideal_clock_calibration(), "clock_final_parameter": str(2*F(interval)),
        "record_order": ["weak_charge", "mass_return", "site1_occupied", "clock_X", "clock_Y"],
        "work_reset_count": 0, "all_outcomes_retained": True,
        "charge_conservation_domain": "completed logical operations after preparation, not every intermediate processor pulse",
        "scratch_failure_policy": "retain every failed scratch/pointer outcome as terminal abort",
        "source_matrices_are_contracted_pulse_words": True,
        "physical_particle_or_clock_identification": False,
    }
    return Result(branches, aborts, prefix_probabilities, resources, specification)


def marginal(result, position, value=0):
    return math.fsum(norm(state) for record, state in result.branches.items() if record[position] == value)


def recover_gap_ratio(result):
    """Infer a finite charged gap/clock ratio using ONLY full record means.

    The mass read occurs at t; clock quadratures occur at2t. The work phase
    is2*Delta*t. Identical copies and particle-hole symmetry allow division
    of the neutral pair frequency by two. Ambiguous phase windows are refused.
    This is an exact-law inverse; finite-shot confidence is not asserted here.
    """
    abort_mass = math.fsum(norm(value) for value in result.aborts.values())
    if abort_mass > 1e-10:
        raise ValueError("unresolved abort mass makes this point inverse inconclusive")
    visibility = result.specification["weak_visibility"]
    contrast = result.specification["clock"]["contrast"]
    mass_cosine = (2*marginal(result, 1)-1)/visibility
    clock_x = (2*marginal(result, 3)-1)/contrast
    clock_y = (2*marginal(result, 4)-1)/contrast
    radius = math.hypot(clock_x, clock_y)
    if abs(radius-1) > 1e-7:
        raise ValueError("clock response differs from the declared calibrated circle")
    if not -1 <= mass_cosine <= 1:
        raise ValueError("record mean outside admitted mass phase window")
    angle = math.atan2(clock_y, clock_x)
    if not 0 < angle < math.pi:
        raise ValueError("clock phase is ambiguous or outside the declared window")
    phase = math.acos(mass_cosine)
    return {"finite_single_charged_gap_over_clock": phase/angle,
            "neutral_pair_phase": phase, "measured_clock_phase": angle,
            "mass_cosine_from_records": mass_cosine,
            "clock_radius_from_records": radius,
            "unresolved_abort_probability": abort_mass,
            "single_gap_contract": "identical cells and particle-hole symmetry",
            "finite_shot_confidence_included": False,
            "infinite_volume_mass_identification": False}


def quantum_distance(left, right):
    """Half trace norm of COMPLETE pure CQ outputs, including references/aborts."""
    def block(a, b):
        p, q = norm(a), norm(b)
        if p == 0:
            return q
        overlap = np.vdot(a, b)
        residual = b - (overlap/p)*a
        return math.sqrt(max(0., (p-q)**2 + 4*p*norm(residual)))
    l = {**left.branches, **left.aborts}
    r = {**right.branches, **right.aborts}
    template = next(iter(l.values()))
    zero = np.zeros_like(template)
    return math.fsum(block(l.get(key, zero), r.get(key, zero)) for key in l.keys() | r.keys())/2
