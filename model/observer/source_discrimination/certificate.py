"""Exact encoded-channel checks and a conditional apparatus error contract.

The channel comparison independently derives floor/ceiling repair. It does
not call the prediction oracle or the native update. Exact permutation and
projector identities concern the mathematical source construction; they do
not certify floating compilation, numerical execution or physical services.
"""
from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations_with_replacement

from . import circuit


APPARATUS_ALLOWANCE = F(1, 64)
ERROR_TERMS = (
    "pulse_operator", "load_channel", "export_channel", "preparation",
    "reset", "coin", "diagonal_randomizer", "record", "storage", "readout",
    "compilation", "numerical_execution",
)
RESOURCE_TERMS = ("native_pulses", "proper_code_load_events", "proper_code_export_events")


def _states():
    """Independent enumeration; no imported carrier or endpoint oracle."""
    states = []
    for first, second in combinations_with_replacement(range(12), 2):
        state = [0] * 12
        state[first] += 1
        state[second] += 1
        states.append(tuple(state))
    return tuple(sorted(states))


def _endpoint_image(state, candidate, coin):
    total = state[0] + state[1]
    low, high = total // 2, (total + 1) // 2
    endpoints = (low, high) if coin == 0 else (high, low)
    output = endpoints + state[2:]
    value = state[0] - endpoints[0] if candidate == "C" else state[0] - state[1]
    return (coin, value), output


def _permutation_image(candidate, coin, index):
    for left, right in circuit.REPAIR_SWAPS[candidate, coin]:
        if index == left:
            index = right
        elif index == right:
            index = left
    return index


def _matrix_unit_image(left, right):
    result = {}
    for coin in (0, 1):
        record_left, target_left = left[coin]
        record_right, target_right = right[coin]
        if record_left == record_right:
            result[record_left, target_left, target_right] = F(1, 2)
    return result


def channel_certificate(candidate):
    """Intertwining on every matrix unit, with all committed records intact.

    For the nine-qubit embedding J, the encoded recorded channel satisfies
    E_code(J rho J*) = (id_record tensor J) E_native(rho) (...)*.
    Matrix-unit equality implies the same identity with an arbitrary
    reference system. It is not just a comparison of basis populations.
    """
    if candidate not in ("C", "E"):
        raise ValueError("candidate C or E required")
    states = _states()
    codes = tuple(circuit.encode(state) for state in states)
    if len(states) != 78 or len(set(codes)) != 78:
        raise AssertionError("the full 78-state carrier must embed injectively")
    if any(circuit.decode(code) != state for code, state in zip(codes, states, strict=True)):
        raise AssertionError("encoding must preserve every native configuration")
    decoders = {"C": {0: 0, 1: 1, 2: -1}, "E": {0: 0, 1: 1, 2: -1, 3: 2, 4: -2}}
    expected, encoded = [], []
    for state, code in zip(states, codes, strict=True):
        native_column, encoded_column = [], []
        for coin in (0, 1):
            record, output = _endpoint_image(state, candidate, coin)
            native_column.append((record, circuit.encode(output)))
            local = _permutation_image(candidate, coin, code % 8)
            value = decoders[candidate].get(local // 8)
            if value is None:
                raise AssertionError("an admitted source column reached a reserved record")
            encoded_column.append(((coin, value), (code // 8) * 8 + local % 8))
        expected.append(tuple(native_column))
        encoded.append(tuple(encoded_column))
    if encoded != expected:
        raise AssertionError("the encoded dilation does not implement floor/ceiling repair")

    # Each joint (coin, record, output) column contains two amplitudes
    # 1/sqrt(2); distinct columns have disjoint support. This proves W*W=I.
    occupied = set()
    for column in encoded:
        support = set(column)
        if len(support) != 2 or occupied.intersection(support):
            raise AssertionError("joint record/source columns are not isometric")
        occupied.update(support)
    units = 0
    for left in range(len(states)):
        for right in range(len(states)):
            actual = _matrix_unit_image(encoded[left], encoded[right])
            target = _matrix_unit_image(expected[left], expected[right])
            if actual != target:
                raise AssertionError("recorded channels disagree on a matrix unit")
            units += 1
    return {
        "candidate": candidate, "native_dimension": len(states),
        "encoded_qubits": circuit.SOURCE_QUBITS, "matrix_units_checked": units,
        "kraus_entries": len(occupied),
        "record_count": len({record for column in encoded for record, _ in column}),
        "column_amplitude_squared": "1/2", "gram_identity_exact": True,
        "floor_ceiling_intertwining_exact": True,
        "arbitrary_reference_extension": True,
        "scope": "exact encoded target channel on the full native carrier, retaining every record; floating compilation and physical identification are separate",
    }


def _encoded_readout_column(address, outcome):
    """Exact U |pointer-sector><pointer-sector| U* column.

    U is the bit-zero Y preparation; two factors 1/sqrt(2) combine to 1/2.
    The pointer flips only at the ALL-zero nine-bit address, not at every
    even address. The inverse basis change is followed by its restoration.
    """
    base = address & ~1
    inverse_column = ((base, 1), (base + 1, 1 if address & 1 else -1))
    result = {}
    for middle, before in inverse_column:
        if (middle == 0) != (outcome == 1):
            continue
        forward_column = ((base, -1 if middle & 1 else 1), (base + 1, 1))
        for target, after in forward_column:
            result[target] = result.get(target, F()) + F(before * after, 2)
    return {target: value for target, value in result.items() if value}


def readout_certificate():
    """Exact full complement, including coherences among all minus states."""
    states = _states()
    codes = tuple(circuit.encode(state) for state in states)
    x = tuple(int(port in (2, 3)) for port in range(12))
    z = tuple(int(port in (1, 2)) for port in range(12))
    if (circuit.encode(x), circuit.encode(z)) != (0, 1):
        raise AssertionError("the fixed witness states must occupy addresses zero and one")
    columns = {}
    for address in range(1 << circuit.SOURCE_QUBITS):
        plus = {0: F(1, 2), 1: F(1, 2)} if address in (0, 1) else {}
        minus = {target: -value for target, value in plus.items()}
        minus[address] = minus.get(address, F()) + 1
        minus = {target: value for target, value in minus.items() if value}
        expected = {1: plus, -1: minus}
        for outcome in (1, -1):
            actual = _encoded_readout_column(address, outcome)
            if actual != expected[outcome]:
                raise AssertionError("the pointer does not implement the full Lueders projector")
            columns[address, outcome] = actual
    # Since each full-register Kraus column equals P+ or I-P+, the matrix
    # unit images are their exact outer products. Check their trace sum too.
    units = 0
    for left in codes:
        for right in codes:
            trace = sum((a * columns[right, outcome].get(target, F())
                         for outcome in (1, -1)
                         for target, a in columns[left, outcome].items()), F())
            if trace != int(left == right):
                raise AssertionError("the complete readout channel is not trace preserving")
            units += 1
    return {
        "encoded_basis_columns_checked": 1 << circuit.SOURCE_QUBITS,
        "native_matrix_units_checked": units, "plus_rank": 1,
        "native_minus_rank": 77, "encoded_minus_rank": (1 << circuit.SOURCE_QUBITS) - 1,
        "both_unnormalized_branches_retained": True,
        "minus_subspace_coherence_preserved": True,
        "old_record_blocks_unchanged": True, "arbitrary_reference_extension": True,
    }


def apparatus_error_contract(resources, error_bounds):
    """Sufficient conditional bound for the unchanged eta=1/64 protocol.

    Resources are maxima over COMPLETE attempted-trial branches, including
    the recorded fair coin and diagonal-control randomizer. Pulse errors
    are operator-norm bounds. Load/export errors are half-diamond bounds
    per invocation, on complete channels including failures. All other
    terms are separately supplied aggregate whole-trial trace/channel
    bounds, covering additional errors not already charged to these words.
    Every term must be explicitly supplied as an exact nonnegative rational.
    This arithmetic never establishes any of those physical premises.
    """
    if not isinstance(resources, dict) or any(name not in resources for name in RESOURCE_TERMS):
        raise ValueError("complete pulse, load and export resource counts required")
    if any(type(resources[name]) is not int or resources[name] < 0 for name in RESOURCE_TERMS):
        raise ValueError("nonnegative integer resource counts required")
    if not isinstance(error_bounds, dict) or set(error_bounds) != set(ERROR_TERMS):
        raise ValueError("every declared error term must be explicitly supplied")
    if any(not isinstance(value, F) or value < 0 for value in error_bounds.values()):
        raise ValueError("nonnegative exact rational error bounds required")
    contributions = {
        "pulse_operator": resources["native_pulses"] * error_bounds["pulse_operator"],
        "load_channel": resources["proper_code_load_events"] * error_bounds["load_channel"],
        "export_channel": resources["proper_code_export_events"] * error_bounds["export_channel"],
    }
    contributions.update({name: error_bounds[name] for name in ERROR_TERMS[3:]})
    additive = sum(contributions.values(), F())
    bound = min(F(1), additive)
    return {
        "contributions": {name: str(value) for name, value in contributions.items()},
        "additive_error_upper": str(additive),
        "complete_instrument_trace_distance_upper_if_premises_hold": str(bound),
        "every_joint_plus_and_abort_probability_error_upper_if_premises_hold": str(bound),
        "frozen_protocol_apparatus_allowance": str(APPARATUS_ALLOWANCE),
        "sufficient_for_frozen_protocol_if_premises_hold": additive <= APPARATUS_ALLOWANCE,
        "physical_error_certificate": False, "certified_numerical_execution": False,
        "scope": "conditional telescoping implication only; all supplied upper bounds require independent justification; all failed attempts remain abort outcomes",
    }
