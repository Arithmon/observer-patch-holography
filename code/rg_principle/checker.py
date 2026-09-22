"""Independent strict circuit execution; this module never imports the compiler."""

import hashlib
import json


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def validate(program):
    if type(program) is not dict or set(program) != {"n", "m", "work", "table", "gates"}:
        raise ValueError("exact circuit schema required")
    n, m, work = (program[key] for key in ("n", "m", "work"))
    if any(type(x) is not int for x in (n, m, work)) or not (0 <= n <= 8 and 1 <= m <= 8 and 0 <= work <= 32):
        raise ValueError("invalid register counts")
    rows = program["table"]
    if type(rows) is not list or len(rows) != 2**n:
        raise ValueError("missing truth table rows")
    for row in rows:
        if type(row) is not list or len(row) != m or any(type(x) is not int or x not in (0, 1) for x in row):
            raise ValueError("invalid binary output row")
    gates = program["gates"]
    if type(gates) is not list or len(gates) > 100000:
        raise ValueError("invalid gate ledger")
    for gate in gates:
        if (type(gate) is not list or not 1 <= len(gate) <= 3 or
                any(type(x) is not int or not 0 <= x < n + m + work for x in gate) or
                len(set(gate)) != len(gate)):
            raise ValueError("invalid NOT/CNOT/Toffoli operands")
    return n, m, work


def step(state, gate):
    result = list(state)
    if all(state[control] == 1 for control in gate[:-1]):
        result[gate[-1]] = 1 - state[gate[-1]]
    return result


def execute(state, gates):
    history = [list(state)]
    for gate in gates:
        history.append(step(history[-1], gate))
    return history


def history_admissible(program, history):
    """Feasibility uses actual transitions, never agreement with the target truth table."""
    n, m, work = validate(program)
    if type(history) is not list or len(history) != len(program["gates"]) + 1:
        return False
    if any(type(state) is not list or len(state) != n + m + work or
           any(type(bit) is not int or bit not in (0, 1) for bit in state) for state in history):
        return False
    if any(history[0][n:]):
        return False
    return all(after == step(before, gate) for before, after, gate in
               zip(history, history[1:], program["gates"]))


def check_program(program, interventions=False):
    n, m, work = validate(program)
    traces, intervention_traces = [], []
    intervention_gates = 0
    for row in range(2**n):
        # Deliberately independent of compiler.bits and its minterm construction.
        source = [int(bool(row & (2**j))) for j in range(n)]
        start = source + [0] * (m + work)
        history = execute(start, program["gates"])
        expected = source + program["table"][row] + [0] * work
        if history[-1] != expected:
            raise ValueError(f"truth/retention/clean-work failure on input {row}")
        if execute(history[-1], list(reversed(program["gates"])))[-1] != start:
            raise ValueError("inverse ledger failed")
        traces.append(history)
        if interventions:
            for boundary, before in enumerate(history):
                for register in range(n + m + work):
                    changed = list(before)
                    changed[register] ^= 1
                    suffix = program["gates"][boundary:]
                    after = execute(changed, suffix)[-1]
                    intervention_gates += len(suffix)
                    # A permutation must distinguish any isolated change of its full state.
                    if after == history[-1]:
                        raise ValueError("localized intervention disappeared from full state")
                    # Dependency closure supplies a separate noninterference control.
                    influenced = {register}
                    for gate in suffix:
                        if influenced.intersection(gate):
                            influenced.add(gate[-1])
                    if any(after[j] != history[-1][j] for j in range(len(after)) if j not in influenced):
                        raise ValueError("influence outside complete gate dependency closure")
                    intervention_traces.append([row, boundary, register, after])
    return {"program_sha256": digest(program), "registers": n + m + work,
            "work_registers": work, "gates_per_execution": len(program["gates"]),
            "inputs": len(traces), "forward_gates": len(traces) * len(program["gates"]),
            "inverse_gates": len(traces) * len(program["gates"]),
            "intervention_gates": intervention_gates,
            "executed_gates": 2 * len(traces) * len(program["gates"]) + intervention_gates,
            "trace_sha256": digest(traces), "interventions": len(intervention_traces),
            "intervention_sha256": digest(intervention_traces)}


def check_boundary_lowering(logical, native):
    """Diagnostic for same-register/same-boundary lowerings, including cached substitution."""
    n, m, work = validate(logical)
    if validate(native) != (n, m, work) or len(native["gates"]) != len(logical["gates"]):
        raise ValueError("diagnostic requires the same register and boundary interface")
    for value in range(2**n):
        start = [int(bool(value & (2**j))) for j in range(n)] + [0] * (m + work)
        left = execute(start, logical["gates"])
        right = execute(start, native["gates"])
        if left != right:
            raise ValueError("baseline boundary mismatch")
        for boundary, state in enumerate(left):
            for register in range(len(state)):
                intervened = list(state)
                intervened[register] ^= 1
                a = execute(intervened, logical["gates"][boundary:])[-1]
                b = execute(intervened, native["gates"][boundary:])[-1]
                if a != b:
                    raise ValueError("localized intermediate intervention mismatch")
    return True


def decode_is_lossless(words, encode, decode):
    """Check the declared complete finite input set, not sampled collisions."""
    images = set()
    for word in words:
        transcript = encode(word)
        if decode(transcript) != word:
            return False
        images.add(transcript)
    return len(images) == len(words)
