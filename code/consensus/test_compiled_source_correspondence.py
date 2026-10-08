"""Exact original-transition controls for the finite settling compiler.

The oracle uses the supplied primitive tables, current local registers and
declared readback maps.  It does not use the compiler's derived tables,
layouts, source slots, extension, or settling potential.
"""
from __future__ import annotations

from dataclasses import replace
from itertools import product
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compiled_lattice_settling_certificate as cert


def _bits(size):
    return tuple(product((0, 1), repeat=size))


def _primitive(name, inputs, registers, outputs, readback, truth, update):
    return cert.Primitive(
        name, tuple(inputs), tuple(registers), tuple(outputs), tuple(readback),
        {ports: tuple(truth(ports)) for ports in _bits(len(inputs))},
        {(ports, state): tuple(update(ports, state))
         for ports in _bits(len(inputs)) for state in _bits(len(registers))},
    )


def _wire(*, padded=False, stateful=False, hidden_cycle=False):
    if hidden_cycle:
        return _primitive("WIRE", ("x",), ("q", "hidden"), ("z",), (0,),
                          lambda p: (p[0],), lambda p, s: (p[0], 1 - s[1]))
    if padded:
        return _primitive("WIRE", ("x",), ("scratch", "q"), ("z",), (1,),
                          lambda p: (p[0],), lambda p, s: (1 - p[0], p[0]))
    return _primitive("WIRE", ("x",), ("q",), ("z",), (0,),
                      lambda p: (p[0],),
                      lambda p, s: (p[0] ^ s[0],) if stateful else (p[0],))


def _one_wire_net():
    return cert.Netlist("one_wire", ("x",),
                        (cert.Instance("W", "WIRE", (("in", "x"),)),), (("answer", "W"),))


def _identity_circuit():
    return cert.GateCircuit("identity", ("x",), (), ("x",), lambda env: {"x": env["x"]})


def _original_states(net, primitives, flat_state):
    result = {}
    cursor = 0
    for instance in net.instances:
        size = len(primitives[instance.kind].registers)
        result[instance.name] = tuple(flat_state[cursor:cursor + size])
        cursor += size
    assert cursor == len(flat_state)
    return result


def _original_step(net, primitives, inputs, state):
    current = _original_states(net, primitives, state)
    external = dict(zip(net.inputs, inputs))
    instances = {instance.name: instance for instance in net.instances}
    next_state = []
    for instance in net.instances:
        ports = []
        for reference in instance.sources:
            if reference[0] == "in":
                ports.append(external[reference[1]])
            else:
                parent = instances[reference[1]]
                primitive = primitives[parent.kind]
                register = primitive.readback[primitive.out_ports.index(reference[2])]
                ports.append(current[parent.name][register])
        primitive = primitives[instance.kind]
        next_state.extend(primitive.update[(tuple(ports), current[instance.name])])
    return tuple(next_state)


def _original_outputs(net, primitives, state):
    current = _original_states(net, primitives, state)
    instances = {instance.name: instance for instance in net.instances}
    result = {}
    for label, name in net.outputs:
        primitive = primitives[instances[name].kind]
        assert len(primitive.out_ports) == 1
        result[label] = current[name][primitive.readback[0]]
    return result


@pytest.mark.parametrize("padded", [False, True])
def test_compiled_step_readback_and_settling_match_the_original_kernel(padded):
    primitives = cert.reference_primitives()
    primitives["WIRE"] = _wire(padded=padded)
    # Deliberately not in topological list order. Each output consumes its own
    # fan-out port; all initial registers, including scratch bits, are varied.
    net = cert.Netlist("branched", ("x",), (
        cert.Instance("right", "WIRE", (("out", "fork", "z1"),)),
        cert.Instance("fork", "FANOUT2", (("in", "x"),)),
        cert.Instance("left", "WIRE", (("out", "fork", "z0"),)),
    ), (("left", "left"), ("right", "right")))
    compiled = cert.compile_net(net, primitives)
    size = sum(len(primitives[item.kind].registers) for item in net.instances)
    all_states = _bits(size)
    for inputs in _bits(1):
        fixed = [state for state in all_states if _original_step(net, primitives, inputs, state) == state]
        assert len(fixed) == 1
        assert cert.extension_state(compiled, inputs) == fixed[0]
        for state in all_states:
            following = _original_step(net, primitives, inputs, state)
            assert cert.realized_step(compiled, inputs, state) == following
            assert cert.output_values(compiled, state) == _original_outputs(net, primitives, state)
            # This two-layer source kernel settles in at most two updates.
            assert _original_step(net, primitives, inputs, following) == fixed[0]
            expected_time = 0 if state == fixed[0] else (1 if following == fixed[0] else 2)
            assert cert.settle_trajectory(compiled, inputs, state, fixed[0])[0] == expected_time


@pytest.mark.parametrize("exhaustive", [False, True])
def test_nondefault_register_encodings_preserve_circuit_semantics(exhaustive):
    primitives = cert.reference_primitives()
    primitives["WIRE"] = _wire(padded=True)
    primitives["NAND"] = _primitive(
        "NAND", ("x", "y"), ("scratch", "q"), ("z",), (1,),
        lambda p: (1 - (p[0] & p[1]),), lambda p, s: (1, 1 - (p[0] & p[1])),
    )
    circuit = cert.GateCircuit("not_via_nand", ("x",), (("not_x", ("x", "x")),),
                               ("not_x",), lambda env: {"not_x": 1 - env["x"]})
    for primitive in primitives.values():
        assert cert.verify_primitive(primitive)["evidence_bundle"]["state_independent_update"]
    report = cert.circuit_block(circuit, primitives, exhaustive=exhaustive)
    assert report["extension_semantics_checks"] == 2
    assert report["settling"]["bound_verified"] is True
    assert report["rank_ladder"]["verified"] is True
    assert report["register_count"] == 6
    assert report["settling"]["worst_settling_time"] == 3


@pytest.mark.parametrize("hidden_cycle", [False, True])
@pytest.mark.parametrize("entrypoint", ["compile", "exhaustive", "family"])
def test_stateful_original_kernel_cannot_be_replaced_by_its_zero_state_row(hidden_cycle, entrypoint):
    primitive = _wire(hidden_cycle=True) if hidden_cycle else _wire(stateful=True)
    primitives = {**cert.reference_primitives(), "WIRE": primitive}
    net = _one_wire_net()
    states = _bits(len(primitive.registers))
    # Independent original-input obstruction: neither kernel has any fixed
    # register state at x=1, irrespective of its declared output truth.
    assert all(_original_step(net, primitives, (1,), state) != state for state in states)
    initial = states[0]
    first = _original_step(net, primitives, (1,), initial)
    if not hidden_cycle:
        assert _original_step(net, primitives, (1,), first) == initial
    with pytest.raises(cert.CertificateError):
        if entrypoint == "compile":
            cert.compile_net(net, primitives)
        else:
            cert.circuit_block(_identity_circuit(), primitives, exhaustive=entrypoint == "exhaustive")


def _bad_domain(kind):
    primitive = _wire()
    truth, update = dict(primitive.truth), dict(primitive.update)
    if kind == "missing_truth":
        del truth[(1,)]
    elif kind == "missing_state":
        del update[((1,), (1,))]
    elif kind == "extra_state":
        update[((1,), (2,))] = (1,)
    elif kind == "extra_port":
        truth[(2,)] = (0,)
    elif kind == "wrong_width":
        update = {key: value + (0,) for key, value in update.items()}
    elif kind == "nonboolean_fixed_state":
        truth = {ports: (2,) for ports in truth}
        update = {key: (2,) for key in update}
        update.update({(ports, (2,)): (2,) for ports in truth})
    else:
        raise AssertionError(kind)
    return replace(primitive, truth=truth, update=update)


@pytest.mark.parametrize("kind", ["missing_truth", "missing_state", "extra_state", "extra_port", "wrong_width", "nonboolean_fixed_state"])
@pytest.mark.parametrize("entrypoint", ["primitive", "compile"])
def test_certificate_requires_the_complete_boolean_transition_domain(kind, entrypoint):
    with pytest.raises(cert.CertificateError):
        primitive = _bad_domain(kind)
        if entrypoint == "primitive":
            cert.verify_primitive(primitive)
        else:
            cert.compile_net(_one_wire_net(), {**cert.reference_primitives(), "WIRE": primitive})


@pytest.mark.parametrize("stateful", [False, True])
def test_rank_certificate_is_bound_to_the_actual_source_library(stateful):
    primitives = cert.reference_primitives()
    compiled = cert.compile_net(_one_wire_net(), primitives)
    replacement = _wire(stateful=True) if stateful else _primitive(
        "WIRE", ("x",), ("q",), ("z",), (0,), lambda p: (0,), lambda p, s: (0,),
    )
    with pytest.raises(cert.CertificateError):
        cert.rank_ladder_certification(compiled, {**primitives, "WIRE": replacement})


@pytest.mark.parametrize("field", ["tables", "depth"])
def test_rank_certificate_rejects_altered_derived_compilation(field):
    primitives = cert.reference_primitives()
    compiled = cert.compile_net(_one_wire_net(), primitives)
    with pytest.raises(cert.CertificateError):
        if field == "tables":
            tables = {name: dict(table) for name, table in compiled.tables.items()}
            key = next(iter(tables["W"]))
            tables["W"][key] = tuple(1 - bit for bit in tables["W"][key])
            changed = replace(compiled, tables=tables)
        else:
            changed = replace(compiled, analysis={**compiled.analysis, "depth": 0})
        cert.rank_ladder_certification(changed, primitives)


@pytest.mark.parametrize("case", ["duplicate_inputs", "duplicate_output_labels"])
def test_distinct_interface_coordinates_cannot_be_overwritten(case):
    if case == "duplicate_inputs":
        net = replace(_one_wire_net(), inputs=("x", "x"))
    else:
        net = cert.Netlist("two_outputs", ("x", "y"), (
            cert.Instance("W1", "WIRE", (("in", "x"),)),
            cert.Instance("W2", "WIRE", (("in", "y"),)),
        ), (("answer", "W1"), ("answer", "W2")))
    with pytest.raises(cert.CertificateError):
        cert.compile_net(net, cert.reference_primitives())


def test_gate_compiler_rejects_duplicate_input_coordinates():
    circuit = replace(_identity_circuit(), inputs=("x", "x"))
    with pytest.raises(cert.CertificateError):
        cert.circuit_block(circuit, cert.reference_primitives(), exhaustive=True)
