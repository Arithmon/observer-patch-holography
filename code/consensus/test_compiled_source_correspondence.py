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


@pytest.mark.parametrize("kernel", range(16))
@pytest.mark.parametrize("truth", range(4))
def test_all_one_bit_kernels_are_classified_by_the_original_transition(kernel, truth):
    # Four independently specified bits give every possible (x, q) -> q'
    # transition. Two further bits give every declared x -> z function.
    update = {(ports, state): ((kernel >> (2 * ports[0] + state[0])) & 1,)
              for ports in _bits(1) for state in _bits(1)}
    declared = {ports: ((truth >> ports[0]) & 1,) for ports in _bits(1)}
    admissible = all(update[ports, state] == declared[ports]
                     for ports in _bits(1) for state in _bits(1))
    # Exactly four of the 64 pairs are memoryless and match their declaration.
    assert admissible == (kernel == ((truth & 1) * 3 + ((truth >> 1) & 1) * 12))
    primitive = cert.Primitive("CELL", ("x",), ("q",), ("z",), (0,), declared, update)
    primitives = {"CELL": primitive}
    net = cert.Netlist("one_cell", ("input",),
                       (cert.Instance("cell", "CELL", (("in", "input"),)),),
                       (("answer", "cell"),))
    if not admissible:
        with pytest.raises(cert.CertificateError):
            cert.compile_net(net, primitives)
        return
    compiled = cert.compile_net(net, primitives)
    assert compiled.analysis["depth"] == 1
    assert cert.rank_ladder_certification(compiled, primitives)["verified"] is True
    for inputs in _bits(1):
        fixed = declared[inputs]
        assert cert.extension_state(compiled, inputs) == fixed
        for state in _bits(1):
            assert cert.realized_step(compiled, inputs, state) == _original_step(net, primitives, inputs, state)
            assert cert.output_values(compiled, state) == _original_outputs(net, primitives, state)
            assert cert.settle_trajectory(compiled, inputs, state, fixed)[0] == int(state != fixed)


def test_compiled_kernel_is_isolated_from_mutable_source_aliases():
    primitive = _wire()
    update, truth = dict(primitive.update), dict(primitive.truth)
    primitive = replace(primitive, update=update, truth=truth)
    library = {"WIRE": primitive}
    compiled = cert.compile_net(_one_wire_net(), library)
    update.update({key: (0,) for key in update})
    truth.update({key: (0,) for key in truth})
    library["WIRE"] = _wire(stateful=True)
    for state in _bits(1):
        assert cert.realized_step(compiled, (1,), state) == (1,)
    assert cert.extension_state(compiled, (1,)) == (1,)
    assert cert.rank_ladder_certification(compiled, {"WIRE": _wire()})["verified"] is True


@pytest.mark.parametrize("field", ["tables", "rank"])
def test_reconstructed_compilation_is_isolated_from_mutable_derived_aliases(field):
    compiled = cert.compile_net(_one_wire_net(), {"WIRE": _wire()})
    if field == "tables":
        tables = {name: dict(table) for name, table in compiled.tables.items()}
        copied = replace(compiled, tables=tables)
        tables["W"].update({key: (0,) for key in tables["W"]})
        assert cert.realized_step(copied, (1,), (0,)) == (1,)
    else:
        rank = dict(compiled.analysis["rank"])
        copied = replace(compiled, analysis={**compiled.analysis, "rank": rank})
        rank["W"] = 23
        assert copied.analysis["rank"]["W"] == 1
    assert cert.rank_ladder_certification(copied, {"WIRE": _wire()})["verified"] is True


@pytest.mark.parametrize("field", ["net", "primitives", "src_slots", "output_slots"])
def test_source_and_execution_coordinates_cannot_be_rebound_to_stale_derivations(field):
    compiled = cert.compile_net(_one_wire_net(), {"WIRE": _wire()})
    if field == "net":
        replacement = cert.Netlist("changed_source", ("x", "y"),
                                  (cert.Instance("W", "WIRE", (("in", "y"),)),),
                                  (("answer", "W"),))
    elif field == "primitives":
        replacement = {"WIRE": _primitive("WIRE", ("x",), ("q",), ("z",), (0,),
                                          lambda p: (0,), lambda p, s: (0,))}
    elif field == "src_slots":
        replacement = {"W": (("r", 0),)}
    else:
        replacement = {"answer": 1}
    with pytest.raises(cert.CertificateError):
        changed = replace(compiled, **{field: replacement})
        cert.rank_ladder_certification(changed, {"WIRE": _wire()})


@pytest.mark.parametrize("readback", [(-1,), (1,), (True,), (), (0, 0)])
def test_readback_must_be_a_complete_valid_register_coordinate(readback):
    with pytest.raises(cert.CertificateError):
        primitive = replace(_wire(), readback=readback)
        cert.compile_net(_one_wire_net(), {"WIRE": primitive})


@pytest.mark.parametrize("field", ["in_ports", "registers", "out_ports"])
def test_primitive_coordinates_must_have_distinct_names(field):
    primitive = _primitive("CELL", ("x", "y"), ("q", "r"), ("z", "w"), (0, 1),
                           lambda p: p, lambda p, s: p)
    with pytest.raises(cert.CertificateError):
        cert.verify_primitive(replace(primitive, **{field: ("same", "same")}))


@pytest.mark.parametrize("case", ["duplicate_gates", "gate_input_collision", "duplicate_outputs"])
def test_gate_circuit_names_cannot_overwrite_distinct_signal_coordinates(case):
    if case == "duplicate_gates":
        circuit = cert.GateCircuit("duplicate_gates", ("a", "b"),
                                   (("n", ("a", "b")), ("n", ("a", "b"))),
                                   ("n",), lambda env: {"n": 1 - (env["a"] & env["b"])})
    elif case == "gate_input_collision":
        circuit = cert.GateCircuit("gate_input_collision", ("a", "b"),
                                   (("a", ("a", "b")),), ("a",),
                                   lambda env: {"a": 1 - (env["a"] & env["b"])})
    else:
        circuit = replace(_identity_circuit(), outputs=("x", "x"))
    with pytest.raises(cert.CertificateError):
        cert.circuit_block(circuit, cert.reference_primitives(), exhaustive=True)


def test_valid_source_names_are_not_restricted_by_generated_patch_prefixes():
    circuit = cert.GateCircuit("prefix_collision", ("g_n", "b"),
                               (("n", ("g_n", "b")),), ("n",),
                               lambda env: {"n": 1 - (env["g_n"] & env["b"])})
    report = cert.circuit_block(circuit, cert.reference_primitives(), exhaustive=True)
    assert report["extension_semantics_checks"] == 4
    assert report["settling"]["worst_settling_time"] == 2
    assert report["rank_ladder"]["verified"] is True


def test_zero_input_constant_cell_has_an_exact_one_round_extension():
    primitive = _primitive("CONST", (), ("q",), ("z",), (0,),
                           lambda p: (1,), lambda p, s: (1,))
    net = cert.Netlist("constant", (), (cert.Instance("C", "CONST", ()),), (("one", "C"),))
    compiled = cert.compile_net(net, {"CONST": primitive})
    assert cert.extension_state(compiled, ()) == (1,)
    for state in _bits(1):
        assert cert.realized_step(compiled, (), state) == _original_step(net, {"CONST": primitive}, (), state)
        assert cert.output_values(compiled, state) == _original_outputs(net, {"CONST": primitive}, state)
        assert cert.settle_trajectory(compiled, (), state, (1,))[0] == int(state != (1,))


def test_aliased_readback_ports_preserve_source_dynamics_and_path_multiplicity():
    # Two distinct output ports may observe the same register. Both port uses
    # count as paths even when they reconverge at the same consuming instance.
    fork = _primitive("ALIAS", ("x",), ("q",), ("a", "b"), (0, 0),
                      lambda p: (p[0], p[0]), lambda p, s: (p[0],))
    primitives = {"ALIAS": fork, "NAND": cert.reference_primitives()["NAND"]}
    net = cert.Netlist("aliased_fork", ("x",), (
        cert.Instance("G", "NAND", (("out", "F", "a"), ("out", "F", "b"))),
        cert.Instance("F", "ALIAS", (("in", "x"),)),
    ), (("answer", "G"),))
    compiled = cert.compile_net(net, primitives)
    # F has its length-zero path and two length-one paths to G.
    assert dict(compiled.analysis["weight"]) == {"F": 3, "G": 1}
    assert dict(compiled.analysis["rank"]) == {"F": 1, "G": 2}
    assert cert.rank_ladder_certification(compiled, primitives)["verified"] is True
    states = _bits(2)
    for inputs in _bits(1):
        fixed = [state for state in states
                 if _original_step(net, primitives, inputs, state) == state]
        assert fixed == [(1 - inputs[0], inputs[0])]
        assert cert.extension_state(compiled, inputs) == fixed[0]
        for state in states:
            following = _original_step(net, primitives, inputs, state)
            assert cert.realized_step(compiled, inputs, state) == following
            assert cert.output_values(compiled, state) == _original_outputs(net, primitives, state)
            assert cert.potential(compiled, state, fixed[0]) == (
                int(state[0] != fixed[0][0]) + 3 * int(state[1] != fixed[0][1])
            )
            assert _original_step(net, primitives, inputs, following) == fixed[0]
            expected_time = 0 if state == fixed[0] else (1 if following == fixed[0] else 2)
            assert cert.settle_trajectory(compiled, inputs, state, fixed[0])[0] == expected_time


def test_settling_reports_the_exact_nonzero_potential_margin():
    primitives = {"WIRE": _wire(), "NAND": cert.reference_primitives()["NAND"]}
    net = cert.Netlist("masked_wire", ("x", "mask"), (
        cert.Instance("W", "WIRE", (("in", "x"),)),
        cert.Instance("G", "NAND", (("out", "W", "z"), ("in", "mask"))),
    ), (("answer", "G"),))
    inputs, initial, fixed = (0, 0), (1, 0), (0, 1)
    # The zero mask makes NAND settle despite its unsettled wire input.
    assert _original_step(net, primitives, inputs, initial) == fixed
    assert _original_step(net, primitives, inputs, fixed) == fixed
    compiled = cert.compile_net(net, primitives)
    assert cert.extension_state(compiled, inputs) == fixed
    assert cert.realized_step(compiled, inputs, initial) == fixed
    assert dict(compiled.analysis["weight"]) == {"W": 2, "G": 1}
    assert cert.potential(compiled, initial, fixed) == 3
    assert cert.potential(compiled, fixed, fixed) == 0
    # Exact margin = potential drop 3 - two initially unsettled instances.
    assert cert.settle_trajectory(compiled, inputs, initial, fixed) == (1, 1)


@pytest.mark.parametrize("has_one_row", [False, True])
def test_incomplete_large_state_declaration_is_rejected_before_enumeration(has_one_row):
    # The supplied data are tiny. The missing 2**65 transition rows must not
    # induce enumeration of the claimed domain just to discover incompleteness.
    state = (0,) * 64
    update = {((0,), state): state} if has_one_row else {}
    primitive = cert.Primitive("CELL", ("x",), tuple(f"q{k}" for k in range(64)),
                               ("z",), (63,), {(0,): (0,), (1,): (1,)}, update)
    net = cert.Netlist("incomplete", ("x",),
                       (cert.Instance("C", "CELL", (("in", "x"),)),), (("answer", "C"),))
    with pytest.raises(cert.CertificateError):
        cert.compile_net(net, {"CELL": primitive})
