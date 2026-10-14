"""Independent finite-matrix and hostile-input checks of protocol composition.

The oracle constructs the declared Fock and charge matrices from their original
parameters. It imports no producer Hamiltonian, preparation, pointer or splitting
helper; comparing full branch maps also tests arbitrary reference coherence.
"""
from dataclasses import replace
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import protocol as p


def _independent_matrices(d, coupling):
    lowering = np.zeros((d, d))
    for occupation in range(1, d):
        lowering[occupation-1, occupation] = math.sqrt(occupation)
    q = (lowering+lowering.T)/math.sqrt(2)
    momentum = (lowering-lowering.T)/(1j*math.sqrt(2))
    # Project Q^2 and P^2 before taking the cutoff, preserving the upper edge.
    boundary = np.diag([0.]*(d-1)+[d/2])
    onsite = ((momentum@momentum).real+boundary
              +(1+9/32)*(q@q+boundary))/2
    charge_x = np.array([[0., 1.], [1., 0.]])
    x0, x1 = np.kron(np.eye(2), charge_x), np.kron(charge_x, np.eye(2))
    space = np.eye(d*d)
    terms = {
        "hopping": np.kron(space, x0+x1),
        "electric": np.kron(space, np.diag([0., 1., 1., 0.])),
        "onsite": np.kron(np.kron(np.eye(d), onsite)
                          +np.kron(onsite, np.eye(d)), np.eye(4)),
        "bond": -(9/32)*np.kron(np.kron(q, q), np.eye(4)),
        "interaction": float(coupling)*math.sqrt(2)
        *np.kron(np.kron(np.eye(d), q), np.diag([1., -1., 1., -1.])),
    }
    coordinate = np.kron(np.kron(np.eye(d), q), np.eye(4))
    return terms, coordinate, np.kron(space, x0)


def _independent_split(terms, duration, steps):
    dt = float(duration)/steps
    hopping = expm(-.25j*dt*terms["hopping"])
    gauge = hopping@expm(-.5j*dt*terms["electric"])@hopping
    onsite = expm(-.5j*dt*terms["onsite"])
    bond = expm(-.5j*dt*terms["bond"])
    step = gauge@onsite@bond@expm(-1j*dt*terms["interaction"])@bond@onsite@gauge
    return np.linalg.matrix_power(step, steps)


def _oracle(program):
    terms, q, feedback_x = _independent_matrices(program.dimension, program.coupling)
    size = len(q)
    minus, plus = expm(-.5j*q), expm(.5j*q)
    ramsey = ((minus-1j*plus)/2, (minus+1j*plus)/2)
    gauge = tuple(np.diag([int(index % 2 == bit) for index in range(size)])
                  for bit in (0, 1))
    branches = {(): np.eye(size, dtype=complex)}
    for operation in program.operations:
        updated = {}
        for record, branch in branches.items():
            if operation.kind == "evolve":
                updated[record] = _independent_split(terms, operation.time, program.steps)@branch
            elif operation.kind == "feedback":
                angle = math.atan2(4, 3) if operation.enabled and record[operation.record] else 0.
                updated[record] = expm(-1j*angle*feedback_x)@branch
            else:
                maps = ramsey if operation.kind == "smooth_read" else gauge
                for bit, operator in enumerate(maps):
                    updated[record+(bit,)] = operator@branch
        branches = updated
    return branches


def test_original_input_branch_oracle_with_entangled_reference_and_changed_parameters():
    program = p.Program(dimension=2, steps=3, coupling=F(-2, 7), operations=(
        p.Operation("smooth_read"), p.Operation("feedback", record=0),
        p.Operation("evolve", time=F(2, 9)), p.Operation("gauge_read"),
        p.Operation("feedback", record=1, enabled=False),
        p.Operation("evolve", time=F(1, 13)), p.Operation("smooth_read")))
    initial = np.array([complex(index % 7-3, index % 5-2) for index in range(32)]).reshape(16, 2)
    initial /= np.linalg.norm(initial)
    assert np.linalg.matrix_rank(initial) == 2
    history = p.execute(program, initial=initial)
    expected = _oracle(program)
    assert {branch.record for branch in history.branches} == set(expected)
    for branch in history.branches:
        np.testing.assert_allclose(branch.kraus, expected[branch.record], atol=3e-12, rtol=0)
        np.testing.assert_allclose(branch.state, expected[branch.record]@initial, atol=3e-12, rtol=0)
    # Vectorization tests the Choi state on a maximally entangled work/reference
    # input and its within-branch coherence. The preceding map checks also retain
    # circuit phases that cancel from each branch's Choi operator.
    actual_choi = sum(np.outer(branch.kraus.ravel(), branch.kraus.ravel().conj())
                      for branch in history.branches)/16
    expected_choi = sum(np.outer(operator.ravel(), operator.ravel().conj())
                        for operator in expected.values())/16
    np.testing.assert_allclose(actual_choi, expected_choi, atol=3e-12, rtol=0)
    assert abs(np.trace(actual_choi)-1) < 3e-12


@pytest.mark.parametrize("changes", [
    {"dimension": True}, {"dimension": np.bool_(False)}, {"dimension": 3},
    {"steps": True}, {"steps": np.bool_(True)}, {"steps": 0}, {"steps": 2.5},
    {"coupling": True}, {"coupling": np.bool_(True)}, {"coupling": math.nan},
    {"coupling": math.inf}, {"coupling": -math.inf}, {"coupling": 1j},
    {"native_current": True}, {"native_current": np.bool_(True)}, {"native_current": 2},
    {"operations": ("smooth_read",)},
])
def test_invalid_program_parameters_fail_closed(changes):
    with pytest.raises((ValueError, TypeError)):
        p.execute(p.Program(**changes), backend="product")


@pytest.mark.parametrize("kind, changes", [
    ("missing", {}), ("evolve", {"time": math.nan}),
    ("evolve", {"time": math.inf}), ("evolve", {"time": -F(1, 3)}),
    ("evolve", {"time": True}), ("evolve", {"record": 0}),
    ("smooth_read", {"time": F(1)}), ("gauge_read", {"enabled": False}),
    ("feedback", {"record": True}), ("feedback", {"record": -1}),
    ("feedback", {"record": 0, "enabled": 1}),
    ("feedback", {"record": 0, "time": F(1)}),
])
def test_invalid_or_silently_ignored_operation_fields_fail_closed(kind, changes):
    with pytest.raises((ValueError, TypeError)):
        p.execute(p.Program(operations=(p.Operation("gauge_read"), p.Operation(kind, **changes))),
                  backend="product")


@pytest.mark.parametrize("initial", [
    np.zeros(16), np.ones(16), np.ones(15)/math.sqrt(15), np.empty((16, 0)),
    np.ones((2, 8))/4, np.ones((16, 1, 1))/4,
    np.array([math.nan]+[0.]*15), np.array([math.inf]+[0.]*15),
    np.array([complex(0., math.inf)]+[0.]*15),
])
def test_invalid_initial_quantum_states_are_not_normalized_or_reshaped_silently(initial):
    with pytest.raises((ValueError, TypeError)):
        p.execute(p.Program(operations=(p.Operation("smooth_read"),)),
                  initial=initial, backend="product")


@pytest.mark.parametrize("initial", [
    [1., False]+[0.]*14,
    [[1.], [np.bool_(False)]]+[[0.]]*14,
    [np.array([1.]), np.array([False])]+[np.array([0.])]*14,
    np.array([1]+[0]*15, dtype=bool),
    [1., None]+[0.]*14,
    np.array([1]+[0]*15, dtype=object),
    np.ma.array([1.]+[0.]*15, mask=[False]*15+[True]),
])
def test_original_amplitude_types_and_missingness_survive_validation(initial):
    with pytest.raises((ValueError, TypeError)):
        p.execute(p.Program(), initial=initial, backend="product")


def test_valid_alternate_initial_state_is_retained_without_vacuum_replacement():
    initial = np.zeros((16, 2), dtype=np.complex128)
    initial[3, 0], initial[12, 1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    result = p.execute(p.Program(), initial=initial, backend="product")
    np.testing.assert_array_equal(result.initial_state, initial)
    np.testing.assert_array_equal(result.branches[0].state, initial)
    assert result.preparation_mode == "supplied_input"
    assert dict(result.preparation_resources)["work_preparations"] == 0


@pytest.fixture(scope="module")
def checkpoint_pair():
    program = p.Program(steps=3, operations=(p.Operation("smooth_read"),
        p.Operation("feedback", record=0), p.Operation("evolve", time=F(1, 7)),
        p.Operation("gauge_read")))
    return program, p.execute(program, backend="product", stop_after=1)


def test_checkpoint_cannot_omit_duplicate_relabel_or_reorder_a_retained_branch(checkpoint_pair):
    program, checkpoint = checkpoint_pair
    a, b = checkpoint.branches
    for branches in ((a,), (a, a), (b, a), (replace(a, record=b.record), b)):
        with pytest.raises((ValueError, TypeError)):
            p.execute(program, backend="product", checkpoint=replace(checkpoint, branches=branches))


def test_checkpoint_seal_binds_quantum_phases_maps_cursor_and_original_preparation(checkpoint_pair):
    program, checkpoint = checkpoint_pair
    branch = checkpoint.branches[0]
    # A phase mutation preserves every saved probability and marginal.
    replacements = [replace(checkpoint, cursor=checkpoint.cursor+1),
                    replace(checkpoint, initial_state=-checkpoint.initial_state)]
    for name in ("state", "kraus"):
        modified = replace(branch, **{name: 1j*getattr(branch, name)})
        replacements.append(replace(checkpoint, branches=(modified, *checkpoint.branches[1:])))
    for altered in replacements:
        with pytest.raises((ValueError, TypeError)):
            p.execute(program, backend="product", checkpoint=altered)


def test_checkpoint_cannot_change_remaining_semantics_or_numerical_backend(checkpoint_pair):
    program, checkpoint = checkpoint_pair
    changed = replace(program, operations=(program.operations[0],
        p.Operation("feedback", record=0, enabled=False), *program.operations[2:]))
    with pytest.raises((ValueError, TypeError)):
        p.execute(changed, backend="product", checkpoint=checkpoint)
    with pytest.raises((ValueError, TypeError)):
        p.execute(program, backend="continuous", checkpoint=checkpoint)


@pytest.mark.parametrize("cursor", [True, np.bool_(True), -1, 99, 1.5])
def test_invalid_pause_positions_do_not_coerce_to_a_different_program(cursor):
    with pytest.raises((ValueError, TypeError)):
        p.execute(p.Program(operations=(p.Operation("smooth_read"),)),
                  backend="product", stop_after=cursor)


@pytest.mark.parametrize("maps", [
    (np.eye(16),),
    (np.eye(16), np.zeros((16, 16)), np.zeros((16, 16))),
    (np.eye(16)/2, np.zeros((16, 16))),
    (np.full((16, 16), math.nan), np.zeros((16, 16))),
    (np.eye(16), np.full((16, 16), math.inf)),
    (np.eye(15), np.zeros((15, 15))),
])
def test_mutated_source_read_cannot_hide_missing_outcomes_or_invalid_matrices(monkeypatch, maps):
    p._operation_maps.cache_clear()
    monkeypatch.setattr(p.smooth, "pointer", lambda dimension: (maps, {}))
    try:
        with pytest.raises((ValueError, ArithmeticError)):
            p.execute(p.Program(operations=(p.Operation("smooth_read"),)))
    finally:
        p._operation_maps.cache_clear()


def test_complete_menu_with_zero_probability_branch_remains_valid(monkeypatch):
    p._operation_maps.cache_clear()
    # The validation boundary admits a complete instrument with a zero branch;
    # this mutation control does not identify it with the real Ramsey detector.
    monkeypatch.setattr(p.smooth, "pointer", lambda dimension:
                        ((np.eye(16), np.zeros((16, 16))), {}))
    try:
        result = p.execute(p.Program(operations=(p.Operation("smooth_read"),)))
        assert result.records == ((0,), (1,))
        assert result.branches[1].probability == 0
        assert abs(result.branches[0].probability-1) < 2e-12
    finally:
        p._operation_maps.cache_clear()


@pytest.mark.parametrize("path, value, reason", [
    (("cursor",), True, "incompatible"),
    (("initial_state", "shape", 0), True, "shape"),
    (("initial_state", "values", 0, 0), False, "amplitudes"),
    (("branches", 0, "record", 0), False, "ordered binary"),
    (("branches", 0, "state", "values", 0, 0), None, "amplitudes"),
    (("preparation_resources", 0, 1), True, "nonnegative integer"),
    (("stages", 0, "kind"), "gauge_read", "stage journal"),
    (("program", "operations", 1, "enabled"), False, "incompatible"),
])
def test_serialized_checkpoint_rejects_original_types_and_changed_semantics(
        checkpoint_pair, path, value, reason):
    program, checkpoint = checkpoint_pair
    payload = json.loads(p.checkpoint_to_json(checkpoint))
    parent = payload
    for key in path[:-1]:
        parent = parent[key]
    parent[path[-1]] = value
    with pytest.raises(ValueError, match=reason):
        p.checkpoint_from_json(program, json.dumps(payload), backend="product")


def test_serialized_checkpoint_rejects_duplicate_keys_nonfinite_values_and_omission(checkpoint_pair):
    program, checkpoint = checkpoint_pair
    encoded = p.checkpoint_to_json(checkpoint)
    duplicate = '{"schema":"duplicate",'+encoded[1:]
    with pytest.raises(ValueError, match="duplicate"):
        p.checkpoint_from_json(program, duplicate, backend="product")
    for value in (math.nan, math.inf, -math.inf):
        payload = json.loads(encoded)
        payload["branches"][0]["state"]["values"][0][0] = value
        with pytest.raises(ValueError, match="nonfinite"):
            p.checkpoint_from_json(program, json.dumps(payload), backend="product")
    payload = json.loads(encoded)
    payload["branches"].pop()
    with pytest.raises(ValueError, match="omitted"):
        p.checkpoint_from_json(program, json.dumps(payload), backend="product")


def test_import_replays_consistently_resealed_but_wrong_instruments(checkpoint_pair):
    program, checkpoint = checkpoint_pair
    # This unitary preserves Gram completeness, every branch probability and
    # state=K*input; a public digest alone cannot detect its consistent forgery.
    swap = np.eye(16)[np.arange(16) ^ 1]
    branches = tuple(replace(branch, state=swap@branch.state, kraus=swap@branch.kraus)
                     for branch in checkpoint.branches)
    forged = replace(checkpoint, branches=branches, seal="")
    forged = replace(forged, seal=p._seal(forged))
    assert abs(forged.probability_mass-checkpoint.probability_mass) < 2e-12
    assert forged.completeness_error < 2e-12
    encoded = p.checkpoint_to_json(forged)
    with pytest.raises(ValueError, match="does not reproduce"):
        p.checkpoint_from_json(program, encoded, backend="product")


def test_valid_reference_checkpoint_roundtrip_preserves_every_complex_component():
    initial = np.zeros((16, 2), complex)
    initial[0, 0], initial[7, 1], initial[12, 0] = .6, .8j/math.sqrt(2), .8/math.sqrt(2)
    program = p.Program(operations=(p.Operation("smooth_read"), p.Operation("gauge_read")))
    history = p.execute(program, initial=initial, backend="product", stop_after=1)
    loaded = p.checkpoint_from_json(program, p.checkpoint_to_json(history), backend="product")
    np.testing.assert_array_equal(loaded.initial_state, initial)
    for actual, saved in zip(loaded.branches, history.branches, strict=True):
        np.testing.assert_array_equal(actual.state, saved.state)
        np.testing.assert_array_equal(actual.kraus, saved.kraus)
    resumed = p.execute(program, backend="product", checkpoint=loaded)
    uninterrupted = p.execute(program, backend="product", initial=initial)
    for actual, expected in zip(resumed.branches, uninterrupted.branches, strict=True):
        np.testing.assert_array_equal(actual.state, expected.state)


@pytest.mark.parametrize("backend", ("source", "product", "continuous"))
@pytest.mark.parametrize("reference_columns", (1, 3))
@pytest.mark.parametrize("cursor", (0, 1))
def test_checkpoint_roundtrip_preserves_signed_zero_bytes(backend, reference_columns, cursor):
    initial = np.zeros((16, reference_columns), complex)
    initial[0, 0] = complex(.6, -0.0)
    initial[12, -1] = .8j
    initial[3, 0] = complex(-0.0, 0.0)
    initial[6, 0] = complex(0.0, -0.0)
    initial[7, 0] = complex(-0.0, -0.0)
    program = p.Program(operations=(p.Operation("smooth_read"),
        p.Operation("feedback", record=0), p.Operation("gauge_read")))
    original = p.execute(program, initial=initial, backend=backend, stop_after=cursor)
    restored = p.checkpoint_from_json(program, p.checkpoint_to_json(original), backend=backend)
    assert restored.seal == original.seal
    assert restored.initial_state.tobytes() == original.initial_state.tobytes()
    for saved, loaded in zip(original.branches, restored.branches, strict=True):
        assert loaded.state.tobytes() == saved.state.tobytes()
        assert loaded.kraus.tobytes() == saved.kraus.tobytes()
    resumed = p.execute(program, checkpoint=restored, backend=backend)
    uninterrupted = p.execute(program, initial=initial, backend=backend)
    for actual, expected in zip(resumed.branches, uninterrupted.branches, strict=True):
        np.testing.assert_array_equal(actual.state, expected.state)
