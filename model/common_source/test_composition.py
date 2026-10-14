"""Independent native attachment, continued-clock and receipt mutation checks."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from itertools import combinations
import math
import sys

import numpy as np
import pytest
from scipy.linalg import expm

from . import composition as c, protocol as p
from .test_protocol_adversarial import _independent_matrices


@pytest.fixture(scope="module")
def programs():
    return {current: p.continuation_program(2, 16, native_current=current)
            for current in (-1, 0, 1)}


@pytest.fixture(scope="module")
def receipt():
    return c.run(16)


def _endpoint_rule(loads, seam, coin):
    """Independent conserved split: coin places the odd unit at one endpoint."""
    a, b = seam
    total = loads[a]+loads[b]
    result = list(loads)
    result[a] = total//2+coin*(total % 2)
    result[b] = total-result[a]
    return tuple(result), (coin, loads[a]-result[a])


def _native_oracle(histories):
    contributions = {name: [] for name in (
        "probability_mass", "mean_final_port0_load",
        "mean_final_port0_load_disconnected_actuator", "last_smooth_plus_probability",
        "clock_bit_one_probability")}
    records = set()
    paths = 0
    for pair in combinations(range(12), 2):
        original = tuple(int(index in pair) for index in range(12))
        for first_coin in (0, 1):
            retained, first = _endpoint_rule(original, (0, 1), first_coin)
            for branch in histories[first[1]].branches:
                plus = int(branch.record[2] == 0)
                for second_coin in (0, 1):
                    final, last = _endpoint_rule(retained, (0, 5) if plus else (0, 1), second_coin)
                    disconnected, _ = _endpoint_rule(retained, (0, 1), second_coin)
                    weight = float(np.vdot(branch.state, branch.state).real)/(66*2*2)
                    values = (1, final[0], disconnected[0], plus, branch.record[3])
                    for name, value in zip(contributions, values, strict=True):
                        contributions[name].append(weight*value)
                    records.add((first, branch.record, last))
                    paths += 1
    return {name: math.fsum(values) for name, values in contributions.items()}, paths, records


def test_native_attachment_matches_independent_integer_rule_and_original_branch_norms(programs):
    actual = c.native_experiment(programs)
    expected, paths, records = _native_oracle({j: c._executed(program) for j, program in programs.items()})
    assert paths == actual["complete_joint_outcome_count"] == 66*2*16*2
    assert len(records) == actual["declared_distinct_public_record_words"]
    for name, value in expected.items():
        assert actual["summary"][name] == pytest.approx(value, abs=2e-12)
    assert actual["fixed_width_public_record_storage_bits"] == 2*math.ceil(math.log2(6))+4
    assert actual["arbitrary_quantum_native_input_claimed"] is False


def test_altered_valid_native_programs_recompute_the_response(programs):
    changed = {j: replace(program, operations=tuple(
        replace(operation, enabled=False) if operation.kind == "feedback" else operation
        for operation in program.operations)) for j, program in programs.items()}
    actual = c.native_experiment(changed)
    expected, _, _ = _native_oracle({j: c._executed(program) for j, program in changed.items()})
    for name, value in expected.items():
        assert actual["summary"][name] == pytest.approx(value, abs=2e-12)
    assert abs(actual["summary"]["clock_bit_one_probability"]
               -c.native_experiment(programs)["summary"]["clock_bit_one_probability"]) > .01


def test_same_bit_count_cannot_relabel_gauge_reads_as_smooth_observations(programs):
    wrong = {j: replace(program, operations=tuple(
        p.Operation("gauge_read") if operation.kind == "smooth_read" else operation
        for operation in program.operations)) for j, program in programs.items()}
    with pytest.raises(ValueError):
        c.native_experiment(wrong)
    with pytest.raises(ValueError):
        c._clock(wrong[0])


def test_native_attachment_rejects_same_count_duplicate_branch_and_mutated_amplitude(programs, monkeypatch):
    original = {j: c._executed(program) for j, program in programs.items()}
    history = original[0]
    wrong_menus = ((history.branches[0], history.branches[0], *history.branches[2:]),
                   (replace(history.branches[0], state=2*history.branches[0].state),
                    *history.branches[1:]))
    for branches in wrong_menus:
        altered = {**original, 0: replace(history, branches=branches)}
        monkeypatch.setattr(c, "_executed", lambda program: altered[program.native_current])
        with pytest.raises(ValueError):
            c.native_experiment(programs)


def _continuous_oracle(program, reset=False):
    """d=2 original-input control, independent of runtime Hamiltonian helpers."""
    assert program.dimension == 2
    terms, q, x = _independent_matrices(2, program.coupling)
    initial = np.eye(16, dtype=complex)[:, :1]
    initial = expm(-1j*program.native_current*math.atan2(4, 3)*x)@initial
    h = sum(terms.values())
    minus, plus = expm(-.5j*q), expm(.5j*q)
    smooth = ((minus-1j*plus)/2, (minus+1j*plus)/2)
    gauge = tuple(np.diag([int(index % 2 == bit) for index in range(16)]) for bit in (0, 1))
    branches, read_count = {(): initial}, 0
    for operation in program.operations:
        updated = {}
        for record, state in branches.items():
            if operation.kind == "evolve":
                updated[record] = expm(-1j*float(operation.time)*h)@state
            elif operation.kind == "feedback":
                angle = math.atan2(4, 3) if operation.enabled and record[operation.record] else 0
                updated[record] = expm(-1j*angle*x)@state
            else:
                for bit, operator in enumerate(smooth if operation.kind == "smooth_read" else gauge):
                    output = operator@state
                    if reset and operation.kind == "smooth_read" and read_count == 0:
                        output = math.sqrt(float(np.vdot(output, output).real))*initial
                    updated[record+(bit,)] = output
        branches = updated
        read_count += operation.kind == "smooth_read"
    return branches


def test_retained_and_reset_experiments_have_independent_complete_state_oracles(programs):
    program = programs[1]
    expected = _continuous_oracle(program)
    reset = _continuous_oracle(program, reset=True)
    for actual, oracle in ((c._direct_continuous(program), expected),
                           (c._direct_continuous(program, reset_after_first_scalar=True), reset)):
        assert actual.keys() == oracle.keys()
        for record in oracle:
            np.testing.assert_allclose(actual[record], oracle[record], atol=3e-12, rtol=0)
    # Every old two-bit probability survives the intentionally destructive reset.
    for prefix in ((0, 0), (0, 1), (1, 0), (1, 1)):
        weights = [math.fsum(float(np.vdot(state, state).real)
                             for record, state in rows.items() if record[:2] == prefix)
                   for rows in (expected, reset)]
        assert weights[0] == pytest.approx(weights[1], abs=2e-12)
    assert c.complete_distance(expected, reset) > .01


def test_source_preparation_error_cannot_hide_inside_the_temporal_budget(monkeypatch):
    original = c.scalar.preparation

    def perturbed_preparation(dimension):
        state, metadata = original(dimension)
        direction = np.zeros_like(state)
        direction[-1] = 1
        direction -= state*np.vdot(state, direction)
        direction /= np.linalg.norm(direction)
        changed = state+1e-6*direction
        return changed/np.linalg.norm(changed), metadata

    c._executed.cache_clear()
    monkeypatch.setattr(c.scalar, 'preparation', perturbed_preparation)
    try:
        program = p.continuation_program(2, 16)
        source, product = c._executed(program), c._executed(program, 'product')
        comparison = p.compare_histories(source, product)
        assert comparison['maximum_kraus_operator_error'] < c.TOLERANCE
        assert comparison['complete_flagged_vector_error'] > 100*c.TOLERANCE
        # The process maps agree and the temporal allowance is much larger
        # than the defect. Independent preparation/output agreement must fail.
        assert c.complete_distance(c._states(source), c._direct_continuous(program)) < float(
            F(c.bounds.program_bound(program)['ideal_temporal_instrument_error_upper']))
        with pytest.raises(AssertionError, match='continued-source composition'):
            c._one(program)
    finally:
        c._executed.cache_clear()


def test_complete_trace_distance_preserves_branch_weights_and_reference_coherence():
    left = {(0,): np.array([math.sqrt(.3), 0.]), (1,): np.array([math.sqrt(.7), 0.])}
    right = {(0,): np.array([0., math.sqrt(.4)]), (1,): np.array([math.sqrt(.6), 0.])}
    assert c.complete_distance(left, right) == pytest.approx(.4, abs=1e-14)
    plus = np.diag([1., 1j])/math.sqrt(2)
    minus = np.diag([1., -1j])/math.sqrt(2)
    assert c.complete_distance({(): plus}, {(): minus}) == pytest.approx(1., abs=1e-14)
    assert c.complete_distance({(): plus}, {(): 1j*plus}) == pytest.approx(0., abs=1e-14)


@pytest.mark.parametrize("left,right", [
    ({}, {}), ({(): np.array([])}, {(): np.array([])}),
    ({(): np.ones((2, 2))/2}, {(): np.ones(4)/2}),
    ({(0,): np.ones(2)}, {(1,): np.ones(2)}),
    ({(): np.array([math.nan])}, {(): np.ones(1)}),
])
def test_complete_distance_refuses_empty_or_incompatible_interfaces(left, right):
    with pytest.raises(ValueError):
        c.complete_distance(left, right)


def test_clock_prediction_uses_the_actual_two_read_prefix_and_original_hamiltonian(programs):
    program = programs[0]
    rows = c._clock(program)
    prefix = c._executed(program, stop_after=len(program.operations)-2)
    assert {len(branch.record) for branch in prefix.branches} == {3}
    h = sum(_independent_matrices(2, program.coupling)[0].values())
    effect = np.diag([index % 2 for index in range(16)])
    for setting in rows["settings"]:
        time = F(setting["model_time"])
        unitary = expm(-1j*float(time)*h)
        expected = math.fsum(float(np.vdot(unitary@branch.state, effect@unitary@branch.state).real)
                             for branch in prefix.branches)
        assert setting["continuous_prediction_on_same_retained_prefix"] == pytest.approx(expected, abs=3e-12)
        assert F(setting["modular_parameter_beta_one"]) == -time
    assert rows["readout_probability_range"] > .1
    assert rows["gibbs_reference_prepared_as_actual_state"] is False


def test_finite_clock_time_budget_uses_its_actual_changed_coupling(programs):
    changed = replace(programs[0], coupling=F(3, 8))
    result = c._clock(changed)
    scale = 10**15
    floor = math.isqrt(6*scale*scale)
    root_upper = F(floor+int(floor*floor < 6*scale*scale), scale)
    norm_sum = 3+(1+2*F(9, 32))*3+abs(changed.coupling)*root_upper
    for row in result["settings"]:
        time = F(row["model_time"])
        expected = min(F(1), norm_sum**3*time**3/(3*changed.steps**2))
        assert F(row["ideal_appended_time_error_upper"]) == expected


def test_realized_path_resources_are_separate_from_expanded_tree_and_missing_services(receipt):
    for level in receipt["levels"]:
        for history in level["histories"]:
            path = history["resources_per_realized_path"]
            tree = history["expanded_tree_resources_not_physical_path_cost"]
            assert path["fresh_pointer_blanks"] == path["committed_classical_bits"] == 4
            assert tree["fresh_pointer_blanks"] == tree["committed_classical_bits"] == 1+2+4+8
            assert path["work_preparations"] == tree["work_preparations"] == 1
            assert history["number_of_retained_quantum_outcomes"] == 16
        services = level["physical_service_coverage"]
        assert services["complete_service_key_coverage"] is False
        assert services["conditional_complete_error_upper"] is None
        assert set(services["missing_service_stage_keys"]) == set(services["required_service_stage_keys"])
        assert "native_preparation_and_first_record" in services["missing_service_stage_keys"]
        assert "native_actuation_and_last_record" in services["missing_service_stage_keys"]
        assert services["ideal_instrument_error_upper"] == level[
            "ideal_error_budget_uniform_over_native_currents"][
                "complete_ideal_product_to_infinite_instrument_error_upper"]


@pytest.mark.parametrize("mutation", ("probability", "branch_omission", "decoder", "scope", "bound"))
def test_receipt_replay_rejects_changed_science_and_complete_history_fields(receipt, mutation):
    changed = deepcopy(receipt)
    history = changed["levels"][0]["histories"][0]
    if mutation == "probability":
        history["branches"][0]["probability"] += .01
    elif mutation == "branch_omission":
        history["branches"].pop()
    elif mutation == "decoder":
        changed["levels"][0]["native_observer"]["actuator_input"] = "second smooth bit without inversion"
    elif mutation == "scope":
        changed["scope"]["model_time_is_physical_duration"] = True
    else:
        history["ideal_finite_temporal_error_upper"] = "0"
    with pytest.raises(ValueError):
        c.compare_receipt(changed, receipt)


def test_replaying_mutated_native_bit_decoder_changes_the_checked_receipt(programs, monkeypatch):
    expected = c.native_experiment(programs)
    original = c.native_join.native_continuation
    monkeypatch.setattr(c.native_join, "native_continuation",
                        lambda loads, bit, coin: original(loads, 1-bit, coin))
    changed = c.native_experiment(programs)
    with pytest.raises(ValueError):
        c.compare_receipt(expected, changed)


def test_cli_check_accepts_fresh_json_and_refuses_tampered_or_duplicate_json(receipt, tmp_path, monkeypatch):
    target = tmp_path/"composition.json"
    monkeypatch.setattr(c, "run", lambda: receipt)
    monkeypatch.setattr(sys, "argv", ["composition", "--check", str(target)])
    target.write_text(c.canonical(receipt), encoding="utf-8")
    c.main()
    changed = deepcopy(receipt)
    changed["levels"][0]["histories"][0]["branches"][0]["probability"] += .01
    target.write_text(c.canonical(changed), encoding="utf-8")
    with pytest.raises(ValueError):
        c.main()
    target.write_text('{"all_pass":false,'+c.canonical(receipt).lstrip()[1:], encoding="utf-8")
    with pytest.raises(ValueError):
        c.main()
