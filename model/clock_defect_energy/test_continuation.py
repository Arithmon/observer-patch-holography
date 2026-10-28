"""Independent endpoint, reference, export and negative controls."""
from collections import defaultdict
from fractions import Fraction as F
import json

import numpy as np
import pytest

from model.completion import observer
from model.completion.family import FamilyParameters
from model.native_modular_reference import history_transport
from . import continuation as c


def arithmetic_step(state, seam, coin, law):
    """Primitive floor/ceiling arithmetic independent of retained map tables."""
    a, b = seam
    result = list(state)
    result[a] = (state[a]+state[b]+coin)//2
    result[b] = (state[a]+state[b]+1-coin)//2
    current = state[a]-result[a]
    label = current if law == "C" else state[a]-state[b]
    return tuple(result), label, current


@pytest.mark.parametrize("law", ["C", "E"])
def test_all_native_outcomes_against_independent_endpoint_arithmetic(law):
    params = FamilyParameters(c_weight=F(int(law == "C")))
    positions = {state: i for i, state in enumerate(c.BASIS)}
    choices = c.branches((), (0, 1), params)
    for branch in choices:
        if branch.record.kind != "repair" or branch.record.law != law:
            continue
        record = branch.record
        independent = []
        for i, state in enumerate(c.BASIS):
            output, label, current = arithmetic_step(state, record.seam, record.coin, law)
            if label == record.label:
                assert current == record.current
                independent.append((i, positions[output]))
        assert branch.pairs == tuple(independent)
    assert len(choices) == 497  # semantic zeros are present too
    assert all(sum(b.weight for b in choices if any(i == x for x, _ in b.pairs)) == 1
               for i in range(66))


def test_all_core_elementary_modular_certificate():
    result = c.elementary_certificate()
    assert result["core_locations"] == 66
    assert result["transported_columns_checked"] == 522720
    assert result["complete_native_input_columns_checked"] == 3960
    assert result["matrix_units_certified_by_column_identity"] > 10_000_000


@pytest.mark.parametrize("c_weight", [F(0), F(1), F(1, 2)])
def test_complex_entangled_reference_and_all_branch_modular_conjugations(c_weight):
    params = FamilyParameters(c_weight=c_weight, repair_probability=F(3, 4))
    core, parameter = (1, 5), .73
    rng = np.random.default_rng(1136)
    vector = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    vector /= np.linalg.norm(vector)
    rho = np.einsum("ia,jb->iajb", vector, vector.conj())
    sigma = np.array(c.reference(core), float)
    phase = np.exp(1j*parameter*np.log(sigma))
    rotated = phase[:, None, None, None]*rho*phase.conj()[None, None, :, None]
    total = np.zeros_like(rho)
    for branch in c.branches((), core, params):
        if not branch.weight or not branch.pairs:
            continue
        inputs, outputs = map(np.array, zip(*branch.pairs))
        block = float(branch.weight)*rho[inputs][:, :, inputs, :]
        wanted = float(branch.weight)*rotated[inputs][:, :, inputs, :]
        transported = np.array(c.reference(branch.core), float)[outputs]
        q = np.exp(1j*parameter*np.log(transported))
        actual = q[:, None, None, None]*block*q.conj()[None, None, :, None]
        assert np.max(abs(actual-wanted)) < 1e-17
        for a, i in enumerate(outputs):
            for b, j in enumerate(outputs):
                total[i, :, j, :] += block[a, :, b, :]
    assert np.max(abs(np.einsum("iaib->ab", total)-np.einsum("iaib->ab", rho))) < 2e-15


@pytest.mark.parametrize("law", ["C", "E"])
def test_twenty_four_attempt_replay_no_reset_and_exact_export_reconstruction(law):
    params = FamilyParameters(c_weight=F(int(law == "C")))
    state = c.BASIS[15]
    tape = []
    for n in range(24):
        edge = observer.SEAMS[(7*n+3) % 30]
        state, label, current = arithmetic_step(state, edge, n % 2, law)
        tape.append(observer.Record("repair", "explore", edge, law, n % 2, label, current))
    result = c.replay(tuple(tape), params)
    assert result["exported"]+result["memory"] == tuple(tape)
    assert len(result["exported"]) == 23
    expected = {}
    for i, initial in enumerate(c.BASIS):
        state = initial
        for record in tape:
            state, label, current = arithmetic_step(state, record.seam, record.coin, law)
            if (label, current) != (record.label, record.current):
                break
        else:
            j = c.BASIS.index(state)
            expected[j, 0, j, 0] = c.reference()[i]/120**24
            assert c.reference(result["core"])[j] == c.reference()[i]
    assert expected and result["density"] == expected
    assert observer.trace(result["density"]) < 1


def test_two_step_full_tagged_tree_equals_core_fifo_marginal_with_reference():
    params = FamilyParameters(c_weight=F(1, 2), repair_probability=F(3, 4))
    # A sparse entangled state keeps full enumeration tractable.
    rho = {(15, 0, 15, 0): F(1, 2), (51, 1, 51, 1): F(1, 2),
           (15, 0, 51, 1): F(1, 2), (51, 1, 15, 0): F(1, 2)}
    expected = {}
    for first in c.branches((), c.INITIAL_CORE, params):
        middle = c.apply_branch(rho, first)
        if not middle:
            continue
        for second in c.branches(first.memory, first.core, params):
            value = c.apply_branch(middle, second)
            if value:
                assert second.exported == first.record
                key = second.core, second.memory
                if key not in expected:
                    expected[key] = defaultdict(F)
                observer.add_into(expected[key], value)
    blocks = {(c.INITIAL_CORE, ()): rho}
    blocks = c.advance_marginal(c.advance_marginal(blocks, params), params)
    assert blocks == expected
    final = observer.summed_density(blocks)
    assert observer.trace(final) == 1
    assert observer.reference_marginal(final) == observer.reference_marginal(rho)
    assert any(key[0] not in c.compiled_cores() for key in blocks)


def test_compressed_kernel_independent_uniform_plus_feedback_arithmetic():
    params = FamilyParameters(repair_probability=F(3, 4))
    denominator, kernel = c._integer_kernel(params)
    positions = {state: i for i, state in enumerate(c.BASIS)}
    for i in (0, 15, 51, 65):
        for core in ((0, 1), (1, 5), (4, 11)):
            for bit in (0, 1):
                expected = defaultdict(F)
                expected[(i*66+c.CORE_INDEX[core])*2] += F(1, 4)
                choices = [(edge, F(1, 60)) for edge in observer.SEAMS]
                choices.append((native_edge(bit), F(1, 2)))
                for edge, selection in choices:
                    for coin in (0, 1):
                        state, _, current = arithmetic_step(c.BASIS[i], edge, coin, "C")
                        a, b = edge
                        after = (tuple(sorted(b if p == a else a if p == b else p for p in core))
                                 if current else core)
                        j = positions[state]
                        target = (j*66+c.CORE_INDEX[after])*2+int(current > 0)
                        expected[target] += F(3, 4)*selection/2
                actual = {j: F(w, denominator) for j, w in kernel[(i*66+c.CORE_INDEX[core])*2+bit]}
                assert actual == expected


def native_edge(bit):
    return (0, 5) if bit else (0, 1)


@pytest.mark.parametrize("law", ["C", "E"])
def test_last_record_and_source_cannot_replace_core_metadata(law):
    result = c.last_record_counterexample(law)
    a, b = result["histories"]
    assert a["history"][-1] == b["history"][-1]
    assert a["core"] != b["core"]
    assert (a["reference_ratio_on_shared_output_unit"], b["reference_ratio_on_shared_output_unit"]) == ("8/25", "5/8")
    assert all(F(row["actual_history_probability"]) > 0 for row in (a, b))


def test_four_predecessor_words_do_not_cover_positive_exploration():
    assert c.compiled_cores() == frozenset(((0, 1), (0, 8), (1, 5), (5, 8)))
    witness = c.unavailable_witness()
    assert witness["core"] == (1, 7)
    assert witness["status"] == "unavailable_service"
    assert F(witness["probability"]) > 0
    assert not witness["clock_action_executed"]


def test_hold_timeout_and_zero_branches_are_explicit():
    delayed = FamilyParameters(repair_probability=F(3, 4))
    hold = observer.Record("hold")
    result = c.replay((hold,)*8, delayed)
    assert result["core"] == c.INITIAL_CORE
    assert observer.trace(result["density"]) == F(1, 4**8)
    assert len(result["exported"]) == 7
    params = FamilyParameters(overflow="halt")
    first = next(b for b in c.branches((), c.INITIAL_CORE, params) if b.weight and b.pairs)
    timeout, = c.branches(first.memory, first.core, params)
    assert timeout.halted and timeout.core == first.core and timeout.weight == 1
    assert timeout.record.kind == "timeout" and timeout.exported is None
    with pytest.raises(ValueError, match="halted"):
        c.replay((first.record, timeout.record, timeout.record), params)
    zero = next(b for b in c.branches((), c.INITIAL_CORE, delayed) if not b.pairs)
    assert c.apply_branch({(0, 0, 0, 0): F(1)}, zero) == {}


@pytest.mark.parametrize("bad", [(1, 0), (0, 0), (0, 12), (False, 1), [0, 1]])
def test_invalid_core_rejected(bad):
    with pytest.raises(ValueError):
        c.reference(bad)


def test_double_occupation_and_coherent_reference_not_admitted():
    branch = c.branches((), c.INITIAL_CORE, FamilyParameters())[0]
    with pytest.raises(ValueError, match="binary66"):
        c.apply_branch({(77, 0, 77, 0): F(1)}, branch)
    # The predecessor's noncommuting control remains an actual obstruction;
    # this metadata join does not admit arbitrary non-diagonal references.
    assert history_transport.noncommuting_reference_control()["domain_commutator_nonzero"]


def test_receipt_replay():
    c.check_receipt(c.OUTPUT, c.run())


@pytest.mark.parametrize("mutation", ["duplicate_key", "boolean_for_integer", "changed_probability"])
def test_receipt_rejects_typed_and_numeric_mutations(tmp_path, mutation):
    original = c.OUTPUT.read_text()
    expected = json.loads(original)
    if mutation == "duplicate_key":
        changed = original.replace('"source_preparations": 1',
                                   '"source_preparations": 2, "source_preparations": 1')
    elif mutation == "boolean_for_integer":
        changed = original.replace('"source_preparations": 1', '"source_preparations": true')
    else:
        changed = original.replace('"probability": "53/21600"', '"probability": "1/2"')
    assert changed != original
    path = tmp_path / "mutated.json"
    path.write_text(changed)
    with pytest.raises(ValueError):
        c.check_receipt(path, expected)
