"""Adversarial admission and independent complete-channel controls."""
from collections import defaultdict
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest

from . import interpreter as p
from ..completion import observer as old
from ..completion.family import FamilyParameters
from ..local_repair_family import continuing
from ..federated_quantum_sector import construction as sector
from ..federated_observer import refinement


def fixed(repeat=1, law='mixture', seam=(0, 1)):
    return {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': list(seam)}, 'law': law, 'repeat': repeat}


def entangled(shell='H66'):
    states = p.local.basis(shell == 'H66')
    i = states.index(tuple(int(k in (0, 4)) for k in range(12)))
    j = states.index(tuple(int(k in (1, 4)) for k in range(12)))
    return {(x, a, y, b): F(1, 2) for x, a in ((i, 0), (j, 1)) for y, b in ((i, 0), (j, 1))}


def summed(blocks):
    result = defaultdict(F)
    for block in blocks:
        for key, value in block.items():
            result[key] += value
    return {key: value for key, value in result.items() if value}


@pytest.mark.parametrize('shell', ['H78', 'H66'])
@pytest.mark.parametrize('law', ['C', 'E', 'mixture'])
def test_every_basis_column_complete_and_reference_untouched(shell, law):
    value = p.program(shell=shell, commands=[fixed(law=law)])
    size = int(shell[1:])
    # One diagonal operator simultaneously tests all independently tagged
    # source columns and every reference label, not a sampled distribution.
    density = {(i, i, i, i): F(1, size) for i in range(size)}
    output = p.execute(value, density)
    assert len(output.branches) == 25
    assert p.reference_marginal(b.density for b in output.branches) == p.reference_marginal([density])
    assert sum(p.trace(b.density) for b in output.branches) == 1
    assert any(not b.density for b in output.branches)


def test_same_existing_biased_observer_complete_quantum_blocks():
    value = p.program(); initial = entangled()
    actual = p.execute(value, initial)
    expected = continuing.branches(p._parameters(value['parameters'], value['carrier'])[0])
    assert len(actual.branches) == len(expected) == 721
    for output, branch in zip(actual.branches, expected, strict=True):
        assert (output.record.seam, output.record.law, output.record.coin, output.record.raw) == (
            branch.record.seam, branch.record.law, branch.record.coin, branch.record.raw)
        assert output.density == continuing._apply(initial, branch)
    assert p.reference_marginal(b.density for b in actual.branches) == {(0, 0): F(1, 2), (1, 1): F(1, 2)}


def test_original_H78_record_only_observer_admitted_and_same_channel():
    value = p.program(shell='H78', biased=False, commands=[{
        'op': 'repair', 'selector': {'kind': 'observer', 'positive': [0, 5], 'otherwise': [0, 1]},
        'law': 'mixture', 'repeat': 1}])
    value['parameters']['repair_probability'] = '3/4'
    value['declared_resources'] = p._native_cost(value)
    initial = old.entangled_preparation()
    actual = p.execute(value, initial)
    expected = [old.apply_branch(initial, branch) for branch in old.branches((), FamilyParameters(repair_probability=F(3, 4)))]
    assert summed(b.density for b in actual.branches) == summed(expected)
    assert len(actual.branches) == 745  # Includes raw padding beyond semantic labels.
    def key(record):
        if record.kind == 'hold':
            return ('hold',)
        return ('repair', 'explore' if record.selector == 'uniform' else record.selector,
                record.seam, record.law, record.coin, record.label, record.current)
    actual_cq, expected_cq = defaultdict(F), defaultdict(F)
    for branch in actual.branches:
        for entry, coefficient in branch.density.items():
            actual_cq[key(branch.record), entry] += coefficient
    for branch in old.branches((), FamilyParameters(repair_probability=F(3, 4))):
        for entry, coefficient in old.apply_branch(initial, branch).items():
            expected_cq[key(branch.record), entry] += coefficient
    assert {k: v for k, v in actual_cq.items() if v} == {k: v for k, v in expected_cq.items() if v}


def test_FIFO_export_timeout_and_evaluator_stop_are_complete():
    initial = entangled()
    exports = p.execute(p.program(commands=[fixed(2)]), initial)
    assert exports.status == 'COMPLETE' and len(exports.branches) == 625
    assert all(len(b.archive) == len(b.memory) == 1 for b in exports.branches)
    assert sum(p.trace(b.density) for b in exports.branches) == 1
    timeout = p.execute(p.program(commands=[fixed(2)], overflow='halt'), initial)
    assert timeout.status == 'FIFO_TIMEOUT'
    assert len(timeout.branches) == 25 and all(b.record.kind == 'FIFO_TIMEOUT' for b in timeout.branches)
    stopped = p.execute(p.program(commands=[fixed(2)]), initial, max_branches=25)
    assert stopped.status == 'RESOURCE_STOP' and stopped.attempted_steps == 1
    assert sum(p.trace(b.density) for b in stopped.branches) == 1
    assert p.reference_marginal(b.density for b in stopped.branches) == p.reference_marginal([initial])


def test_feedback_reads_only_the_retained_current():
    value = p.program(commands=[fixed(), {'op': 'repair', 'selector': {
        'kind': 'last_current', 'positive': [0, 5], 'otherwise': [0, 1]}, 'law': 'mixture', 'repeat': 1}])
    output = p.execute(value, entangled())
    for branch in output.branches:
        if branch.record.kind == 'repair':
            expected = (0, 5) if (branch.archive[-1].current or 0) > 0 else (0, 1)
            assert branch.record.seam == expected


@pytest.mark.parametrize('repeat', [2, 3, 65, 100])
def test_intrinsic_FIFO_timeout_ends_execution_before_later_resource_limits(repeat):
    initial = entangled()
    value = p.program(commands=[{'op': 'hold', 'repeat': repeat}, fixed()], overflow='halt')
    output = p.execute(value, initial, max_steps=2)
    assert output.status == 'FIFO_TIMEOUT'
    assert output.attempted_steps == output.resources['attempts_executed'] == 2
    assert len(output.branches) == 1
    branch = output.branches[0]
    assert branch.record.kind == 'FIFO_TIMEOUT' and branch.halted
    assert branch.density == initial
    assert branch.memory == (p.Record('hold'),) and branch.archive == ()


@pytest.mark.parametrize('mutation', ['hidden_selector', 'future_operation', 'underreport', 'boolean_cost', 'overflow'])
def test_admission_rejects_forged_controls_or_resources(mutation):
    value = p.program(commands=[fixed()])
    if mutation == 'hidden_selector':
        value['commands'][0]['selector'] = {'kind': 'hidden_load', 'port': 0}
    elif mutation == 'future_operation':
        value['commands'][0]['op'] = 'target_hamiltonian'
    elif mutation == 'underreport':
        value['declared_resources']['native_repairs_upper'] = 0
    elif mutation == 'boolean_cost':
        value['declared_resources']['attempts_upper'] = True
    else:
        value['memory']['overflow'] = 'discard_failure'
    with pytest.raises(ValueError):
        p.admit(value)


def test_independent_native_product_preserves_entanglement_and_other_factors():
    command = {**fixed(), 'cell': 1}
    value = p.program(cells=2, commands=[command])
    state = entangled()
    # First native cell has a fixed superposition entangled with R. The
    # second receives the repair; its outcome cannot silently measure first.
    initial = {(i+66*3, a, j+66*3, b): v for (i, a, j, b), v in state.items()}
    output = p.execute(value, initial)
    reduced = defaultdict(F)
    for block in output.branches:
        for (i, a, j, b), v in block.density.items():
            if i//66 == j//66:
                reduced[i % 66, a, j % 66, b] += v
    assert dict(reduced) == state


@pytest.mark.parametrize('law', ['C', 'E'])
def test_intercell_global_sector_matches_existing_native_unit(law):
    graph = p.federation_source.federation(0)
    edge = graph.inter_seams[0]
    value = p.program(federation={'level': 0, 'total': 1}, commands=[fixed(law=law, seam=edge)])
    for i in edge:
        for j in edge:
            initial = {(1 << i, 0, 1 << j, 1): F(1)}
            actual = p.execute(value, initial)
            expected = {(a, 0, b, 1): coefficient for (a, b), coefficient in
                        sector.native_unit(law, 1 << i, 1 << j, edge).items()}
            assert summed(b.density for b in actual.branches) == expected
    population = p.execute(value, {(1 << edge[0], 0, 1 << edge[0], 0): F(1)})
    assert sum(v for b in population.branches for (i, a, j, r), v in b.density.items()
               if i == j == 1 << edge[1]) == F(1, 2)
    assert edge[0]//12 != edge[1]//12


def test_growth_is_admitted_with_explicit_unevaluated_stop():
    value = p.program(federation={'level': 3, 'total': 1})
    initial = {(1, 0, 1, 0): F(1)}
    output = p.execute(value, initial)
    assert output.status == 'RESOURCE_STOP_UNEVALUATED_CARRIER'
    assert output.branches[0].density == initial
    assert p.carrier_interface(value['carrier'])['local_total_two_assumed'] is False


def test_quantum_restriction_is_partial_trace_not_coordinate_mean():
    density = {(1, 0, 1, 0): F(1, 2), (2, 1, 2, 1): F(1, 2),
               (1, 0, 2, 1): F(1, 2), (2, 1, 1, 0): F(1, 2)}
    restricted = p.occupation_partial_trace(density, 3, (0,))
    assert restricted == {(1, 0, 1, 0): F(1, 2), (0, 1, 0, 1): F(1, 2)}
    assert p.reference_marginal([restricted]) == p.reference_marginal([density])
    assert refinement.static_certificate(0, 1)['restriction_times_lift_identity']


def test_source_word_all_complement_aborts_and_reference_blocks():
    tape = [p.compiler.pulse('X', .125, [0]), p.compiler.pulse('J', .25, [0, 1]),
            p.compiler.pulse('Z', 0., [1])]
    word = p.source_word(tape, 2)
    state = np.zeros(8, complex); state[0] = 1/np.sqrt(2); state[7] = 1j/np.sqrt(2)
    density = np.outer(state, state.conj())
    output = p.execute_word(word, density, reference_dimension=2)
    assert len(output['blocks']) == 11
    expected = np.kron(p.compiler.execute(tape, 2), np.eye(2))
    np.testing.assert_allclose(output['blocks'][('SUCCESS',)], expected@density@expected.conj().T, atol=1e-14)
    total = sum(output['blocks'].values())
    reduced = np.einsum('iaib->ab', total.reshape(4, 2, 4, 2))
    np.testing.assert_allclose(reduced, np.eye(2)/2, atol=1e-14)
    assert abs(np.trace(total)-1) < 1e-14
    assert word['declared_resources']['zero_angle_requests'] == 1


@pytest.mark.parametrize('mutation', ['future', 'angle', 'underreport', 'matrix'])
def test_rejected_source_word_retains_work_and_executes_nothing(mutation):
    word = p.source_word([p.compiler.pulse('X', .125, [0])], 1)
    if mutation == 'future':
        word['word'][0]['kind'] = 'NEW_GATE'
    elif mutation == 'angle':
        word['word'][0]['angle'] = '1/3'
    elif mutation == 'matrix':
        word['word'][0]['matrix'] = [[1, 0], [0, 1]]
    else:
        word['declared_resources']['proper_code_load_events'] = 0
    rho = np.array([[.5, .5j], [-.5j, .5]])
    output = p.submit_word(word, rho)
    assert output['status'] == 'REJECTED_PROGRAM' and output['pulses_executed'] == 0
    np.testing.assert_array_equal(output['blocks'][('REJECTED_PROGRAM',)], rho)


def test_long_literal_service_is_counted_and_explicitly_resource_stopped():
    word = p.source_word([p.compiler.pulse('Y', .125, [0])], 1, repeat=2**20)
    rho = np.eye(2)/2
    output = p.execute_word(word, rho)
    assert output['status'] == 'RESOURCE_STOP'
    np.testing.assert_array_equal(output['blocks'][('RESOURCE_STOP',)], rho)
    assert output['resources']['native_pulses'] == 2**20


def test_admission_detaches_mutable_descriptions():
    value = p.program(commands=[fixed()]); admitted = p.admit(value)
    value['commands'][0]['selector']['seam'][0] = 11
    assert admitted['commands'][0]['selector']['seam'] == [0, 1]


def test_actual_record_clock_literal_membership_and_terminal_custody():
    from . import members
    admitted = members.record_clock_member()
    assert admitted['status'] == 'ADMITTED'
    assert admitted['resources']['total_native_attempts'] == 31458883
    assert admitted['resources']['complete_supplied_realization_pulses_upper'] == 15938607203238
    assert admitted['composition_custody']['training_exports'] == 31457279
    assert admitted['composition_custody']['native_terminal_state'] == list(
        p.local.basis()[admitted['composition_custody']['native_terminal_index']])
    assert set(admitted['source_words']) == {'probe', 'analyzer', 'half_clock', 'clock'}


def test_record_clock_resource_underreport_is_rejected(monkeypatch):
    from . import members
    retained = members.checked_record_receipt()
    retained['resources']['complete_supplied_realization_pulses_upper'] -= 1
    monkeypatch.setattr(members, 'checked_record_receipt', lambda: retained)
    with pytest.raises(ArithmeticError, match='literal schedule/service counts'):
        members.record_clock_member()
