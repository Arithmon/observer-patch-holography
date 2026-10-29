"""Complete exact Choi comparisons and adversarial host-domain controls."""
from fractions import Fraction as F

import pytest

from .compiler import apply_history, compile_history
from ..native_program import interpreter as native
from ..federated_observer import source as federation


def fixed(repeat=1, *, cell=None):
    return {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]},
            'law': 'mixture', 'repeat': repeat,
            **({'cell': cell} if cell is not None else {})}


def compare(compiled, operator):
    import json
    expected = native.execute(json.loads(compiled.program_json), operator,
        max_branches=compiled.effective_max_branches, max_steps=compiled.max_steps)
    actual = apply_history(compiled, operator)
    assert actual == expected
    return actual


@pytest.mark.parametrize('shell', ['H66', 'H78'])
def test_full_choi_with_record_feedback(shell):
    commands = [fixed(), {'op': 'repair', 'selector': {'kind': 'last_current',
        'positive': [1, 7], 'otherwise': [0, 1]}, 'law': 'mixture', 'repeat': 1}]
    program = native.program(shell=shell, commands=commands)
    compiled = compile_history(program)
    dimension = int(shell[1:])
    choi = {(i, i, j, j): F(1, dimension)
            for i in range(dimension) for j in range(dimension)}
    actual = compare(compiled, choi)
    assert len(compiled.branches) == 625
    assert any(not branch.source_map for branch in compiled.branches)
    assert compiled.map_entries < sum(len(branch.density) for branch in actual.branches)
    if shell == 'H66':
        assert compiled.map_entries == 1056
        assert sum(len(branch.density) for branch in actual.branches) == 38788


def test_complex_reference_parts_and_holds_are_preserved():
    commands = [fixed(), {'op': 'hold', 'repeat': 1}, fixed()]
    compiled = compile_history(native.program(commands=commands, capacity=1))
    # A Bell state with relative phase i, split into exact rational parts.
    real = {(0, 0, 0, 0): F(1, 2), (1, 1, 1, 1): F(1, 2)}
    imaginary = {(0, 0, 1, 1): F(-1, 2), (1, 1, 0, 0): F(1, 2)}
    for operator in (real, imaginary):
        actual = compare(compiled, operator)
        assert native.reference_marginal(b.density for b in actual.branches) == native.reference_marginal([operator])
    assert all(len(branch.archive) == 2 for branch in compiled.branches)


@pytest.mark.parametrize('options,status,attempts', [
    ({'max_steps': 0}, 'RESOURCE_STOP', 0),
    ({'max_branches': 24}, 'RESOURCE_STOP', 0),
    ({'max_branches': 25}, 'RESOURCE_STOP', 1),
    ({'max_map_entries': 66*24}, 'RESOURCE_STOP', 0),
])
def test_host_stops_reuse_the_shared_identity_semantics(options, status, attempts):
    compiled = compile_history(native.program(commands=[fixed(3)]), **options)
    actual = compare(compiled, {(0, 2, 1, 4): F(-3, 7)})
    assert actual.status == status and actual.attempted_steps == attempts
    assert compiled.map_entries <= compiled.map_entry_budget


def test_intrinsic_timeout_finishes_before_later_literal_commands():
    compiled = compile_history(native.program(commands=[fixed(100)], overflow='halt'))
    actual = compare(compiled, {(0, 2, 1, 4): F(-3, 7)})
    assert actual.status == 'FIFO_TIMEOUT' and actual.attempted_steps == 2
    assert all(branch.record.kind == 'FIFO_TIMEOUT' for branch in actual.branches)


def test_product_entangled_spectator():
    compiled = compile_history(native.program(cells=2, commands=[fixed(cell=1)]))
    first, second = 0+66*1, 2+66*3
    operator = {(i, a, j, b): F(1, 2)
        for i, a in ((first, 0), (second, 1)) for j, b in ((first, 0), (second, 1))}
    actual = compare(compiled, operator)
    assert compiled.source_columns == 66**2
    assert native.reference_marginal(b.density for b in actual.branches) == native.reference_marginal([operator])


@pytest.mark.parametrize('total', [0, 1, 239, 240])
def test_federation_complete_sector_and_intercell_seam(total):
    edge = federation.federation(0).inter_seams[0]
    command = fixed()
    command['selector']['seam'] = list(edge)
    compiled = compile_history(native.program(federation={'level': 0, 'total': total}, commands=[command]))
    full = (1 << 240)-1
    addresses = [0] if total == 0 else [full] if total == 240 else (
        [1 << edge[0], 1 << edge[1]] if total == 1 else
        [full ^ (1 << edge[0]), full ^ (1 << edge[1])])
    operator = {(i, a, j, b): F(1, len(addresses))
        for a, i in enumerate(addresses) for b, j in enumerate(addresses)}
    compare(compiled, operator)
    assert compiled.source_columns == (1 if total in (0, 240) else 240)


@pytest.mark.parametrize('program', [
    native.program(cells=3),
    native.program(federation={'level': 0, 'total': 2}),
    native.program(federation={'level': 3, 'total': 0}),
])
def test_compilation_bound_does_not_change_native_admission(program):
    assert native.admit(program) == program
    with pytest.raises(ValueError):
        compile_history(program)


@pytest.mark.parametrize('options', [
    {'max_columns': True}, {'max_columns': 65}, {'max_map_entries': 65},
    {'max_map_entries': 2**24+1}, {'max_steps': 65}, {'max_branches': 0},
])
def test_host_limits_are_validated_before_compilation(options):
    with pytest.raises(ValueError):
        compile_history(native.program(), **options)


def test_program_and_returned_mutations_do_not_change_compiled_results():
    program = native.program(commands=[fixed()])
    compiled = compile_history(program)
    operator = {(0, 0, 0, 0): F(1)}
    expected = apply_history(compiled, operator)
    actual = apply_history(compiled, operator)
    program['commands'][0]['repeat'] = 100
    actual.resources.clear()
    actual.branches[0].density.clear()
    # Even an explicit dataclass bypass affects only the detached output.
    object.__setattr__(actual.branches[0].record, 'kind', 'forged')
    assert apply_history(compiled, operator) == expected
    assert type(compiled.branches[0].source_map) is tuple
    assert type(compiled.branches[0].weight) is tuple


def test_apply_retains_native_operator_domain_validation():
    compiled = compile_history(native.program(commands=[fixed()]))
    for operator in ({(-1, 0, -1, 0): F(1)}, {(66, 0, 66, 0): F(0)},
                     {(0, True, 0, 0): F(1)}, {(0, 0, 0, 0): float('nan')}):
        with pytest.raises(ValueError):
            apply_history(compiled, operator)
