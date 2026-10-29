"""Compile native histories to exact scalar-weighted partial injections.

This is a host evaluator optimization, not a source primitive or a way to
derive quantum action from the diagonal of an arbitrary channel. The native
raw C/E branches are positive scalars times partial injections. Record-only
control, holds and terminal identities preserve that form under composition.
The shared interpreter evaluates tagged diagonal matrix units to recover the
map and squared scalar for every complete history, including zero histories.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
import json

from ..native_program import interpreter as native


MAX_COMPILED_COLUMNS = 8192
MAX_MAP_ENTRY_BUDGET = 2**24
DEFAULT_MAP_ENTRY_BUDGET = 2**20


@dataclass(frozen=True, slots=True)
class CompiledBranch:
    record: native.Record
    memory: tuple[native.Record, ...]
    archive: tuple[native.Record, ...]
    halted: bool
    source_map: tuple[tuple[int, int], ...]
    # Keep rational payloads in immutable integer pairs, not caller-mutable
    # Fraction slots. No compiled object is placed in an authoritative cache.
    weight: tuple[int, int]


@dataclass(frozen=True, slots=True)
class CompiledHistory:
    program_json: str
    program_digest: str
    branches: tuple[CompiledBranch, ...]
    attempted_steps: int
    status: str
    resources: tuple[tuple[str, int], ...]
    source_columns: int
    map_entries: int
    effective_max_branches: int
    max_steps: int
    map_entry_budget: int


def _source_addresses(carrier, max_columns):
    """Enumerate the whole declared carrier only after a bounded size check."""
    kind = carrier['kind']
    if kind in ('local', 'independent_native_product'):
        dimension = int(carrier['shell'][1:])
        cells = carrier.get('cells', 1)
        size = 1
        # Stop before taking an enormous power or iterating an enormous cells
        # count. Every factor has dimension at least 66.
        for _ in range(cells):
            if size > max_columns//dimension:
                raise ValueError('carrier exceeds the host compilation column limit')
            size *= dimension
        return tuple(range(size))
    if carrier['level'] > native.MAX_EVALUATED_LEVEL:
        raise ValueError('host compilation requires an evaluated federation level')
    ports = 240 << (2*carrier['level'])
    total = carrier['total']
    smaller = min(total, ports-total)
    size = 1
    for j in range(1, smaller+1):
        size = size*(ports-j+1)//j
        if size > max_columns:
            raise ValueError('carrier exceeds the host compilation column limit')
    # Near full occupation, enumerate the missing ports instead. This avoids
    # materializing thousands of occupied labels per basis mask.
    full = (1 << ports)-1
    holes = total > ports//2
    return tuple((full ^ mask) if holes else mask for mask in (
        sum(1 << port for port in subset)
        for subset in combinations(range(ports), smaller)))


def compile_history(program, *, max_columns=MAX_COMPILED_COLUMNS,
                    max_map_entries=DEFAULT_MAP_ENTRY_BUDGET,
                    max_branches=native.MAX_EVALUATED_BRANCHES,
                    max_steps=native.MAX_EVALUATED_STEPS):
    """Compile all retained native histories within explicit host limits.

    Oversized carriers raise ValueError without changing grammar admission.
    The worst-case map budget lowers the interpreter's branch cap to
    min(max_branches, max_map_entries // source_columns). Its retained
    RESOURCE_STOP is therefore the ordinary shared evaluator stop, and the
    returned effective cap makes a comparison with direct execution explicit.
    This does not count compilation as a physical source execution.
    """
    value = native.admit(program)
    native._integer(max_columns, 1, MAX_COMPILED_COLUMNS)
    native._integer(max_map_entries, 1, MAX_MAP_ENTRY_BUDGET)
    native._integer(max_branches, 1, native.MAX_EVALUATED_BRANCHES)
    native._integer(max_steps, 0, native.MAX_EVALUATED_STEPS)
    addresses = _source_addresses(value['carrier'], max_columns)
    if len(addresses) > max_map_entries:
        raise ValueError('host map budget cannot hold one complete carrier map')
    effective = min(max_branches, max_map_entries//len(addresses))
    # Each independent reference tag identifies one input source column. This
    # unnormalized operator is algebraic bookkeeping, not a prepared density.
    tagged = {(i, i, i, i): F(1) for i in addresses}
    history = native.execute(value, tagged, max_branches=effective, max_steps=max_steps)
    branches = []
    for branch in history.branches:
        partial, weights = [], set()
        for (output, source, output_again, source_again), weight in branch.density.items():
            if output != output_again or source != source_again:
                raise ArithmeticError('native branch changed its independent source tag')
            partial.append((source, output))
            weights.add(weight)
        if (len(weights) > 1 or len({source for source, _ in partial}) != len(partial)
                or len({out for _, out in partial}) != len(partial)):
            raise ArithmeticError('native history is not a scalar-weighted partial injection')
        weight = next(iter(weights), F(0))
        if weight < 0:
            raise ArithmeticError('native history has a negative channel weight')
        branches.append(CompiledBranch(deepcopy(branch.record), deepcopy(branch.memory),
            deepcopy(branch.archive), branch.halted, tuple(sorted(partial)),
            (weight.numerator, weight.denominator)))
    return CompiledHistory(native.canonical(value), history.program_digest, tuple(branches),
        history.attempted_steps, history.status, tuple(sorted(history.resources.items())),
        len(addresses), sum(len(branch.source_map) for branch in branches),
        effective, max_steps, max_map_entries)


def apply_history(compiled, operator):
    """Apply a compiled native history to any exact source/reference operator.

    As in native.execute, complex operators use separate real and imaginary
    calls. The returned History owns fresh density, resource and record data.
    There is no shared compilation cache or serialized admission authority.
    """
    if type(compiled) is not CompiledHistory:
        raise ValueError('a compiled native history required')
    program = native.admit(json.loads(compiled.program_json))
    density = native._validate_density(operator, program['carrier'])
    branches = []
    for branch in compiled.branches:
        partial = dict(branch.source_map)
        weight = F(*branch.weight)
        output = {(partial[i], a, partial[j], b): weight*coefficient
                  for (i, a, j, b), coefficient in density.items()
                  if weight and i in partial and j in partial}
        branches.append(native.Outcome(deepcopy(branch.record), output,
            deepcopy(branch.memory), deepcopy(branch.archive), branch.halted))
    return native.History(tuple(branches), compiled.attempted_steps,
                          compiled.program_digest, compiled.status, dict(compiled.resources))
