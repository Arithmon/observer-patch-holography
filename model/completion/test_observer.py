"""Independent instrument, coherent reference, memory and feedback checks."""
from collections import defaultdict
from fractions import Fraction as F
import json

import pytest

from .family import FamilyParameters
from . import observer as o


def params(**changes):
    return FamilyParameters(**dict(c_weight=F(1, 2), repair_probability=F(3, 4),
                                    exploration_probability=F(1, 2), memory_capacity=1,
                                    overflow="export_oldest") | changes)


@pytest.mark.parametrize("law", ("C", "E"))
def test_all_raw_maps_against_independent_floor_ceiling(law):
    """All source columns, all seams, both coins; off-diagonal maps follow
    from the same verified partial injections, not a diagonal-only test.
    """
    index = {state: i for i, state in enumerate(o.STATES)}
    for seam in o.SEAMS:
        expected = defaultdict(list)
        a, b = seam
        for i, state in enumerate(o.STATES):
            total = state[a]+state[b]
            for coin in (0, 1):
                target = list(state)
                low, high = total//2, (total+1)//2
                target[a], target[b] = (low, high) if coin == 0 else (high, low)
                current = state[a]-target[a]
                label = current if law == "C" else state[a]-state[b]
                expected[coin, label, current].append((i, index[tuple(target)]))
        for coin, label, current, pairs in o._raw_maps(law, seam):
            assert pairs == tuple(expected.get((coin, label, current), ()))
            assert len({i for i, _ in pairs}) == len(pairs)
            assert len({j for _, j in pairs}) == len(pairs)


@pytest.mark.parametrize("c", (F(0), F(1, 2), F(1)))
@pytest.mark.parametrize("r,e", ((F(1, 16), F(1, 16)), (F(1), F(1))))
def test_complete_instrument_at_family_boundaries(c, r, e):
    p = params(c_weight=c, repair_probability=r, exploration_probability=e)
    assert o.completeness((), p)["sum_effects"] == "I_78"
    full = (o.Record("hold"),)
    assert o.completeness(full, p)["sum_effects"] == "I_78"


def test_reference_conditioning_is_retained_but_no_unconditional_signalling():
    density = o.entangled_preparation()
    p = params()
    outputs = [(b, o.apply_branch(density, b)) for b in o.branches((), p)]
    marginal = defaultdict(F)
    for _, rho in outputs:
        o.add_into(marginal, rho)
    assert o.trace(marginal) == 1
    assert o.reference_marginal(marginal) == o.reference_marginal(density)
    branch, selected = next((b, rho) for b, rho in outputs if o.trace(rho)
                           and b.record.seam == o.native.E and b.record.current == 1)
    assert o.reference_marginal(selected) == {(0, 0): o.trace(selected)}
    assert branch.record.current == 1


def test_c_e_keep_distinct_coherence_laws_and_complex_reference_entries():
    seam = o.native.E
    # coin1 waits for contrast1 and for contrast0. C keeps this cross term;
    # E separates their contrast records and therefore destroys it.
    x = next(i for i, s in enumerate(o.STATES) if s[0] == s[4] == 1)
    y = next(i for i, s in enumerate(o.STATES) if s[4] == s[5] == 1)
    value = 1j
    operator = {(x, 7, y, 11): value}
    outputs = {}
    for law in ("C", "E"):
        result = defaultdict(complex)
        for coin, label, current, pairs in o._raw_maps(law, seam):
            if coin == 1:
                branch = o.Branch(F(1, 2), pairs, o.Record("repair"), (), None)
                o.add_into(result, o.apply_branch(operator, branch))
        outputs[law] = dict(result)
    assert outputs["C"] == {(x, 7, y, 11): value/2}
    assert outputs["E"] == {}


def test_continuing_fifo_exports_exact_history_and_never_peeks_at_source():
    p = params(memory_capacity=2)
    density = o.entangled_preparation()
    word = (o.Record("hold"),)*7
    output, memory, exported = o.history(density, word, p)
    assert o.trace(output) == F(1, 4)**7
    assert len(memory) == 2 and len(exported) == 5
    assert exported+memory == word
    assert o.proposal(memory) == o.native.E
    positive = o.Record("repair", "explore", o.native.F, "C", 0, 1, 1)
    assert o.proposal((positive,)) == o.native.F


def test_halt_on_full_reports_timeout_and_keeps_untouched_source():
    p = params(overflow="halt")
    memory = (o.Record("hold"),)
    branch, = o.branches(memory, p)
    assert branch.halted and branch.record.kind == "timeout"
    assert branch.memory == memory and branch.exported is None
    density = o.entangled_preparation()
    assert o.apply_branch(density, branch) == density
    output, after, exports = o.history(density, (o.Record("hold"), o.Record("timeout")), p)
    assert o.trace(output) == F(1, 4)
    assert after == memory and exports == ()


def test_intervention_changes_only_proposal_and_keeps_exploration_support():
    p = params()
    intervention = o.SEAMS[-1]
    ordinary = o.branches((), p)
    changed = o.branches((), p, preferred=intervention)
    assert [b for b in ordinary if b.record.selector == "explore"] == [
        b for b in changed if b.record.selector == "explore"]
    assert {b.record.seam for b in changed if b.record.selector == "feedback"} == {intervention}
    assert o.completeness((), p)["sum_effects"] == "I_78"


@pytest.mark.parametrize("second", (o.Record("hold"),
    o.Record("repair", "feedback", o.native.F, "C", 0, 0, 0)))
def test_two_step_selected_marginal_equals_explicit_prefix_sum(second):
    """Enumerate every first record for the specified second record. The
    feedback-F outcome exists only after positive first currents, directly
    exercising the grouping by the memory-selected proposal.
    """
    p = params()
    initial = o.entangled_preparation()
    blocks = o.advance_marginal({(): initial}, p)
    final = o.advance_marginal(blocks, p)
    selected_density = final[(second,)]
    direct = defaultdict(F)
    count = supported_paths = 0
    for first in o.branches((), p):
        output, _, _ = o.history(initial, (first.record, second), p)
        o.add_into(direct, output)
        if o.trace(output):
            supported_paths += 1
            if second.kind == "repair":
                assert first.record.current > 0
        count += 1
    assert count == 497
    assert 0 < supported_paths < count
    assert dict(direct) == selected_density


def test_overlapping_regional_views_preserve_shared_event_identity():
    records = (o.Record("repair", "feedback", (0, 1), "C", 0, 1, 1),
               o.Record("repair", "feedback", (0, 5), "C", 0, 0, 0),
               o.Record("hold"),
               o.Record("repair", "feedback", (0, 1), "E", 1, 1, 0))
    a, b = {0, 1, 5}, {0, 1, 11}
    left = tuple(item for item in o.regional_view(records, a) if set(item[1].seam) <= b)
    right = tuple(item for item in o.regional_view(records, b) if set(item[1].seam) <= a)
    assert left == right == ((0, records[0]), (3, records[3]))


def test_invalid_memory_and_intervention_fail_closed():
    p = params()
    for memory in ((o.Record("repair", current=99),), (o.Record("hold"),)*2,
                   [o.Record("hold")], (o.Record("hold", coin=1),)):
        with pytest.raises(ValueError):
            o.branches(memory, p)
    with pytest.raises(ValueError):
        o.branches((), p, preferred=(0, 0))


def test_retained_report_replays_exactly():
    from tools import strict_json
    assert o.report() == strict_json.load(o.HERE/"observer_results.json")
