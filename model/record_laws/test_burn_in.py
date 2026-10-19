"""Independent native branch checks and premise controls for uniform burn-in."""
from collections import Counter, defaultdict
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import pytest

from . import burn_in as b
from . import family as f


@pytest.fixture(scope="module")
def result():
    return b.run()


def test_receipt_replays(result):
    retained = json.loads(Path(__file__).with_name("burn_in_results.json").read_text())
    assert json.dumps(retained, sort_keys=True) == json.dumps(result, sort_keys=True)


def _floor_images(state, seam):
    a, c = seam
    counts = Counter()
    for coin in (0, 1):
        output = list(state)
        output[a] = (state[a]+state[c]+coin)//2
        output[c] = state[a]+state[c]-output[a]
        counts[state[a]-output[a], tuple(output)] += 1
    return {j: (state, F(n, 2)) for (j, state), n in counts.items()}


def _root(p):
    return f.S(1) if p == 1 else f.S(0, F(1, 2))


def _phase(state, seam, complex_law):
    if not complex_law:
        return f.S(1)
    key = tuple(state[i] for i in seam)
    return {(0, 0): f.S(1), (1, 1): f.S(0, 0, 1),
            (1, 0): f.S(0, 0, -1), (0, 1): f.S(0, 0, -1),
            (2, 0): f.S(-1), (0, 2): f.S(-1)}[key]


def _independent_adaptive():
    # Rank-one current Kraus amplitudes are reconstructed directly from
    # floor/ceiling and local phases, without multipliers(), branch(), policy(),
    # _images(), or any sector/reduction/bound function.
    initial = b.adaptive_input()
    histories = {(): initial}
    seams = b.mixing.seams()
    for step in range(2):
        following = {}
        for prefix, state in histories.items():
            complex_law = bool(prefix and prefix[-1][1] > 0)
            for edge, seam in enumerate(seams):
                weight = F(1, 30) if step == 0 else (F(31, 60) if edge == 0 else F(1, 60))
                branches = defaultdict(lambda: defaultdict(f.S))
                for (x, y, r, s), value in state.items():
                    sx, sy = f.STATES[x], f.STATES[y]
                    left, right = _floor_images(sx, seam), _floor_images(sy, seam)
                    for j in left.keys() & right.keys():
                        ox, px = left[j]
                        oy, py = right[j]
                        ax = _root(px)*_phase(sx, seam, complex_law)
                        ay = _root(py)*_phase(sy, seam, complex_law)
                        branches[j][f.INDEX[ox], f.INDEX[oy], r, s] += weight*value*ax*ay.conj()
                for j, output in branches.items():
                    following[prefix+((edge, j),)] = dict(output)
        histories = following
    return histories


def test_every_adaptive_entangled_output_matches_independent_native_kraus_oracle(result):
    actual, independent = b.adaptive_outputs(), _independent_adaptive()
    assert actual == independent
    assert len(actual) == result["adaptive_reference_execution"]["complete_record_blocks"]
    assert b.reference_marginal(b.forget(actual)) == {(0, 0): F(1, 2), (1, 1): F(1, 2)}
    assert any(value.c or value.d for rho in actual.values() for value in rho.values())
    for rho in actual.values():
        for (x, y, r, s), value in rho.items():
            assert rho[y, x, s, r] == value.conj()


def test_adaptive_phase_feedback_is_dominated_pathwise_by_c0():
    actual = b.adaptive_outputs()
    dominant = b.execute(b.adaptive_input(), lambda *_: b.examples()["C0"], b.adaptive_schedules())
    assert actual.keys() == dominant.keys()
    differences = 0
    for prefix, rho in actual.items():
        target = dominant[prefix]
        assert rho.keys() == target.keys()
        for (x, y, r, s), value in rho.items():
            upper = target[x, y, r, s]
            assert (upper*upper.conj()-value*value.conj()).sign() >= 0
            if x == y:
                assert value == upper
            differences += value != upper
    assert differences > 0


def test_exact_dominating_channel_and_phase_examples(result):
    structural = result["structural"]
    assert structural["full_matrix_units"] == 6084
    assert structural["binary_matrix_units"] == 4356
    assert structural["maximum_offdiagonal_row_and_column_sum"] == "43/60+2/15*sqrt2"
    assert structural["attaining_signed_difference_sectors"] == 60
    assert structural["collision_survival"] == "5/6"
    assert result["coefficient_domination"]["admissible_examples"] == 7


def test_exact_budgets_cover_full_carrier_without_stationary_preparation(result):
    for tolerance, row in result["finite_budgets"].items():
        target = F(tolerance)
        n = row["binary_attempts"]
        assert b.binary_bound(n) <= F(row["binary_error_upper"]) <= target
        assert b.binary_bound(n-1) > target
        actual = b.full_bound(row["full_attempts"], row["full_even_split"])
        assert actual["error"] <= F(row["full_error_upper"]) <= target
        assert b.full_bound(row["full_attempts"]-1)["error"] > target
    assert b.binary_bound(0) == b.full_bound(0)["error"] == 1
    assert b.full_bound(300, 100)["error"] > b.full_bound(400, 100)["error"]
    assert b.binary_bound(300, F(1, 2)) > b.binary_bound(300)
    assert b.full_bound(500, coverage=F(1, 2))["error"] > b.full_bound(500)["error"]


@pytest.mark.parametrize("value", [True, -1, 1.0, "1", None])
def test_attempt_count_rejects_ambiguous_values(value):
    with pytest.raises(ValueError):
        b.binary_bound(value)
    with pytest.raises(ValueError):
        b.full_bound(value)


@pytest.mark.parametrize("value", [True, 0, -1, F(3, 2), 0.5, "1"])
def test_coverage_rejects_unsupported_guarantee(value):
    with pytest.raises(ValueError):
        b.binary_bound(10, value)


@pytest.mark.parametrize("split", [True, 3, -2, 12, 2.0])
def test_collision_split_is_exact_even_and_within_history(split):
    with pytest.raises(ValueError):
        b.full_bound(10, split)


def test_schedule_floor_is_checked_and_independence_is_not_inferred():
    assert list(map(b.coverage_of, b.adaptive_schedules())) == [1, F(1, 2)]
    assert b.coverage_of((1,)+(0,)*29) == 0
    for wrong in ((F(1, 30),)*29, (True,)+(0,)*29, (0.5, 0.5)+(0,)*28,
                  (F(-1, 30),)+(F(1, 30),)*29):
        with pytest.raises(ValueError):
            b.validate_schedule(wrong)


def test_retaining_records_and_record_selected_schedules_are_real_obstructions(result):
    control = result["controls"]
    assert control["last_current_plus_and_endpoint_occupied"] == "0"
    assert control["incorrect_source_tape_product_prediction"] == "5/396"
    assert F(control["record_selected_seam_changes_stationary_population_TV"]) > 0


def test_disconnected_schedule_keeps_coherence_and_fails_positive_coverage():
    x, y = f.representative((1, 1)), f.representative((0, 0))
    rho = {(x, y, 0, 1): f.S(1)}
    edge = b.mixing.seams().index((0, 1))
    schedule = tuple(int(i == edge) for i in range(30))
    out = b.execute(rho, lambda *_: b.examples()["C0"], [schedule]*4)
    assert b.forget(out) == rho
    assert b.coverage_of(schedule) == 0


def test_non_cp_adaptive_selection_is_rejected():
    with pytest.raises(ValueError, match="non-CP"):
        b.execute(b.adaptive_input(), lambda *_: f.generic_tables()["outside"],
                  [b.adaptive_schedules()[0]])


def test_runtime_source_closure_is_pinned():
    # A fresh interpreter avoids unrelated pytest imports entering this set.
    script = """import json,sys
from pathlib import Path
from model.record_laws import burn_in as b
b.mixing.seams()
root=b.ROOT.resolve()
paths=set()
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        p=Path(name).resolve()
        if p.is_relative_to(root) and p.suffix=='.py':
            paths.add(p.relative_to(root).as_posix())
print(json.dumps(sorted(paths-set(b.SOURCE_PATHS))))
"""
    output = subprocess.check_output([sys.executable, "-c", script], cwd=b.ROOT, text=True)
    assert json.loads(output) == []
