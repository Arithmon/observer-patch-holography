"""Independent finite-history and error-transfer checks for the scalar adapter."""
from fractions import Fraction as Q
from itertools import combinations, product
import inspect
import json
from pathlib import Path

import pytest

from . import experiment as subject


def endpoint(loads, seam, coin):
    a, b = seam
    total = loads[a]+loads[b]
    out = list(loads)
    out[a] = (total+1)//2 if coin else total//2
    out[b] = total-out[a]
    return tuple(out), (coin, loads[a]-out[a])


def independent_histories(response):
    rows = {}
    for occupied in combinations(range(12), 2):
        initial = tuple(int(k in occupied) for k in range(12))
        for c1, y, c2 in product((0, 1), (-1, 1), (0, 1)):
            first, r1 = endpoint(initial, (0, 1), c1)
            final, r2 = endpoint(first, (0, 5) if y == 1 else (0, 1), c2)
            rows[initial, r1, y, r2] = ((Q(1, 2)+y*r1[1]*response)/264, final)
    return rows


def test_all_paths_and_exact_native_mean_without_producer_rule():
    for r in (Q(-1, 2), Q(-1, 5), Q(0), Q(1, 3), Q(1, 2)):
        expected = independent_histories(r)
        assert len(expected) == 528
        actual = {}
        for a, b, initial, _, final, record in subject.affine_histories():
            actual[initial, record.first, record.detector, record.second] = (a+b*r, final)
        assert actual == expected
        assert sum(w for w, _ in actual.values()) == 1
        assert min(w for w, _ in actual.values()) >= 0
        assert sum(w*key[-1][1] for key, (w, _) in actual.items()) == Q(5, 132)*r
        assert sum(w*final[0] for w, final in actual.values()) == Q(1, 6)-Q(5, 132)*r


def test_complete_history_tv_exact_not_only_marginal():
    for left, right in ((Q(-1, 2), Q(1, 2)), (Q(-1, 5), Q(-1, 6)), (Q(0), Q(0))):
        a, b = independent_histories(left), independent_histories(right)
        variation = sum(abs(a[key][0]-b[key][0]) for key in a)/2
        assert variation == Q(5, 33)*abs(left-right)


def test_committed_records_actuator_access_and_full_tape_halt():
    assert tuple(inspect.signature(subject.decoder_seam).parameters) == ("detector_bit",)
    for initial, y, coin in product(subject.native.states(), (-1, 1), (0, 1)):
        record = subject.Record((0, -1), y)
        final, retained = subject.feedback(initial, record, coin)
        expected, second = endpoint(initial, (0, 5) if y == 1 else (0, 1), coin)
        assert (final, retained.second) == (expected, second)
        assert retained.first == record.first and retained.detector == y
        assert subject.feedback(final, retained, 1-coin) == (final, retained)


def test_exact_pinned_bound_transfer_and_signed_prediction():
    parent = subject.sources()
    bounds = subject.transferred_bounds(parent)
    for output, original in zip(bounds, parent["levels"]):
        epsilon = Q(original["graph_vs_compact_continuum_error_upper"][1])
        assert Q(output["history_total_variation_upper"]) == 5*epsilon/33
        assert Q(output["native_current_error_upper"]) == 5*epsilon/132
    finest = bounds[-1]
    assert Q(finest["native_current_interval"][1]) < 0
    assert float(Q(finest["native_current_magnitude_lower"])) > .00521
    assert float(Q(finest["native_current_error_upper"])) < .0019


def test_no_empirical_input_or_field_receipt_in_execution():
    # The field module assembles/evolves the declared source independently.
    from . import field
    text = inspect.getsource(field)
    assert "read_text(" not in text and "json.load" not in text
    assert set(subject.PINS) == {subject.PARENT, "model/observer/native.py"}


def test_retained_run_predictions_controls_and_every_record():
    result = json.loads(Path(__file__).with_name("results.json").read_text())
    assert result["all_pass"] and all(result["checks"].values())
    assert len(result["executions"]) == 4
    for execution in result["executions"]:
        sampled = execution["sample"]
        assert sum(row["count"] for row in sampled["record_counts"]) == sampled["trials"]
        assert len(sampled["record_counts"]) == 24
        for value in sampled["statistics"].values():
            assert abs(value["mean"]-value["forecast"]) <= 6*value["standard_error"]+1/sampled["trials"]
        assert sampled["statistics"]["current"]["forecast"] < 0
        assert sampled["statistics"]["reversed"]["forecast"] > 0
        assert sampled["statistics"]["independent_donor"]["forecast"] == 0
        assert sampled["statistics"]["disconnected"]["forecast"] == 0
        assert len(execution["terminal_field_branches"]) == 6


def test_replay_compares_exact_claims_and_bounds_but_tolerates_float_roundoff():
    subject.compare_receipt({"q": 13, "bound": "5/33", "p": 0.4},
                            {"q": 13, "bound": "5/33", "p": 0.4+1e-12})
    for changed in ({"q": 14, "bound": "5/33", "p": .4},
                    {"q": 13, "bound": "6/33", "p": .4},
                    {"q": 13, "bound": "5/33", "p": float("nan")},
                    {"q": 13, "bound": "5/33", "p": .400001}):
        with pytest.raises(ValueError):
            subject.compare_receipt({"q": 13, "bound": "5/33", "p": .4}, changed)


def test_small_seeded_run_replays_and_zero_impulse_forecast():
    a = subject.sample(0, trials=2000, seed=91)
    assert a == subject.sample(0, trials=2000, seed=91)
    assert a["statistics"]["current"] == a["statistics"]["reversed"]
    assert a["statistics"]["current"]["forecast"] == 0
    with pytest.raises(ValueError):
        subject.sample(float("nan"))
