"""Reject drift between the two issue packages before claiming a comparison."""
from copy import deepcopy
import numpy as np
import pytest

import support_comparison as comparison


@pytest.mark.parametrize("level", [3, 4, 5])
def test_support_captures_are_identical(level):
    path = (comparison.ROOT/"code/source_routing/support_w12_l3.json" if level == 3
            else comparison.HERE/f"support_w12_l{level}.json")
    packet = comparison.diagnostic.read_json(path)
    with np.load(comparison.diagnostic.HERE/f"geometry/geometry_l{level}.npz") as geometry:
        comparison.compare_support(packet, geometry, level)


@pytest.mark.parametrize("field", ["faces", "intra_carrier_seams", "glued_pairs"])
def test_changed_support_is_not_the_same_experiment(field):
    packet = comparison.diagnostic.read_json(comparison.HERE/"support_w12_l4.json")
    packet[field][0][0] += 1
    with np.load(comparison.diagnostic.HERE/"geometry/geometry_l4.npz") as geometry:
        with pytest.raises(ValueError, match="776/777 support differs"):
            comparison.compare_support(packet, geometry, 4)


@pytest.fixture
def reference():
    with np.load(comparison.diagnostic.DEFAULT/"q13_reads.npz") as data:
        return dict(data)


def test_matching_logical_projection(reference):
    comparison.compare_projection_arrays(reference, reference["indptr"], reference["indices"], reference["values"])


def test_projection_report_keeps_diagnostics_off_json_stdout(reference, capsys):
    path = comparison.ROOT/"evidence/source_net_causal_poset/routed_read_law/q13_baseline.json"
    info = comparison.diagnostic.read_json(path)["inputs"]
    result = comparison.check_q13_projection(info, reference["indptr"], reference["indices"], reference["values"])
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Verified q13 d=3" in captured.err
    assert result["logical_events_compared"] == 10985
    assert result["q13_interval_interior_flags"] == [True, True, True, False]
    assert result["canonical_only_routing_derived"] is False


@pytest.mark.parametrize("mutation, message", [
    ("offset", "read menus differ"), ("neighbour", "read menus differ"),
    ("value", "logical values differ"), ("fractional", "exact logical value format"),
    ("missing_commit", "exact logical value format"),
])
def test_projection_drift_fails(reference, mutation, message):
    offsets, indices, logical = (deepcopy(reference[name]) for name in ("indptr", "indices", "values"))
    if mutation == "offset":
        offsets[1] += 1
    elif mutation == "neighbour":
        indices[0] += 1
    elif mutation == "value":
        logical[-1, 0] += 1
    elif mutation == "fractional":
        logical = logical.astype(float)
        logical[-1, 0] += .5
    else:
        logical = logical[:-1]
    with pytest.raises(ValueError, match=message):
        comparison.compare_projection_arrays(reference, offsets, indices, logical)
