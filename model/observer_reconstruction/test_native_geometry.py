"""Independent observation-only geometry identification and its boundaries."""
from copy import deepcopy
from fractions import Fraction as F
import json

import numpy as np
import pytest

from . import native_geometry as g


@pytest.fixture(scope="module")
def records():
    return g.execute_calibration()


@pytest.fixture(scope="module")
def recovered(records):
    return g.recover(records)


def test_records_contain_only_preparation_and_committed_observations(records):
    assert len(records) == 12*30*2
    assert all(set(row) == g.RECORD_FIELDS for row in records)
    for port in range(12):
        assert sum(F(row["probability"]) for row in records if row["preparation_port"] == port) == 1


def test_estimator_cannot_read_native_state_incidence_or_coordinates(records, recovered, monkeypatch):
    from ..common_source import slow_geometry
    from ..observer import native

    def forbidden(*args, **kwargs):
        raise AssertionError("estimator attempted hidden source access")

    monkeypatch.setattr(native, "native", forbidden)
    monkeypatch.setattr(slow_geometry, "source", forbidden)
    monkeypatch.setattr(slow_geometry, "raw_frame", forbidden)
    assert g.recover(records) == recovered


def test_exact_record_mean_matches_independent_actual_load_audit(recovered):
    from ..common_source import slow_geometry
    from ..observer import native

    expected = [[F(0) for _ in range(12)] for _ in range(12)]
    seams = slow_geometry.source().seams()
    for preparation in range(12):
        initial = tuple(2*int(i == preparation) for i in range(12))
        for seam in seams:
            for coin in (0, 1):
                loads, _ = native.native(initial, seam, coin)
                for i in range(12):
                    expected[i][preparation] += F(loads[i], 120)
    assert recovered["transition_matrix_exact"] == [[str(value) for value in row] for row in expected]
    assert recovered["rank_three_unit_gram_recovered"]
    assert all(recovered["band_certificate"].values())


def test_record_geometry_matches_independent_source_gram_and_distances(recovered):
    from ..common_source import slow_geometry

    expected = slow_geometry.source().intrinsic_gram_exact()
    assert recovered["gram_qsqrt5"] == [[[str(a), str(b)] for a, b in row] for row in expected]
    gram = np.asarray(recovered["gram_numeric"])
    assert np.linalg.matrix_rank(gram, tol=1e-12) == 3
    np.testing.assert_allclose(np.diag(gram), np.ones(12), atol=0, rtol=0)
    np.testing.assert_allclose(gram @ np.ones(12), 0, atol=1e-15, rtol=0)
    assert np.linalg.eigvalsh(gram)[0] > -1e-14
    assert np.min(recovered["squared_distances_numeric"]) == 0


def test_interface_relabeling_transports_estimate_and_isometry_is_unobservable(recovered):
    permutation = (3, 1, 11, 2, 9, 8, 7, 4, 5, 10, 0, 6)
    changed = g.recover(g.execute_calibration(permutation=permutation))
    for key in ("transition_matrix_exact", "laplacian_exact", "gram_qsqrt5"):
        for i in range(12):
            for j in range(12):
                assert changed[key][permutation[i]][permutation[j]] == recovered[key][i][j]
    gram = np.asarray(recovered["gram_numeric"])
    values, vectors = np.linalg.eigh(gram)
    coordinates = vectors[:, -3:] * np.sqrt(values[-3:])
    reflection = np.diag((-1, 1, 1))
    np.testing.assert_allclose((coordinates@reflection)@(coordinates@reflection).T, gram, atol=2e-15)


def test_changed_scheduler_is_not_the_same_response_geometry(records, recovered):
    seams = sorted({tuple(row["seam"]) for row in records})
    biased = {seam: F(1, 31) for seam in seams}
    biased[seams[0]] = F(2, 31)
    for weights in (biased, {seams[0]: F(1)}):
        result = g.recover(g.execute_calibration(weights))
        assert result["transition_matrix_exact"] != recovered["transition_matrix_exact"]
        assert not result["rank_three_unit_gram_recovered"]
        assert result["gram_qsqrt5"] is None


def test_plausible_record_tampering_changes_estimate_instead_of_reading_target(records, recovered):
    erased = [{**row, "current": 0} for row in records]
    result = g.recover(erased)
    assert result["transition_matrix_exact"] != recovered["transition_matrix_exact"]
    assert not result["rank_three_unit_gram_recovered"]
    assert result["laplacian_exact"] == [["0"]*12 for _ in range(12)]


@pytest.mark.parametrize("mutation", (
    lambda rows: rows.pop(),
    lambda rows: rows.append(deepcopy(rows[0])),
    lambda rows: rows[0].update(probability=1/60),
    lambda rows: rows[0].update(probability="-1/60"),
    lambda rows: rows[0].update(probability="1/0"),
    lambda rows: rows[0].update(preparation_port=True),
    lambda rows: rows[0].update(seam=[0, 0]),
    lambda rows: rows[0].update(seam=[0, 12]),
    lambda rows: rows[0].update(coin=True),
    lambda rows: rows[0].update(current=2),
    lambda rows: rows[0].update(hidden_loads=[0]*12),
))
def test_malformed_or_incomplete_record_laws_are_rejected(records, mutation):
    rows = deepcopy(records)
    mutation(rows)
    with pytest.raises(ValueError):
        g.recover(rows)


def test_input_dependent_scheduler_is_rejected(records):
    rows = deepcopy(records)
    indices = [i for i, row in enumerate(rows) if row["preparation_port"] == 0]
    for i in indices[:2]:
        rows[i]["probability"] = "1/120"
    for i in indices[2:4]:
        rows[i]["probability"] = "1/40"
    with pytest.raises(ValueError, match="same calibration scheduler"):
        g.recover(rows)


def test_report_has_complete_controls_honest_native_costs_and_json():
    result = g.report()
    assert result["all_pass"]
    assert result["source_execution"]["native_branch_calls_main"] == 720
    assert result["source_execution"]["native_branch_calls_including_controls"] == 2184
    assert result["source_execution"]["native_quantum_compiler_words_executed"] == 0
    assert not result["scope"]["geometry_is_field_grid_geometry"]
    assert not result["scope"]["physical_length_or_clock_identified"]
    json.dumps(result, allow_nan=False)
