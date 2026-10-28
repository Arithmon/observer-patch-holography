"""Independent complete-record and corrected-cycle controls."""
from collections import Counter
from fractions import Fraction as F

import pytest

from ..exchange_sector import sector
from ..observer import record_recovery as records
from ..particle_screening import screening
from ..screened_recovery import protocol
from . import transport


def test_direct_current_counts_and_mode_second_moment_for_both_laws():
    columns = screening.sharp_code((0, 1))
    profile = transport.PROFILE
    for law in ("C", "E"):
        for pairs in columns:
            nonzero = F()
            squared_mode_change = F()
            for pair in pairs:
                x = tuple(int(port in pair) for port in range(12))
                old_mode = sum(profile[p]*x[p] for p in range(12))
                for seam in records.seams():
                    for coin in (0, 1):
                        y, record = records.recorded_step(law, x, seam, coin)
                        current = records.record_current(law, record)
                        assert tuple(y[p]-x[p] for p in range(12)) == tuple(
                            current*(int(p == seam[1])-int(p == seam[0])) for p in range(12))
                        nonzero += F(bool(current), len(pairs)*60)
                        new_mode = sum(profile[p]*y[p] for p in range(12))
                        squared_mode_change += F((new_mode-old_mode)**2, len(pairs)*60)
            assert nonzero == F(1, 6)
            assert squared_mode_change == F(13, 60)


def test_mode_capacity_depends_on_logical_input_while_current_form_does_not():
    columns = screening.sharp_code((0, 1))
    values = [[sum(transport.PROFILE[p] for p in pair) for pair in col]
              for col in columns]
    assert values == [[0, 0, 1, -1], [2, -2, 1, -1]]
    variance = [F(sum(x*x for x in row), len(row)) for row in values]
    assert variance == [F(1, 2), F(5, 2)]
    report = transport.current_response()
    assert report["mode_current_form"] == "13/120"
    assert [row["current_form_over_variance"] for row in report["mode_logical_columns"]] == [
        "13/60", "13/300"]
    assert report["current_form_rank"] == 11
    assert report["nonzero_current_probability_per_uniform_native_attempt"] == "1/6"
    assert not report["prior_profile_has_logical_input_independent_wave_ratio"]


def test_pre_reset_core_distribution_from_source_columns_without_protocol_probabilities():
    code = screening.sharp_code((0, 1))
    expected = {(0, 1): F(14, 15), (0, 5): F(1, 120), (0, 7): F(1, 120),
                (0, 8): F(1, 120), (0, 9): F(1, 120), (1, 5): F(1, 120),
                (1, 7): F(1, 120), (1, 10): F(1, 120), (1, 11): F(1, 120)}
    for law in ("C", "E"):
        for col in code:
            distribution = Counter()
            for pair in col:
                x = tuple(int(port in pair) for port in range(12))
                for seam in records.seams():
                    for coin in (0, 1):
                        _, record = records.recorded_step(law, x, seam, coin)
                        distribution[sector.moved_core((0, 1), law, record)] += F(1, len(col)*60)
            assert dict(distribution) == expected
    cycle = transport.corrected_cycle()
    assert all(row["ideal_corrected_logical_superoperator_eigenvalues"] == ["1"]*4
               and row["ideal_corrected_logical_generator_eigenvalues_per_attempt"] == ["0"]*4
               and row["post_reset_core_distribution"] == {"(0, 1)": "1"}
               for row in cycle["laws"].values())
    assert cycle["unreset_fixed_code_q1_probability"] == "29/30"
    assert not cycle["common_logical_phase_or_rest_gap_fixed_by_channel"]
    assert cycle["laws"]["C"]["same_code_uniform_seam_expected_source_pulses_from_recovery_receipt"] == "37499681/60"
    assert cycle["laws"]["E"]["same_code_uniform_seam_expected_source_pulses_from_recovery_receipt"] == "36870437/60"


def test_all_complete_flag_weights_are_normalized_and_source_positive_counts_match():
    for law, count in (("C", 90), ("E", 118)):
        positive = 0
        for seam in records.seams():
            weights = [protocol.branch_weight(law, seam, coin, raw)
                       for coin in (0, 1)
                       for raw in range(1 << (2 if law == "C" else 3))]
            assert sum(weights) == 1
            positive += sum(p > 0 for p in weights)
        assert positive == count
    source = transport.corrected_cycle()
    assert [source["laws"][law]["positive_complete_raw_flags"] for law in ("C", "E")] == [90, 118]


def test_campaign_fails_when_owning_recovery_source_drift_is_detected(monkeypatch):
    real = transport.digest
    monkeypatch.setattr(transport, "digest", lambda path:
                        "0"*64 if path == "model/screened_recovery/protocol.py" else real(path))
    with pytest.raises(ValueError, match="predecessor custody changed"):
        transport.campaign_binding()
