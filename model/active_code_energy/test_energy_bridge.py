"""Independent exact controls for the active-code energy bridge."""
from fractions import Fraction as F

import pytest

from ..observer import record_recovery as records
from ..clock_defect_energy import defect_history
from . import energy_bridge as bridge


def test_sector_reweighting_shifts_modular_contrast_without_changing_currents():
    result = bridge.sector_reweighting()
    assert result["connected_charge_sector_sizes"] == {"0": 2970, "1": 1320, "2": 66}
    rows = result["reweight_q1_relative_to_q0_by_factor_two"]
    assert rows["1"]["stationary_sector_masses"] == {"0": "2/5", "1": "5/9", "2": "2/45"}
    assert rows["2"]["stationary_sector_masses"] == {"0": "9/35", "1": "5/7", "2": "1/35"}
    assert rows["1"]["fixed_q1_conditional_state_weight"] == "1/1320"
    assert rows["2"]["fixed_q1_conditional_state_weight"] == "1/1320"
    assert all(row["conditional_q1_source_marginal_each_state"] == "1/66" for row in rows.values())
    assert all(row["nonconstant_q_potential_native_gradient_form"] == "0" for row in rows.values())
    assert result["relative_modular_gap_shift_q1_minus_q0"] == "-log(2)"


def test_native_activity_from_incidence_independent_of_bridge_kernel():
    pairs, edges, *_ = defect_history.geometry()
    degree = [0] * len(pairs)
    for i, j in edges:
        degree[i] += 1
        degree[j] += 1
    # Each eligible edge has one moving coin among 60 uniform seam/coin
    # attempts; all other complete flags are retained as holds.
    code = ((0, 2), (1, 3), (0, 3), (1, 2))
    assert [degree[pairs.index(pair)] for pair in code] == [10] * 4
    screened = (((0, 2), (1, 3), (0, 4), (1, 6)),
                ((0, 3), (1, 2), (0, 6), (1, 4)))
    assert [[degree[pairs.index(pair)] for pair in logical] for logical in screened] == [[10] * 4] * 2
    assert F(sum(degree), 66 * 60) == F(5, 33)
    q1 = [i for i, pair in enumerate(pairs) if len(set(pair) & {0, 1}) == 1]
    assert len(q1) == 20
    assert F(sum(degree[i] for i in q1), len(q1) * 60) == F(23, 150)
    work = bridge.active_code_work()
    assert work["prior_source_executed_code_nonzero_current_per_uniform_seam_attempt"] == "1/6"
    assert work["screened_successor_code"]["any_logical_input_nonzero_current_per_uniform_seam_attempt"] == "1/6"
    assert work["screened_successor_code"]["source_pulse_recovery_executed"] is False
    assert work["excess_over_stationary_q1_background_per_attempt"] == "1/66"
    assert work["excess_over_fixed_core_q1_conditional_background_per_attempt"] == "1/75"


def test_complete_record_current_and_recovery_cost_are_distinct():
    # Reconstruct the moving count directly from native seam/coin outcomes,
    # retaining the zero-current events in the denominator.
    all_screened_pairs = {pair for logical in bridge.SCREENED_SUPPORTS for pair in logical}
    for law in ("C", "E"):
        for pair in all_screened_pairs:
            state = tuple(int(p in pair) for p in range(12))
            moving = 0
            for seam in records.seams():
                for coin in (0, 1):
                    _, record = records.recorded_step(law, state, seam, coin)
                    moving += bool(records.record_current(law, record))
            assert moving == 10
    laws = bridge.active_code_work()["prior_source_executed_code_recovery_laws"]
    assert laws["C"]["mean_recovery_native_pulses_per_attempt"] == "52670"
    assert laws["E"]["mean_recovery_native_pulses_per_attempt"] == "89999/2"
    assert laws["C"]["fixed_routed_repair_native_pulses_per_attempt"] == 386739
    assert laws["E"]["fixed_routed_repair_native_pulses_per_attempt"] == 386742
    assert laws["C"]["mean_repair_plus_recovery_native_pulses_per_attempt"] == "439409"
    assert laws["E"]["mean_repair_plus_recovery_native_pulses_per_attempt"] == "863483/2"
    for row in laws.values():
        assert row["raw_abort_probability_on_code"] == "0"
        assert row["core_reset_value_change_probability"] == "1/4"
        assert sum(F(branch["probability"]) for branch in row["positive_raw_branches"]) == 1


def test_custody_rejects_a_changed_predecessor_digest(monkeypatch):
    real = bridge.digest

    def changed(path):
        return "0" * 64 if path == "model/exchange_sector/sector.py" else real(path)

    monkeypatch.setattr(bridge, "digest", changed)
    with pytest.raises(ValueError, match="predecessor source-bound input changed"):
        bridge.campaign_binding()
