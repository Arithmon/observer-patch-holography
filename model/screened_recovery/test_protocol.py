"""Independent gates for the screened active-service packet."""
from __future__ import annotations

import numpy as np
import pytest

from . import protocol


def test_complete_exact_service_and_emitted_layout():
    science = protocol.exact_service_certificate()
    assert science["raw_slots_checked"] == 720
    assert science["positive_slots_by_law"] == {"C": 90, "E": 118}
    assert science["uniform_seam_core_reset_probability_by_law"] == {
        "C": "1/15", "E": "1/15"}
    assert science["uniform_seam_q1_probability_after_recovery_without_label_reset"] == {
        "C": "29/30", "E": "29/30"}
    layout = protocol.compiled_layout_certificate()
    assert layout["positive_branch_inverse_schedules_checked"] == 208
    assert layout["native_branch_image_slots_checked"] == 960
    costs = protocol.resource_report()
    assert costs["encoder"]["native_pulses"] == 65484
    assert costs["laws"]["C"]["uniform_seam_expected_native_plus_recovery_pulses"] == "37499681/60"
    assert costs["laws"]["E"]["uniform_seam_expected_native_plus_recovery_pulses"] == "36870437/60"
    assert costs["laws"]["C"]["fixed_format_total_classical_bits_per_uniform_attempt"] == 23
    assert costs["laws"]["E"]["fixed_format_total_classical_bits_per_uniform_attempt"] == 24


def test_all_positive_source_recoveries_and_continuing_CE_path():
    replay = protocol.source_replay_report()
    assert replay["all_positive_raw_branch_source_recoveries"]["C"]["positive_source_recoveries_executed"] == 90
    assert replay["all_positive_raw_branch_source_recoveries"]["E"]["positive_source_recoveries_executed"] == 118
    assert replay["maximum_source_replay_error"] < 1e-8
    path = replay["corrected_two_attempt_CE_continuation"]
    assert [step["law"] for step in path] == ["C", "E"]
    assert path[1]["core_archive"]["core_after_native"] == [1, 5]
    assert path[1]["core_archive"]["core_after_recovery_and_reset"] == [0, 1]
    assert path[1]["cumulative_branch_weight_exact"] == "1/16"
    bound = protocol.source_composition_report()
    assert all(row["encoded_t_attempt_operator_error_upper"]["100"] < 1e-6
               for row in bound["laws"].values())


def test_arbitrary_reference_columns_and_zero_abort_slots():
    logical = np.array([[1, 1j, .25], [1j, -1, .75]], complex)
    decoded = protocol.decoding(protocol.encoding(logical))
    assert np.linalg.norm(decoded[:2]-logical) < 1e-9
    assert np.linalg.norm(decoded[2:]) < 1e-9
    source = np.zeros((512, 2), complex)
    source[9, 0] = 1
    source[21, 1] = 1j
    # C, coin zero, raw two at core seam has zero code probability, but is an
    # admitted raw record and receives the declared identity extension.
    assert protocol.branch_weight("C", (0, 1), 0, 2) == 0
    assert np.array_equal(protocol.apply_recovery("C", (0, 1), 0, 2, source), source)
    with pytest.raises(ValueError):
        protocol.apply_recovery("C", (0, 1), 0, 3, source)
    with pytest.raises(ValueError):
        protocol.encoding(np.zeros((3, 2)))
