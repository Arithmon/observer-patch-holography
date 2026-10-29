from dataclasses import replace
from fractions import Fraction as F

import pytest

from .certificate import (clock_history_control, comparison_graph,
                          field_reference_control, native_archive_certificate,
                          native_sector_certificate, reconstruct_offsets)


def test_disconnected_offsets_are_not_fixed_by_choosing_component_zeros():
    result = reconstruct_offsets(4, [(0, 1, F(3, 2)), (2, 3, F(-7, 4))])
    assert result["components"] == [[0, 1], [2, 3]]
    assert result["offsets_with_one_zero_per_component"] == (0, F(3, 2), 0, F(-7, 4))
    assert result["relative_offsets_unidentified"] == 1
    # One admitted cross-component comparison fixes that relative zero.
    joined = reconstruct_offsets(4, [(0, 1, F(3, 2)), (2, 3, F(-7, 4)),
                                     (1, 2, F(5, 2))])
    assert joined["offsets_with_one_zero_per_component"] == (0, F(3, 2), 4, F(9, 4))
    assert joined["relative_offsets_unidentified"] == 0


def test_cycle_consistency_is_a_new_constraint_not_an_extra_offset():
    good = [(0, 1, 2), (1, 2, 3), (2, 0, -5)]
    assert reconstruct_offsets(3, good)["incidence_rank"] == 2
    with pytest.raises(ValueError, match="inconsistent energy-offset cycle"):
        reconstruct_offsets(3, [*good[:-1], (2, 0, -4)])
    with pytest.raises(ValueError, match="inconsistent energy-offset cycle"):
        reconstruct_offsets(2, [(0, 1, 1), (1, 0, 1)])


@pytest.mark.parametrize("vertices,edges", [(True, []), (0, []),
                                           (2, [(0, 0)]), (2, [(0, 2)]),
                                           (2, [(False, 1)])])
def test_graph_rejects_invalid_declared_access(vertices, edges):
    with pytest.raises(ValueError):
        comparison_graph(vertices, edges)


def test_duplicate_access_edges_do_not_inflate_rank():
    graph = comparison_graph(64, [(0, i) for i in range(1, 64)]*2)
    assert graph["edges"] == graph["incidence_rank"] == 63
    assert graph["offset_gauge_dimension"] == 1
    assert comparison_graph(64, [])["relative_offsets_unidentified"] == 63


@pytest.mark.parametrize("difference", [True, 0.5, float("nan"), "1/2"])
def test_offset_certificate_requires_exact_differences(difference):
    with pytest.raises(ValueError, match="exact rational"):
        reconstruct_offsets(2, [(0, 1, difference)])


def test_separate_clock_keeps_gauge_but_coherent_sector_read_exposes_it():
    control = clock_history_control()
    assert control["complete_CQ_states_exactly_equal"]
    assert control["retained_histories"] == 8
    assert control["Ramsey_plus_probability_without_offset"] == "1"
    assert control["Ramsey_plus_probability_with_pi_over_2_phase"] == "1/2"


def test_actual_native_branches_obey_operator_sector_condition():
    result = native_sector_certificate()
    assert result["source_core_dimension"] == 4356
    assert result["matrix_entry_checks"] == {"C": 261360, "E": 261360}
    assert not result["all_source_hardware_interventions_are_sector_preserving"]


def test_one_charge_mixing_raw_entry_invalidates_the_claim(monkeypatch):
    from ..completion import observer
    from ..exchange_sector import sector

    x, y = sector.BASIS[:2]
    source, target = observer.STATES.index(x), observer.STATES.index(y)
    monkeypatch.setattr(observer, "_raw_maps", lambda law, seam: [(0, 0, 0, [(source, target)])])
    with pytest.raises(ArithmeticError, match="mixes charge sectors"):
        native_sector_certificate()


def test_zero_modular_energy_extends_to_full_native_archive_reference():
    result = native_archive_certificate()
    assert result["complete_raw_slots"] == 720
    assert result["initial_off_diagonal_difference_nonzero"]
    assert result["all_lengths_source_and_complete_archive_modular_energy_excess"] == "0"
    assert not result["reference_on_growing_archive_is_stationary"]


def test_a_coherent_merger_breaks_the_diagonal_history_induction(monkeypatch):
    from ..observer import adaptive_record_recovery as adaptive

    branches = adaptive.branch_maps(binary=True)
    branch = next(b for b in branches if len(b.inputs) >= 2)
    broken = replace(branch, outputs=(branch.outputs[0],)*len(branch.outputs))
    monkeypatch.setattr(adaptive, "branch_maps", lambda **kwargs: (broken,))
    with pytest.raises(ArithmeticError, match="partial injection"):
        native_archive_certificate()


def test_all_actual_field_fibers_need_relative_normalization_comparisons():
    result = field_reference_control()
    assert result["field_dimension"] == 8128
    assert result["momentum_fibers"] == 64
    assert result["additional_comparison_rank_for_one_global_offset"] == 63
    assert result["partition_log_spread_float"] > .01
    assert result["tree_reconstruction_residual_float"] < 1e-12
    assert not result["tree_matrix_units_source_accessibility_certified"]
