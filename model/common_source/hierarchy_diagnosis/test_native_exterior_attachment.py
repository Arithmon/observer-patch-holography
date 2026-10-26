"""Native 15-carrier attachment tested against the actual source transfers."""
import pytest
import sympy as s

from . import native_exterior_attachment as exterior


def test_native_exterior_carrier_has_the_right_size_but_wrong_direct_attachment():
    report = exterior.report()
    assert report["same_carrier_exterior_candidate"] == {
        "projector_rank": 15, "wedge_five_rank": 10,
        "constant_wedge_five_rank": 5, "all_projector_identities_exact": True}
    charge = report["conditional_additive_one_body_charge"]
    assert charge["native_exterior_charge_multiplicities"] != charge["decoded_M15_charge_multiplicities"]
    assert charge["single_constant_line_charge_shift_needed_for_three_and_two"] == [4, -6]
    assert not charge["one_additive_one_body_charge_action_recovers_M15"]
    dynamics = report["direct_native_transfer_attachment"]
    assert dynamics["mixed_cap_one_round_accepted_weight_exact"] == "27/32"
    assert dynamics["mixed_cap_projector_commutator"]["frobenius_norm_squared_exact"] == "961/8192"
    assert dynamics["complete_C_E_accepted_quantum_leakage_from_Q15_over_15_exact"] == {
        "C": "379/2560", "E": "703/2560"}
    assert dynamics["uniform_all_seam_control_projector_commutator"]["frobenius_norm_squared_exact"] == "29/600"
    assert not dynamics["fifteen_dimensional_candidate_invariant_under_either_transfer"]
    signed = report["declared_signed_exterior_control"]
    assert signed["fifteen_dimensional_carrier_invariant_exact"]
    assert not signed["this_CAR_lift_selected_by_recorded_native_repair"]
    rank = report["fixed_work_carrier_invariant_rank_test"]
    assert rank["finite_selector_gap_lower_exact"] == "1/800"
    assert rank["C_and_E_complete_work_channels_bistochastic_exact"]
    assert 15 not in rank["all_fixed_projector_ranks_exact"]
    assert not rank["rank_fifteen_noiseless_invariant_work_subspace_possible"]


def test_projected_core_is_an_actual_choice_not_an_automatic_two_plane(monkeypatch):
    # Antipodal ports coincide in the even five-band. The declared adjacent
    # protected core gives a two-plane; its dimension is not a free label.
    p5 = exterior._band(2)
    assert p5[0, 3] == p5[0, 0]
    monkeypatch.setattr(exterior.local_preparation, "CORE", (0, 3))
    with pytest.raises(ArithmeticError, match="two-plane"):
        exterior._candidate()


def test_direct_transfer_leaks_from_the_native_fifteen_carrier():
    _, _, _, _, q15, _, _, _ = exterior._candidate()
    selector = exterior.local_preparation.selectors(1)
    selected = exterior._accepted_transfer(selector["completion_accepted_seam_weights"])
    assert q15*q15 == q15 and s.trace(q15) == 15
    leakage = (s.eye(66)-q15)*selected*q15
    assert leakage != s.zeros(66)
    assert s.trace(leakage.T*leakage) == s.Rational(961, 16384)
