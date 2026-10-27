"""Exact native-object persistence, actual group action and counterfeit controls."""
from copy import deepcopy
import json

import pytest
import sympy as s

from . import native_gauge_object as n


@pytest.fixture(scope="module")
def packet():
    return n.report()


@pytest.mark.parametrize("choice", (0, 1))
def test_source_color_classes_are_disjoint_perfect_matchings(choice):
    matchings = n.equivariant_matchings(choice)
    assert len(matchings) == 6
    assert {pair for matching in matchings for pair in matching} == set(n.geometry()["nonedges"])
    assert sum(map(len, matchings)) == 36
    for matching in matchings:
        assert sorted(i for pair in matching for i in pair) == list(range(12))
    W = n.color_matrix(choice)
    assert W.T*W == 6*s.eye(6)
    for port in range(12):
        N = s.diag(*[int(port in pair) for pair in n.geometry()["pairs"]])
        assert W.T*N*W == s.eye(6)


def test_color_covariance_is_the_actual_port_permutation_not_a_copied_matrix():
    for choice in (0, 1):
        W = n.color_matrix(choice)
        for g in n.geometry()["rotations"]:
            U = n.permutation_matrix(n.pair_action(g))
            R = n.permutation_matrix(n.axis_action(g))
            assert U*W == W*R


def test_both_orientation_choices_share_only_the_uniform_state(packet):
    V0, V1 = n.matching_isometry(0), n.matching_isometry(1)
    Q = n.five_projector()
    assert n._zero(V0.T*V1-s.ones(6)/6)
    assert n._zero(Q*V0.T*V1*Q)
    assert packet["protected_object"]["orientation_relation"]["actual_inversion_exchanges_colorings_with_identity_axis_map"]
    assert not packet["protected_object"]["orientation_relation"]["orientation_exchange_is_a_physical_parity_or_CP_claim"]


def test_all_raw_CE_records_are_checked_without_success_postselection(packet):
    r = packet["protected_object"]["record_protection"]
    assert r["complete_raw_branches"] == 720
    assert r["input_coin_addresses"] == 2*30*2*66
    assert r["identically_zero_raw_branches"] > 0
    assert r["probability_counts_by_law"]["C"] == {"0": 120, "1/12": 60, "5/12": 60}
    assert r["probability_counts_by_law"]["E"] == {"0": 300, "1/12": 120, "1/3": 60}
    assert not r["unread_native_work_channel_preserves_the_code"]


def test_replacing_one_source_pair_breaks_record_protection_and_covariance():
    W = s.MutableDenseMatrix(n.color_matrix())
    pairs = n.geometry()["pairs"]
    original = n.equivariant_matchings()[0][0]
    replacement = n.geometry()["seams"][0]
    W[pairs.index(original), 0] = 0
    W[pairs.index(replacement), 0] = 1
    assert any(W.T*s.diag(*[int(port in pair) for pair in pairs])*W != s.eye(6) for port in range(12))
    assert any(n.permutation_matrix(n.pair_action(g))*W != W*n.permutation_matrix(n.axis_action(g))
               for g in n.geometry()["rotations"])


def test_two_object_geometric_sector_retains_actual_carrier_dimensions(packet):
    p = n.geometric_projectors()
    for name, rank in (("three", 3), ("three_prime", 3), ("six", 6)):
        assert n._zero(p[name]*p[name]-p[name])
        assert p[name].T == p[name]
        assert s.trace(p[name]) == rank
    assert n._zero(p["three"]*p["three_prime"])
    assert n._zero(p["six"]-p["three"]-p["three_prime"])
    q = s.kronecker_product(n.five_projector(), n.five_projector())
    assert s.trace(q) == 25
    assert n._zero(q*p["six"]-p["six"])
    assert packet["two_object_geometry"]["selected_six_is_antisymmetric_and_sum_zero_on_each_color_factor"]


def test_actual_simultaneous_rotations_intertwine_both_geometric_irreps():
    frames = n.geometric_intertwiners()
    for name, record in frames.items():
        W, metric = record["frame"], record["metric"]
        assert metric == s.Rational(1, 90)
        assert n._zero(W.T*W-metric*s.eye(3))
        assert n._zero(W*W.T/metric-n.geometric_projectors()[name])
    assert n._zero(frames["three"]["frame"].T*frames["three_prime"]["frame"])


def test_pair_carrier_is_not_silently_replaced_by_one_six_dimensional_native_block(packet):
    result = packet["two_object_geometry"]
    assert result["protected_product_dimension"] == 36
    assert result["five_product_dimension"] == 25
    assert not result["single_carrier_rank_six_limit_applies_to_this_two_carrier_subcode"]
    assert result["independent_recorded_errors_have_scalar_Gram_on_selected_six"]


def test_one_carrier_failure_follows_actual_branch_effects_and_nonzero_variances(packet):
    boundary = packet["one_carrier_boundary"]
    assert boundary["effect_span"]["all_12_port_and_30_edge_identities_exact"]
    assert boundary["effect_span"]["rank_one_native_edge_compression_vanishes_for_correctable_rank_at_least_two"]
    for name, sign in (("three", -1), ("three_prime", 1)):
        witness = boundary["vector_irreps"][name]
        assert witness["unique_nonedge_copy"]
        value = s.sympify(witness["positive_scalar_compression_defect_squared"])
        assert value == s.Rational(7, 240)+sign*s.sqrt(5)/80
        assert value > 0
    assert not boundary["single_carrier_A5_invariant_correctable_code_can_contain_geometric3_or3prime"]
    assert not boundary["multiple_corrected_carriers_excluded"]


def test_finite_gauge_response_match_does_not_import_its_generators(packet):
    attachment = packet["finite_gauge_theorem_attachment"]
    assert attachment["actual_geometric_character_matches_retained_M6"]
    assert attachment["order_character_counts"] == [
        {"order": 1, "character": "6", "count": 1},
        {"order": 2, "character": "-2", "count": 15},
        {"order": 3, "character": "0", "count": 20},
        {"order": 5, "character": "1", "count": 24}]
    assert not attachment["retained_infinitesimal_response_generators_inserted_into_native_code"]
    assert not attachment["same_continuous_port_current_law_derived_here"]
    assert not attachment["spin_half_lift_or_3_plus_2_internal_module_derived_here"]


@pytest.mark.parametrize("bad", (True, -1, 2, 0.0, "0", None))
def test_coloring_domain_rejects_unregistered_or_coerced_requests(bad):
    with pytest.raises(ValueError):
        n.equivariant_matchings(bad)


def test_permutation_and_rotation_domains_are_explicit():
    with pytest.raises(ValueError):
        n.permutation_matrix((0, 0))
    with pytest.raises(ValueError):
        n.axis_action(tuple(reversed(range(12))))


def test_scope_preserves_persistence_services_and_unproved_physical_identifications(packet):
    scope = packet["scope"]
    assert scope["same_native_persistent_object_and_geometric_M6_action_constructed"]
    assert not scope["autonomous_native_defect_creation_or_topological_protection_derived"]
    assert not scope["geometric_code5_is_the_spin_lifted_3_plus_2_by_dimension"]
    assert not scope["antisymmetric_carrier_tensor_sector_establishes_particle_statistics"]
    assert not scope["physical_spin_chirality_gauge_dynamics_or_mass_spectrum_recovered"]


def test_exact_receipt_replay_and_tamper_refusal(tmp_path, packet):
    data = deepcopy(packet); data["source_bindings"] = n.source_bindings()
    path = tmp_path/"native_object.json"
    path.write_text(json.dumps(data))
    assert n.check(path) == packet
    data["finite_gauge_theorem_attachment"]["same_continuous_port_current_law_derived_here"] = True
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="exact native"):
        n.check(path)


def test_duplicate_json_keys_refuse_before_replaying(tmp_path):
    path = tmp_path/"duplicate.json"
    path.write_text('{"schema":"oph.native_gauge_object.v1","schema":"oph.native_gauge_object.v1"}')
    with pytest.raises(ValueError, match="duplicate"):
        n.check(path)


def test_exact_receipt_comparison_rejects_Boolean_integer_substitution(packet):
    altered = deepcopy(packet)
    altered["scope"]["same_native_persistent_object_and_geometric_M6_action_constructed"] = 1
    with pytest.raises(ValueError, match="JSON type"):
        n._same(packet, altered)
