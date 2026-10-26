"""Independent native-primitive, full-reference and stopped-history checks."""
from fractions import Fraction as F
from dataclasses import replace
import hashlib
import json
import math

import numpy as np
import pytest

from ...observer import native
from ...native_modular_reference import local_preparation
from . import quantum_response as q


@pytest.mark.parametrize("law", ["C", "E"])
def test_accepted_quantum_channel_matches_independent_native_rule_on_entangled_reference(law):
    states, index = native.states(binary=True), q.INDEX
    rng = np.random.default_rng(26026)
    initial = rng.normal(size=(66, 2))+1j*rng.normal(size=(66, 2))
    initial /= np.linalg.norm(initial)
    density = np.outer(initial.ravel(), initial.ravel().conj())
    weights, timeout, survival = q.selector(1)
    expected = np.zeros((132, 132), complex)
    for seam, weight in weights.items():
        for coin in (0, 1):
            groups = {}
            for i, state in enumerate(states):
                output, (_, current) = native.native(state, seam, coin)
                label = current if law == "C" else state[seam[0]]-state[seam[1]]
                groups.setdefault(label, []).append((i, index[output]))
            for mapping in groups.values():
                output = np.zeros_like(initial)
                for i, target in mapping:
                    output[target] = math.sqrt(float(weight/2))*initial[i]
                expected += np.outer(output.ravel(), output.ravel().conj())
    actual = q.apply_channel(density, law, selector_rounds=1)
    np.testing.assert_allclose(actual, expected, atol=2e-16, rtol=0)
    assert np.trace(actual).real == pytest.approx(float(survival), abs=3e-15)
    # Both native field and arbitrary initially entangled reference are retained.
    marginal = np.trace(actual.reshape(66, 2, 66, 2), axis1=0, axis2=2)
    original = np.trace(density.reshape(66, 2, 66, 2), axis1=0, axis2=2)
    np.testing.assert_allclose(marginal, float(survival)*original, atol=2e-15)
    assert survival+timeout == 1


@pytest.mark.parametrize("law", ["C", "E"])
def test_all4356_units_have_positive_small_block_transfer_with_protected_quantum_memory(law):
    blocks, lower, _ = q.reconstruction(law)
    assert sum(len(block.indices) for block in blocks) == 4356
    assert max(len(block.indices) for block in blocks) == 66
    assert lower > F(1, 10)
    assert sum(np.count_nonzero(abs(block.values-1) < 1e-10) for block in blocks) == 6
    assert min(block.values.min() for block in blocks) > 0
    sigma = np.diag(list(map(float, local_preparation.reference())))
    np.testing.assert_allclose(q.apply_channel(sigma, law, conditional=True), sigma, atol=1e-16)
    core0, core1 = local_preparation.CORE
    address = lambda occupied: q.INDEX[tuple(int(i in occupied) for i in range(12))]
    memory = np.zeros((66, 66), complex)
    for spectator in range(12):
        if spectator not in local_preparation.CORE:
            memory[address((core0, spectator)), address((core1, spectator))] = 1
    np.testing.assert_allclose(q.apply_channel(memory, law, conditional=True), memory, atol=3e-16)


@pytest.mark.parametrize("law", ["C", "E"])
def test_relative_charge_probe_damps_without_charged_support_propagation(law):
    probe = q.coherence_probe(law)
    x, y = probe["native_matrix_unit"]
    unit = np.zeros((66, 66), complex); unit[x, y] = 1
    expected = float(F(probe["conditional_native_eigenvalue_exact"]))*unit
    actual = q.apply_channel(unit, law, conditional=True)
    np.testing.assert_allclose(actual, expected, atol=2e-16, rtol=0)
    assert np.count_nonzero(actual) == 1
    assert np.linalg.norm(actual) < np.linalg.norm(unit)
    assert probe["native_coherent_frequency_per_work_proposal"] == 0
    assert probe["relative_core_charge_change"] == 1 and probe["total_native_charge_change"] == 0
    assert probe["operator_is_even_global_transfer_not_odd_CAR"]


def test_actual_finite_selector_timeouts_terminate_and_all_quantum_outcomes_remain():
    _, timeout, survival = q.selector(1)
    assert timeout > F(1, 1000)  # Large enough to detect accidental recycling.
    for law, accepted_events in (("C", 80), ("E", 120)):
        report = q.history_diagnostic(law, steps=2, selector_rounds=1)
        assert report["paths_by_terminal_flag"] == {
            "accepted": 2*accepted_events**2, "selector_timeout": 2*(1+accepted_events)}
        assert report["probability_by_terminal_flag"]["accepted"] == pytest.approx(float(survival**2), abs=3e-15)
        assert report["probability_by_terminal_flag"]["selector_timeout"] == pytest.approx(float(1-survival**2), abs=3e-15)
        assert report["zero_probability_paths_retained"] > 0
        assert report["complete_response_and_reference_residual"] < 3e-14
        assert not report["terminal_timeouts_recycled"]
        assert not report["postselection"]
        assert report["endpoint_state_shape"] == [66, 2]
        assert report["work_resets"] == 0
    # The first terminal timeout has no fabricated later native work record.
    first = next(q.iter_history_outputs("C", steps=2, selector_rounds=1))
    assert first.records == (("selector_timeout",),)
    assert first.accepted_proposals == 0 and first.terminal_flag == "selector_timeout"


def test_exact_record_digest_is_independent_of_floating_endpoint_roundoff(monkeypatch):
    baseline = q.history_diagnostic("C", steps=1, selector_rounds=1)
    digest = hashlib.sha256()
    for path in q.iter_history_outputs("C", steps=1, selector_rounds=1):
        digest.update(json.dumps([path.records, path.terminal_flag, path.accepted_proposals, path.pointer],
                                 separators=(",", ":")).encode("utf-8"))
    assert baseline["endpoint_record_stream_sha256"] == digest.hexdigest()
    original = q.iter_history_outputs
    def rounded(*args, **kwargs):
        for path in original(*args, **kwargs):
            yield replace(path, state=path.state*np.nextafter(1., 2.))
    monkeypatch.setattr(q, "iter_history_outputs", rounded)
    altered = q.history_diagnostic("C", steps=1, selector_rounds=1)
    assert altered["endpoint_record_stream_sha256"] == baseline["endpoint_record_stream_sha256"]
    assert altered["floating_endpoint_states_recomputed_not_byte_hashed"]
    assert altered["complete_response_and_reference_residual"] < 3e-14
    assert "endpoint_state_stream_sha256" not in altered


def test_native_branch_invariant_and_fixed_load_obstruct_bare_CAR_field_identification():
    report = q.transport_obstruction()
    assert report["partial_injection_branches_checked"] == 300
    assert report["every_committed_branch_preserves_signed_difference"]
    assert not report["positive_and_negative_difference_sites_propagate"]
    assert not report["classical_record_memory_or_feedback_evades_this_invariant"]
    assert report["annihilation_addresses_leaving_carrier"] == 132
    assert report["bare_odd_port_annihilation_restriction"] == "P_N=2*c_port*P_N=2=0"
    assert not report["full12_mode_unital_CAR_representation_on_this_carrier_possible"]
    assert not report["larger_network_or_composite_field_constructions_excluded"]


def test_same_reference_modular_clock_has_only_seven_frequency_classes_and_no_transport():
    structure = q.modular_structure()
    assert structure["all_matrix_units"] == 4356
    assert structure["exact_weight_ratio_multiplicities"] == {
        "1/5": 45, "8/25": 900, "5/8": 20, "1": 2426,
        "8/5": 20, "25/8": 900, "5": 45}
    assert 3.4 < structure["positive_frequency_max_over_min"] < 3.5
    assert not structure["a_calibrated_clock_supplies_charged_field_transport"]
    sigma = np.array(list(map(float, local_preparation.reference())))
    phase = np.diag(np.exp(.37j*np.log(sigma)))
    x, y = q.coherence_probe("C")["native_matrix_unit"]
    source = np.zeros((66, 66), complex); source[x, y] = 1+.23j
    for law in ("C", "E"):
        left = q.apply_channel(phase@source@phase.conj().T, law, conditional=True)
        right = phase@q.apply_channel(source, law, conditional=True)@phase.conj().T
        np.testing.assert_allclose(left, right, atol=3e-16, rtol=0)


def test_positive_transfer_two_point_form_does_not_assert_physical_unitary_evolution():
    difference = next(key for key, indices in q._sectors().items() if len(indices) == 10)
    matrix = np.asarray(q.exact_block("E", difference), float)
    rng = np.random.default_rng(1212)
    vectors = rng.normal(size=(3, 10))+1j*rng.normal(size=(3, 10))
    gram = np.array([[np.vdot(vectors[i], np.linalg.matrix_power(matrix, i+j+1)@vectors[j])
                      for j in range(3)] for i in range(3)])
    assert np.linalg.eigvalsh(gram).min() > -1e-12
    witness = q.coherence_probe("E")
    assert 0 < F(witness["conditional_native_eigenvalue_exact"]) < 1
    assert witness["native_coherent_frequency_per_work_proposal"] == 0


def test_finite_incident_clock_retains_remote_memory_and_both_kinds_of_terminal_timeout():
    certificate = q.seam_positivity_certificate()
    assert all(row["all_seams"] == 30 for row in certificate["laws"].values())
    for law in ("C", "E"):
        for port in (5, 3):
            result = q.incident_clock_response(law, selector_rounds=1, outside_port=port, deadline=3)
            j, r = map(F, (result["incident_coherence_multiplier_exact"], result["remote_coherence_multiplier_exact"]))
            p_j, p_r = map(F, (result["incident_probability_per_proposal_exact"], result["remote_probability_per_proposal_exact"]))
            success = sum((F(row["probability_exact"]) for row in result["accepted_tick_by_remote_count"]), F())
            timeout = sum((F(row["probability_exact"]) for row in result["selector_timeout_by_remote_count"]), F())
            timeout += F(result["deadline_timeout"]["probability_exact"])
            assert success+timeout == 1
            assert F(result["conditional_tick_coherence_eigenvalue_exact"]) == (
                j*sum((r**k for k in range(3)), F())/(p_j*sum((p_r**k for k in range(3)), F())))
            assert result["new_autonomous_oscillation_frequency"] == 0
            assert not result["source_clock_reconstruction_has_Lorentzian_particle_poles"]
    # Port3 is touched by remote repairs, so waiting changes its coherence.
    response = q.incident_clock_response("E", selector_rounds=1, outside_port=3)
    assert F(response["remote_coherence_multiplier_exact"]) < F(response["remote_probability_per_proposal_exact"])


def test_complex_response_separates_actual_native_attenuation_and_the_same_modular_frequency():
    sigma = np.array(list(map(float, local_preparation.reference())))
    modular = np.diag(np.exp(1j*np.log(sigma)))
    frequencies, damping = [], []
    for law in ("C", "E"):
        result = q.complex_pole_decomposition(law, outside_port=3)
        x, y = q.coherence_probe(law, outside_port=3)["native_matrix_unit"]
        unit = np.zeros((66, 66), complex); unit[x, y] = 1
        observed = modular@q.apply_channel(unit, law, conditional=True)@modular.conj().T
        row = result["rows"][0]
        pole = complex(*row["complex_multiplier"])
        np.testing.assert_allclose(observed, pole*unit, atol=3e-16, rtol=0)
        frequencies.append([entry["modular_frequency"] for entry in result["rows"]])
        damping.append(row["damping_minus_log_modulus"])
        assert result["rows"][0]["attenuation_exact"] != result["rows"][1]["attenuation_exact"]
        assert all(abs(entry["phase_per_unit_modular_step"]-math.log(8/5)) < 2e-16 for entry in result["rows"])
    assert frequencies[0] == frequencies[1]
    assert damping[0] != damping[1]
    joint = q.joint_modular_spectrum("C")
    assert len(joint["rows"]) == 7 and sum(row["operator_modes"] for row in joint["rows"]) == 4356


def test_exact_finite_operator_log_spread_certificate_covers_every_nonfixed_mode():
    certificate = q.operator_log_spread_certificate()
    assert sorted(row["source_sector_dimension"] for row in certificate["six_fixed_sectors"]) == [1, 10, 10, 10, 10, 45]
    assert sum(row["dimension"] for row in certificate["six_fixed_sectors"]) == 80
    for row in certificate["six_fixed_sectors"]:
        assert row["exact_LDL_reconstruction_checked"]
        assert all(F(pivot) > 0 for pivot in row["positive_pivots_exact"])
    assert certificate["fixed_operator_dimension_exact"] == 6
    for law in ("C", "E"):
        row = certificate["other_sectors"][law]
        assert row["other_matrix_units"] == 4270
        assert F(row["maximum_other_row_sum_exact"]) < F(99, 100)
        assert F(row["positive_transfer_floor_exact"]) > F(7, 50)
    assert certificate["log_rate_spread_upper_strict"] == "200"
    assert F(certificate["exact_exp_two_Taylor_lower"]) > F(50, 7)
    assert not certificate["imported_exclusion_gap_theorem_required"]
    assert not certificate["existing_1_over600_source_bound_modified"]
    assert not certificate["these_operator_logs_are_physical_masses"]


@pytest.mark.parametrize("law", ["C", "E"])
def test_all_quantum_incident_clock_modes_have_positive_similarity_without_recycling_timeout(law):
    result = q.full_incident_clock_spectrum(law)
    assert result["operator_modes_retained"] == 4356
    assert result["largest_block_inverse"] == 66
    assert result["positive_similarity_residual"] < 2e-12
    assert result["fixed_modes_at_numerical_tolerance"] == 6
    assert result["finite_selector_timeout_is_terminal"]
    assert not result["remote_wait_has_finite_deadline"]
    assert not result["Lorentzian_particle_mass_identification_established"]
    assert result["exact_zero_modes_from_identically_zero_incident_blocks"] <= result["zero_transfer_modes_at_numerical_tolerance"]
    # Check the full clock spectrum includes the independently derived scalar
    # selected-coherence pole in the infinite remote-wait limit.
    probe = q.incident_clock_response(law, outside_port=5)
    j, r = F(probe["incident_coherence_multiplier_exact"]), F(probe["remote_coherence_multiplier_exact"])
    p_j = F(probe["incident_probability_per_proposal_exact"])
    _, timeout, _ = q.selector()
    scalar = float(j/(1-r)/(p_j/(p_j+timeout)))
    assert result["minimum_resolved_positive_transfer_eigenvalue"] < scalar < result["maximum_resolved_nonfixed_transfer_eigenvalue"]


def test_invalid_domains_are_rejected_without_overriding_source_protocol():
    with pytest.raises(ValueError):
        q.selector(0)
    with pytest.raises(ValueError):
        q.coherence_probe("X")
    with pytest.raises(ValueError):
        q.coherence_probe("C", outside_port=0)
    with pytest.raises(ValueError):
        q.apply_channel(np.eye(66), "C", conditional=1)
    with pytest.raises(ValueError):
        q.exact_block("C", (False,)*12)
    with pytest.raises(ValueError):
        list(q.iter_history_outputs("C", steps=3))
