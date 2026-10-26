import numpy as np
import pytest
import sympy as s

from . import coherent_arc as c

# The compiled reflection contains a long sequence of floating point pulses.
# This numerical regression tolerance is not a service-error certificate.
SOURCE_TOLERANCE = 3e-11


def test_all_native_arcs_reverse_and_reproduce_the_exact_population_kernel():
    data = c.construction()
    reversal = np.array(data["reversal"])
    assert np.array_equal(reversal[reversal], np.arange(66*64))
    basis = c.embed(np.eye(66))
    observed = c.compress(basis[reversal])
    assert np.max(abs(observed-np.array(data["transition"], float))) < 1e-14
    assert c.construction()["transition"]*s.ones(66, 1) == s.ones(66, 1)


def test_source_coin_word_is_checked_on_every_column_including_padding():
    preparation, reflection, resources = c.source_coin()
    coin = c.ideal_coin()
    expected = 2*np.outer(coin, coin.conj())-np.eye(64)
    assert np.max(abs(preparation[:, 0]-coin)) < SOURCE_TOLERANCE
    assert np.max(abs(reflection-expected)) < SOURCE_TOLERANCE
    assert np.max(abs(reflection.conj().T@reflection-np.eye(64))) < SOURCE_TOLERANCE
    assert resources["all_coin_input_columns_executed"] == 64


def test_walk_retains_reference_and_padding_inputs_without_a_projection():
    rng = np.random.default_rng(271828)
    vector = rng.normal(size=(66*64, 2))+1j*rng.normal(size=(66*64, 2))
    vector /= np.linalg.norm(vector)
    ideal, source = c.advance(vector), c.advance(vector, source=True)
    assert np.linalg.norm(ideal-source) < SOURCE_TOLERANCE
    assert abs(np.linalg.norm(source)-1) < SOURCE_TOLERANCE
    padding = np.zeros((66*64, 1), complex)
    padding[63] = 1
    assert np.array_equal(c.advance(padding), -padding)
    assert np.linalg.norm(c.advance(padding, source=True)+padding) < SOURCE_TOLERANCE


def test_full_walk_obeys_second_order_memory_not_repeated_classical_compression():
    native = np.eye(66)[:, :3]
    initial = c.embed(native)
    first, other = c.advance(initial), c.advance(initial, phase=1)
    assert np.linalg.norm(first-other) < 1e-14
    second = c.compress(c.advance(first))
    t = np.array(c.construction()["transition"], float)
    assert np.max(abs(second-(2*t@t-np.eye(66))@native)) < 1e-13
    assert np.linalg.norm(second-t@t@native) > .01
    assert np.linalg.norm(c.advance(first)-c.advance(other, phase=1)) > .01


def test_elementary_phase_bound_accounts_for_all_workspace_modes():
    spectrum = c.spectrum()
    assert spectrum["fixed_component_sizes"] == [1, 10, 10, 45]
    total = (spectrum["cyclic_zero_phase_multiplicity"]
             + spectrum["cyclic_positive_phase_count_including_multiplicity"]
             + spectrum["cyclic_negative_phase_count_including_multiplicity"]
             + spectrum["complement_zero_phase_multiplicity"]
             + spectrum["complement_pi_phase_multiplicity"])
    assert total == spectrum["total_executed_arc_dimension"] == 4224
    assert spectrum["full_nonzero_principal_phase_spread_diagnostic"] < 396/7
    assert spectrum["cyclic_principal_phase_spread_diagnostic"] < 10


def test_complete_ramsey_outcomes_and_scope():
    report = c.report()
    for row in report["continuing_arc_diagnostics"]:
        assert abs(sum(row["Ramsey_complete_outcome_probabilities"])-1) < 2e-11
        assert row["source_vs_ideal_state_residual"] < 2e-11
        assert row["ideal_compression_recurrence_residual"] < 2e-12
    assert not report["interpretation"]["Standard_Model_charged_field_identification"]
    assert not report["interpretation"]["continuation_selected_by_committed_native_records"]
    assert not report["source_execution"]["complete_native_record_reversal_pulse_compilation"]


def test_invalid_arc_domains_rejected():
    with pytest.raises(ValueError):
        c.advance(np.zeros((2640, 1)))
    with pytest.raises(ValueError):
        c.advance(np.zeros((4224, 1)), phase=True)
    with pytest.raises(ValueError):
        c.embed(np.zeros((65, 1)))
