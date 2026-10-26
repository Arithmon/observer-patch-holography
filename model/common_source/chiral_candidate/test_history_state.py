"""Independent STATE-level join of full Gaussian field and literal Fock factor.

This reconstructs conditional density operators from every retained clock
and purification output. Weight checks remain separate: conditional-state
normalization in this comparison does not filter the continuing history.
"""
import math

import numpy as np
import pytest

from . import gaussian_source, history


@pytest.fixture(scope="module")
def complete_history():
    return history.execute()


def _six_mode_commutators():
    """Dense CAR operators reconstructed independently from occupation masks."""
    majoranas = []
    for mode in range(6):
        annihilator = np.zeros((64, 64), complex)
        for mask in range(64):
            if mask & (1 << mode):
                remaining = mask ^ (1 << mode)
                sign = (-1)**((mask & ((1 << mode)-1)).bit_count())
                annihilator[remaining, mask] = sign
        majoranas.extend((annihilator+annihilator.conj().T,
                          -1j*(annihilator-annihilator.conj().T)))
    return np.array([[.5j*(left@right-right@left) for right in majoranas]
                     for left in majoranas])


def test_all_sixteen_full_field_conditionals_match_literal_reduced_quantum_states(complete_history):
    report, arrays = complete_history
    assert arrays["gaussian_covariances"].shape == (16, 180, 180)
    assert arrays["probe_outputs"].shape == (64, 256, 4)
    assert set(arrays["field_records"]) == {format(i, "04b") for i in range(16)}
    commutators = _six_mode_commutators()
    odd = np.array([mask.bit_count() % 2 for mask in range(64)], bool)
    for record, weight, full_gamma in zip(arrays["field_records"], arrays["gaussian_weights"],
                                          arrays["gaussian_covariances"], strict=True):
        density = np.zeros((64, 64), complex)
        outcomes = 0
        for final, output in zip(arrays["records"], arrays["probe_outputs"], strict=True):
            if str(final).startswith(str(record)):
                # Source bits0..5 are the field factor, bits6..7 its clocks;
                # the four untouched purifier columns remain quantum outputs.
                field_columns = output.reshape(4, 64, 4).transpose(1, 0, 2).reshape(64, 16)
                density += field_columns@field_columns.conj().T
                outcomes += 1
        assert outcomes == 4  # Both outcomes of each clock, without selection.
        literal_weight = float(np.trace(density).real)
        assert weight > 0 and literal_weight > 0
        assert abs(literal_weight-weight) < 2e-12
        actual_gamma = np.einsum("ij,abji->ab", density/literal_weight, commutators).real
        np.testing.assert_allclose(actual_gamma, full_gamma[:12, :12], atol=2e-10, rtol=0)
        # Compare with the independently propagated Gaussian branch weight too.
        weighted_gamma = np.einsum("ij,abji->ab", density/weight, commutators).real
        np.testing.assert_allclose(weighted_gamma, full_gamma[:12, :12], atol=2e-10, rtol=0)
        assert float(np.trace(density[np.ix_(odd, odd)]).real) < 2e-26
        # The six neutrino zero modes remain in the FULL Gaussian state.
        np.testing.assert_allclose(full_gamma[-12:, -12:], gaussian_source.vacuum_covariance(6),
                                   atol=2e-13, rtol=0)
        np.testing.assert_allclose(full_gamma[-12:, :-12], 0., atol=2e-13, rtol=0)
    assert report["work_preparations"] == 1 and report["work_resets"] == 0


def test_both_quantum_clock_purifiers_survive_the_complete_record_tree(complete_history):
    report, arrays = complete_history
    assert set(arrays["records"]) == {format(i, "06b") for i in range(64)}
    reference = sum(output.conj().T@output for output in arrays["probe_outputs"])
    np.testing.assert_allclose(reference, np.diag([1, 5, 5, 25])/36, atol=2e-11, rtol=0)
    assert abs(math.fsum(arrays["gaussian_weights"])-1) < 2e-11
    probabilities = report["source_record_probabilities"]
    for record, weight in zip(arrays["field_records"], arrays["gaussian_weights"], strict=True):
        marginal = math.fsum(value for key, value in probabilities.items() if key.startswith(str(record)))
        assert abs(marginal-weight) < 2e-12
