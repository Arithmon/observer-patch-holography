"""Same full interacting action, noncommuting variation, and regional limits."""
import pytest
from .entropy import run


@pytest.fixture(scope='module')
def entropy_packet():
    return run()


def test_quantum_first_law_is_not_a_commuting_or_zero_response_test(entropy_packet):
    row = entropy_packet
    assert row['reference_dimension'] == 1431
    assert row['reference_minimum_eigenvalue'] > 0
    assert row['full_Gibbs_modular_max_entry_residual'] < 2e-11
    assert row['variation_reference_commutator_Hilbert_Schmidt'] > 1e-6
    assert row['energy_derivative'] > .9
    assert row['entropy_derivative'] > 1
    assert row['first_law_derivative_residual'] < 2e-11
    assert row['quantum_relative_entropy'] > .01
    assert row['exact_identity_numerical_residual'] < 2e-11
    assert row['quantum_relative_entropy'] < row['analytic_relative_entropy_upper_evaluated_numerically']


def test_regional_modular_operator_is_not_affine_bare_energy(entropy_packet):
    row = entropy_packet
    assert row['regional_dimension'] == 11
    assert [len(x) for x in row['regional_basis']] == [0]+[1]*4+[2]*6
    assert row['regional_trace_residual'] < 2e-11
    assert row['regional_minimum_eigenvalue'] > 0
    assert row['regional_modular_bare_energy_commutator_Hilbert_Schmidt'] > .0007
    simple,extended = row['bare_regional_modular_comparisons']
    assert simple['Hilbert_Schmidt_residual'] > .48
    assert extended['Hilbert_Schmidt_residual'] > .0058
    assert not row['density_preparation_source_executed']
    assert not row['physical_clock_or_temperature_calibrated']
    assert not row['outward_certified_floating_error']
