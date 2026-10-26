"""Independent finite thermal, charge, quantum, source and record controls."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import reciprocal_vacuum as r


@pytest.fixture(scope='module')
def scalar():
    return r.finite_model()


@pytest.fixture(scope='module')
def evidence():
    return r.build()


def test_same_actual_previous_hamiltonians_and_noncommuting_state_family(scalar):
    for exact, geometry in zip(scalar.exact_hamiltonians, r.GEOMETRIES, strict=True):
        np.testing.assert_allclose(np.asarray(exact, float), r.vacuum_repair.hamiltonian(2, 2, geometry), rtol=0, atol=2e-13)
    assert scalar.exact_hamiltonians[0]*scalar.exact_hamiltonians[1] != scalar.exact_hamiltonians[1]*scalar.exact_hamiltonians[0]
    assert all(row['minimum_resolved_sector_gap'] > 1e-6 for row in scalar.spectral_diagnostics)


def test_exact_pair_theorem_and_full_quantum_gibbs_not_only_populations(scalar):
    assert r.exact_balance_identity()['collision_energy_commutator_exact']
    expected = expm(-float(scalar.beta)*scalar.hamiltonian)
    expected /= np.trace(expected)
    np.testing.assert_allclose(scalar.gibbs, expected, rtol=0, atol=2e-13)
    assert np.linalg.norm(scalar.gibbs-np.diag(np.diag(scalar.gibbs)), 2) > .01
    np.testing.assert_allclose(r.apply_channel(scalar, scalar.gibbs), scalar.gibbs, rtol=0, atol=2e-13)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_native_proposals_rejections_and_bath_energy_outputs(scalar, law):
    rows, weights, heat, native = r.records(scalar, law)
    assert rows.shape == (512, 8) and native.shape == (2, 16, 4, 3)
    assert set(rows[:, 7]) == {0, 1}
    for i in scalar.active:
        assert abs(weights[rows[:, 0] == i].sum()-1) < 2e-14
    assert abs(sum(scalar.probabilities[row[0]]*weight*q for row, weight, q in zip(rows, weights, heat))) < 2e-13
    # No-op proposals and failed bath activations remain part of the law.
    assert np.any((rows[:, 2] == rows[:, 0]//8) & (rows[:, 7] == 0))
    assert np.any((rows[:, 2] != rows[:, 0]//8) & (rows[:, 7] == 0))
    for row, q in zip(rows, heat, strict=True):
        assert q == scalar.energies[row[6]]-scalar.energies[row[0]]


def test_degenerate_gap_uses_one_bath_orientation_and_reverse_joint_flows():
    # Exact degeneracies must not make both transition directions consume
    # bath1. The same fixed bath mode/unitary serves either initial state.
    h = s.diag(*range(8))
    model = r.finite_model((h, h))
    rows, weights, _, _ = r.records(model)
    flows = {(int(row[0]), int(row[6]), int(row[4]), int(row[5])):
             model.probabilities[row[0]]*weight
             for row, weight in zip(rows, weights, strict=True) if row[7]}
    assert flows
    for (i, j, before, after), probability in flows.items():
        assert (j, i, after, before) in flows
        assert probability == pytest.approx(flows[j, i, after, before], abs=2e-15)


def test_nonuniform_geometry_and_derived_free_energy_response(scalar):
    result = r.thermodynamics(scalar)
    p = result['geometry_probabilities']
    assert .57 < p[0] < .58 and p[0] > p[1]
    assert result['uniform_geometry_conditional_Gibbs_state_stationarity_defect_L1'] > .05
    assert result['isothermal_free_energy_cost_of_forcing_that_uniform_state'] > .04
    assert abs(result['geometry_free_energy_first_variation']) < 1e-12
    assert result['geometry_free_energy_second_variation'] > 0
    assert abs(result['entropy_chain_rule_residual']) < 1e-12
    assert result['population_gap_diagnostic_on_admitted_sector'] > 0
    step = s.Rational(1, 100000)
    plus = r.finite_model((scalar.exact_hamiltonians[0], scalar.exact_hamiltonians[1]+step*s.eye(8)))
    minus = r.finite_model((scalar.exact_hamiltonians[0], scalar.exact_hamiltonians[1]-step*s.eye(8)))
    measured = (sum(plus.probabilities[8:])-sum(minus.probabilities[8:]))/(2*float(step))
    assert measured == pytest.approx(result['geometry_probability_response_to_added_branch1_energy'], abs=2e-10)
    # Reciprocity of the same two conjugate thermal responses, rather than
    # separately fitting geometry and field update rates.
    beta_step = s.Rational(1, 100000)
    plus_beta = r.finite_model(beta=scalar.beta+beta_step)
    minus_beta = r.finite_model(beta=scalar.beta-beta_step)
    dp_dbeta = (sum(plus_beta.probabilities[8:])-sum(minus_beta.probabilities[8:]))/(2*float(beta_step))
    de_dpotential = (plus.probabilities@scalar.energies-minus.probabilities@scalar.energies)/(2*float(step))
    assert dp_dbeta == pytest.approx(result['geometry_probability_beta_response'], abs=2e-10)
    assert de_dpotential == pytest.approx(float(scalar.beta)*dp_dbeta, abs=5e-10)


def test_joint_probability_entropy_variation_is_thermodynamic_not_bare_area(scalar):
    p = sum(scalar.probabilities[8:]); beta = float(scalar.beta)
    conditionals = [scalar.gibbs[:8, :8]/(1-p), scalar.gibbs[8:, 8:]/p]
    def free(x):
        rho = np.zeros((16, 16), complex)
        rho[:8, :8], rho[8:, 8:] = (1-x)*conditionals[0], x*conditionals[1]
        return float(np.trace(scalar.hamiltonian@rho).real)-r._entropy(rho)/beta
    step = 1e-4
    first = (free(p+step)-free(p-step))/(2*step)
    second = (free(p+step)-2*free(p)+free(p-step))/step**2
    assert abs(first) < 2e-7
    assert second == pytest.approx((1/p+1/(1-p))/beta, abs=5e-6)


def test_full_quantum_channel_decreases_free_energy_from_coherent_input(scalar):
    v = np.zeros(16, complex); v[0], v[9] = 1/math.sqrt(2), 1j/math.sqrt(2)
    rho = np.outer(v, v.conj())
    after = r.apply_channel(scalar, rho)
    heat = float(np.trace(scalar.hamiltonian@(after-rho)).real)
    entropy_change = r._entropy(after)-r._entropy(rho)
    assert entropy_change-float(scalar.beta)*heat > 0
    assert abs(np.trace(after)-1) < 2e-14


def test_source_full_collision_and_purified_reference_outputs(evidence):
    report, arrays = evidence
    assert arrays['collision_unitary'].shape == (32, 32)
    assert arrays['purified_collision_isometry'].shape == (64, 16)
    assert arrays['energy_readout_isometry'].shape == (256, 16)
    assert report['source_collision']['source_operator_error'] < r.SOURCE_ADMISSION
    assert report['source_collision']['literal_endpoint_total_energy_defect'] < 100*r.SOURCE_ADMISSION
    np.testing.assert_allclose(arrays['purified_collision_isometry'].conj().T@arrays['purified_collision_isometry'], np.eye(16), rtol=0, atol=2e-10)
    assert not report['source_collision']['physical_control_work_or_error_certified']
    assert report['source_energy_readout']['source_all_column_operator_error'] < r.SOURCE_ADMISSION
    assert report['finite_bath_contract']['fresh_thermal_inputs_or_rethermalization_required_for_repeated_Markov_law']


def test_fixed_charge_matter_case_preserves_particles_and_is_not_grandcanonical(evidence):
    report, arrays = evidence
    model = r.fermion_model()
    assert len(model.active) == 6 and {model.charges[i] for i in model.active} == {1}
    np.testing.assert_allclose(r.apply_channel(model, model.gibbs), model.gibbs, rtol=0, atol=2e-13)
    rows = arrays['fermion_C_record_labels']
    assert all(model.charges[a] == model.charges[b] == 1 for a, b in rows[:, [0, 6]])
    matter = report['native_CAR_matter_control']
    assert matter['source_encoded_charge_commutator'] < r.SOURCE_ADMISSION
    assert not matter['bath_particle_exchange_or_chemical_potential_assumed']
    assert not matter['unique_grandcanonical_state_from_arbitrary_charge_weights_claimed']
    empty = np.zeros((16, 16), complex); empty[0, 0] = 1
    with pytest.raises(ValueError, match='sector'):
        r.apply_channel(model, empty)


def test_generic_small_case_and_degeneracy_domain_are_explicit():
    z = s.diag(-1, 1)
    m = r.finite_model((z, 2*z), beta=s.Integer(1))
    np.testing.assert_allclose(r.apply_channel(m, m.gibbs), m.gibbs, atol=1e-14)
    with pytest.raises(ValueError, match='fixed8-level'):
        r.records(m)
    with pytest.raises(ValueError, match='16-level'):
        r.source_collision(m)
    # A non-diagonal unresolved degenerate eigenspace requires finer exact
    # symmetry labels; an arbitrary LAPACK rotation is never source authority.
    v = s.Matrix([[1, 1, 0, 0], [1, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])/s.sqrt(2)
    v[2, 2] = v[3, 3] = 1
    h = v*s.diag(0, 1, 1, 2)*v.T
    with pytest.raises(ValueError, match='degenerate'):
        r.finite_model((h, h))


def test_invalid_unnormalized_and_nonsymmetric_inputs_refuse(scalar):
    with pytest.raises(ValueError, match='normalized'):
        r.apply_channel(scalar, 2*scalar.gibbs)
    with pytest.raises(ValueError, match='masked'):
        r.apply_channel(scalar, np.ma.array(scalar.gibbs, mask=False))
    with pytest.raises(ValueError):
        r.finite_model((s.Matrix([[0, 1], [0, 0]]), s.eye(2)))
    with pytest.raises(ValueError):
        r.finite_model(beta=s.Float(1.))
    with pytest.raises(ValueError):
        r.finite_model(beta=True)
    with pytest.raises(ValueError):
        r.proposal_records('wrong')


def test_receipt_cannot_claim_physical_autonomy_or_change_fixed_source_gate(evidence):
    report, _ = evidence
    for field in ('scope', 'gate'):
        changed = deepcopy(report)
        if field == 'scope':
            changed['scope']['isolated_autonomous_cooling_or_reservoir_preparation_derived'] = True
        else:
            changed['source_collision']['source_all_input_operator_admission'] += 1e-10
        with pytest.raises(ValueError):
            r._compare(report, changed)
