"""Actual-port, protection, symmetry and complete-output controls."""
from copy import deepcopy
from fractions import Fraction as F
import itertools

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from . import native_virtual_exchange as v
from . import record_code as code


def test_complete_native_carrier_and_all_twelve_primitive_swaps():
    masks, work, permutations, _, _ = v.geometry()
    assert len(masks) == 10626 == __import__('math').comb(24, 4)
    assert work.shape == (36, 36)
    for port, permutation in enumerate(permutations):
        assert np.array_equal(permutation[permutation], np.arange(len(masks)))
        for before, after in enumerate(permutation):
            first, second = masks[before], masks[after]
            assert first.bit_count() == second.bit_count() == 4
            assert (first ^ second) & ~((1 << port) | (1 << (port+12))) == 0
            assert ((first >> port) & 1) == ((second >> (port+12)) & 1)
    for first, second in itertools.combinations(permutations, 2):
        assert np.array_equal(first[second], second[first])


def test_exact_matching_projectors_and_first_order_are_computed_from_native_ports():
    native, u, pa, pb = v.source_matrices()
    assert np.array_equal((u.T@u).toarray(), 36*np.eye(36, dtype=int))
    assert (pa@pa-6*pa).nnz == 0
    assert (pb@pb-6*pb).nnz == 0
    assert (pa@pb-pb@pa).nnz == 0
    assert np.array_equal((u.T@native@u).toarray(), 120*np.eye(36, dtype=int))
    b, a = v.exact_interaction()
    assert (pa@a).nnz == (pb@a).nnz == 0
    assert s.Matrix((a.T@a).toarray())/324 == b
    assert b != s.eye(36)*b[0, 0]


def test_load_only_protection_leaves_both_same_shell_and_interload_syndromes():
    b, a = v.exact_interaction()
    del b
    masks = v.geometry()[0]
    changed = np.array([(mask & 4095).bit_count() != 2 for mask in masks])
    assert a[changed].nnz > 0
    assert a[~changed].nnz > 0
    # The original per-carrier shell is66-dimensional, not the rank6 work code.
    assert sum(not flag for flag in changed) == 66*66


def test_complete_native_penalty_spectrum_has_the_declared_controlled_error():
    result, arrays = v.finite_penalty_spectrum()
    cert = result['protection_certificate']
    assert cert['admitted']
    actual, approximate = arrays['complete_low_energies'], arrays['second_order_low_energies']
    assert actual.shape == approximate.shape == (36,)
    assert np.all(actual >= -v.NUMERICAL_GATE)
    assert max(actual) <= float(F(cert['low_band_energy_upper']))+v.NUMERICAL_GATE
    assert max(abs(actual-approximate)) <= float(F(cert['ordered_second_order_energy_error_upper']))
    assert not result['spectra_are_physical_particle_masses']
    assert not cert['continuous_penalty_is_the_recorded_CE_cycle']


def test_protection_gap_domain_refuses_unseparated_band():
    assert not v.protection_bound(F(1), F(1))['admitted']
    assert v.protection_bound(F(1), F(1))['ordered_second_order_energy_error_upper'] is None
    cert = v.protection_bound(F(1, 100), F(1))
    assert cert['admitted']
    assert F(cert['virtual_Gram_operator_norm_upper']) == F(80, 9)


@pytest.mark.parametrize('bad', (True, 0, -1, 0.1, float('nan'), '1/3'))
def test_original_exact_parameter_domain(bad):
    with pytest.raises(ValueError):
        v.protection_bound(bad, F(1))
    with pytest.raises(ValueError):
        v.commutator_word(bad)


def test_native_port_word_retains_all_outputs_and_has_the_correct_generator_sign():
    coarse, first = v.commutator_word(F(1, 128))
    fine, second = v.commutator_word(F(1, 256))
    for report, arrays in ((coarse, first), (fine, second)):
        output, target = arrays['port_word_output'], arrays['port_word_target']
        assert output.shape == target.shape == (10626, 36)
        assert np.linalg.norm(output.conj().T@output-np.eye(36), 2) < v.NUMERICAL_GATE
        assert np.linalg.norm(output-target, 2) < float(F(report['complete_operator_error_upper']))
        assert not report['code_projection_or_postselection_applied']
        assert not report['continuous_penalty_required_for_this_word']
    _, u, _, _ = v.source_matrices()
    frame = u/6
    b, _ = v.exact_interaction()
    observed = (frame.T@second['port_word_output']-np.eye(36))/(2j*float(F(1, 256)**2))
    assert np.linalg.norm(observed-np.asarray(b, float), 2) < 0.1
    error_first = np.linalg.norm(first['port_word_output']-first['port_word_target'], 2)
    error_second = np.linalg.norm(second['port_word_output']-second['port_word_target'], 2)
    assert 7 < error_first/error_second < 9


def test_finite_word_angle_domain():
    with pytest.raises(ValueError):
        v.commutator_word(F(1, 4))


def test_source_compiler_executes_all_local_port_macro_columns():
    report, arrays = v.port_macro_controls()
    assert len(arrays) == 4
    assert report['word_port_macro_source_pulses'] == 25056
    assert not report['whole_word_joint_code_phases_compiled_here']
    assert not report['local_macro_admission_is_a_physical_service_certificate']
    for index, row in enumerate(report['local_source_macros']):
        angle = float(F(row['angle_exact']))
        actual = arrays[f'port_macro_{index}']
        symmetric = np.array([0, 1, 1, 0])/np.sqrt(2)
        antisymmetric = np.array([0, 1, -1, 0])/np.sqrt(2)
        assert np.linalg.norm(actual@symmetric-symmetric) < v.NUMERICAL_GATE
        assert np.linalg.norm(actual@antisymmetric-np.exp(-2j*angle)*antisymmetric) < v.NUMERICAL_GATE
        assert np.linalg.norm(actual.conj().T@actual-np.eye(4), 2) < v.NUMERICAL_GATE


def test_existing_logical_charge_and_parity_are_not_assumed_conserved():
    result, _ = v.symmetry_contract()
    assert result['raw_occupation_commutator_nonzeros'] == 162
    assert result['raw_fermion_parity_commutator_nonzeros'] == 104
    assert not result['logical_occupation_identified_with_electric_charge']
    assert not result['occupation_average_is_an_SM_gauge_attachment']


def test_gauge_projection_uses_actual_generators_and_is_an_orthogonal_projection():
    b, _ = v.exact_interaction()
    projected, _, _ = v.gauge_projection(b)
    repeated, _, _ = v.gauge_projection(projected)
    assert repeated == projected
    assert s.trace(b) == s.trace(projected)
    assert s.trace(projected.T*(b-projected)) == 0
    assert v.gauge_projection(s.eye(36))[0] == s.eye(36)
    for d in v.source_generators():
        current = s.kronecker_product(d, s.eye(6))+s.kronecker_product(s.eye(6), d)
        assert any((b*current-current*b).applyfunc(s.expand))
        assert not any((projected*current-current*projected).applyfunc(s.expand))
    assert projected.eigenvals() == {
        s.Rational(32, 9): 15, s.Rational(44, 9): 9, s.Rational(47, 9): 6,
        s.Rational(154, 27): 1, s.Rational(692, 135): 5}


def test_actual_port_stabilizer_and_symmetry_closed_span():
    result, arrays = v.persistent_port_symmetry()
    assert result['native_incidence_automorphisms'] == 120
    assert result['actual_code_stabilizer_order'] == 8
    assert len(result['distinct_logical_stabilizer_actions']) == 4
    assert result['exact_symmetry_closed_span_dimension'] == 20
    assert result['distinct_matching_orbit_vectors'] == 66
    projector = arrays['symmetry_orbit_projector']
    assert np.linalg.norm(projector@projector-projector, 2) < 1e-13
    assert np.trace(projector) == pytest.approx(20)
    assert not result['symmetry_closed_span_is_correctable_for_complete_recorded_menu']
    assert not result['full_response_M6_action_on_this_object_derived']


def test_orbit_correctability_witness_uses_the_actual_retained_C_record():
    result, _ = v.persistent_port_symmetry()
    witness = result['correctability_witness']
    seam = tuple(witness['seam'])
    assert seam in code.SEAMS
    vectors = []
    effects = []
    for label in ('first_matching_pairs', 'second_matching_pairs'):
        support = {tuple(pair) for pair in witness[label]}
        vectors.append(support)
    shared = vectors[0] & vectors[1]
    for pair in shared:
        x = [int(port in pair) for port in range(12)]
        left = (x[seam[0]]+x[seam[1]])//2  # actual coin0 floor repair
        effects.append(F(int(x[seam[0]]-left == 1), 2))
    overlap = F(len(shared), 6)
    branch = sum(effects, F())/6
    defect = branch-F(1, 12)*overlap
    assert str(overlap) == witness['normalized_vector_overlap']
    assert str(branch) == witness['normalized_branch_off_diagonal']
    assert defect == F(witness['Knill_Laflamme_defect_exact']) == F(5, 72)
    # Both individual matching vectors retain the original1/12 probability.
    for support in vectors:
        prob = sum(F(int(seam[0] in pair and seam[1] not in pair), 2) for pair in support)/6
        assert prob == F(1, 12)


def test_receipt_claims_and_fixed_admissions_cannot_be_forged():
    report, _ = v.build()
    assert report['empirical_inputs_to_prediction'] == []
    assert not report['scope']['observed_particle_spectrum_or_physical_mass_scale_derived']
    forged = deepcopy(report)
    forged['scope']['observed_particle_spectrum_or_physical_mass_scale_derived'] = True
    with pytest.raises(ValueError):
        v._compare(report, forged)
    forged = deepcopy(report)
    forged['numerical_admission'] *= 2
    with pytest.raises(ValueError):
        v._compare(report, forged)
    forged = deepcopy(report)
    forged['source_executed_local_port_macros']['numerical_operator_admission'] *= 1.01
    with pytest.raises(ValueError):
        v._compare(report, forged)
    forged = deepcopy(report)
    forged['complete_penalty_spectrum']['complete_low_energies'][0] = float('nan')
    with pytest.raises(ValueError):
        v._compare(report, forged)
