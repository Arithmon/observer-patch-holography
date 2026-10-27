"""Native source words, complete records, reference domains and composition."""
import numpy as np
import pytest
from scipy.linalg import expm

from . import native_repair_mass_feedback as f


@pytest.fixture(scope='module')
def execution(): return f.execute()


def test_actual_port_rays_attach_the_existing_kernel_without_mass_eigenvectors():
    d = f.geometry()
    assert len(set(d['ports'])) == 5
    np.testing.assert_allclose(d['frame'].conj().T@d['frame'], np.eye(3), atol=1e-14)
    for ray in d['rays']:
        assert abs(ray@ray-1) < 1e-14
        np.testing.assert_allclose(d['frame']@d['frame'].conj().T@ray, ray, atol=1e-14)


def test_source_geometry_does_not_read_averaged_kernel_or_its_spectrum(monkeypatch):
    original = f.geometry()['rays'].copy()
    law = dict(f.feedback.architecture())
    for key in ('Y', 'e2', 'determinant', 'stationary', 'transition'): law.pop(key)
    monkeypatch.setattr(f.feedback, 'architecture', lambda: law)
    f.geometry.cache_clear()
    np.testing.assert_array_equal(f.geometry()['rays'], original)


def test_literal_full128_feedback_words_and_reference_columns(execution):
    report, arrays = execution
    X = np.asarray([[0, 1], [1, 0]])
    for ray, U in zip(arrays['feedback_port_rays'], arrays['feedback_source_unitaries'], strict=True):
        target = expm(-1j*f.STEP*np.kron(X, np.outer(ray, ray)))
        np.testing.assert_allclose(U, target, atol=2e-11)
        np.testing.assert_allclose(U.conj().T@U, np.eye(128), atol=2e-11)
    m = report['execution']
    assert m['six_unknown_input_feedback_error_diagnostic'] < f.SOURCE_ADMISSION
    assert m['complete128_native_feedback_Gram_defect_diagnostic'] < f.SOURCE_ADMISSION
    assert m['factor_proved_prepared_dilation_error_diagnostic'] < m['prepared_six_input_dilation_allowance']
    assert not report['composition']['target_projectors_are_applied_instead_of_source_words']


def test_complete_memory_records_and_syndrome_outputs(execution):
    report, arrays = execution
    labels = arrays['feedback_branch_labels']
    assert len(set(map(tuple, labels))) == 640
    for cc, cr, ec, er, old, new in labels:
        assert new == ((old+1) % 5 if (cc, cr) == (0, 1) else 0)
    d = report['complete_domains']
    assert d['old_decoder_syndrome_quantum_dimension'] == 4096
    assert not d['old_decoder_syndromes_erased_or_conditioned_away']
    assert arrays['family_C'].shape == (8, 512, 8) and arrays['family_E'].shape == (16, 512, 8)
    assert not report['scope']['effective_tiny_time_readout_intervals_transferred_to_this_source_history']
    m = report['execution']
    assert m['retained_low_cycle_CP_norm_upper'] > 1
    assert m['sixty_four_cycle_geometric_amplification_sum'] > 64
    assert m['sixty_four_low_cycle_telescope_channel_allowance'] > 64*m['low_register_cycle_channel_allowance']
    assert m['prepared_sixty_four_cycle_final_native_output_channel_allowance'] > m['sixty_four_low_cycle_telescope_channel_allowance']
    assert report['composition']['additional_final_reencoder_executed_once']


def test_complete_gram_contraction_against_literal_tensor_maps(execution):
    _, source = execution
    a = {key: value.copy() for key, value in source.items()}
    a['family_C'][:] = 0; a['family_E'][:] = 0
    # Two raw records and two quantum syndromes on each side, including a
    # nonzero padded-input operation. The independent sum uses literal128
    # Kronecker maps rather than the channel-index contraction.
    c_terms, e_terms = [], []
    for name, terms, records in (('C', c_terms, (1, 2)), ('E', e_terms, (0, 3))):
        for n, raw in enumerate(records):
            for h in (0, 7):
                V = np.roll(np.eye(8), n+(h != 0), axis=0)/2
                a['family_'+name][raw, 8*h:8*h+8] = V
                terms.append((raw, V))
    a['feedback_encoder_C'][0, 0] += .02
    encC, encE = a['feedback_encoder_C'], a['feedback_encoder_E']
    W = np.kron(np.eye(2), np.kron(encC.conj().T@encC, encE.conj().T@encE))
    computed = f.complete_grams(a)
    for old in (0, 4):
        literal = np.zeros((128, 128), complex)
        for raw, C in c_terms:
            cc, cr = a['family_C_labels'][raw]
            new = (old+1) % 5 if (cc, cr) == (0, 1) else 0
            for _, E in e_terms:
                K = a['feedback_source_unitaries'][new]@np.kron(np.eye(2), np.kron(C, E))
                literal += K.conj().T@W@K
        np.testing.assert_allclose(computed[old], literal, atol=2e-14)


@pytest.mark.parametrize('mutation', ['missing', 'dtype', 'shape', 'nan', 'memory_label'])
def test_complete_factor_schema_fails_closed(execution, mutation):
    _, source = execution; a = {key: value.copy() for key, value in source.items()}
    if mutation == 'missing': del a['feedback_encoder_E']
    if mutation == 'dtype': a['feedback_source_unitaries'] = a['feedback_source_unitaries'].astype(np.complex64)
    if mutation == 'shape': a['feedback_source_preparation'] = a['feedback_source_preparation'][:63]
    if mutation == 'nan': a['family_E'][0, 511, 7] = np.nan
    if mutation == 'memory_label': a['feedback_branch_labels'][-1, -1] = 4
    with pytest.raises(ValueError): f.validate_arrays(a)


def test_padded_word_column_and_old_syndrome_tampering_are_detected(execution):
    _, source = execution
    for key, index in (('feedback_source_unitaries', (4, 127, 127)), ('family_E', (0, 511, 7))):
        a = {name: value.copy() for name, value in source.items()}
        a[key][index] += .01
        with pytest.raises(ValueError): f.metrics(a)


def test_literal_sources_do_not_import_empirical_mass_material():
    source = (f.ROOT/'model/common_source/native_repair_mass_feedback.py').read_text(encoding='utf-8')
    assert 'sm_calibration' not in source and 'parameter_packet' not in source
    assert 'chiral_candidate' not in source and 'np.linalg.eigh(' not in source
