"""Complete native family selection, records, source factors and references."""
import numpy as np
import pytest
import sympy as s

from . import native_family_selection as f


@pytest.fixture(scope='module')
def execution(): return f.execute()


def test_incidence_spectrum_computes_band_rank_and_retains_every_sector():
    d = f.geometry()
    x = s.Symbol('x')
    assert d['polynomial'] == x**20*(x-12)**2*(x+3)**8*(x*x-4*x-16)**3
    assert d['rank'] == 3 and s.trace(d['Q']) == 33
    assert s.trace(d['fixed']) == 2
    assert f._zero(sum(d['projectors'].values(), s.zeros(36))-s.eye(36))
    assert d['seed'] == 8 and s.simplify(d['scale']-5-s.sqrt(5)) == 0


def test_full_twelve_coefficient_interface_and_extra24_levels():
    d = f.geometry(); W = d['W']
    assert f._zero(W.T*W-s.eye(12))
    assert 36-int(s.trace(f._clean(W.T*W))) == 24
    assert f._zero(d['C']*W-W*d['C12'])
    assert sorted(d['coefficient_roots'].values()) == [1, 3, 3, 5]


def test_complete433_logical_records_and_terminal_padding_have_tp_extension():
    K, labels = f.cooling_maps()
    assert K.shape == (605, 64, 64) and len(set(map(tuple, labels))) == 605
    assert np.count_nonzero(labels[:, 0] == 0) == 1
    assert np.count_nonzero(labels[:, 0] == 1) == 432
    assert np.count_nonzero(labels[:, 0] == 2) == 28
    assert np.count_nonzero(labels[:, 0] == 3) == 144
    gram = np.einsum('rij,rik->jk', K.conj(), K, optimize=True)
    np.testing.assert_allclose(gram, np.eye(64), atol=2e-14)


def test_record_labels_are_covariant_under_every_native_rotation():
    d = f.geometry(); P = np.asarray(d['P'], float); A = np.asarray(d['frame'], float)
    for g, action in zip(d['rotations'], d['pair_actions'], strict=True):
        R = np.asarray(f.objects.permutation_matrix(action), float)
        G = np.asarray(f.objects.permutation_matrix(g), float)
        np.testing.assert_allclose(R@P, P@R, atol=1e-14)
        np.testing.assert_allclose(R@A, A@G, atol=1e-14)
        # This is the full factor identity for every i,j reset map; Q rows
        # transform by the same ordinary permutation as the 36 color labels.
        np.testing.assert_allclose(R@(np.eye(36)-P), (np.eye(36)-P)@R, atol=1e-14)


def test_fixed_state_space_and_full_complement_reset_without_postselection():
    d = f.geometry(); P = np.asarray(d['P'], complex)
    rng = np.random.default_rng(18)
    value = rng.normal(size=(36, 7))+1j*rng.normal(size=(36, 7))
    rho = value@value.conj().T; rho /= np.trace(rho)
    output = f.apply_cooling(rho)
    assert np.trace(output).real == pytest.approx(1)
    np.testing.assert_allclose(P@output@P, output, atol=1e-14)
    np.testing.assert_allclose(f.apply_cooling(output), output, atol=1e-14)
    sigma = P@rho@P; sigma /= np.trace(sigma)
    np.testing.assert_allclose(f.apply_cooling(sigma), sigma, atol=1e-14)


def test_heisenberg_fixed_operator_system_is_not_full_ordinary_algebra():
    d = f.geometry(); A = np.asarray(d['frame'], complex); P = np.asarray(d['P'], complex)
    vector = 2*A[:, 0]; B = np.outer(vector, vector.conj())
    fixed = f.adjoint_cooling(B)
    np.testing.assert_allclose(f.adjoint_cooling(fixed), fixed, atol=1e-14)
    assert np.linalg.norm(fixed@fixed-f.adjoint_cooling(B@B)) > .1
    np.testing.assert_allclose(fixed, B+np.trace(B)*(np.eye(36)-P)/3, atol=1e-14)


def test_source_compilation_retains_full_unknown_input_domain(execution):
    report, arrays = execution
    m = report['actual_native_composition']
    assert m['unknown_native_logical_inputs'] == 36 and m['full_encoded_register_inputs'] == 64
    assert m['all_inputs_and_all_raw_records_in_composition']
    assert m['full_register_native_plus_cooling_records'] == 77440
    assert m['native_plus_cooling_error_bound_input_domain'].startswith('all36 protected')
    assert not m['full64_native_Gram_check_proves_identity_extension_closeness']
    assert m['source_cooling']['all64_unknown_input_and_reference_columns_executed']
    assert m['source_cooling']['source_cooling_dilation_error_diagnostic'] < f.SOURCE_ADMISSION
    assert m['complete_native_plus_cooling_Gram_defect_diagnostic'] < 1e-9
    assert arrays['family_ray_preparations'].shape == (12, 64, 64)
    for U in arrays['family_ray_preparations']:
        np.testing.assert_allclose(U.conj().T@U, np.eye(64), atol=2e-11)


def test_compiled_coordinate_and_flip_matrices_are_source_outputs(execution):
    _, arrays = execution
    k = arrays['family_coordinate_local']
    np.testing.assert_allclose(k[0].conj().T@k[0]+k[1].conj().T@k[1], np.eye(2), atol=1e-13)
    for bit in (0, 1):
        ideal = np.zeros((2, 2)); ideal[bit, bit] = 1
        np.testing.assert_allclose(k[bit], ideal, atol=1e-13)
    np.testing.assert_allclose(arrays['family_flip_local'], [[0, 1], [1, 0]], atol=1e-13)
    assert abs(np.linalg.norm(arrays['family_selector_outputs'])-1) < 1e-12


def test_raw_reference_contraction_matches_literal_kronecker(execution):
    _, arrays = execution
    # An independent four-syndrome subinstrument check of the same contraction
    # used for the complete source/reference Gram; no sample is a physics test.
    L = arrays['family_C'][0].reshape(64, 8, 8)[:2]
    R = arrays['family_E'][0].reshape(64, 8, 8)[:2]
    G = f.cooling_gram(arrays)
    literal = sum((np.kron(a, b).conj().T@G@np.kron(a, b) for a in L for b in R), np.zeros((64, 64), complex))
    a = np.einsum('hia,hjc->ijac', L.conj(), L)
    b = np.einsum('hkb,hld->klbd', R.conj(), R)
    factored = np.einsum('ijac,ikjl,klbd->abcd', a, G.reshape(8, 8, 8, 8), b).reshape(64, 64)
    np.testing.assert_allclose(factored, literal, atol=1e-14)


@pytest.mark.parametrize('kind', ['missing', 'extra', 'dtype', 'shape', 'nonfinite', 'labels'])
def test_source_menu_and_all_arrays_fail_closed(execution, kind):
    _, original = execution; arrays = {k: v.copy() for k, v in original.items()}
    if kind == 'missing': del arrays['family_selector_outputs']
    if kind == 'extra': arrays['untracked'] = np.zeros(1)
    if kind == 'dtype': arrays['family_cooling_labels'] = arrays['family_cooling_labels'].astype(float)
    if kind == 'shape': arrays['family_C'] = arrays['family_C'][:, :, :6]
    if kind == 'nonfinite': arrays['family_ray_preparations'][0, 0, 0] = np.nan
    if kind == 'labels': arrays['family_cooling_labels'][-1, 1] = 0
    with pytest.raises(ValueError): f.validate_arrays(arrays)


def test_source_selector_failure_and_ray_complement_cannot_be_dropped(execution):
    _, original = execution
    arrays = {k: v.copy() for k, v in original.items()}
    arrays['family_selector_outputs'][12, 0] += 1e-4
    with pytest.raises(ValueError): f.cooling_metrics(arrays)
    arrays = {k: v.copy() for k, v in original.items()}
    arrays['family_ray_preparations'][0, :, 9] += 1e-4
    with pytest.raises(ValueError): f.cooling_metrics(arrays)
