"""Full-register clock, intersector pulse and spectral detector certificates."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import expm

from . import prepared_clock_source as m


@pytest.fixture(scope='module')
def source():
    return m.compiled_source()


def test_spectrum_comes_from_exact_native_seams_with_full_constant_and_padding():
    row = m.native_spectral_data()
    independent = sp.zeros(12)
    for a, b in row['seams']:
        vector = sp.zeros(12, 1)
        vector[a], vector[b] = 1, -1
        independent += vector*vector.T
    assert independent == row['laplacian']
    assert independent.eigenvals() == {0: 1, 5-sp.sqrt(5): 3, 6: 5, 5+sp.sqrt(5): 3}
    assert tuple(row['basis_exact'][:12, 0]) == (1/sp.sqrt(12),)*12
    for band, addresses in enumerate(row['band_addresses']):
        for address in addresses:
            vector = row['basis_exact'][:12, address]
            assert all(sp.simplify(x) == 0
                       for x in independent*vector-row['band_eigenvalues'][band]*vector)
    np.testing.assert_array_equal(row['basis'][12:, 12:], np.eye(4))
    np.testing.assert_array_equal(row['basis'][:12, 12:], np.zeros((12, 4)))


def test_reference_is_unnormalized_sampler_polynomial_before_global_normalization(source):
    row, spectrum = source['reference'], source['spectrum']
    lazy = sp.eye(12)-spectrum['laplacian']/120
    independent = sum((sp.Rational(w)*lazy**(2*n) for n, w in enumerate(row['weights'])),
                      sp.zeros(12))
    assert row['polynomial'] == independent
    assert row['density_exact'] == independent/sp.trace(independent)
    assert row['m0'] == sum(row['weights'])
    assert all(x > 0 for x in np.linalg.eigvalsh(np.array(row['density_exact'], float)))
    # Normalizing each sampled component before mixing changes the reference.
    wrong = sum((float(w)*np.array(lazy**(2*n), float)/float(sp.trace(lazy**(2*n)))
                 for n, w in enumerate(row['weights'])), np.zeros((12, 12)))
    wrong /= float(sum(row['weights']))
    assert np.linalg.norm(wrong-row['density'][:12, :12], 'fro') > 1e-5


def test_raw_complete_clocks_match_independent_matrix_functions_and_sign(source):
    laplacian = np.array(source['spectrum']['laplacian'], float)
    values, vectors = np.linalg.eigh(laplacian)
    values = np.maximum(values, 0)
    weights = source['reference']['weights']
    masses = sum(float(w)*(1-values/120)**(2*n) for n, w in enumerate(weights))
    hamiltonian = (vectors*np.log(float(sum(weights))/masses))@vectors.T
    square_root = (vectors*np.sqrt(values))@vectors.T
    for name, generator in (('clock1', hamiltonian), ('sqrt_clock1', square_root)):
        expected = np.eye(16, dtype=complex)
        expected[:12, :12] = expm(-1j*generator)
        np.testing.assert_allclose(source['raw'][name], expected,
                                   atol=float(source['guards'][name])+2e-8, rtol=0)
        assert np.linalg.norm(source['raw'][name]-expected.conj().T, 'fro') > .1
        np.testing.assert_allclose(source['raw'][name][12:, 12:], np.eye(4), atol=1e-9, rtol=0)
    assert np.linalg.norm(source['raw']['clock1']-source['raw']['sqrt_clock1'], 'fro') > 1
    state = source['reference']['density']
    actual = source['raw']['clock1']
    assert np.linalg.norm(actual@state@actual.conj().T-state, 'fro') < 1e-9
    # On the faithful support, sigma_1(A)=alpha_-1(A), not alpha_+1(A).
    density_eigenvalues, density_basis = np.linalg.eigh(state[:12, :12])
    modular = (density_basis*np.exp(1j*np.log(density_eigenvalues)))@density_basis.conj().T
    observable = np.zeros((12, 12)); observable[0, 1] = observable[1, 0] = 1
    lhs = modular@observable@modular.conj().T
    rhs = actual[:12, :12]@observable@actual[:12, :12].conj().T
    np.testing.assert_allclose(lhs, rhs, atol=1e-9, rtol=0)


@pytest.mark.parametrize('band', (1, 2, 3))
def test_same_register_kicks_and_complete_detectors_on_entangled_columns(source, band):
    basis = source['spectrum']['basis']
    address = source['pair_addresses'][band]
    x = np.zeros((16, 16), complex); x[0, address] = x[address, 0] = 1
    y = np.zeros((16, 16), complex); y[0, address], y[address, 0] = -1j, 1j
    for sign, name in ((1, f'kick_plus_{band}'), (-1, f'kick_minus_{band}')):
        expected = basis@expm(sign*1j*np.pi*x/4)@basis.T
        np.testing.assert_allclose(source['raw'][name], expected, atol=1e-9, rtol=0)
    zero_minus_pair = np.zeros((16, 16)); zero_minus_pair[0, 0] = 1
    zero_minus_pair[address, address] = -1
    for axis, observable in (('x', x), ('y', y)):
        detector = source['raw'][f'detector_{axis}_{band}']
        np.testing.assert_allclose(detector.conj().T@zero_minus_pair@detector,
                                   basis@observable@basis.T, atol=1e-9, rtol=0)
    # Every system/reference column continues linearly, including actual padding.
    rng = np.random.default_rng(180+band)
    columns = rng.normal(size=(16, 12))+1j*rng.normal(size=(16, 12))
    columns /= np.linalg.norm(columns)
    detector = source['raw'][f'detector_x_{band}']
    output = detector@source['raw']['clock1']@source['raw'][f'kick_plus_{band}']@columns
    expected = (source['targets'][f'detector_x_{band}']@source['targets']['clock1']
                @source['targets'][f'kick_plus_{band}']@columns)
    assert np.linalg.norm(output-expected, 'fro') < 1e-9
    probabilities = np.sum(abs(output)**2, axis=1)
    assert len(probabilities) == 16 and min(probabilities) > 0
    assert abs(sum(probabilities)-1) < 1e-9
    assert sum(probabilities[12:]) > .1


def test_wrong_phase_sign_and_padding_are_rejected_by_full_source_guards(source):
    name = 'modular_phase'
    wrong_sign = m.compiler.inverse(source['words'][name])
    with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
        m._execute(wrong_sign, source['targets'][name], source['guards'][name], 'wrong_sign')
    wrong_padding = source['raw']['clock1'].copy()
    wrong_padding[12, 12] *= -1
    with pytest.raises(ArithmeticError, match='complete-register numerical guard'):
        m._check_matrix(wrong_padding, source['targets']['clock1'],
                        source['guards']['clock1'], 'wrong_padding')
    dropped_padding = source['raw']['detector_x_2'][:12, :12]
    with pytest.raises(ArithmeticError, match='complete finite register'):
        m._check_matrix(dropped_padding, source['targets']['detector_x_2'],
                        source['guards']['detector_x_2'], 'dropped_padding')


def test_literal_clock_word_and_raw_routine_reuse_agree_without_polar_correction(source, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError('polar correction must not replace executed source maps')
    monkeypatch.setattr(m.compiler, 'stabilize_executed_unitary', forbidden)
    columns = np.eye(16, dtype=complex)[:, [0, 1, 12, 15]]
    literal = m.compiler.execute(m.clock_word(), 4, columns)
    np.testing.assert_allclose(literal, source['raw']['clock1']@columns, atol=1e-10, rtol=0)
    assert source['raw_source_routines_used_without_polar_correction']


def test_public_results_are_detached_and_api_rejects_bool_and_unselected_settings(source):
    changed = m.compiled_source()
    changed['raw']['clock1'][0, 0] = 17
    changed['words']['clock1'][0]['angle'] = 17
    fresh = m.compiled_source()
    np.testing.assert_array_equal(fresh['raw']['clock1'], source['raw']['clock1'])
    assert fresh['words']['clock1'][0] == source['words']['clock1'][0]
    for band in (True, 0, 4, 1.):
        with pytest.raises(ValueError):
            m.detector_word(band)
    with pytest.raises(ValueError):
        m.pair_word(1, True)
    with pytest.raises(ValueError):
        m.clock_word('fitted')
    with pytest.raises(ValueError):
        m.detector_word(1, 'Z')
    for value in (np.str_('X'), ['X'], 1, True):
        with pytest.raises(ValueError):
            m.detector_word(1, value)
        with pytest.raises(ValueError):
            m._pair_spectral_word(1, value, F(1, 4))
    with pytest.raises(ValueError):
        m.clock_word(np.str_('modular'))
    for wrong in (np.eye(16, dtype=bool), np.eye(16, dtype=object),
                  np.full((16, 16), '0')):
        with pytest.raises(ArithmeticError, match='complete finite register'):
            m._check_matrix(wrong, np.eye(16), F(1, 100), 'wrong_dtype')


def test_exact_basis_source_words_do_not_depend_on_platform_libm_angle_paths(source, monkeypatch):
    def forbidden(*_args):
        raise AssertionError('basis synthesis used platform libm angle arithmetic')
    for name in ('atan2', 'cos', 'sin', 'acos'):
        monkeypatch.setattr(m.math, name, forbidden)
    assert m._basis_word(source['spectrum']['basis_exact']) == source['words']['basis']
    certificate = source['spectrum']['exact_basis_certificate']
    assert certificate['symmetric_source_laplacian']
    assert certificate['all_twelve_eigenvector_equations_checked_exactly']
    assert certificate['within_each_band_orthogonality_and_positive_norms_checked_exactly']
