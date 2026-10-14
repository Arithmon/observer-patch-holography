"""Original-input controls for fixed-arithmetic source-word generation."""
import math

import mpmath
import numpy as np
import pytest
from scipy.linalg import expm

from model.observer.source_scalar import compiler as c
from . import canonical_source as canonical
from . import executable_source as source


@pytest.mark.parametrize('dimension,alpha', [(2, 0.), (4, -.37), (8, .23), (16, -.11)])
def test_canonical_displacement_executes_independent_original_input(dimension, alpha):
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    target = expm(alpha * (ladder.T - ladder))
    word = canonical.displacement_word(dimension, alpha)
    actual = c.execute(word, dimension.bit_length() - 1)
    assert np.linalg.norm(actual - target, 2) < 2e-11
    assert c.resources(word)['proper_code_service_events'] == 2 * len(word)
    assert c.resources(word)['native_pulses'] == {2: 5, 4: 82, 8: 1046, 16: 11638}[dimension]


@pytest.mark.parametrize('dimension', [2, 4, 8, 16])
def test_canonical_coordinate_basis_retains_full_phase_and_inverse(dimension):
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    values, target = np.linalg.eigh((ladder + ladder.T) / 2)
    for column in target.T:
        if column[np.flatnonzero(abs(column) > 1e-12)[0]] < 0:
            column *= -1
    word = canonical.coordinate_word(dimension)
    actual = c.execute(word, dimension.bit_length() - 1)
    backward = c.execute(c.inverse(word), dimension.bit_length() - 1)
    assert np.linalg.norm(actual - target, 2) < 2e-11
    assert np.linalg.norm(backward - target.T, 2) < 2e-11
    assert np.linalg.norm(actual @ np.diag(values) @ actual.conj().T
                          - (ladder + ladder.T) / 2, 2) < 2e-11


def test_ambient_precision_cannot_change_canonical_words_or_context():
    old_precision = mpmath.mp.prec
    words = []
    try:
        for digits in (15, 37, 113):
            mpmath.mp.dps = digits
            context_before = mpmath.mp.prec
            canonical._coordinate_data.cache_clear()
            words.append((canonical.displacement_word(4, -.3125),
                          canonical.coordinate_word(4),
                          canonical.quadrature_words(4, .2, -.7, math.pi / 2)))
            assert mpmath.mp.prec == context_before
    finally:
        mpmath.mp.prec = old_precision
        canonical._coordinate_data.cache_clear()
    assert words[0] == words[1] == words[2]
    assert len({source._word_hash(word_group) for word_group in words}) == 1


def test_word_generation_does_not_consume_scipy_or_numpy_matrix_targets(monkeypatch):
    expected = (canonical.displacement_word(4, .23), canonical.coordinate_word(4),
                canonical.quadrature_words(4, -.2, .4, -.17))
    canonical._coordinate_data.cache_clear()

    def unavailable(*args, **kwargs):
        raise AssertionError('platform numerical target entered the word producer')

    monkeypatch.setattr(np.linalg, 'eigh', unavailable)
    monkeypatch.setattr(source, 'expm', unavailable)
    assert (canonical.displacement_word(4, .23), canonical.coordinate_word(4),
            canonical.quadrature_words(4, -.2, .4, -.17)) == expected


def test_distinct_original_binary64_angles_are_not_rounded_to_common_input():
    first = .23
    second = np.nextafter(first, math.inf)
    assert canonical.displacement_word(2, first) != canonical.displacement_word(2, second)


def test_fixed_source_word_hashes_bind_all_canonical_generators():
    _, _, rotation, detector = canonical.quadrature_words(8, .2, -.7, .17)
    words = [canonical.displacement_word(8, -.3125), canonical.coordinate_word(8),
             rotation, detector]
    assert [source._word_hash(word) for word in words] == [
        '5246a16e563e422351944215eb2f056fb27d44490e16942ffed5bcfdd813d596',
        '4131e588ae539e5de0487bc3fda08f75ff6eb4cdf79192211e10567f807d0428',
        '03b8e2b283fc545c43da235f41b6096b7aa20255569407d065b0261f0724befd',
        '1a5af757ef33e325bb41a04b9a3ec6cb73753dd5457ddb692b746f7b224a9046',
    ]


@pytest.mark.parametrize('u,v,theta', [(0., -1., 0.), (-0., -1., 0.),
                                      (.2, -.7, .3), (-.5, .4, -math.pi / 2)])
def test_canonical_setting_words_match_independent_phase_targets(u, v, theta):
    dimension = 4
    ladder = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    values = np.linalg.eigvalsh((ladder + ladder.T) / 2)
    phi, radius, rotation, detector = canonical.quadrature_words(dimension, u, v, theta)
    expected_phi = math.atan2(-u, v)
    expected_radius = math.hypot(u, v)
    expected_angles = expected_radius * values + theta / 2
    assert radius == pytest.approx(expected_radius, abs=1e-15)
    assert np.linalg.norm(c.execute(rotation, 2) - np.diag(
        np.exp(-1j * expected_phi * np.arange(dimension))), 2) < 2e-12
    assert np.linalg.norm(c.execute(detector, 3) - np.diag(
        np.exp(1j * np.r_[-expected_angles, expected_angles])), 2) < 2e-12
    if u == 0 and v < 0:
        assert phi == math.pi
        assert canonical.quadrature_words(dimension, 0., v, theta) == (
            canonical.quadrature_words(dimension, -0., v, theta))


@pytest.mark.parametrize('dimension', [True, 1, 3, 128, 4.])
def test_canonical_dimension_validation(dimension):
    with pytest.raises(ValueError):
        canonical.displacement_word(dimension, .1)
    with pytest.raises(ValueError):
        canonical.coordinate_word(dimension)
    with pytest.raises(ValueError):
        canonical.quadrature_words(dimension, .1, .2, .3)


@pytest.mark.parametrize('value', [True, np.bool_(False), math.nan, math.inf, '0.2'])
def test_canonical_original_parameter_validation(value):
    with pytest.raises(ValueError):
        canonical.displacement_word(2, value)
    for index in range(3):
        arguments = [.1, .2, .3]
        arguments[index] = value
        with pytest.raises(ValueError):
            canonical.quadrature_words(2, *arguments)


def test_independent_target_cannot_be_replaced_by_canonical_word(monkeypatch):
    source._displacement.cache_clear()
    monkeypatch.setattr(source, 'expm', lambda value: np.eye(len(value)))
    try:
        with pytest.raises(ArithmeticError, match='independent target'):
            source._displacement(2, .234567)
    finally:
        source._displacement.cache_clear()


@pytest.mark.parametrize('dimension', [2, 64])
def test_wrong_canonical_word_cannot_replace_source_execution(monkeypatch, dimension):
    source._displacement.cache_clear()
    monkeypatch.setattr(canonical, 'displacement_word',
                        lambda dimension, alpha: [c.pulse('I', math.pi / 2, [0])])
    try:
        with pytest.raises(ArithmeticError, match='independent target'):
            source._displacement(dimension, .234567)
    finally:
        source._displacement.cache_clear()


def test_canonical_policy_does_not_claim_interval_or_physical_precision():
    first = canonical.policy()
    first['decimal_digits'] = 1
    assert canonical.policy()['decimal_digits'] == 80
    assert not canonical.policy()['source_primitive_or_execution_changed']
    assert not canonical.policy()['outward_roundoff_or_physical_error_certificate']
