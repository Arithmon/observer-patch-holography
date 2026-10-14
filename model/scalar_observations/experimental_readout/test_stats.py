"""Independent interval, holdout and sampling-contract checks."""
from copy import deepcopy
from fractions import Fraction as F
import itertools
import math

import pytest

from model.scalar_observations.experimental_readout import stats as s


def test_hoeffding_counts_real_components_and_uses_minus_one_to_one_range():
    assert s.hoeffding_width(42) >= math.sqrt(2*math.log(84/.025)/250)
    assert s.hoeffding_width(42) == pytest.approx(math.sqrt(2*math.log(84/.025)/250))
    assert s.hoeffding_width(200) > s.hoeffding_width(42)


@pytest.mark.parametrize('statistics', [42, 200])
def test_width_has_independent_exact_positive_taylor_certificate(statistics):
    width = s.hoeffding_width(statistics)
    certificate = s.concentration_certificate(statistics, 250, width)
    argument = 250*F(width)**2/2
    assert F(certificate['width_exact']) == F(width)
    assert F(certificate['exponential_argument']) == argument
    assert certificate['taylor_degree'] == 128
    # Construct each term separately, independently of the implementation's
    # multiplicative recurrence, and check its recorded downward rounding.
    total = sum((argument**k/math.factorial(k) for k in range(129)), F(0))
    lower = F(certificate['exp_lower_partial_sum_rounded_down'])
    assert 0 <= total-lower < F(1, 2**64)
    assert F(certificate['familywise_error_upper']) == 2*statistics/lower
    assert F(certificate['familywise_error_upper']) <= F(1, 40)


def test_float_width_estimate_is_only_a_starting_guess(monkeypatch):
    original_sqrt = math.sqrt
    # Deliberately perturb the untrusted initial guess. Exact checks must
    # restore coverage independently of the logarithm/square-root estimate.
    monkeypatch.setattr(s.math, 'sqrt', lambda value: .9*original_sqrt(value))
    width = s.hoeffding_width(42)
    certificate = s.concentration_certificate(42, 250, width)
    assert F(certificate['familywise_error_upper']) <= F(1, 40)


def test_extreme_failure_budget_terminates_with_certified_conservative_width():
    width = s.hoeffding_width(2, delta=1e-300)
    certificate = s.concentration_certificate(2, 250, width)
    assert F(certificate['familywise_error_upper']) <= F('1e-300')
    assert s.mean_box(0j, width) == {'real': (-1., 1.), 'imag': (-1., 1.)}


@pytest.mark.parametrize('width', [-1, float('nan'), float('inf'), True])
def test_concentration_certificate_rejects_invalid_width(width):
    with pytest.raises(ValueError):
        s.concentration_certificate(42, 250, width)


@pytest.mark.parametrize('statistics,samples,delta', [
    (0, 250, .025), (True, 250, .025), (1., 250, .025),
    (2, 0, .025), (2, -1, .025), (2, True, .025), (2, 250., .025),
    (2, 250, 0), (2, 250, 1), (2, 250, -1),
    (2, 250, float('nan')), (2, 250, float('inf')), (2, 250, True)])
def test_invalid_sampling_contract_is_rejected(statistics, samples, delta):
    with pytest.raises(ValueError):
        s.hoeffding_width(statistics, samples, delta)


@pytest.mark.parametrize('mean', [1.0001, -1.0001j, float('nan'), complex(0, float('inf')), True, '0'])
def test_invalid_observed_component_mean_is_rejected(mean):
    with pytest.raises(ValueError):
        s.calibrate_axes({(1, 0): mean})


def test_independent_component_means_need_not_lie_in_unit_disk():
    # Different real/imaginary samples can put their noisy combined mean
    # outside the unit disk; component range is the applicable data check.
    calibration = s.calibrate_axes({(1, 0): .9+.9j, (0, 1): .8+.8j})
    assert s.predict(calibration, [(1, 1)])[(1, 1)]['mean'] == pytest.approx(1.44j)


def test_box_product_contains_all_complex_vertex_products_with_mixed_signs():
    left = {'real': (-.7, .2), 'imag': (-.1, .9)}
    right = {'real': (-.8, .6), 'imag': (-.4, .3)}
    box = s.box_product(left, right)
    products = [complex(a, b)*complex(c, d)
                for a, b, c, d in itertools.product(left['real'], left['imag'],
                                                   right['real'], right['imag'])]
    for z in products:
        assert box['real'][0] <= z.real <= box['real'][1]
        assert box['imag'][0] <= z.imag <= box['imag'][1]
    assert box['real'] == pytest.approx((min(z.real for z in products), max(z.real for z in products)))
    assert box['imag'] == pytest.approx((min(z.imag for z in products), max(z.imag for z in products)))


def test_box_conjugation_changes_phase_with_correct_interval_order():
    box = {'real': (-.2, .7), 'imag': (-.3, .8)}
    assert s.conjugate_box(box) == {'real': (-.2, .7), 'imag': (-.8, .3)}
    assert s.conjugate_box(s.conjugate_box(box)) == box


def test_mean_intervals_respect_component_range_without_clipping_mean():
    assert s.mean_box(1-1j, .25) == {'real': (pytest.approx(.75), 1.),
                                    'imag': (-1., pytest.approx(-.75))}
    assert s.mean_box(.2+.3j, 0) == {'real': (.2, .2), 'imag': (.3, .3)}


@pytest.mark.parametrize('box', [
    {'real': (1, 0), 'imag': (0, 1)},
    {'real': (0, float('inf')), 'imag': (0, 1)},
    {'real': (0, 1)}, {'real': (0, 1), 'imag': (0, 1), 'extra': 0}])
def test_malformed_boxes_rejected(box):
    with pytest.raises(ValueError):
        s.box_product(box, {'real': (0, 1), 'imag': (0, 1)})


def test_calibration_copy_and_off_axis_rejection_protect_withheld_data():
    source = {(0, 0): .99+.01j, (1, 0): .7+.02j, (0, 1): .3-.6j}
    calibration = s.calibrate_axes(source)
    frozen = deepcopy(calibration)
    before = s.predict(calibration, [(1, 1)])
    validation = {(1, 1): .2-.4j}
    s.score(before, validation)
    validation[(1, 1)] = -1+1j
    s.score(before, validation)
    source[(1, 0)] = -1
    assert calibration == frozen
    assert s.predict(calibration, [(1, 1)]) == before
    with pytest.raises(ValueError, match='only nonnegative coordinate axes'):
        s.calibrate_axes({**source, **validation})


def test_21_axis_records_count_42_real_statistics_including_origin():
    means = {(0, 0): .99}
    means.update({(i, 0): .5 for i in range(1, 11)})
    means.update({(0, j): .5 for j in range(1, 11)})
    calibration = s.calibrate_axes(means)
    assert calibration['real_statistics'] == 42
    predicted = s.predict(calibration, [(i, j) for i in range(1, 11) for j in range(1, 11)])
    assert len(predicted) == 100
    result = s.score(predicted, {key: .25 for key in predicted})
    assert result['real_statistics'] == 200
    assert result['complex_rmse'] == result['complex_mae'] == 0
    assert result['coverage'] == 'conditional_on_independent_bounded_raw_observations'
    assert result['empirical_confirmation'] is False
    assert result['calibration_failure_budget_included'] is False


def test_predict_all_quadrants_uses_conjugation_without_additional_data():
    a, b = .6+.1j, .4-.5j
    calibration = s.calibrate_axes({(1, 0): a, (0, 2): b})
    predicted = s.predict(calibration, [(1, 2), (-1, 2), (1, -2), (-1, -2)])
    assert predicted[(1, 2)]['mean'] == a*b
    assert predicted[(-1, 2)]['mean'] == a.conjugate()*b
    assert predicted[(1, -2)]['mean'] == a*b.conjugate()
    assert predicted[(-1, -2)]['mean'] == (a*b).conjugate()


def test_origin_is_not_used_for_response_renormalization():
    calibration = s.calibrate_axes({(0, 0): .2+.1j, (1, 0): .6, (0, 1): -.2j})
    prediction = s.predict(calibration, [(0, 0), (1, 0), (1, 1)])
    assert prediction[(0, 0)]['mean'] == 1
    assert prediction[(1, 0)]['mean'] == .6
    assert prediction[(1, 1)]['mean'] == -.12j
    assert calibration['axis_means'][(0, 0)] == .2+.1j


@pytest.mark.parametrize('keys', [[(1, 1), (1, 1)], [(2, 1)], [], [(True, 1)], [(1., 1)]])
def test_missing_or_ambiguous_prediction_coordinates_rejected(keys):
    calibration = s.calibrate_axes({(1, 0): .6, (0, 1): .4})
    with pytest.raises(ValueError):
        s.predict(calibration, keys)


def test_score_uses_complex_distance_and_rejects_disjoint_simultaneous_boxes():
    predictions = {(1, 1): {'mean': 0j, 'box': {'real': (0., 0.), 'imag': (0., 0.)}},
                   (1, 2): {'mean': .2j, 'box': {'real': (0., 0.), 'imag': (.2, .2)}}}
    result = s.score(predictions, {(1, 1): .3+.4j, (1, 2): .2j}, samples=10_000)
    assert result['complex_rmse'] == pytest.approx(math.sqrt(.25/2))
    assert result['complex_mae'] == pytest.approx(.25)
    assert result['box_violations'] == [(1, 1)]
    assert result['rows'][(1, 1)]['disjoint_components'] == ['real', 'imag']
    assert result['verdict'] == 'conditional_box_violation'
    assert result['empirical_confirmation'] is False
    assert result['provenance_verified'] is False


def test_score_requires_exactly_the_withheld_coordinate_set():
    predictions = {(1, 1): {'mean': 0j, 'box': {'real': (0., 0.), 'imag': (0., 0.)}}}
    with pytest.raises(ValueError, match='identical coordinates'):
        s.score(predictions, {(1, 2): 0j})


def test_nominal_gaussian_has_specified_phase_factorization_and_normalization():
    u, v = .7, 1.1
    z = s.nominal_gaussian(u, v)
    assert z == pytest.approx(complex(math.cos(-math.sqrt(2)*1.5*v),
                                     math.sin(-math.sqrt(2)*1.5*v))*math.exp(-.54*(u*u+v*v)))
    assert z == pytest.approx(s.nominal_gaussian(u, 0)*s.nominal_gaussian(0, v))
    assert s.nominal_gaussian(-u, -v) == pytest.approx(z.conjugate())
    assert s.nominal_gaussian(0, 0) == 1
    assert s.nominal_gaussian(0, .1).imag < 0


@pytest.mark.parametrize('kwargs', [{'u': True}, {'v': float('nan')},
                                  {'thermal': -.1}, {'amplitude': float('inf')}])
def test_nominal_gaussian_rejects_invalid_physical_parameters(kwargs):
    with pytest.raises(ValueError):
        s.nominal_gaussian(**{'u': 0, 'v': 0, **kwargs})
