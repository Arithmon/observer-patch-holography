"""Joint-parameter, outward-envelope and single-state regression controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
import math

import mpmath as mp
import pytest

from . import shared_calibration as shared
from .data import load_observations
from .executable_calibration import AXIS_KEYS


def _means(parameters, keys, extent=4):
    alpha, thermal = parameters['alpha_real'], parameters['thermal_occupation']
    values = {}
    for i, j in keys:
        u, v = extent*i/10, extent*j/10
        radius, phase = math.exp(-(thermal+.5)*(u*u+v*v)), 2*alpha*v
        values[i, j] = radius*complex(math.cos(phase), math.sin(phase))
    return values


def test_exact_characteristic_boxes_cover_independent_complex_values():
    ctx = mp.mp.clone()
    ctx.dps = 90
    cells = (((F(2, 25), F(1, 8)), (F(-22, 7), F(22, 7))),
             ((F(3, 32), F(7, 64)), (F(3, 4), F(4, 5))),
             ((F(1, 8), F(1, 8)), (F(-3), F(-3))))
    exact = lambda value: ctx.mpf(value.numerator)/value.denominator
    for cell in cells:
        for key in ((0, 0), (1, 1), (3, 7), (10, 10)):
            box = shared.characteristic_box(cell, key)
            for damping in (cell[0][0], sum(cell[0])/2, cell[0][1]):
                for theta in (cell[1][0], sum(cell[1])/2, cell[1][1]):
                    i, j = key
                    value = ctx.exp(-exact(damping)*(i*i+j*j)+1j*exact(theta)*j)
                    for part, number in (('real', value.real), ('imag', value.imag)):
                        assert exact(box[part][0]) <= number <= exact(box[part][1])


def test_axis_cover_retains_true_state_without_having_validation_inputs():
    parameters = {'alpha_real': .75, 'thermal_occupation': .125}
    axes = _means(parameters, AXIS_KEYS)
    cells = shared.axis_cover(axes, extent=4, width=F(1, 100), depth=6)
    a, theta = (F(parameters['thermal_occupation'])+F(1, 2))*F(4, 10)**2, 2*F(parameters['alpha_real'])*F(4, 10)
    assert any(cell[0][0] <= a <= cell[0][1] and cell[1][0] <= theta <= cell[1][1] for cell in cells)
    with pytest.raises(ValueError, match='declared observations'):
        shared.axis_cover({**axes, (1, 1): 0j}, width=F(1, 100))


def test_simultaneous_intersections_cannot_be_replaced_by_marginal_envelopes():
    # Each proposed observation matches one of two entire states. The states
    # disagree, so there is no single witness for their incompatible mixture.
    cells = (((F(1, 10), F(1, 10)), (F(3, 4), F(3, 4))),
             ((F(1, 10), F(1, 10)), (F(-3, 4), F(-3, 4))))
    keys = ((1, 1), (1, 2))
    envelopes = shared.forward_envelopes(cells, keys)
    observed = {}
    for index, key in enumerate(keys):
        exact = shared.characteristic_box(cells[index], key)
        observed[key] = {part: (lo-F(1, 1000), hi+F(1, 1000)) for part, (lo, hi) in exact.items()}
        assert all(envelopes[key][part][0] <= (lo+hi)/2 <= envelopes[key][part][1]
                   for part, (lo, hi) in observed[key].items())
    assert not any(shared._compatible(cell, observed) for cell in cells)


def test_one_axis_selected_witness_is_checked_on_every_quadrature():
    parameters = {'alpha_real': .75, 'thermal_occupation': .125}
    axes, validation = (_means(parameters, keys) for keys in (AXIS_KEYS, shared.VALIDATION_KEYS))
    result = shared.common_witness(parameters, axes, validation, extent=4,
                                    axis_width=F(1, 1000), validation_width=F(1, 1000))
    assert result['all_exact_prediction_intervals_inside_simultaneous_boxes']
    bad = dict(validation)
    bad[3, 7] += .1j
    result = shared.common_witness(parameters, axes, bad, extent=4,
                                    axis_width=F(1, 1000), validation_width=F(1, 1000))
    assert not result['all_exact_prediction_intervals_inside_simultaneous_boxes']
    assert result['parameters'] == parameters


def test_empty_and_unbounded_calibrations_have_distinct_results():
    parameters = {'alpha_real': .75, 'thermal_occupation': .125}
    axes = _means(parameters, AXIS_KEYS)
    bad = dict(axes)
    bad[0, 0] = .5+0j
    assert shared.axis_cover(bad, width=F(1, 100), depth=2) == ()
    with pytest.raises(ValueError, match='do not bound'):
        shared.axis_cover(axes, width=F(2), depth=2)


@pytest.mark.parametrize('bad', [True, 3.0, 3, F(4)])
def test_coordinate_convention_requires_original_published_extent(bad):
    axes = _means({'alpha_real': .75, 'thermal_occupation': .125}, AXIS_KEYS)
    with pytest.raises(ValueError, match='published extent'):
        shared.axis_cover(axes, extent=bad)


def test_source_link_rejects_a_different_preparation_or_dropped_outcome():
    observations = load_observations()
    receipt = json.loads((shared.FOLDER/'executable_results.json').read_text())
    original = receipt['runs'][0]
    parameters = original['calibration']['parameters']
    good = shared.source_witness_link(original, parameters, observations['validation'], F(3, 10))
    assert good['same_binary_preparation_and_fixed_detector_dictionary']
    for mutation in ('thermal', 'detector', 'coordinate', 'phase'):
        bad = deepcopy(original)
        if mutation == 'thermal':
            bad['source_execution']['preparation']['thermal']['mean_thermal_occupation_parameter'] += .1
        elif mutation == 'detector':
            bad['source_execution']['records'][0]['both_detector_outcomes_retained'] = False
        elif mutation == 'coordinate':
            bad['predictions'][0]['u'] += .1
        else:
            bad['source_execution']['records'][0]['theta'] += .1
        with pytest.raises(ValueError, match='dictionary|grid'):
            shared.source_witness_link(bad, parameters, observations['validation'], F(3, 10))


def test_original_axis_confidence_event_is_shared_by_both_extent_hypotheses():
    for statistics in (42, 200):
        width, certificate = shared.sampling_contract(statistics)
        assert F(certificate['familywise_error_upper']) <= F(1, 40)
        assert F(certificate['width_exact']) == width


def test_sampling_radius_does_not_use_platform_libm_transcendentals(monkeypatch):
    expected = shared.sampling_contract(200)
    def forbidden(*args, **kwargs):
        raise AssertionError('platform libm selected a receipt radius')
    monkeypatch.setattr(math, 'log', forbidden)
    monkeypatch.setattr(math, 'sqrt', forbidden)
    assert shared.sampling_contract(200) == expected


def test_retained_shared_comparison_is_exact_and_has_actual_common_witnesses():
    actual = shared.build()
    retained = json.loads(shared.OUTPUT.read_text())
    assert retained == actual
    assert actual['dataset']['validation_complex_points'] == 100
    assert len(actual['runs']) == 2
    for run in actual['runs']:
        assert run['axis_cover_cells'] >= run['cells_not_excluded_by_all_validation_boxes'] > 0
        assert len(run['prediction']['envelopes']) == 100
        assert run['verdict'] == 'one_axis_selected_state_is_simultaneously_compatible'
        assert run['shared_state_witness']['all_exact_prediction_intervals_inside_simultaneous_boxes']
        assert run['same_state_native_source_link']['all_retained_source_means_inside_validation_boxes']
    assert actual['scope']['empirical_confirmation'] is False
    assert actual['scope']['physical_clock_or_interacting_action_tested'] is False
    assert sum((F(row['familywise_error_upper']) for row in actual['sampling'].values()), F()) <= F(1, 20)
