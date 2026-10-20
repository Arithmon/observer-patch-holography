"""Exact inverse moments, dependent-tape coverage and complete readout controls."""
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import random
import subprocess
import sys

import pytest

from . import inverse_clock as m


@pytest.fixture(scope='module')
def moments():
    return m.source_moments()


@pytest.fixture(scope='module')
def exact_fit(moments):
    mean, joint = moments
    return m.identify((mean, mean), tuple((x, x) for x in joint[:4]))


def enclosed(value, box):
    assert box[0] <= value <= box[1]


def test_native_inverse_recovers_closed_form_spectrum_without_using_forward_fit(moments, exact_fit):
    mean, joint = moments
    assert mean == F(5, 33)
    c = tuple(x-mean**2 for x in joint)
    assert c[:4] == (F(1, 3630), F(79, 326700), F(521, 2450250), F(3059, 16335000))
    assert c[0] != mean*(1-mean)  # Consecutive-event covariance, not variance.
    assert exact_fit['Hankel_determinant_interval'] == (F(1, 9702990000),)*2
    enclosed(F(5, 3), exact_fit['recurrence_trace_interval'])
    enclosed(F(311, 450), exact_fit['recurrence_product_interval'])
    with localcontext() as context:
        context.prec = 80
        root = Decimal(3).sqrt()
        roots = ((25+root)/30, (25-root)/30)
        weights = (Decimal(1)/7260+root/16335, Decimal(1)/7260-root/16335)
        for values, boxes in ((roots, exact_fit['transfer_root_intervals']),
                              (weights, exact_fit['positive_weight_intervals']),
                              (tuple(-x.ln() for x in roots), exact_fit['energy_intervals'])):
            for value, (lo, hi) in zip(values, boxes, strict=True):
                assert Decimal(lo.numerator)/Decimal(lo.denominator) <= value
                assert value <= Decimal(hi.numerator)/Decimal(hi.denominator)
    assert not exact_fit['existence_of_model_in_probability_box_certified']


def test_interval_inversion_covers_shared_mean_and_asymmetric_probability_errors():
    rng = random.Random(1029)
    successes = 0
    for _ in range(80):
        roots = (F(rng.randrange(82, 97), 100), F(rng.randrange(53, 76), 100))
        weights = (F(rng.randrange(1, 10), 1000), F(rng.randrange(1, 10), 1000))
        mean = F(rng.randrange(25, 65), 100)
        joints = tuple(mean**2+sum((w*r**n for w, r in zip(weights, roots)), F()) for n in range(4))
        def uncertain(value):
            return value-F(rng.randrange(1, 10), 10**13), value+F(rng.randrange(1, 10), 10**13)
        mean_box = uncertain(mean)
        fit = m.identify(mean_box, tuple(map(uncertain, joints)))
        if not fit['certified']:
            continue
        successes += 1
        for values, boxes in ((roots, fit['transfer_root_intervals']),
                              (weights, fit['positive_weight_intervals'])):
            for value, box in zip(values, boxes, strict=True):
                enclosed(value, box)
        for lag in (4, 7):
            actual = mean**2+sum((w*r**lag for w, r in zip(weights, roots)), F())
            enclosed(actual, m.predict_joint(mean_box, fit, lag))
    assert successes == 80


@pytest.mark.parametrize('roots,weights', [
    ((F(4, 5), F(4, 5)), (F(1, 100), F(1, 100))),
    ((F(4, 5), F(3, 5)), (F(1, 100), F(0))),
    ((F(4, 5), F(3, 5)), (F(1, 100), F(-1, 1000))),
    ((F(4, 5), F(2, 5)), (F(1, 100), F(1, 100))),
    ((F(11, 10), F(3, 5)), (F(1, 100), F(1, 100))),
])
def test_degenerate_nonpositive_and_out_of_regime_models_are_inconclusive(roots, weights):
    joint = tuple(F(1, 4)+sum((w*r**n for w, r in zip(weights, roots)), F()) for n in range(4))
    result = m.identify((F(1, 2),)*2, tuple((x, x) for x in joint))
    assert not result['certified']
    assert result['reason']


def test_broad_tubes_are_inconclusive_and_single_mode_matches_only_first_two(moments):
    mean, joint = moments
    assert not m.identify((mean-F(1, 100), mean+F(1, 100)),
                         tuple((x-F(1, 100), x+F(1, 100)) for x in joint[:4]))['certified']
    c = tuple(x-mean**2 for x in joint)
    root = c[1]/c[0]
    assert c[0]*root == c[1]
    assert c[2]-c[0]*root**2 == F(1, 2673000)


def test_frequency_and_quadratures_cannot_be_refit_independently(exact_fit):
    original = m.predict_readout(exact_fit, 1, 4, 'y')
    altered = deepcopy(exact_fit)
    altered['energy_intervals'] = ((F(10), F(10)),)*2
    assert m.predict_readout(altered, 1, 4, 'y') == original
    e1, e2 = (sum(box)/2 for box in exact_fit['energy_intervals'])
    for common_scale in (F(1, 100), F(7), F(10**18)):
        assert (common_scale*e2)/(common_scale*e1) == e2/e1
    assert 2*e2/e1 != e2/e1


def test_raw_array_bounds_preserve_padding_and_do_not_renormalize(exact_fit):
    result = m.predict_readout(exact_fit, 1, 0, raw_error=F(1, 1000))
    assert result['ideal_complete_probability_intervals'] == ((F(1), F(1)), (F(0), F(0)),
                                                              (F(0), F(0)), (F(0), F(0)))
    assert result['raw_complete_weight_intervals'][0][1] == F(501, 500)
    assert result['raw_complete_weight_intervals'][3] == (F(0), F(1, 500))
    assert result['event_weight_error'] == F(1, 500)


def test_complete_tape_counts_keep_zero_current_and_common_window_count():
    tape = (1, 0, 1, 1, 0, 1, 0, 0)
    result = m.summarize_tape(tape)
    assert result == {'windows': 3, 'mean_ones': 2, 'joint_ones': (1, 1, 2, 0, 1)}
    confidence = m.tape_confidence(result, F(1, 10))
    assert confidence['coverage_failure_upper'] == 1
    assert not confidence['overlapping_windows_treated_as_independent']
    assert confidence['stationary_union_variance_coefficient'] == 39


def test_marked_chain_window_covariances_obey_dependence_bound():
    reduction = m.source.clock.lumping(m.source.transfer.native_counts())
    def event_probability(indices):
        row = reduction['stationary']
        for i in range(max(indices)+1):
            matrix = reduction['marked'][1] if i in indices else reduction['matrix']
            row = tuple(sum(row[a]*matrix[a][b] for a in range(3)) for b in range(3))
        return sum(row)
    nonzero = False
    for length in range(1, 7):
        selected = {0, length-1}
        mean = event_probability(selected)
        covariances = []
        for lag in range(12):
            covariance = event_probability(selected | {x+lag for x in selected})-mean**2
            bound = F(1, 4)*(m.RATE_UPPER**(lag-length) if lag >= length else 1)
            assert abs(covariance) <= bound
            covariances.append(covariance)
        nonzero |= covariances[1] != 0
        count = 12
        variance = (count*covariances[0]+2*sum((count-lag)*covariances[lag]
                                              for lag in range(1, count)))/count**2
        assert 0 <= variance <= m.variance_coefficient(length)/count
    assert nonzero


def test_confidence_budget_adds_warmup_once_and_is_only_sufficient():
    radius = F(1, 10**12)
    count = 4333333333333333333333333334
    result = m.tape_confidence({'windows': count, 'mean_ones': 0, 'joint_ones': (0,)*5}, radius)
    assert result['coverage_failure_upper'] == F(39)/(count*radius**2)+m.burn_in.binary_bound(221)
    assert result['coverage_failure_upper'] <= F(1, 100)
    assert m.tape_confidence({'windows': count, 'mean_ones': 0, 'joint_ones': (0,)*5},
                             radius, warmup=0)['coverage_failure_upper'] == 1


def test_input_contract_rejects_bools_float_aliases_and_malformed_boxes(exact_fit):
    broken = F(1, 2)
    broken._numerator = 2
    for bad in (True, 0.2, broken):
        with pytest.raises(ValueError):
            m.interval((bad, F(1)))
    for bad in ([F(0), F(1)], (F(1), F(0)), (F(0),)):
        with pytest.raises(ValueError):
            m.interval(bad)
    for bad in ((0,)*5, (True,)*6, [0]*6):
        with pytest.raises(ValueError):
            m.summarize_tape(bad)
    for bad in (True, 0, 1.2, 33):
        with pytest.raises(ValueError):
            m.variance_coefficient(bad)
    for bad in (None, {}, {'certified': 1}, {**exact_fit, 'transfer_root_intervals': ((F(1),)*2,)*2}):
        with pytest.raises(ValueError):
            m.predict_readout(bad, 1)
    for arguments in ((True, 4, 'x'), (1, True, 'x'), (1, 4, 'z')):
        with pytest.raises(ValueError):
            m.predict_readout(exact_fit, *arguments)
    with pytest.raises(ValueError):
        m.tape_confidence({'windows': True, 'mean_ones': 0, 'joint_ones': (0,)*5}, F(1, 10))
    with pytest.raises(ValueError):
        m.tape_confidence({'windows': 2, 'mean_ones': 0, 'joint_ones': (0,)*5}, F(0))


def test_replay_complete_coherent_predictions_and_scope():
    result = m.run()
    assert len(result['coherent_readout_forecasts']) == 4
    assert all(row['all_padding_outcomes_retained'] for row in result['coherent_readout_forecasts'])
    assert result['declared_probability_tubes']['synthetic_source_centered_example']
    assert result['unused_lag_validation']['shift_by_1_over_1000_rejected']
    assert not result['scope']['matched_physical_record_dataset_available']
    assert not result['scope']['native_C_E_instrument_identified_from_classical_bits']


@pytest.mark.parametrize('parameter,heralds', [
    (F(0), 1), (F(1, 30), 1), (F(1, 30), 2), (F(1, 30), 100),
    (F(1), 30), (F(1), 1000), (F(7), 210),
])
def test_lazy_filter_gibbs_limit_bounds_state_and_free_entropy(parameter, heralds):
    import numpy as np
    bound = m.edge_gibbs_bound(parameter, heralds)
    eigenvalues = np.array([0., 5-np.sqrt(5), 6., 5+np.sqrt(5)])
    dimensions = np.array([1, 3, 5, 3])
    t = float(parameter)
    log_finite = 2*heralds*np.log1p(-t*eigenvalues/(2*heralds))
    weights, target = dimensions*np.exp(log_finite), dimensions*np.exp(-t*eigenvalues)
    success = weights.sum()/12
    state, gibbs = weights/weights.sum(), target/target.sum()
    entropy = -sum(state*np.log(state/dimensions))
    energy = state@eigenvalues
    deficit = np.log(target.sum())-(entropy-t*energy)
    relative_entropy = state@np.log(state/gibbs)
    assert success >= float(bound['all_success_probability_lower'])
    assert 0 <= np.max(-log_finite-t*eigenvalues) <= float(bound['spectral_log_error_upper'])+1e-13
    gap_errors = -log_finite-t*eigenvalues
    assert np.ptp(gap_errors) <= float(bound['spectral_log_error_upper'])+1e-13
    assert sum(abs(state-gibbs))/2 <= float(bound['trace_distance_to_edge_Gibbs_upper'])+1e-13
    assert abs(deficit-relative_entropy) <= 1e-13
    assert -1e-13 <= deficit <= float(bound['free_entropy_deficit_upper'])+1e-13
    assert not bound['physical_seconds_or_modular_time']


def test_prepared_edge_reference_is_coherent_and_keeps_every_herald_failure():
    row = m.prepared_edge_gibbs()
    assert row['executed_full_register_input_columns'] == 12
    assert row['sector_dimensions'] == [1, 3, 5, 3]
    assert row['representation_algebra_dimension'] == 44
    outcomes = row['retained_first_failure_and_success_probabilities']
    assert all(0 < p < 1 for p in outcomes)
    assert sum(outcomes) == 1
    assert row['normalized_success_matrix_element_0_1'] > 0
    assert row['normalized_success_purity'] > F(1, 12)
    readout = row['coherence_readout']
    assert len(readout['ideal_unnormalized_address_probabilities']) == 16
    assert sum(readout['ideal_unnormalized_address_probabilities']) == outcomes[-1]
    assert readout['normalized_zero_minus_one_contrast'] == F(260343, 43825510)
    assert readout['dephased_input_contrast'] == 0
    assert not readout['sector_projector_detector_compiled']
    assert row['one_filter_resources']['maximum_work_qubits'] == 15
    assert row['one_filter_resources']['service_failure_outcomes_retained']
    assert not row['event_or_scratch_failure_outputs_discarded_from_total_probability']
    assert not row['native_passive_equilibrium_prepares_this_state']
    assert not row['edge_and_wave_modular_generators_identified']
    assert not row['full_group_regular_or_compact_group_sewing_established']


def test_gibbs_family_rejects_nonprobabilistic_and_nonexact_preparations():
    for args in ((True, 1), (0.1, 1), (F(-1), 30), (F(1), 29),
                 (F(1, 30), True), (F(1, 30), 0)):
        with pytest.raises(ValueError):
            m.edge_gibbs_bound(*args)


def test_edge_sector_thermal_ratios_do_not_share_the_wave_square_root_clock():
    import numpy as np
    bands = np.array([5-np.sqrt(5), 6., 5+np.sqrt(5)])
    # A common scale would require lambda_j = b*sqrt(lambda_j) for all j.
    assert np.ptp(bands/np.sqrt(bands)) > .5
    # Multiplicity-free trace weights differ from regular-group d_R^2 weights.
    dimensions = np.array([1, 3, 5, 3])
    assert dimensions.sum() == 12
    assert (dimensions**2).sum() == 44 < 60
    assert not np.allclose(dimensions/12, dimensions**2/44)


def test_edge_gibbs_is_central_on_group_image_but_not_on_full_register_algebra():
    import numpy as np
    from ..common_source.slow_geometry import source
    laplacian = np.zeros((12, 12), dtype=int)
    for a, b in m.continuation.quantum_mixing.seams():
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    for rotation in source().rotations():
        matrix = np.eye(12, dtype=int)[:, rotation]
        assert np.array_equal(laplacian@matrix, matrix@laplacian)
    observable = np.zeros((12, 12), dtype=int)
    observable[0, 1] = observable[1, 0] = 1
    assert np.linalg.norm(laplacian@observable-observable@laplacian) > 0


def test_all_repo_local_runtime_imports_are_pinned():
    script = '''
import json, sys
from pathlib import Path
from model.predictions import inverse_clock as m
m.run()
missing=[]
for module in list(sys.modules.values()):
    path=getattr(module, '__file__', None)
    if path is None: continue
    try: relative=Path(path).resolve().relative_to(m.ROOT).as_posix()
    except ValueError: continue
    if relative.startswith('.'): continue
    if relative.endswith('.py') and relative not in m.SOURCE_PATHS: missing.append(relative)
print(json.dumps(sorted(set(missing))))
'''
    result = subprocess.run([sys.executable, '-c', script], cwd=m.ROOT, capture_output=True,
                            text=True, check=True)
    assert json.loads(result.stdout) == []
