from fractions import Fraction as F

import pytest

from . import continuation_bounds, open_bounds as bounds, open_environment as env, protocol


def schedule():
    return env.Environment(tuple(env.PhaseNoise(after, 0, F(1, 64), 8)
                                 for after in (1, 4, 6, 9)))


@pytest.mark.parametrize('d', [2, 4])
def test_full_noise_budget_has_one_charge_per_stage_and_distinct_targets(d):
    program = protocol.continuation_program(d, 512)
    ordinary = continuation_bounds.program_bound(program)
    value = bounds.budget(program, schedule(), prepared_noiseless_error=F(1, 5))
    assert value['matched_phase_noise_occupation_certificate'] == ordinary['occupation_certificate']
    assert value['number_moment_growth_and_cutoff_leakage_from_noise_exact'] == '0'
    temporal = F(ordinary['ideal_temporal_instrument_error_upper'])
    quadrature = F(value['environment_budget']['ideal_finite_to_gaussian_half_diamond_upper'])
    assert F(value['ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper']) == temporal+quadrature
    reference = value['prepared_noiseless_reference_comparison']
    assert F(reference['ideal_noisy_finite_to_NOISELESS_infinite_error_upper']) == F(1, 5)+F((d-1)**2, 64)
    assert reference['target_contains_phase_noise'] is False
    service = value['physical_service_coverage']
    assert service['conditional_complete_error_upper'] is None
    assert len(service['missing_service_stage_keys']) == len(program.operations)+2+4
    assert value['finite_comparison_includes_arbitrary_work_reference_inputs']
    assert value['occupation_comparison_requires_declared_preparation']


def test_analytic_gaussian_refinement_decreases_without_claiming_source_execution():
    rows = [bounds.gaussian_refinement(d, (F(1, 64),)*4) for d in (2, 4, 16, 64)]
    errors = [F(row['ideal_exact_parameter_quadrature_error_upper']) for row in rows]
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert all(row['orders'][0] >= row['dimension'] for row in rows)
    assert [row['source_execution_available'] for row in rows] == [True, True, False, False]


@pytest.mark.parametrize('error', [True, -.1, F(-1), F(2), .2, None])
def test_supplied_prepared_certificate_contract(error):
    if error is None:
        assert bounds.budget(protocol.continuation_program(), env.Environment())[
            'prepared_noiseless_reference_comparison'] is None
    else:
        with pytest.raises(ValueError):
            bounds.budget(protocol.continuation_program(), env.Environment(), prepared_noiseless_error=error)


def test_zero_environment_retains_existing_temporal_and_occupation_errors():
    program = protocol.continuation_program()
    old = continuation_bounds.program_bound(program)
    new = bounds.budget(program, env.Environment())
    assert new['ideal_finite_source_to_continuous_GAUSSIAN_history_error_upper'] == old['ideal_temporal_instrument_error_upper']
    assert new['ideal_finite_source_to_infinite_GAUSSIAN_history_error_upper'] == old['complete_ideal_product_to_infinite_instrument_error_upper']


def test_sharp_detector_does_not_acquire_smooth_occupation_bound():
    with pytest.raises(ValueError):
        bounds.budget(protocol.sharp_program(), env.Environment())
