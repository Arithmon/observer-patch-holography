"""Retained same-observable comparisons for an interacting fixed mode band."""
from fractions import Fraction as F
import math

from ..observer.source_scalar import compiler
from . import bounds, protocol, source
from .lattice import GalerkinParameters, SevenModeLattice


GRIDS = (5, 9, 17, 33, 65)
DIMENSIONS = (2, 4)
STEPS = 64
REFINEMENTS = (16, 32, STEPS)
SETTINGS = protocol.Settings()
NUMERICAL_GUARD = 2e-8


def _histories(lattice, **kwargs):
    return [protocol.history(lattice, settings=SETTINGS, sign=sign, **kwargs)
            for sign in (-1, 1)]


def _difference(left, right):
    return max(protocol.tv(a, b) for a, b in zip(left, right))


def _mass(histories):
    return [math.fsum(row['probability'] for row in history) for history in histories]


def _local_checks(lattice):
    p, d = lattice.parameters, lattice.dimensions[0]
    dt = float(SETTINGS.interval)/STEPS
    return {
        'onsite': {name: source.onsite(d, mode, dt/2, p)[1]
                   for name, mode in (('zero', 0), ('first_harmonic', 1))},
        'pair': {str(coefficient): source.pair(d, d, coefficient, dt, p.omega)[1]
                 for coefficient in sorted({lattice.pair_coefficient(i, j)
                     for i in range(7) for j in range(i+1, 7)})},
        'pointer': source.instrument(d, SETTINGS.eta, SETTINGS.theta, p.omega)[1],
        'signed_preparations': {str(sign): source.feedback(d, sign*SETTINGS.initial_kick, p.omega)[1]
                                for sign in (-1, 1)},
        'feedback': source.feedback(d, SETTINGS.feedback_kick, p.omega)[1],
    }


def occupation_diagnostics():
    """Independent larger-code references; no source or limiting-error claim.

    The stronger setting exposed a cutoff problem during protocol design.
    Retain it alongside the selected weaker setting, without treating either
    finite difference as an error bound to unrestricted oscillators.
    """
    strong = protocol.Settings(eta=2., initial_kick=1., feedback_kick=1.)
    rows = []
    for name, settings in (('principal', SETTINGS), ('strong_detector_stress', strong)):
        codes = []
        previous = None
        base = None
        for dimensions in ((4,)*7, (4, 8, 4, 4, 4, 4, 4), (6,)*7):
            lattice = SevenModeLattice(dimensions)
            histories = [protocol.history(lattice, settings=settings, sign=sign,
                         source_execution=False) for sign in (-1, 1)]
            if base is None:
                base = histories
            codes.append({'dimensions': list(dimensions), 'tensor_dimension': lattice.dimension,
                'reference_histories_by_native_sign': histories,
                'tv_from_uniform_four': _difference(histories, base),
                'tv_from_previous': None if previous is None else _difference(histories, previous),
                'source_executed': False})
            previous = histories
        rows.append({'setting': name, 'eta': settings.eta,
            'initial_kick': settings.initial_kick, 'feedback_kick': settings.feedback_kick,
            'codes': codes})
    return {'spatial_gradient': 'exact within the same seven-mode band', 'settings': rows,
        'selection_used_internal_numerical_exploration': True,
        'selection_used_experimental_data': False,
        'finite_differences_are_bounds_to_unrestricted_oscillators': False,
        'all_occupation_dimensions_source_executed': False}


def run():
    levels = []
    previous_target = None
    times = (SETTINGS.interval, SETTINGS.interval)
    for d in DIMENSIONS:
        dimensions = (d,)*7
        target_lattice = SevenModeLattice(dimensions)
        target = _histories(target_lattice, source_execution=False)
        product = bounds.splitting_history_bound(dimensions, times, STEPS)
        grids = []
        for grid in GRIDS:
            lattice = SevenModeLattice(dimensions, GalerkinParameters(grid=grid))
            reference = _histories(lattice, source_execution=False)
            actual = _histories(lattice, steps=STEPS)
            spatial = bounds.spatial_history_bound(dimensions, grid, sum(times))
            total = min(F(1), spatial+product)
            split_difference = _difference(actual, reference)
            spatial_difference = _difference(reference, target)
            total_difference = _difference(actual, target)
            grids.append({
                'grid_per_axis': grid, 'kappa': lattice.parameters.kappa,
                'steps_per_interval': STEPS,
                'source_histories_by_native_sign': actual,
                'reference_histories_by_native_sign': reference,
                'probability_sums': _mass(actual),
                'source_reference_tv': split_difference,
                'spatial_reference_tv': spatial_difference,
                'source_to_same_occupation_exact_gradient_tv': total_difference,
                'ideal_product_history_tv_upper_exact': str(product),
                'ideal_spatial_history_tv_upper_exact': str(spatial),
                'combined_ideal_history_tv_upper_exact': str(total),
                'combined_ideal_history_tv_upper_float': float(total),
                'source_error_within_budget_plus_numerical_guard': split_difference <= float(product)+NUMERICAL_GUARD,
                'spatial_error_within_budget_plus_numerical_guard': spatial_difference <= float(spatial)+NUMERICAL_GUARD,
                'combined_error_within_budget_plus_numerical_guard': total_difference <= float(total)+NUMERICAL_GUARD,
            })
        finest = SevenModeLattice(dimensions, GalerkinParameters(grid=GRIDS[-1]))
        best = grids[-1]['source_histories_by_native_sign']
        refinements = []
        for steps in REFINEMENTS:
            actual = best if steps == STEPS else _histories(finest, steps=steps)
            budget = bounds.splitting_history_bound(dimensions, times, steps)
            difference = _difference(actual, grids[-1]['reference_histories_by_native_sign'])
            refinements.append({'steps_per_interval': steps,
                'source_histories_by_native_sign': actual,
                'reference_tv': difference, 'ideal_history_tv_upper_exact': str(budget),
                'within_budget_plus_numerical_guard': difference <= float(budget)+NUMERICAL_GUARD})
        free_lattice = SevenModeLattice(dimensions, GalerkinParameters(grid=GRIDS[-1], coupling=0))
        free = _histories(free_lattice, steps=STEPS)
        no_feedback = _histories(finest, steps=STEPS, feedback=False)
        reset = _histories(finest, steps=STEPS, reset=True)
        interaction_difference = _difference(best, free)
        comparison_budget = product+bounds.splitting_history_bound(dimensions, times, STEPS, coupling=0)
        records = protocol.native_records(best)
        detached = protocol.native_records(best, feedback=False)
        levels.append({
            'occupation_per_mode': d, 'tensor_dimension': finest.dimension,
            'field_and_pointer_qubits': 7*(d.bit_length()-1)+1,
            'exact_gradient_reference_histories_by_native_sign': target,
            'spatial_refinements': grids,
            'temporal_refinements_at_finest_grid': refinements,
            'calibration_consistency': protocol.detector_calibration(finest),
            'local_source_checks_at_finest_grid': _local_checks(finest),
            'worst_branch_quantum_compiler_resources': protocol.source_cost(finest, STEPS, SETTINGS),
            'additional_native_repairs_per_history': 2,
            'classical_routing_control_overhead_costed': False,
            'interventions': {
                'interaction_deleted_histories': free,
                'feedback_deleted_histories': no_feedback,
                'fresh_state_reset_histories': reset,
                'interaction_deletion_tv': interaction_difference,
                'interaction_comparison_product_budget_exact': str(comparison_budget),
                'interaction_contrast_minus_product_budget_and_numerical_guard': interaction_difference-float(comparison_budget)-NUMERICAL_GUARD,
                'feedback_deletion_tv': _difference(best, no_feedback),
                'fresh_state_reset_tv': _difference(best, reset),
                'first_record_marginal_change_under_feedback_deletion': max(
                    abs(sum(row['probability'] for row in a if row['bits'][0] == bit)
                       -sum(row['probability'] for row in b if row['bits'][0] == bit))
                    for a, b in zip(best, no_feedback) for bit in (0, 1)),
                'floating_intervention_contrast_is_interval_certified': False,
            },
            'complete_native_history_law': records,
            'complete_native_history_count': len(records),
            'complete_probability_sum': math.fsum(row['probability'] for row in records),
            'native_feedback_port_zero_difference': math.fsum(
                a['probability']*a['final_port_zero']-b['probability']*b['final_port_zero']
                for a, b in zip(records, detached)),
            'previous_occupation_reference_tv': None if previous_target is None else _difference(target, previous_target),
            'complete_occupation_error_bound': None,
            'source_to_unrestricted_history_error_bound': None,
        })
        previous_target = target
    checks = {
        'complete_history_probability_conserved': all(abs(row['complete_probability_sum']-1) < NUMERICAL_GUARD for row in levels),
        'every_native_branch_retained': all(row['complete_native_history_count'] == 1056 for row in levels),
        'all_spatial_and_temporal_comparisons_within_bounds': all(
            grid['source_error_within_budget_plus_numerical_guard']
            and grid['spatial_error_within_budget_plus_numerical_guard']
            and grid['combined_error_within_budget_plus_numerical_guard']
            for row in levels for grid in row['spatial_refinements']),
        'temporal_refinements_within_bounds': all(x['within_budget_plus_numerical_guard']
            for row in levels for x in row['temporal_refinements_at_finest_grid']),
        'interaction_contrast_resolved_beyond_product_budget': all(
            row['interventions']['interaction_contrast_minus_product_budget_and_numerical_guard'] > 0
            for row in levels),
        'feedback_changes_observations': all(row['interventions']['feedback_deletion_tv'] > 1e-4 for row in levels),
        'reset_changes_observations': all(row['interventions']['fresh_state_reset_tv'] > .01 for row in levels),
        'later_feedback_preserves_earlier_record_law': all(
            row['interventions']['first_record_marginal_change_under_feedback_deletion'] < NUMERICAL_GUARD
            for row in levels),
    }
    result = {
        'schema': 'oph-fixed-band-scalar-observations-v1',
        'scope': 'Common smeared-field observations for a supplied interacting seven-mode action; complete finite-code spatial and temporal bounds, and a separate unrestricted fixed-band spatial theorem. No occupation-to-unrestricted certificate, band removal, Lorentz locality or empirical agreement.',
        'physical_meaning': {
            'torus_side': '2*pi', 'speed': '1', 'hbar': '1', 'retained_wavevectors': '0, +/-e_x, +/-e_y, +/-e_z',
            'bare_mass': '1', 'coupling_lambda_over_volume': '1', 'physical_lambda': '(2*pi)^3',
            'preparation': 'Product Omega=2 blank with exp(+/-i Q_cos_x/4) kick; current>0 chooses plus, zero and negative current choose minus; not an interacting vacuum',
            'measured_field': 'Q_cos_x = integral sqrt(2/V) cos(x) phi(x) d^3x',
            'instrument': 'K0=cos((Q_cos_x+pi/2)/2), K1=sin((Q_cos_x+pi/2)/2)',
            'feedback': 'First pointer bit one applies exp(i Q_cos_x/4); last bit selects native seam F rather than E',
            'finite_implementation': 'Every source preparation, instrument and feedback uses Q_d=P_d Q P_d in the declared occupation code, not the projected unrestricted operation. Their difference has no complete-history certificate here.',
            'intervals_in_hamiltonian_model_time': [str(t) for t in times],
            'calibrated_physical_seconds': False, 'experimental_data': None,
            'parameters_fitted_to_validation_records': [],
            'the_action_and_spatial_band_are_supplied': True,
        },
        'source': compiler.source_summary(),
        'levels': levels,
        'unrestricted_fixed_band_spatial_certificates': bounds.serializable([
            bounds.unrestricted_spatial_history(grid, times, eta=1,
                feedback_kick=F(1, 4), initial_kick=F(1, 4)) for grid in GRIDS]),
        'occupation_diagnostics': occupation_diagnostics(),
        'numerical_guard': NUMERICAL_GUARD,
        'floating_results_are_interval_certified': False,
        'physical_gate_error_certified': False,
        'checks_scope': 'Software and stated finite-model contracts; not physical completion or empirical validation',
        'checks': checks, 'all_pass': all(checks.values()),
    }
    if not result['all_pass']:
        raise ArithmeticError(f'fixed-band scalar observation checks failed: {checks}')
    return result
