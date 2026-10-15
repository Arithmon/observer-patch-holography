"""Canonical observer access on the same native, mixed and recorded process.

The report retains the original action, preparation, coherent source services,
noise and subsystem access assumptions. It consumes no experimental data.
Numerical identities attach existing finite mathematics; they are not physical
error certificates or a construction of a calibrated clock or spatial metric.
"""
from __future__ import annotations

from fractions import Fraction
import math

import numpy as np

from tools import strict_json
from ..common_source import open_environment, protocol
from . import algebra, common, runtime


DIMENSION = 2
STEPS = 8
TOLERANCE = 2e-9
FACTORS = (2, 2, 4)  # scalar1, scalar0, complete gauge/matter/link sector


def canonical_environment():
    return open_environment.Environment((
        open_environment.PhaseNoise(1, 0, Fraction(1, 20)),))


def reference_input():
    """Complex nonmaximal logical/reference entanglement; work is independent."""
    value = np.array([[1., 2j, 0.], [1-1j, 0., 3.]], complex)
    return value/np.linalg.norm(value)


def _residual(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ArithmeticError(f'{name} cannot be Boolean')
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ArithmeticError(f'{name} must be a finite nonnegative numerical diagnostic')
    return value


def same_history_report(history):
    """Attach two regional observers, public records and a forecast to this history."""
    history = runtime.validate_history(history)
    if history.work_dimension != math.prod(FACTORS) or history.reference_dimension != 1:
        raise ValueError('the canonical independently prepared two-scalar work state is required')
    if history.cursor != len(history.program.events):
        raise ValueError('the complete canonical four-read history is required')
    canonical = protocol.continuation_program(DIMENSION, STEPS)
    declaration = strict_json.loads(history.program.binding)
    if (declaration.get('canonical_program') != canonical.specification()
            or declaration.get('canonical_environment') != canonical_environment().specification()):
        raise ValueError('the declared canonical current-zero work history and environment are required')
    expected_program = common.compile_program(canonical, canonical_environment())
    if (history.program.digest != expected_program.digest
            or history.preparation_mode != 'declared'):
        raise ValueError('the executed canonical program and its declared preparation are required')
    complete = algebra.observer_view(history)
    regional = algebra.observer_view(history, FACTORS, (1,))
    if len(complete.records) != 16 or any(len(record) != 4 for record in complete.records):
        raise ValueError('the complete four-bit record menu is required')
    pauli_x = np.array([[0., 1.], [1., 0.]])
    pauli_z = np.diag([1., -1.])
    probe = pauli_x+.3*pauli_z
    embedded = algebra.embed_operator(probe, FACTORS, (1,))
    pairing = _residual(abs(algebra.expectation(regional, probe)
                            -algebra.expectation(complete, embedded)), 'regional pairing')
    mass_error = _residual(max(abs(complete.probabilities[r]-regional.probabilities[r])
                               for r in complete.records), 'record restriction')
    inclusion, local_gns, common_gns = algebra.gns_embedding(complete, FACTORS, (1,))
    isometry = _residual(np.linalg.norm(inclusion.conj().T@inclusion
                                       -np.eye(local_gns.dimension), 2), 'GNS isometry')
    lift_error = _residual(np.linalg.norm(inclusion@local_gns.lift(probe)
                                         -common_gns.lift(embedded)), 'GNS observable lift')
    public = {record: 1-2*record[-1] for record in complete.records}
    public_mean = float(algebra.public_expectation(complete, public).real)

    # Every regional block is retained; no favorable branch is selected to
    # make modular flow available. The implementation rejects nonfaithfulness.
    masses = tuple(regional.probabilities.values())
    if any(mass <= 0 for mass in masses):
        raise ArithmeticError('full record-aware modular access requires all blocks to be faithful')
    least_eigenvalue = min(float(np.linalg.eigvalsh(block/mass)[0])
                          for block, mass in zip(regional.blocks, masses, strict=True))
    parameter = .37
    flowed = algebra.modular_flow(regional, probe, parameter)
    modular_motion = _residual(max(np.linalg.norm(matrix-probe, 2) for matrix in flowed.values()),
                               'modular observable motion')
    modular_pairing = _residual(abs(algebra.expectation(regional, flowed)
                                    -algebra.expectation(regional, probe)), 'modular pairing')
    public_matrices = {record: value*np.eye(2) for record, value in public.items()}
    public_flow = algebra.modular_flow(regional, public_matrices, parameter)
    public_fixed = _residual(max(np.linalg.norm(public_flow[record]-matrix, 2)
                                 for record, matrix in public_matrices.items()), 'modular public record')

    # Forecast only the declared final gauge instrument. Its supported region
    # is checked against every complete source Kraus map before restriction.
    last = history.program.events[-1]
    if not last.visible or last.control is not None or len(last.instruments) != 1:
        raise ValueError('the final source read must be an uncontrolled complete instrument')
    instrument = last.instruments[0]
    prefix = runtime.execute(history.program, stop_after=history.cursor-1)
    restored = runtime.checkpoint_from_json(history.program, runtime.checkpoint_to_json(prefix))
    resumed = runtime.execute(history.program, checkpoint=restored)
    replay_error = _residual(runtime.complete_distance(history, resumed), 'checkpoint continuation')
    prefix_view = algebra.observer_view(prefix)
    full_forecast = algebra.forecast(prefix_view, instrument.kraus)
    regional_kraus = tuple(tuple(algebra.partial_trace(matrix, FACTORS, (2,))/4
                                  for matrix in group) for group in instrument.kraus)
    support_error = _residual(max(np.linalg.norm(matrix-algebra.embed_operator(local, FACTORS, (2,)), 2)
                                  for group, locals_ in zip(instrument.kraus, regional_kraus, strict=True)
                                  for matrix, local in zip(group, locals_, strict=True)), 'forecast instrument support')
    if support_error > TOLERANCE:
        raise ArithmeticError('regional forecast cannot discard factors used by its instrument')
    gauge_prefix = algebra.observer_view(prefix, FACTORS, (2,))
    regional_forecast = algebra.forecast(gauge_prefix, regional_kraus)
    forecast_error = _residual(max(abs(full_forecast[r]-complete.probabilities[r])
                                  for r in complete.records), 'complete source forecast')
    regional_forecast_error = _residual(max(abs(regional_forecast[r]-complete.probabilities[r])
                                           for r in complete.records), 'supported regional forecast')
    old_record_error = _residual(max(abs(mass-math.fsum(probability for r, probability
                                                      in complete.probabilities.items() if r[:-1] == record))
                                    for record, mass in prefix_view.probabilities.items()), 'prefix preservation')

    # The expectation onto the accessible region is CP, but cannot fill a
    # consensus tower's multiplicative algebra-refinement field.
    omitted_x = algebra.embed_operator(pauli_x, FACTORS, (0,))
    restricted_x = algebra.conditional_expectation(omitted_x, FACTORS, (1,))
    nonmultiplicative = _residual(np.linalg.norm(
        algebra.conditional_expectation(omitted_x@omitted_x, FACTORS, (1,))
        -restricted_x@restricted_x, 2), 'nonmultiplicativity control')
    checks = {
        'same_state_regional_pairing': pairing < TOLERANCE,
        'all_record_weights_preserved': mass_error < TOLERANCE,
        'regional_GNS_isometry': isometry < TOLERANCE,
        'regional_GNS_observable_lift': lift_error < TOLERANCE,
        'all_regional_blocks_faithful': math.isfinite(least_eigenvalue) and least_eigenvalue > 1e-12,
        'modular_access_moves_quantum_probe': modular_motion > 1e-4,
        'modular_state_pairing_preserved': modular_pairing < TOLERANCE,
        'modular_public_records_fixed': public_fixed < TOLERANCE,
        'checkpoint_replay_preserves_complete_output': replay_error < TOLERANCE,
        'full_work_forecast_matches_source_continuation': forecast_error < TOLERANCE,
        'region_supported_forecast_matches_source_continuation': regional_forecast_error < TOLERANCE,
        'old_record_marginals_preserved': old_record_error < TOLERANCE,
        'CP_restriction_is_not_claimed_multiplicative': nonmultiplicative > .99,
    }
    if not all(checks.values()) or not math.isfinite(public_mean):
        raise ArithmeticError('same-history observer attachment failed')
    return {
        'work_history_canonical_declaration_digest': canonical.digest,
        'visible_record_count': len(complete.records),
        'record_law': [{'record': list(record), 'probability': mass}
                       for record, mass in complete.probabilities.items()],
        'probability_mass': history.probability_mass,
        'final_public_gauge_bit_mean': public_mean,
        'regional_pairing_residual': pairing, 'regional_record_mass_residual': mass_error,
        'regional_GNS_dimension': local_gns.dimension, 'common_GNS_dimension': common_gns.dimension,
        'regional_GNS_isometry_residual': isometry, 'regional_GNS_lift_residual': lift_error,
        'smallest_normalized_regional_block_eigenvalue_diagnostic': least_eigenvalue,
        'modular_parameter_s': parameter,
        'maximum_modular_probe_motion': modular_motion,
        'modular_state_pairing_residual': modular_pairing,
        'modular_public_record_residual': public_fixed,
        'checkpoint_complete_output_distance': replay_error,
        'forecast_instrument_region_support_residual': support_error,
        'full_work_forecast_probability_residual': forecast_error,
        'regional_forecast_probability_residual': regional_forecast_error,
        'old_record_marginal_residual': old_record_error,
        'conditional_expectation_multiplicativity_failure': nonmultiplicative,
        'forecast_scope': 'complete state or the explicitly checked final gauge-supported instrument',
        'regional_autonomous_interacting_future_inferred': False,
        'modular_parameter_identified_with_physical_time': False,
        'supplied_factor_order': ['scalar_site1', 'scalar_site0', 'gauge_matter_link'],
        'numerical_rank_diagnostics_are_eigenvalue_certificates': False,
        'checks': checks, 'all_pass': all(checks.values()),
    }


def _native_reference_report(history, initial):
    from . import native
    reference = np.zeros((initial.shape[1], initial.shape[1]), complex)
    mass_terms, abort_terms = [], []
    records = []
    for branch, (record, amplitude, density) in zip(history.branches, native.output(history, initial), strict=True):
        if record != branch.record:
            raise ArithmeticError('native stream changed its complete record ordering')
        weight = float(np.trace(density).real)
        mass = weight*float(np.vdot(amplitude, amplitude).real)
        mass_terms.append(mass)
        abort_terms.append(mass if branch.abort else 0.)
        reference += weight*(amplitude.T@amplitude.conj())
        records.append(record)
    probability_mass = math.fsum(mass_terms)
    reference_error = _residual(np.linalg.norm(reference-initial.T@initial.conj(), 2), 'native spectator reference')
    if (not math.isfinite(probability_mass) or abs(probability_mass-1) > TOLERANCE
            or reference_error > TOLERANCE or len(records) != len(set(records))):
        raise ArithmeticError('native joint output lost a record or reference marginal')
    return {'reference_dimension': initial.shape[1], 'complete_raw_outcomes': len(records),
            'joint_probability_mass': probability_mass,
            'joint_abort_probability': math.fsum(abort_terms),
            'reference_marginal_residual': reference_error,
            'native_work_initially_independent': True,
            'reference_is_not_available_to_feedback': True}


def run():
    """Produce the fixed finite observer join, with independent compatibility oracles."""
    from . import native
    environment = canonical_environment()
    initial = reference_input()
    rows, histories, memory = [], {}, []
    for law in ('C', 'E'):
        for mode in ('protected_qubit', 'passive_qubit'):
            history = native.execute(environment, law, DIMENSION, STEPS, native_mode=mode)
            histories[law, mode] = history
            rows.append({'law': law, 'mode': mode,
                         'joint_reference_output': _native_reference_report(history, initial),
                         'native_channel_audit': native.audit(history),
                         'source_resources': native.resources(history),
                         'independent_original_wrapper_comparison': native.compatibility_report(history, initial=initial)})
        memory.append({'law': law, **native.memory_witness(histories[law, 'protected_qubit'],
                                                          histories[law, 'passive_qubit'])})
    # This is the SAME current-zero work history used by the active native
    # wrapper, with the same disturbance, records and source operations.
    work = histories['C', 'protected_qubit'].work_histories[1]
    access = same_history_report(work)
    canonical = protocol.continuation_program(DIMENSION, STEPS)
    compatibility = common.compatibility_report(canonical, environment)
    noiseless = common.execute(canonical)
    noise_distance = _residual(runtime.complete_distance(work, noiseless), 'unobserved environment effect')
    checks = {
        'same_history_observer_access': access['all_pass'],
        'complete_native_reference_outputs': all(abs(row['joint_reference_output']['joint_probability_mass']-1)
                                               < TOLERANCE for row in rows),
        'environment_has_nonzero_complete_output_effect': noise_distance > 1e-5,
        'both_laws_and_both_memory_controls_executed': len(rows) == 4,
        'independent_native_wrappers_agree': all(
            _residual(row['independent_original_wrapper_comparison']['complete_reference_output_trace_distance'],
                      'independent native reference comparison') < TOLERANCE
            and _residual(row['independent_original_wrapper_comparison']['maximum_work_block_trace_distance'],
                          'independent native work comparison') < TOLERANCE
            and row['independent_original_wrapper_comparison']['every_native_map_and_raw_record_agrees'] is True
            and row['independent_original_wrapper_comparison']['source_resources_agree'] is True
            for row in rows),
        'independent_work_interpreter_agrees': _residual(
            compatibility['complete_cq_trace_distance'], 'independent work comparison') < TOLERANCE,
    }
    if not all(checks.values()):
        raise ArithmeticError('primitive observer same-history join failed')
    return {
        'schema': 'oph.primitive_observer.observer_report.v1',
        'scope': 'finite native/source/mixed-state observer with compatible access on the same histories',
        'dimension': DIMENSION, 'product_steps_per_interval': STEPS,
        'environment': environment.specification(),
        'native_and_reference_reports': rows, 'active_passive_memory_comparisons': memory,
        'same_native_work_history_access': access,
        'original_work_interpreter_comparison': compatibility,
        'noise_deletion_complete_output_distance': noise_distance,
        'supplied': ['action and coupling', 'vacuum preparation and coherent code services',
                     'environment law, variance and independent draw schedule',
                     'observer access and tensor factors', 'dimensionless Hamiltonian model coordinate'],
        'empirical_data_consumed': False, 'physical_service_error_certified': False,
        'spatial_metric_or_physical_clock_recovered': False,
        'full_consensus_tower_or_OperationalObserver_instantiated': False,
        'checks': checks, 'all_pass': all(checks.values()),
    }
