"""Complete source-executed records of a supplied canonical scalar regulator."""
from fractions import Fraction as F
import math

import numpy as np

from ..observer import native
from ..observer.source_scalar import compiler
from . import bounds, exact, exact_history, source
from .lattice import ScalarLattice, ScalarParameters, free_dispersion

INTERVAL = F(1, 128)
WITNESS_TIME = F(1, 4096)
LEVELS = ((8,)+(2,)*7, (8,)+(4,)*7)
STEPS = (8, 32)


def projectors(d=8, sign=-1):
    """Independent local algebraic instrument, without compiled gate matrices."""
    if d < 5 or sign not in (-1, 1):
        raise ValueError('five or more levels and sign +/-1 required')
    e = np.zeros(d, complex)
    e[0], e[4] = 1/5, sign*1j*math.sqrt(24)/5
    yes = np.outer(e, e.conj())
    return np.eye(d)-yes, yes


def evolve(lattice, state, model_time, steps):
    """Symmetric onsite/bond steps; every local matrix is source-executed."""
    if type(steps) is not int or steps < 1:
        raise ValueError('positive integer symmetric steps required')
    dt = float(model_time)/steps
    onsite = [source.onsite(d, int(lattice.degree[i]), dt/2, lattice.parameters)[0]
              for i, d in enumerate(lattice.dimensions)]
    bonds = [(edge, source.bond(lattice.dimensions[edge[0]], lattice.dimensions[edge[1]],
                                dt, lattice.parameters)[0]) for edge in lattice.edges]
    out = np.asarray(state, dtype=complex).copy()
    for _ in range(steps):
        for i, matrix in enumerate(onsite):
            out = lattice.apply_local(out, [i], matrix)
        # All bond generators are products of the same local q matrices and
        # commute, including bonds sharing a vertex. Their order is fixed.
        for edge, matrix in bonds:
            out = lattice.apply_local(out, edge, matrix)
        for i, matrix in reversed(list(enumerate(onsite))):
            out = lattice.apply_local(out, [i], matrix)
    return out


def history(lattice, steps=16, sign=-1, connected=True, source_execution=True):
    """Retain both reads and their disturbed states; no reset or postselection."""
    propagate = ((lambda v: evolve(lattice, v, INTERVAL, steps)) if source_execution
                 else (lambda v: lattice.evolve(v, float(INTERVAL))))
    state = propagate(lattice.vacuum_state())
    kraus = (source.instrument(lattice.dimensions[0], sign)[0] if source_execution
             else projectors(lattice.dimensions[0], sign))
    branches = np.column_stack([lattice.apply_local(state, [0], k) for k in kraus])
    if connected:
        if source_execution:
            feedback = source.feedback(lattice.dimensions[0], math.pi/2)[0]
        else:
            feedback = np.eye(lattice.dimensions[0], dtype=complex)
            feedback[4, 4] = 1j
        branches[:, 1] = lattice.apply_local(branches[:, 1], [0], feedback)
    branches = propagate(branches)
    rows = []
    for first in (0, 1):
        for last, k in enumerate(kraus):
            final = lattice.apply_local(branches[:, first], [0], k)
            rows.append({'bits': [first, last], 'probability': float(np.vdot(final, final).real)})
    return rows


def tv(a, b):
    if [row['bits'] for row in a] != [row['bits'] for row in b]:
        raise ValueError('record alphabets differ')
    return math.fsum(abs(x['probability']-y['probability']) for x, y in zip(a, b))/2


def native_records(histories, connected=True):
    """Classical 66-state attachment; first current selects pointer phase."""
    rows = []
    for initial, loads in enumerate(native.states(binary=True)):
        for first_coin in (0, 1):
            middle, (_, first_current) = native.native(loads, native.E, first_coin)
            phase = int(first_current > 0)
            for quantum in histories[phase]:
                first, last = quantum['bits']
                seam = native.F if connected and last else native.E
                for last_coin in (0, 1):
                    final, (_, last_current) = native.native(middle, seam, last_coin)
                    rows.append({'initial': initial,
                        'records': [first_coin, first_current, first, last, last_coin, last_current],
                        'final_port_zero': final[0], 'probability': quantum['probability']/264})
    return rows


def source_cost(lattice, steps, model_time=INTERVAL, intervals=2, reads=2):
    dt = float(model_time)/steps
    requests = []
    for i, d in enumerate(lattice.dimensions):
        requests.append((2*steps*intervals,
            source.onsite(d, int(lattice.degree[i]), dt/2, lattice.parameters)[1]))
    for i, j in lattice.edges:
        requests.append((steps*intervals,
            source.bond(lattice.dimensions[i], lattice.dimensions[j], dt, lattice.parameters)[1]))
    requests += [(reads, source.instrument(lattice.dimensions[0], -1)[1]),
                 (int(reads > 1), source.feedback(lattice.dimensions[0], math.pi/2)[1])]
    return {key: sum(count*info['resources'][key] for count, info in requests)
            for key in source.COUNT_KEYS}


def local_checks(lattice, steps):
    dt = float(INTERVAL)/steps
    return {
        'onsite': {str(d): source.onsite(d, 3, dt/2, lattice.parameters)[1]
                   for d in sorted(set(lattice.dimensions))},
        'bond': {f'{dl}x{dr}': source.bond(dl, dr, dt, lattice.parameters)[1]
                 for dl, dr in sorted({(lattice.dimensions[i], lattice.dimensions[j])
                                       for i, j in lattice.edges})},
        'pointer': {str(sign): source.instrument(lattice.dimensions[0], sign)[1]
                    for sign in (-1, 1)},
        'feedback': source.feedback(lattice.dimensions[0], math.pi/2)[1],
    }


def witness(dimensions, steps=16):
    predictions = []
    for coupling in (0, 1):
        lattice = ScalarLattice(dimensions, ScalarParameters(coupling=coupling))
        state = evolve(lattice, lattice.vacuum_state(), WITNESS_TIME, steps)
        yes = lattice.apply_local(state, [0], source.instrument(dimensions[0], -1)[0][1])
        predictions.append(float(np.vdot(yes, yes).real))
    finite = exact.interaction_witness_certificate(WITNESS_TIME, dimensions)
    error = sum((bounds.product_error(dimensions, WITNESS_TIME, steps, coupling)
                 for coupling in (0, 1)), F(0))
    lower = finite['contrast_lower']-error
    return {'model_time_exact': str(WITNESS_TIME), 'steps': steps,
            'free_probability': predictions[0], 'interacting_probability': predictions[1],
            'contrast': predictions[1]-predictions[0],
            'finite_continuous_certificate': bounds.serializable(finite),
            'ideal_product_probability_error_exact': str(error),
            'ideal_executed_contrast_lower_exact': str(lower),
            'ideal_executed_contrast_lower_float': float(lower),
            'strictly_positive_certified_margin': lower > 0}


def attach_certificates(result):
    """Bind each executed finite history to the same eight infinite oscillators.

    Exact polynomial coefficients and outward rational root bounds establish
    the certificates. Floating reference probabilities only check the binding
    to the numerical executor; they are not inputs to the proofs.
    """
    intervals_contain_reference = True
    for row in result['levels']:
        dimensions = tuple(row['local_dimensions'])
        comparisons = []
        for phase, sign in enumerate((-1, 1)):
            comparison = exact_history.comparison_certificate(
                dimensions, time=INTERVAL, degree=2, instrument_sign=sign)
            finite = comparison['finite']
            for certificate, observation in zip(finite['branches'], row['reference_histories_by_native_phase'][phase]):
                low, high = certificate['probability_interval']
                intervals_contain_reference &= (certificate['records'] == observation['bits']
                    and float(low)-2e-12 <= observation['probability'] <= float(high)+2e-12)
            comparisons.append(comparison)
        cutoff = max(item['complete_history_tv_upper'] for item in comparisons)
        product = F(row['evolution_refinements'][-1]['ideal_product_history_tv_upper_exact'])
        combined = min(F(1), cutoff+product)
        row['complete_occupation_comparisons_by_native_phase'] = bounds.serializable(comparisons)
        row['complete_occupation_history_tv_upper_exact'] = str(cutoff)
        row['ideal_source_to_untruncated_history_tv_upper_exact'] = str(combined)
        row['ideal_source_to_untruncated_history_tv_upper_float'] = float(combined)
    result['checks']['exact_finite_intervals_contain_independent_reference'] = bool(intervals_contain_reference)
    result['checks']['complete_untruncated_history_bound_nonvacuous'] = all(
        row['ideal_source_to_untruncated_history_tv_upper_float'] < 0.01 for row in result['levels'])
    result['all_pass'] = all(result['checks'].values())
    if not result['all_pass']:
        raise ArithmeticError(f'canonical scalar certificate checks failed: {result["checks"]}')
    return result


def run():
    levels = []
    previous = None
    for dimensions in LEVELS:
        lattice = ScalarLattice(dimensions)
        reference = [history(lattice, sign=sign, source_execution=False) for sign in (-1, 1)]
        refinements = []
        for steps in STEPS:
            actual = [history(lattice, steps, sign) for sign in (-1, 1)]
            budget = min(F(1), 2*bounds.product_error(dimensions, INTERVAL, steps))
            error = max(tv(a, b) for a, b in zip(actual, reference))
            refinements.append({'steps': steps, 'histories_by_native_phase': actual,
                'source_reference_history_tv': error, 'ideal_product_history_tv_upper_exact': str(budget),
                'within_bound_plus_numerical_guard': error <= float(budget)+2e-8,
                'worst_branch_source_resources': source_cost(lattice, steps)})
        best = refinements[-1]['histories_by_native_phase']
        disconnected = [history(lattice, STEPS[-1], sign, connected=False) for sign in (-1, 1)]
        free_lattice = ScalarLattice(dimensions, ScalarParameters(coupling=0))
        free = [history(free_lattice, STEPS[-1], sign) for sign in (-1, 1)]
        records = native_records(best)
        detached = native_records(best, connected=False)
        levels.append({'spatial_dimensions': 3, 'sites': 8, 'bonds': 12,
            'local_dimensions': list(dimensions), 'tensor_dimension': lattice.dimension,
            'lattice_and_pointer_qubits': sum(d.bit_length()-1 for d in dimensions)+1,
            'reference_histories_by_native_phase': reference, 'evolution_refinements': refinements,
            'local_source_checks_at_finest_step': local_checks(lattice, STEPS[-1]),
            'interaction_deletion_history_tv': max(tv(a, b) for a, b in zip(best, free)),
            'feedback_deletion_history_tv': max(tv(a, b) for a, b in zip(best, disconnected)),
            'previous_occupation_reference_tv': None if previous is None else max(
                tv(a, b) for a, b in zip(reference, previous)),
            'native_feedback_port_zero_difference': math.fsum(a['probability']*a['final_port_zero']
                -b['probability']*b['final_port_zero'] for a, b in zip(records, detached)),
            'complete_native_quantum_history_count': len(records),
            'complete_probability_sum': math.fsum(row['probability'] for row in records),
            'native_history_law': records,
            'interaction_witness': witness(dimensions)})
        previous = reference
    windows = [bounds.serializable(bounds.free_window(a, F(1), F(1)))
               for a in (F(1), F(1, 2), F(1, 4), F(1, 8))]
    checks = {
        'complete_history_probability_conserved': all(abs(row['complete_probability_sum']-1) < 2e-8 for row in levels),
        'independent_evolution_within_bound': all(x['within_bound_plus_numerical_guard']
            for row in levels for x in row['evolution_refinements']),
        'interaction_witness_has_analytic_margin': all(row['interaction_witness']['strictly_positive_certified_margin'] for row in levels),
        'executed_witness_matches_margin': all(row['interaction_witness']['contrast'] >=
            row['interaction_witness']['ideal_executed_contrast_lower_float']-2e-8 for row in levels),
        'feedback_changes_records': all(row['feedback_deletion_history_tv'] > 1e-3 for row in levels),
    }
    result = {'schema': 'oph-canonical-scalar-regulator-v1',
        'scope': 'Source-executed eight-site canonical phi4 regulator; full-oscillator finite-volume gap and bounded interaction witness; free Klein--Gordon spatial comparison only. No interacting continuum or physical clock identification.',
        'parameters_exact': {'mass': '1', 'speed': '1', 'spacing': '1', 'coupling': '1',
            'oscillator_basis_frequency': '2', 'hbar': '1', 'interval_model_time': str(INTERVAL)},
        'preparation': 'Supplied product Omega=2 oscillator blanks; not the interacting vacuum',
        'source': compiler.source_summary(), 'gap': bounds.gap_statement(),
        'free_spatial_windows': windows,
        'untruncated_interaction_witness': bounds.serializable(exact.interaction_witness_certificate(WITNESS_TIME)),
        'levels': levels, 'numerical_guard': 2e-8, 'floating_results_are_interval_certified': False,
        'checks': checks, 'all_pass': all(checks.values())}
    if not result['all_pass']:
        raise ArithmeticError(f'canonical scalar checks failed: {checks}')
    return attach_certificates(result)
