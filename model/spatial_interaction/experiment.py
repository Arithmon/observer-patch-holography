"""One supplied spatial action, source operations, and complete observer laws."""
from fractions import Fraction as Q
import math

import numpy as np

from ..observer import native
from ..observer.source_scalar import compiler
from . import bounds, source
from .lattice import CubicBosonLattice

HOPPING = Q(1, 16)
INTERACTION = Q(1)
MASS = Q(1)
INTERVAL = Q(1, 16)
RADII = (1, 2, 3)
STEPS = (8, 16)


def apparatus(lattice):
    return (lattice.site_index[(0, 0, 0)], lattice.site_index[(1, 0, 0)])


def prepare(lattice, phase=0, source_execution=True):
    edge = apparatus(lattice)
    state = np.zeros(lattice.dimension, complex)
    state[lattice.index[(edge[0], edge[0])]] = 1
    if source_execution:
        state = lattice.apply_local(state, edge, source.preparation()[0])
        return lattice.apply_local(state, edge, source.feedback(math.pi*phase)[0])
    state = lattice.initial_state()
    state[lattice.index[tuple(sorted(edge))]] *= (-1)**phase
    return state


def evolve(lattice, state, model_time, steps):
    """Symmetric six-matching / onsite product of executed local source words."""
    if type(steps) is not int or steps <= 0:
        raise ValueError('positive integer evolution steps required')
    delta = float(model_time)/steps
    bond = source.bond(delta*lattice.hopping/2)[0]
    onsite = source.onsite(delta*lattice.interaction)[0]
    diagonal = np.diag(onsite)
    sites = len(lattice.sites)
    # This is the product over ALL sites, including vacuum local requests.
    split_factor = diagonal[0]**(sites-2)*diagonal[1]**2
    double_factor = diagonal[0]**(sites-1)*diagonal[2]
    phase = np.array([double_factor if i == j else split_factor for i, j in lattice.basis])
    out = state.copy()
    for _ in range(steps):
        for group in lattice.edge_groups:
            out = lattice.apply_matching(out, group, bond)
        out *= phase if out.ndim == 1 else phase[:, None]
        for group in reversed(lattice.edge_groups):
            out = lattice.apply_matching(out, group, bond)
    return out


def history(lattice, steps, phase=0, connected=True, source_execution=True):
    """Two retained binary reads with feedback; no reset or postselection."""
    state = prepare(lattice, phase, source_execution)
    propagate = (lambda v: evolve(lattice, v, INTERVAL, steps)) if source_execution else (
        lambda v: lattice.reference_evolve(v, float(INTERVAL)))
    edge = apparatus(lattice)
    if source_execution:
        kraus = source.instrument('ramsey')[0]
        read = lambda v: tuple(lattice.apply_local(v, edge, k) for k in kraus)
    else:
        # Independent algebraic projectors, without the source gate matrices.
        read = lattice.ramsey_branches
    state = propagate(state)
    out = []
    for first, branch in enumerate(read(state)):
        angle = math.pi/2 if connected and first else 0.
        if source_execution:
            branch = lattice.apply_local(branch, edge, source.feedback(angle)[0])
        else:
            branch = branch.copy()
            branch[lattice.index[tuple(sorted(edge))]] *= np.exp(1j*angle)
        branch = propagate(branch)
        for last, final in enumerate(read(branch)):
            out.append({'bits': [first, last], 'probability': float(np.vdot(final, final).real)})
    return out


def probabilities(rows):
    return np.array([row['probability'] for row in rows])


def total_variation(a, b):
    return float(np.sum(abs(probabilities(a)-probabilities(b)))/2)


def native_records(histories, connected=True):
    """Complete classical native attachment, including both fair repair coins.

    The first native current controls the local preparation phase. The last
    field bit controls E/F routing. Quantum coherence across native inputs is
    not claimed by this attachment.
    """
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
                    rows.append({
                        'initial': initial, 'records': [first_coin, first_current,
                            first, last, last_coin, last_current],
                        'final_port_zero': final[0],
                        'probability': quantum['probability']/264,
                    })
    return rows


def source_cost(lattice, steps, model_time=INTERVAL, intervals=2):
    """Worst branch pulse/service counts, including every idle-site request."""
    dt = float(model_time)/steps
    edge_count = sum(map(len, lattice.edge_groups))
    terms = [
        (2*edge_count*steps*intervals, source.bond(dt*lattice.hopping/2)[1]),
        (len(lattice.sites)*steps*intervals, source.onsite(dt*lattice.interaction)[1]),
        (1, source.preparation()[1]), (2, source.feedback(math.pi/2)[1]),
        (2, source.instrument('ramsey')[1]),
    ]
    return {key: sum(count*info['resources'][key] for count, info in terms)
            for key in source.COUNT_KEYS}


def witness(radius, steps=16):
    certificate = bounds.interaction_witness_bound(HOPPING, INTERACTION)
    # Kept explicit so the execution protocol is fixed independently of data.
    model_time = Q(1, 32)/bounds.generator_norm_bound(HOPPING, INTERACTION)
    predictions = []
    for interaction in (Q(0), INTERACTION):
        lattice = CubicBosonLattice(radius, float(HOPPING), float(interaction))
        state = evolve(lattice, prepare(lattice), model_time, steps)
        yes = lattice.apply_local(state, apparatus(lattice), source.instrument('ramsey')[0][1])
        predictions.append(float(np.vdot(yes, yes).real))
    error = sum(bounds.generator_norm_bound(HOPPING, u)**3*model_time**3/(3*steps**2)
                for u in (Q(0), INTERACTION))
    lower = 11*INTERACTION*model_time/30-error
    return {'radius': radius, 'model_time_exact': str(model_time), 'steps': steps,
            'free_probability': predictions[0], 'interacting_probability': predictions[1],
            'interaction_contrast': predictions[1]-predictions[0],
            'analytic_continuous_certificate': certificate,
            'ideal_source_product_error_exact': str(error),
            'ideal_source_product_contrast_lower_exact': str(lower),
            'ideal_source_product_contrast_lower_float': float(lower),
            'strictly_positive_certified_margin': lower > 0}


def run():
    levels = []
    previous = None
    for radius in RADII:
        lattice = CubicBosonLattice(radius, float(HOPPING), float(INTERACTION))
        reference = [history(lattice, 1, phase, source_execution=False) for phase in (0, 1)]
        approximations = []
        for steps in STEPS:
            actual = [history(lattice, steps, phase) for phase in (0, 1)]
            error = max(total_variation(a, b) for a, b in zip(actual, reference))
            budget = 2*bounds.generator_norm_bound(HOPPING, INTERACTION)**3*INTERVAL**3/(3*steps**2)
            approximations.append({'steps': steps, 'histories_by_native_phase': actual,
                'maximum_reference_history_tv': error,
                'history_error_certificate': bounds.history_bound(
                    [(INTERVAL, steps)]*2, HOPPING, INTERACTION),
                'within_ideal_bound_plus_numerical_guard': error <= float(budget)+2e-8,
                'worst_branch_source_resources': source_cost(lattice, steps)})
        best = approximations[-1]['histories_by_native_phase']
        disconnected = [history(lattice, STEPS[-1], phase, connected=False) for phase in (0, 1)]
        free = CubicBosonLattice(radius, float(HOPPING), 0.)
        free_histories = [history(free, STEPS[-1], phase) for phase in (0, 1)]
        records = native_records(best)
        detached = native_records(best, connected=False)
        row = {'radius': radius, 'side': 2*radius+1, 'spatial_dimensions': 3,
               'sites': len(lattice.sites), 'two_particle_dimension': lattice.dimension,
               'bonds': sum(map(len, lattice.edge_groups)),
               'declared_lattice_and_pointer_qubits': 2*len(lattice.sites)+1,
               'reference_histories_by_native_phase': reference,
               'evolution_refinements': approximations,
               'infinite_volume_history_certificate': bounds.boundary_bound(
                   radius, 2*INTERVAL, HOPPING, INTERACTION),
               'previous_box_reference_history_tv': None if previous is None else
                   max(total_variation(a, b) for a, b in zip(reference, previous)),
               'interaction_deletion_history_tv': max(total_variation(a, b) for a, b in zip(best, free_histories)),
               'field_feedback_deletion_history_tv': max(total_variation(a, b) for a, b in zip(best, disconnected)),
               'native_feedback_port_zero_difference': sum(a['probability']*a['final_port_zero']
                   -b['probability']*b['final_port_zero'] for a, b in zip(records, detached)),
               'complete_native_quantum_history_count': len(records),
               'complete_native_quantum_probability_sum': math.fsum(r['probability'] for r in records),
               'native_history_law': records,
               'interaction_witness': witness(radius)}
        combined = min(Q(1), Q(row['infinite_volume_history_certificate']['complete_history_tv_upper_exact'])
            + Q(approximations[-1]['history_error_certificate']['complete_history_tv_upper_exact']))
        row['finest_ideal_source_to_infinite_history_tv_upper_exact'] = str(combined)
        row['finest_ideal_source_to_infinite_history_tv_upper_float'] = float(combined)
        levels.append(row)
        previous = reference
    gate_info = {
        'preparation': source.preparation()[1], 'feedback': source.feedback(math.pi/2)[1],
        'ramsey_instrument': source.instrument('ramsey')[1],
        'bond_at_finest_step': source.bond(float(INTERVAL)*float(HOPPING)/(2*STEPS[-1]))[1],
        'onsite_at_finest_step': source.onsite(float(INTERVAL)*float(INTERACTION)/STEPS[-1])[1],
    }
    checks = {
        'every_evolution_refinement_within_bound': all(a['within_ideal_bound_plus_numerical_guard']
            for row in levels for a in row['evolution_refinements']),
        'complete_history_probability_conserved': all(abs(row['complete_native_quantum_probability_sum']-1)<2e-8 for row in levels),
        'interaction_witness_has_strict_analytic_margin': all(row['interaction_witness']['strictly_positive_certified_margin'] for row in levels),
        'executed_witness_matches_certificate': all(row['interaction_witness']['interaction_contrast'] >=
            row['interaction_witness']['ideal_source_product_contrast_lower_float']-2e-8 for row in levels),
        'interaction_changes_full_history': all(row['interaction_deletion_history_tv']>1e-5 for row in levels),
        'feedback_changes_full_history': all(row['field_feedback_deletion_history_tv']>1e-3 for row in levels),
    }
    result = {'schema': 'oph-spatial-interaction-v1',
        'scope': 'Supplied nonrelativistic Bose-Hubbard action on growing 3D boxes; fixed spacing; exact N=2 histories; ideal source-word comparison bounds; no relativistic continuum or physical clock identification.',
        'parameters_exact': {'hopping': str(HOPPING), 'interaction': str(INTERACTION),
            'mass': str(MASS), 'interval_model_time': str(INTERVAL)},
        'gap': bounds.gap_bound(MASS, HOPPING, INTERACTION),
        'source': compiler.source_summary(), 'local_source_checks': gate_info,
        'numerical_guard': 2e-8, 'floating_results_are_interval_certified': False,
        'levels': levels, 'checks': checks, 'all_pass': all(checks.values())}
    if not result['all_pass']:
        raise ArithmeticError(f'spatial observer checks failed: {checks}')
    return result
