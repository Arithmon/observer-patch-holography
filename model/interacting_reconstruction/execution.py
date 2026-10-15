"""Hidden supplied field backend for withheld fully coupled executions.

Only this backend knows the permutation, action parameters and physical cube.
The predictor receives probabilities and protocol declarations instead. The
isolated calibration controls and full dynamics are two supplied schedules
of the same local action, not autonomously generated observer apparatus.
"""
from __future__ import annotations

from numbers import Real

import numpy as np

from ..relativistic_scalar import source as local_source
from ..relativistic_scalar.experiment import evolve
from ..relativistic_scalar.lattice import ScalarLattice, ScalarParameters, local_operators
from .source import PERMUTATION


def coupled_history(dimensions, duration, *, read_site=0, sign=-1, steps=32,
                    feedback_phase=np.pi/2, parameters=None):
    """Run complete source words on the eight-site field and both pointer reads."""
    if (not isinstance(dimensions, (tuple, list)) or len(dimensions) != 8
            or any(type(d) is not int or d not in (2, 4, 8) for d in dimensions)
            or type(read_site) is not int or not 0 <= read_site < 8
            or dimensions[read_site] != 8):
        raise ValueError('eight opaque sites and an eight-level read site required')
    if (isinstance(duration, (bool, np.bool_)) or not isinstance(duration, Real)
            or not np.isfinite(float(duration)) or duration <= 0
            or type(steps) is not int or steps < 1):
        raise ValueError('positive model-coordinate duration and integer step count required')
    physical_dimensions = [0]*8
    for label, physical in enumerate(PERMUTATION):
        physical_dimensions[physical] = dimensions[label]
    lattice = ScalarLattice(physical_dimensions, parameters)
    site = PERMUTATION[read_site]
    pointer, pointer_info = local_source.instrument(8, sign)
    feedback, feedback_info = local_source.feedback(8, feedback_phase)
    state = evolve(lattice, lattice.vacuum_state(), duration, steps)
    branches = np.column_stack([lattice.apply_local(state, [site], k) for k in pointer])
    branches[:, 1] = lattice.apply_local(branches[:, 1], [site], feedback)
    branches = evolve(lattice, branches, duration, steps)
    rows = []
    for first in (0, 1):
        for last, k in enumerate(pointer):
            final = lattice.apply_local(branches[:, first], [site], k)
            rows.append({'bits': [first, last], 'probability': float(np.vdot(final, final).real)})
    dt = float(duration)/steps
    requests = [(2, pointer_info), (1, feedback_info)]
    for i, d in enumerate(lattice.dimensions):
        requests.append((4*steps, local_source.onsite(d, 3, dt/2, lattice.parameters)[1]))
    for i, j in lattice.edges:
        requests.append((2*steps, local_source.bond(lattice.dimensions[i], lattice.dimensions[j],
                                                   dt, lattice.parameters)[1]))
    return rows, {
        'worst_branch_resources': {key: sum(n*meta['resources'][key] for n, meta in requests)
                                   for key in local_source.COUNT_KEYS},
        'maximum_local_target_residual': max(meta['full_phase_target_residual'] for _, meta in requests),
        'steps_per_interval': steps, 'intervals': 2, 'reads': 2,
        'all_outcomes_retained': True, 'reset_between_reads': False,
        'physical_gate_error_certified': False,
    }


def target_audit():
    """Open only after prediction: independent hidden graph/parameter diagnostics."""
    parameters = ScalarParameters()
    lattice = ScalarLattice((2,)*8, parameters)
    permutation = np.asarray(PERMUTATION)
    adjacency = (lattice.stiffness < 0).astype(int)[np.ix_(permutation, permutation)]
    op = local_operators(4, parameters.omega)
    onsite = op['p2']/2+(parameters.mass**2+3*parameters.kappa)*op['q2']/2+parameters.g*op['q4']
    pair = -parameters.kappa*np.kron(lattice.operators[0]['q'][:2, :2],
                                     lattice.operators[0]['q'][:2, :2])
    return {'mass_squared': parameters.mass**2, 'quartic_g': parameters.g,
            'bond_kappa': parameters.kappa, 'adjacency': adjacency,
            'onsite_norm': float(np.linalg.norm(onsite, 2)),
            'pair_norm': float(np.linalg.norm(pair, 2))}
