"""Exact density closure and controlled propagation on the pinned federation.

Native attempts, Poissonized edge time and reconstructed unitary/modular
parameters are separate. The cycle control concerns conserved density
modes; it is not a no-go for other observer sectors or continuum theories.
"""
from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, product
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, expm_multiply

from ..observer import native
from . import source


SOURCE_PATHS = tuple(sorted(set((*source.SOURCE_PATHS,
    'model/federated_observer/scaling.py',
    'Lean/Thermodynamics/NativeRepairFourCycle.lean',
    'Lean/Thermodynamics/GraphDiffusion.lean',
    'Lean/ObserverPatchHolography/Locality/DependencyCone.lean',
    'Lean/ObserverPatchHolography/Locality/AdaptiveScheduler.lean',
    'Lean/SourceLaw/TransferMatrixLift.lean',
    'Lean/Dynamics/PrivateInner.lean'))))


def local_density_identity():
    """Exact basis proof of the edge identity for every binary configuration.

    Arbitrary real/complex coefficients follow by linearity. Summing edges
    proves G F_f=F_(L f) on every fixed-number sector, with no independence
    or product-state assumption.
    """
    rows = []
    for loads in product((0, 1), repeat=2):
        outputs = [native.native(loads, (0, 1), coin)[0] for coin in (0, 1)]
        for endpoint in (0, 1):
            actual = 2*loads[endpoint]-sum(out[endpoint] for out in outputs)
            expected = loads[endpoint]-loads[1-endpoint]
            if actual != expected or any(sum(out) != sum(loads) for out in outputs):
                raise ArithmeticError('native endpoint update lost the exact density closure')
            rows.append((loads, endpoint, actual))
    return tuple(rows)


def triangle_gram():
    """Exact8×8 Gram certificate for three binary endpoint averages."""
    states = tuple(product((0, 1), repeat=3))
    index = {x: i for i, x in enumerate(states)}
    b = np.zeros((8, 8), dtype=object)
    for a, d in combinations(range(3), 2):
        for i, state in enumerate(states):
            changed = list(state); changed[a], changed[d] = changed[d], changed[a]
            b[i, index[tuple(changed)]] += 1
    # A=(3I+B)/2, so A−3I/2=B/2=BᵀB/6 is positive.
    if not np.array_equal(b.T@b, 3*b):
        raise ArithmeticError('binary triangle positivity certificate failed')
    return tuple(tuple(int(value) for value in row) for row in b)


def transfer_floor(level):
    """Full binary configuration-transfer floor, uniform across number sectors."""
    graph = source.federation(level)
    triangle_gram()
    seams = set(graph.seams)
    template = tuple((a, b) for a, b in graph.seams if a < 12 and b < 12)
    faces = tuple((a, b, d) for a, b, d in combinations(range(12), 3)
                  if {(a, b), (a, d), (b, d)} <= set(template))
    counts = Counter(edge for face in faces for edge in combinations(face, 2))
    if len(template) != 30 or len(faces) != 20 or set(counts.values()) != {2}:
        raise ArithmeticError('source cell lacks the declared triangle cover')
    for cell in range(graph.cells):
        if any((12*cell+a, 12*cell+b) not in seams for a, b in template):
            raise ArithmeticError('a federation member changed its elementary source rule')
    lower = F(15*graph.cells, len(graph.seams))
    if not 0 < lower <= 1:
        raise ArithmeticError('transfer positivity floor is outside its domain')
    return lower


def log_remainder(level, order):
    if type(order) is not int or not 0 <= order <= 128:
        raise ValueError('integer logarithm order from zero through128 required')
    floor = transfer_floor(level)
    return (1-floor)**(order+1)/(F(order+1)*floor)


def laplacian(level):
    """Fresh sparse integer matrix from all declared source seams."""
    graph = source.federation(level)
    rows, cols, values = [], [], []
    for a, b in graph.seams:
        rows.extend((a, b, a, b)); cols.extend((a, b, b, a)); values.extend((1, 1, -1, -1))
    return sparse.coo_matrix((values, (rows, cols)), shape=(graph.ports, graph.ports), dtype=np.int64).tocsr()


def spectral_diagnostic(level):
    """Low density-mode poles, not certified continuum eigenvalues."""
    graph = source.federation(level)
    operator = laplacian(level).astype(float)
    start = np.linspace(1, 2, graph.ports)
    values, vectors = eigsh(operator, k=9, which='SM', v0=start, tol=2e-12, maxiter=30*graph.ports)
    order = np.argsort(values); values, vectors = values[order], vectors[:, order]
    residual = float(np.linalg.norm(operator@vectors-vectors*values, ord='fro'))
    if (not np.isfinite(values).all() or not math.isfinite(residual) or residual > 2e-9
            or abs(values[0]) > 2e-10 or values[1] <= 0):
        raise ArithmeticError('low density spectrum failed its numerical diagnostic')
    positive = values[1:]
    energies = -np.log1p(-positive/(2*len(graph.seams)))
    return {'level': level, 'ports': graph.ports, 'cells': graph.cells,
        'graph_laplacian_positive_eigenvalues_diagnostic': positive.tolist(),
        'native_attempt_transfer_energies_diagnostic': energies.tolist(),
        'per_edge_sweep_rescaled_energies_diagnostic': (len(graph.seams)*energies).tolist(),
        'eigenpair_residual_diagnostic': residual,
        'exact_density_relation': 'K F_f=(1-lambda/(2|E|)) F_f; H F_f=-log(1-lambda/(2|E|)) F_f',
        'full_binary_transfer_lower_bound': str(transfer_floor(level)),
        'order16_full_transfer_log_operator_remainder': str(log_remainder(level, 16)),
        'continuum_dispersion_certified': False}


def cycle_bounds(size):
    """Rational bounds on E(2p)/E(p); common positive clock scale cancels.

    For a cycle, x_k=(1-cos(2pi k/N))/N and E_k=-log(1-x_k).
    x <= -log(1-x) <= x/(1-x), pi²<10 and sin(x)<=x give
    the exported bounds. They converge to4, excluding the massless linear
    comparison2 for every declared size.
    """
    if type(size) is not int or not 8 <= size <= 65536:
        raise ValueError('integer cycle size from8 through65536 required')
    lower = 4*(1-F(10, size**2))*(1-F(20, size**3))
    upper = 4/(1-F(80, size**3))
    return lower, upper


def cycle_diagnostic(size):
    lower, upper = cycle_bounds(size)
    energies = [-math.log1p(-(2*math.sin(math.pi*k/size)**2)/size) for k in (1, 2)]
    ratio = energies[1]/energies[0]
    if not float(lower) <= ratio <= float(upper) or lower <= 2:
        raise ArithmeticError('exact cycle dispersion interval failed its control')
    return {'control_cycle_ports': size, 'ratio_lower': str(lower), 'ratio_upper': str(upper),
        'energy_ratio_diagnostic': ratio, 'gapless_linear_comparison': '2',
        'ratio_limit': '4', 'clock_rescaling_changes_ratio': False,
        'this_cycle_is_the_declared_federation': False}


def locality_bound(max_degree, time, radius, reads=1):
    """Exact upper bound for occupation reads under edge Poissonization.

    Each edge attempts at rateone and swaps with probability1/2. A backward
    occupation label jumps at rate<=max_degree/2. To cross radiusR needs
    at leastR jumps. The Poisson factorial-moment bound and a union bound
    give reads*(max_degree*time/2)^R/R!, capped byone.
    Visible native coins are summed out; their joint record law is not the
    state-independent swap-coin law used by this graphical coupling.
    """
    if (type(max_degree) is not int or max_degree < 0 or type(time) is not F or time < 0
            or type(radius) is not int or not 1 <= radius <= 128
            or type(reads) is not int or not 1 <= reads <= 1000000):
        raise ValueError('nonnegative degree/rational time, positive bounded radius/read count required')
    return min(F(1), reads*(max_degree*time/2)**radius/math.factorial(radius))


def distances(level, origin=0):
    graph = source.federation(level)
    if type(origin) is not int or not 0 <= origin < graph.ports:
        raise ValueError('declared federation origin required')
    result = [-1]*graph.ports; result[origin] = 0
    todo = deque([origin])
    while todo:
        a = todo.popleft()
        for b in graph.neighbours[a]:
            if result[b] == -1:
                result[b] = result[a]+1; todo.append(b)
    if any(d < 0 for d in result):
        raise ArithmeticError('declared federation is disconnected')
    return tuple(result)


def propagation_diagnostic(level):
    graph = source.federation(level)
    time = F(1, 4)
    initial = np.zeros(graph.ports); initial[0] = 1
    output = expm_multiply(-float(time)*laplacian(level)/2, initial)
    if not np.isfinite(output).all() or output.min() < -2e-14 or abs(output.sum()-1) > 2e-13:
        raise ArithmeticError('Poissonized density diagnostic lost probability')
    degree = max(map(len, graph.neighbours)); distance = np.asarray(distances(level))
    rows = []
    for radius in (2, 4, 8):
        bound = locality_bound(degree, time, radius)
        outside = float(output[distance >= radius].sum())
        if not -2e-14 <= outside <= float(bound)+2e-13:
            raise ArithmeticError('density propagation exceeded its exact tail allowance')
        rows.append({'radius': radius, 'probability_at_distance_at_least_radius_diagnostic': outside,
                     'uniform_exact_probability_bound': str(bound)})
    return {'level': level, 'maximum_port_degree': degree,
        'declared_per_edge_attempt_Poisson_parameter': str(time), 'origin': 0,
        'density_semigroup': 'exp(-time L/2)', 'tail_checks': rows,
        'fixed_word_dependency_ball_sizes': {str(n): int(np.count_nonzero(distance <= n)) for n in (0, 1, 2, 4)},
        'fixed_word_cone_parameter': 'number of committed endpoint attempts; no duration assigned',
        'native_CP_coherent_propagation_claimed': False,
        'visible_native_coin_and_current_history_bound_claimed': False,
        'physical_light_cone_claimed': False}


def run():
    local = local_density_identity()
    spectra = [spectral_diagnostic(level) for level in (0, 1, 2)]
    ratios = []
    for coarse, fine in zip(spectra, spectra[1:]):
        a = coarse['per_edge_sweep_rescaled_energies_diagnostic'][0]
        b = fine['per_edge_sweep_rescaled_energies_diagnostic'][0]
        ratios.append({'levels': [coarse['level'], fine['level']],
                      'lowest_density_energy_ratio_diagnostic': a/b,
                      'twofold_screen_resolution_exponent_diagnostic': math.log(a/b)/math.log(2),
                      'scaling_limit_inferred_from_three_levels': False})
    return {'exact_local_density_basis_cases': len(local),
        'arbitrary_binary_fixed_number_sector_closure': True,
        'triangle_gram_dimension': len(triangle_gram()),
        'spectra': spectra, 'pinned_family_gap_ratios': ratios,
        'cycle_controls': [cycle_diagnostic(n) for n in (8, 16, 32, 64, 128)],
        'propagation': [propagation_diagnostic(level) for level in (0, 1, 2)],
        'scope': {'spectra_and_propagation_are_numerical_diagnostics': True,
            'local_density_identity_triangle_positivity_and_tail_bounds_are_exact': True,
            'native_attempts_equal_modular_or_physical_time': False,
            'native_CP_channel_equals_transfer_unitary': False,
            'Poisson_bound_preserves_joint_visible_native_coin_records': False,
            'pinned_refinement_local_Cauchy_or_Lorentz_limit_proved': False}}
