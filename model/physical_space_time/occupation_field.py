"""One-history native occupation fields and quantitative Green-limit controls.

The observation-count CLT is distinct from spatial refinement. No coherent
Hamiltonian, independent replica or physical time unit enters either proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import math
from pathlib import Path
from types import MappingProxyType

import numpy as np
import sympy as sp

from . import construction as passage
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*passage.SOURCE_PATHS,
    'model/physical_space_time/occupation_field.py',
    'model/observer/scalar_response/experiment.py', 'tools/strict_json.py',
    'model/__init__.py', 'model/common_source/__init__.py',
    'model/common_source/experiment.py', 'model/common_source/gauge_matter.py',
    'model/common_source/native_join.py', 'model/common_source/records.py',
    'model/common_source/scalar_coupling.py', 'model/observer/__init__.py',
    'model/observer/p_response/__init__.py', 'model/observer/p_response/experiment.py',
    'model/observer/scalar_response/__init__.py',
    'model/observer/source_scalar/__init__.py', 'model/observer/source_scalar/bounds.py',
    'model/observer/source_scalar/compiler.py', 'model/observer/source_scalar/cutoff.py',
    'model/observer/source_scalar/experiment.py', 'model/observer/source_scalar/field.py',
    'model/observer_reconstruction/__init__.py', 'model/physical_space_time/__init__.py',
    'model/relative_observer_geometry/__init__.py'))))
CHI = sp.Rational(5, 33)
P0 = sp.ImmutableMatrix(sp.eye(12)-sp.ones(12)/12)
MODES = ((1, 0), (0, 1), (1, 1), (2, 0))


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original bounded integer is required')
    return value


def _rows(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def _laplacian(size, edges):
    value = sp.zeros(size)
    for a, b in edges:
        value[a, a] += 1; value[b, b] += 1
        value[a, b] -= 1; value[b, a] -= 1
    return sp.ImmutableMatrix(value)


@lru_cache(maxsize=1)
def native_counts():
    """Execute all 3960 native branches and decode every field by currents."""
    states = tuple(native.states(binary=True))
    index = {state: i for i, state in enumerate(states)}
    seams, laplacian = passage.recovered_carrier()
    for bound in (laplacian-2*P0, 10*P0-laplacian):
        _, certificate = bound[:11, :11].LDLdecomposition(hermitian=False)
        if any(certificate[i, i] <= 0 for i in range(11)):
            raise ArithmeticError('the exact source spectral bounds changed')
    a, b = sp.zeros(66), sp.zeros(66)
    checked = 0
    for column, before in enumerate(states):
        for edge in seams:
            for coin in (0, 1):
                after, (_, current) = native.native(before, edge, coin)
                decoded = tuple(before[i]+current*(int(i == edge[1])-int(i == edge[0]))
                                for i in range(12))
                if decoded != after:
                    raise ArithmeticError('native current disclosure lost occupation readback')
                (b if 0 in edge else a)[index[after], column] += 1
                checked += 1
    z = sp.Matrix(12, 66, lambda i, j: sp.Rational(states[j][i])-sp.Rational(1, 6))
    if z*z.T/66 != CHI*P0:
        raise ArithmeticError('canonical fixed-total covariance changed')
    if a != a.T or b != b.T or sp.ones(1, 66)*a != 50*sp.ones(1, 66):
        raise ArithmeticError('remote native counts lost symmetry or exposure')
    if sp.ones(1, 66)*b != 10*sp.ones(1, 66):
        raise ArithmeticError('incident exposure changed')
    if z*(a+b)/60 != (sp.eye(12)-laplacian/60)*z:
        raise ArithmeticError('native density closure failed')
    return states, sp.ImmutableMatrix(a), sp.ImmutableMatrix(b), sp.ImmutableMatrix(z), checked


@lru_cache(maxsize=2, typed=True)
def field(clock='attempt'):
    """Exact 12-dimensional density map and CLT covariance in the named clock."""
    if type(clock) is not str or clock not in ('attempt', 'incident'):
        raise ValueError('attempt or incident clock required')
    seams, laplacian = passage.recovered_carrier()
    green = (laplacian+sp.ones(12)/12).inv()-sp.ones(12)/12
    incident = seams if clock == 'attempt' else tuple(e for e in seams if 0 in e)
    m = len(incident)
    remote = _laplacian(12, tuple(e for e in seams if e not in incident))
    d = (2*m*sp.eye(12)+remote)/60
    b = sp.Rational(m, 30)*sp.eye(12)-(laplacian-remote)/60
    k = b*d.inv()
    covariance = CHI*(4*m*green+remote*green+green*remote-P0)
    inverse = (sp.eye(12)-k+sp.ones(12)/12).inv()-sp.ones(12)/12
    if covariance != CHI*(inverse+inverse.T-P0):
        raise ArithmeticError('stopped density Green sum disagrees with the clock correction')
    if inverse*(CHI*(P0-k*k.T))*inverse.T != covariance:
        raise ArithmeticError('the martingale Poisson covariance differs from the Green sum')
    # Principal minors certify positivity on the constant-free space exactly.
    _, pivots = covariance[:11, :11].LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(11)):
        raise ArithmeticError('limiting field covariance is not strictly positive')
    precision = (covariance+sp.ones(12)/12).inv()-sp.ones(12)/12
    if covariance*precision != P0:
        raise ArithmeticError('the field action is not the inverse covariance')
    return MappingProxyType({'clock': clock, 'm': m, 'laplacian': laplacian,
            'remote_laplacian': remote, 'green': sp.ImmutableMatrix(green),
            'kernel': sp.ImmutableMatrix(k), 'covariance': sp.ImmutableMatrix(covariance),
            'precision': sp.ImmutableMatrix(precision),
            'positive_LDL_pivots': tuple(pivots[i, i] for i in range(11))})


@lru_cache(maxsize=2, typed=True)
def configuration_kernel(clock='attempt'):
    """Independent complete 66-state kernel, retaining all native alternatives."""
    _, a, b, z, _ = native_counts()
    if clock == 'attempt':
        k = (a+b)/60
    elif clock == 'incident':
        k = b*(60*sp.eye(66)-a).inv(method='DM')
    else:
        raise ValueError('attempt or incident clock required')
    if z*k != field(clock)['kernel']*z or sp.ones(1, 66)*k != sp.ones(1, 66):
        raise ArithmeticError('the complete state law and reduced field law differ')
    if k*sp.ones(66, 1) != sp.ones(66, 1) or any(x < 0 for x in k):
        raise ArithmeticError('the uniform stopped state is not stationary')
    return sp.ImmutableMatrix(k)


def covariance_at_count(clock, count):
    """Covariance of count**(-1/2) sum_{j=0}^{count-1} eta_j; no CLT used."""
    _integer(count, 1, 128)
    k = field(clock)['kernel']
    answer, power = count*P0, P0
    for lag in range(1, count):
        power = k*power
        answer += (count-lag)*(power+power.T)
    return sp.ImmutableMatrix(CHI*answer/count)


def covariance_error_bound(clock, count):
    """Rational Euclidean operator bound, valid for every positive count."""
    _integer(count, 1, 10**12)
    field(clock)
    # L>=2 P0 and L<=10 P0; for incident D between (1/6)I and (1/3)I.
    # Similarity condition number <=sqrt(2)<2; full-attempt D=I.
    r, condition = (F(29, 30), F(1)) if clock == 'attempt' else (F(9, 10), F(2))
    return 2*F(CHI)*condition*r/(count*(1-r)**2)


def scalar_moments(clock, count):
    """Exact first four moments of the SAME history's port-zero sum.

    Integer score 6*n_0-1 is converted back to centered occupation. The
    recursion retains every terminal configuration and all counter moments.
    """
    _integer(count, 1, 64)
    states, _, _, _, _ = native_counts()
    kernel = configuration_kernel(clock)
    denominator = int(sp.ilcm(*(x.q for x in kernel)))
    counts = [[int(x*denominator) for x in row] for row in kernel.tolist()]
    destinations = [[(i, counts[i][j]) for i in range(66) if counts[i][j]] for j in range(66)]
    score = [6*state[0]-1 for state in states]
    values = [[s**order for s in score] for order in range(5)]
    for _ in range(1, count):
        propagated = [[0]*66 for _ in range(5)]
        for j in range(66):
            for i, weight in destinations[j]:
                for order in range(5):
                    propagated[order][i] += weight*values[order][j]
        values = [[sum(math.comb(order, old)*score[i]**(order-old)*propagated[old][i]
                       for old in range(order+1)) for i in range(66)] for order in range(5)]
    normalization = 66*denominator**(count-1)
    moments = tuple(F(sum(values[order]), normalization*6**order) for order in range(5))
    if moments[0] != 1 or moments[1] != 0:
        raise ArithmeticError('the full counter law lost normalization or centering')
    if moments[2]/count != covariance_at_count(clock, count)[0, 0]:
        raise ArithmeticError('independent complete-history second moment disagrees')
    return {'clock': clock, 'count': count, 'sum_moments': [str(x) for x in moments],
            'normalized_fourth_cumulant': str((moments[4]-3*moments[2]**2)/count**2),
            'every_terminal_configuration_retained': True}


def incomplete_tape_control():
    """Same initial state/local incident tape, different remote field values."""
    before = tuple(int(i in (1, 2)) for i in range(12))
    seams, _ = passage.recovered_carrier()
    if (1, 5) not in seams or (0, 10) not in seams:
        raise ArithmeticError('the native access counterexample changed its seams')
    outputs, tapes = [], []
    for coin in (0, 1):
        after, _ = native.native(before, (1, 5), coin)
        after, record = native.native(after, (0, 10), 0)
        outputs.append(after); tapes.append(record)
    if tapes[0] != tapes[1] or outputs[0] == outputs[1] or outputs[0][0] != outputs[1][0]:
        raise ArithmeticError('the incomplete-tape control no longer separates access')
    return {'initial_occupied_ports': [1, 2], 'remote_seam': [1, 5],
            'first_incident_seam': [0, 10], 'identical_incident_record': list(tapes[0]),
            'terminal_occupied_ports': [[i for i, x in enumerate(out) if x] for out in outputs],
            'port_only_tape_determines_the_full_field': False}


def torus_edges(size):
    _integer(size, 4, 256)
    return tuple(sorted({tuple(sorted((size*x+y, size*((x+dx) % size)+(y+dy) % size)))
                         for x, y in product(range(size), repeat=2) for dx, dy in ((1, 0), (0, 1))}))


def torus_native_control(size=8):
    """Execute the unchanged native rule; topology is explicitly supplied."""
    _integer(size, 4, 8)
    edges = torus_edges(size)
    ports, m = size**2, len(edges)
    # Single-particle enumeration audits the complete density matrix. The local
    # four-state identity proves its extension to every fixed-number sector.
    matrix = sp.zeros(ports)
    cases = 0
    for before in product((0, 1), repeat=2):
        outputs = [native.native(before, (0, 1), coin)[0] for coin in (0, 1)]
        if sorted(outputs) != sorted((before, before[::-1])):
            raise ArithmeticError('native binary exchange identity failed')
    for source in range(ports):
        before = tuple(int(i == source) for i in range(ports))
        for edge in edges:
            for coin in (0, 1):
                after, _ = native.native(before, edge, coin)
                matrix[after.index(1), source] += 1
                cases += 1
    if matrix != 2*m*sp.eye(ports)-_laplacian(ports, edges):
        raise ArithmeticError('declared torus changed the actual native density law')
    return {'supplied_torus_side': size, 'ports': ports, 'seams': m,
            'native_branch_checks': cases, 'arbitrary_binary_sector_identity_cases': 4,
            'this_is_the_existing_spherical_source_screen': False}


def torus_bounds(size):
    """Certified fixed-band Green and incident-clock error, using pi**2<10.

    Each listed integer wavevector supplies real cos and sin coordinates.
    The spatial scale is a=N**2/(4*pi**2). This is a declared torus chart.
    """
    _integer(size, 8, 65536)
    errors, positive = [], []
    for k in MODES:
        energy = sum(x*x for x in k)
        loss = F(10*sum(x**4 for x in k), 3*size**2)
        lower = energy-loss
        if lower <= 0:
            raise ArithmeticError('fixed Fourier band exceeds the certified positive domain')
        errors.append(loss/(energy*lower)); positive.append(lower)
    green = max(errors)
    # ||W^T L_inc W|| <= 16*pi**2 sum_real_modes |k|**2 / N**4.
    # W is an exact orthonormal sampling of the selected real Fourier modes.
    total_energy = 2*sum(sum(x*x for x in k) for k in MODES)
    incident = F(10, 4*size**2)+F(8*10*total_energy, 4*size**4)/min(positive)
    attempt = F(5, size**4)
    return {'size': size, 'real_mode_count': 2*len(MODES),
            'canonical_total': size**2//3,
            'susceptibility': str(F((size**2//3)*(size**2-size**2//3), size**2*(size**2-1))),
            'green_operator_error_bound': str(green),
            'attempt_clock_correction_bound': str(attempt),
            'incident_clock_correction_bound': str(incident),
            'attempt_total_covariance_error_bound': str(green+attempt),
            'incident_total_covariance_error_bound': str(green+incident),
            'clock_to_clock_covariance_error_bound': str(attempt+incident)}


def torus_diagnostic(size):
    """Independent floating band-matrix check; rational bounds certify errors."""
    bounds = torus_bounds(size)
    points = np.asarray(tuple(product(range(size), repeat=2)), dtype=float)
    columns, energies = [], []
    for k in MODES:
        phase = 2*np.pi*(points@np.asarray(k))/size
        columns.extend((np.sqrt(2)/size*np.cos(phase), np.sqrt(2)/size*np.sin(phase)))
        energy = (size/np.pi)**2*sum(math.sin(math.pi*x/size)**2 for x in k)
        energies.extend((energy, energy))
    w = np.asarray(columns).T
    if np.linalg.norm(w.T@w-np.eye(len(columns))) > 1e-12 or np.linalg.norm(w.sum(axis=0)) > 1e-11:
        raise ArithmeticError('the declared Fourier observation map lost orthonormality')
    target = np.diag([1/sum(x*x for x in k) for k in MODES for _ in range(2)])
    green = np.diag(1/np.asarray(energies))
    local = np.zeros_like(green)
    for a, b in torus_edges(size):
        if 0 in (a, b):
            delta = w[a]-w[b]
            local += np.outer(delta, delta)
    a_scale = size**2/(4*np.pi**2)
    attempt = green-np.eye(len(columns))/(4*(2*size**2)*a_scale)
    incident = green+np.eye(len(columns))/(16*a_scale)-(local@green+green@local)/16
    errors = {'attempt': np.linalg.norm(attempt-target, 2),
              'incident': np.linalg.norm(incident-target, 2)}
    for clock, error in errors.items():
        if error > float(F(bounds[clock+'_total_covariance_error_bound']))+1e-12:
            raise ArithmeticError('Fourier covariance diagnostic exceeded its analytic allowance')
    return {**bounds, 'attempt_error_diagnostic': float(errors['attempt']),
            'incident_error_diagnostic': float(errors['incident']),
            'incident_cross_mode_diagnostic': float(incident[1, 5])}


def run():
    _, _, _, z, checks = native_counts()
    rows = []
    for clock in ('attempt', 'incident'):
        data = field(clock)
        complete = configuration_kernel(clock)
        # Independent full-state lag laws, including noncommuting incident map.
        power = sp.eye(66)
        for lag in range(4):
            if z*power*z.T/66 != data['kernel']**lag*(CHI*P0):
                raise ArithmeticError('complete native lag covariance disagrees with the field')
            power = complete*power
        rows.append({'clock': clock, 'incident_seams': data['m'],
                     'covariance': _rows(data['covariance']),
                     'action_precision': _rows(data['precision']),
                     'positive_LDL_pivots': [str(x) for x in data['positive_LDL_pivots']],
                     'finite_count_covariance_error_bounds': {
                         str(n): str(covariance_error_bound(clock, n)) for n in (128, 4096, 65536)}})
    attempt, incident = field('attempt'), field('incident')
    correction = incident['covariance']/5-attempt['covariance']/30
    if correction == sp.zeros(12):
        raise ArithmeticError('the negative control lost its random-clock sampling defect')
    return {'schema': 'oph.physical_space_time.occupation_field.v1',
            'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            'contract': {'stationary_binary_total': 2, 'ports': 12,
                         'stationary_configurations': 66, 'native_branch_checks': checks,
                         'clock_names': ['native attempts', 'port-zero incident attempts'],
                         'initial_occupation_disclosure_and_all_seam_currents_supplied': True,
                         'same_history_no_independent_replicas': True,
                         'physical_seconds_or_Hamiltonian_supplied': False},
            'field_laws': rows,
            'clock_negative_control': {'mean_attempts_per_incident_tick': 6,
                'equal_mean_duration_rescaling_makes_fields_equal': False,
                'scaled_covariance_difference_00': str(correction[0, 0]),
                'scaled_covariance_difference_11': str(correction[1, 1]),
                'scaled_covariance_difference_01': str(correction[0, 1])},
            'observer_access_negative_control': incomplete_tape_control(),
            'finite_history_controls': [scalar_moments(clock, n) for clock, ns in
                (('attempt', (1, 4, 16, 64)), ('incident', (1, 4, 16))) for n in ns],
            'finite_cap_control': {'requested_incident_ticks': 64,
                'attempt_cap_per_interval': 256,
                'complete_output_timeout_bound': str(64*F(5, 6)**256),
                'retaining_timeouts_without_success_renormalization_required': True,
                'full_cap_256_history_distribution_executed': False},
            'declared_torus_native_control': torus_native_control(),
            'declared_torus_green_convergence': [torus_diagnostic(n) for n in (8, 16, 32, 64)],
            'acceptance': {'one_history_Gaussian_field_limit_and_action_proved': True,
                'actual_incident_clock_correction_derived': True,
                'conditional_graph_Green_convergence_criterion_proved': True,
                'supplied_torus_fixed_band_iterated_limit_proved': True,
                'uniform_finite_count_Gaussian_error_bound_proved': False,
                'existing_spherical_source_continuum_criterion_verified': False,
                'instantaneous_native_occupations_are_Gaussian': False,
                'interacting_CFT_or_Lorentzian_spacetime_derived': False,
                'observed_Universe_identified': False},
            'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--output', type=Path)
    modes.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Native occupation field, clock correction and Green-limit controls replayed.')
    elif args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
