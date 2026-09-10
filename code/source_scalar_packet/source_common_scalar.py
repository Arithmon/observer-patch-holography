"""Certified massive free scalar readouts on prepared golden source addresses.

The tensor action, protected address register, Dirichlet boundary, quantization
and model time are declared. The certificate includes every orthogonal lattice
mode through resolvent bounds; it does not execute a full causal read trace.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from source_scalar_intervals import I, PI, SCALE, sin, exp, dot, nonnegative

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'source_common_scalar_receipt.json'
LEVELS = (7, 9, 11, 13)
MODES = 9
GAMMA = 4
PINS = (
    'code/source_scalar_packet/source_scalar_intervals.py',
    'code/source_scalar_packet/source_common_scalar.py',
    'code/source_scalar_packet/verify_source_common_scalar.py',
    'code/source_scalar_packet/test_source_common_scalar.py',
    'paper/tex_fragments/SOURCE_COMMON_SCALAR_PACKET.tex',
    'code/causal_refinement/source_net_causet.py',
    'Lean/Screen/PrimitivePortFrameQuotient.lean',
    'Lean/Screen/PortFrameGram.lean',
)
SCOPE = {
    'actual_golden_source_addresses': True,
    'prepared_addresses_protected_separately_from_field': True,
    'tensor_action_axes_dirichlet_boundary_quantization_supplied': True,
    'old_causal_radius_contains_all_field_edges': True,
    'same_as_previous_radial_kernel_action': False,
    'full_orthogonal_stiffness_and_frequency_leakage_bounded': True,
    'full_fock_coherent_states_and_bounded_effects': True,
    'compact_disjoint_preparation_detector_slabs': True,
    'all_causal_read_write_events_executed': False,
    'count_clock_custody_for_field_execution_established': False,
    'native_action_or_physical_time_selected': False,
    'interacting_or_charged_quantum_continuum': False,
    'floating_simulation_used_as_certificate': False,
}


def canonical(x):
    return (json.dumps(x, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True, allow_nan=False)+'\n').encode('ascii')


def fibonacci(n):
    if type(n) is not int or n not in LEVELS:
        raise ValueError('only the four declared finite refinement levels are accepted')
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a, b


def grid(n):
    q, p = fibonacci(n)
    phi = (1+I.of(5).sqrt())/2
    # b*p mod q is the exact golden-orbit order for these odd levels.
    indices = sorted(range(q), key=lambda b: b*p % q)
    pairs = [(-((b+isqrt(5*b*b))//2), b) for b in indices]
    points = [a+b*phi for a, b in pairs]+[I.of(1)]
    gaps = [b-a for a, b in zip(points, points[1:])]
    if min(v.lo for v in gaps) <= 0:
        raise ValueError('golden sites are not strictly ordered')
    if max((q*h*h).hi for h in gaps) >= SCALE:
        raise ValueError('a field edge exceeds the declared causal radius')
    weights = [(a+b)/2 for a, b in zip(gaps, gaps[1:])]
    return q, pairs, points, gaps, weights


def solve(diagonal, off, rhs):
    """Outward Thomas elimination for a symmetric positive tridiagonal form."""
    d, r = list(diagonal), list(rhs)
    if len(off)+1 != len(d) or len(r) != len(d):
        raise ValueError('invalid tridiagonal dimensions')
    for i in range(1, len(d)):
        if d[i-1].lo <= 0:
            raise ValueError('uncertified positive elimination pivot')
        a = off[i-1]/d[i-1]
        d[i], r[i] = d[i]-a*off[i-1], r[i]-a*r[i-1]
    if d[-1].lo <= 0:
        raise ValueError('uncertified final pivot')
    out = [I.of(0)]*len(d)
    out[-1] = r[-1]/d[-1]
    for i in range(len(d)-2, -1, -1):
        out[i] = (r[i]-off[i]*out[i+1])/d[i]
    return out


def coefficients():
    a = [F(0)]*(MODES+1)
    a[1] = F(1)
    for _ in range(MODES-1):
        b = list(a)
        for j in range(1, len(a)):
            if j+1 < len(a):
                b[j+1] += a[j]/2
            if j-1 > 0:
                b[j-1] += a[j]/2
        a = b
    a = a[1:]
    norm = sum(v*v for v in a)
    return a, norm, [I.of(v)/I.of(norm).sqrt() for v in a]


def continuum_tail(coeff):
    lower = F(9, 20)
    total = I.of(0)
    for k, a in enumerate(coeff, 1):
        for l, b in enumerate(coeff, 1):
            if k == l:
                integ = 1-lower+sin(2*k*PI*lower)/(2*k*PI)
            else:
                integ = -sin((k-l)*PI*lower)/((k-l)*PI)+sin((k+l)*PI*lower)/((k+l)*PI)
            total += a*b*integ
    return nonnegative(total)


def probability_error(b, alpha, db, da):
    """Gaussian characteristic function Lipschitz bound, all arguments >=0."""
    return db*(2*b+db)/4+db*(alpha+da)+b*da


def above(point, cut):
    boundary = I.of(cut)
    if point.lo > boundary.hi:
        return True
    if point.hi < boundary.lo:
        return False
    raise ValueError('slab boundary classification is not certified')


def continuum():
    raw, norm, coeff = coefficients()
    omega = [(1+PI*PI*(k*k+2)).sqrt() for k in range(1, MODES+1)]
    bnorm = nonnegative(sum((a*a/w for a, w in zip(coeff, omega)), I.of(0))/2).sqrt()
    tail = continuum_tail(coeff)
    # Cut the actual preparation and detector to disjoint slabs; only their
    # smooth comparison vectors have a finite Fourier expansion.
    db = (tail/2).sqrt().upper()
    da = GAMMA*db
    cut_error = probability_error(bnorm.upper(), GAMMA*bnorm.upper(), db, da).upper()
    samples = []
    # 100 closed cells cover the ENTIRE window, not just sampled times.
    for j in range(100):
        t = I.of(F(19, 20)+F(j, 2000))
        t = I(t.lo, I.of(F(19, 20)+F(j+1, 2000)).hi)
        mean = GAMMA*sum((a*a*((-1)**(k+1))*sin(w*t)/w
                         for k, (a, w) in enumerate(zip(coeff, omega), 1)), I.of(0))
        response = exp(-bnorm*bnorm/2)*sin(mean)/2
        samples.append(response)
    response = I(min(v.lo for v in samples), max(v.hi for v in samples))
    return raw, norm, coeff, omega, bnorm, tail, cut_error, response


def build_level(n, context):
    q, pairs, points, gaps, weights = grid(n)
    raw, norm, coeff, omega, bnorm, tail, cut_error, response = context
    x = points[1:-1]
    f = [[I.of(2).sqrt()*sin(k*PI*z) for z in x] for k in range(1, MODES+1)]
    diagonal = [1/a+1/b for a, b in zip(gaps, gaps[1:])]
    off = [-1/h for h in gaps[1:-1]]
    residual = []
    for k, row in enumerate(f, 1):
        vals = []
        for i, v in enumerate(row):
            kv = diagonal[i]*v
            if i:
                kv += off[i-1]*row[i-1]
            if i+1 < len(row):
                kv += off[i]*row[i+1]
            vals.append(kv/weights[i]-(k*PI)**2*v)
        residual.append(vals)
    gram = [[dot(weights, [a*b for a, b in zip(u, v)]) for v in f] for u in f]
    g1 = gram[0][0]
    G = [[g1*g1*a for a in row] for row in gram]
    delta = I.of(0)
    for k, row in enumerate(G):
        row_sum = sum(((a-int(k == l)).abs() for l, a in enumerate(row)), I.of(0))
        delta = I(max(delta.lo, row_sum.lo), max(delta.hi, row_sum.hi))
    delta = delta.upper()
    if delta.hi >= SCALE:
        raise ValueError('sampled Gram matrix has no certified inverse')
    norms = [nonnegative(gram[k][k]).sqrt() for k in range(MODES)]
    energies, columns = [], []
    for k, w in enumerate(omega):
        d = [a+(1+w*w)*b for a, b in zip(diagonal, weights)]
        qs = []
        for j in (k, 0):
            rhs = [a*b for a, b in zip(weights, residual[j])]
            y = solve(d, off, rhs)
            qs.append(nonnegative(dot(rhs, y)))
        energies.append([v.encode() for v in qs])
        columns.append(qs[0].sqrt()*norms[0]**2+2*qs[1].sqrt()*norms[k]*norms[0])
    sampling = ((1+delta).sqrt()*(1/(1-delta).sqrt()-1)).upper()
    root_columns = (sum((v*v for v in columns), I.of(0))).sqrt()
    frequency = (root_columns/(1-delta).sqrt()
                 +2*omega[-1]*sampling).upper()
    # A coefficient-sensitive FULL residual bound; it remains valid after
    # componentwise diagonal contraction, which the resolvent proof consumes.
    correction = root_columns*(1/(1-delta).sqrt()-1)+2*omega[-1]*sampling
    kappa_f = (sum((v*a.abs() for v, a in zip(columns, coeff)), I.of(0))+correction).upper()
    kappa_v = (sum((v*a.abs()/w.sqrt() for v, a, w in zip(columns, coeff, omega)), I.of(0))
               +correction*I.of(2).sqrt()*bnorm).upper()
    # m=1: resolvent integration gives the negative-square-root bound kappa_f/2.
    db = ((kappa_f/2+sampling)/I.of(2).sqrt()).upper()
    alpha = GAMMA*bnorm.upper()
    da = (GAMMA*db+GAMMA*kappa_v/I.of(2).sqrt()).upper()  # uniform |t|<=1
    band_error = probability_error(bnorm.upper(), alpha, db, da).upper()
    # Independently localise the actual GRAPH smearing; no continuous-cell
    # tail bound is silently substituted for this sampled tail.
    fvalues = [sum((a*row[i] for a, row in zip(coeff, f)), I.of(0)) for i in range(len(x))]
    gvalues = [sum((a*((-1)**k)*row[i] for k, (a, row) in enumerate(zip(coeff, f))), I.of(0))
               for i in range(len(x))]
    ft = g1*g1*sum((weights[i]*v*v for i, v in enumerate(fvalues)
                   if above(x[i], F(9, 20))), I.of(0))
    gt = g1*g1*sum((weights[i]*v*v for i, v in enumerate(gvalues)
                   if not above(x[i], F(11, 20))), I.of(0))
    energy = 8*g1*g1*sum((weights[i]*v*v for i, v in enumerate(fvalues)
                          if not above(x[i], F(9, 20))), I.of(0))
    fcut, gcut = (nonnegative(ft)/2).sqrt().upper(), (nonnegative(gt)/2).sqrt().upper()
    graph_cut_error = probability_error((bnorm+db).upper(), (alpha+GAMMA*db).upper(),
                                       gcut, GAMMA*fcut).upper()
    # Comparison actual graph compact -> band continuum -> actual continuum compact.
    total = (band_error+graph_cut_error+cut_error).upper()
    graph_signal_error = (band_error+graph_cut_error).upper()
    resolved = -response.hi-graph_signal_error.hi
    return {
        'fibonacci_level': n, 'q': q, 'dynamic_oscillators': (q-1)**3,
        'tensor_algorithm_largest_matrix': q-1,
        'ordered_source_coordinates_Qphi': pairs,
        'dimensionless_max_gap': I(max(h.lo for h in gaps), max(h.hi for h in gaps)).encode(),
        'all_field_edges_within_old_radius': True,
        'gram_row_sum_error_upper': delta.encode(),
        'full_resolvent_column_energies': energies,
        'full_frequency_residual_upper': frequency.encode(),
        'full_preparation_frequency_residual_upper': kappa_f.encode(),
        'full_coherent_frequency_residual_upper': kappa_v.encode(),
        'sampling_isometry_error_upper': sampling.encode(),
        'sampled_preparation_tail_mass': nonnegative(ft).encode(),
        'sampled_detector_tail_mass': nonnegative(gt).encode(),
        'compact_preparation_energy': nonnegative(energy).encode(),
        'band_probability_error_upper': band_error.encode(),
        'graph_localisation_error_upper': graph_cut_error.encode(),
        'graph_vs_compact_continuum_error_upper': I(min(SCALE, total.lo), min(SCALE, total.hi)).encode(),
        'graph_compact_induced_response': I(max(-SCALE//2, response.lo-graph_signal_error.hi),
                                            min(SCALE//2, response.hi+graph_signal_error.hi)).encode(),
        'resolved_graph_signal_lower': I(max(0, resolved), max(0, resolved)).encode()[0],
        'error_smaller_than_resolved_signal': total.hi < resolved,
    }


def build():
    context = continuum()
    raw, norm, coeff, omega, bnorm, tail, cut_error, response = context
    pins = {name: {'sha256': sha256((ROOT/name).read_bytes()).hexdigest(),
                   'bytes': (ROOT/name).stat().st_size} for name in PINS}
    return {
        'schema': 'oph.source-common-massive-scalar.v1', 'scope': SCOPE,
        'units': {'source_length_L': '2/sqrt(phi+2)', 'coordinates': 'x/L',
                  'time': 'c*t/L', 'hbar': '1', 'dimensionless_mass': '1',
                  'causal_radius': 'L/sqrt(q)', 'causal_layer_duration': 'L/(c*sqrt(q))'},
        'preparation': {'sine_coefficients': [str(v) for v in raw],
                        'coefficient_squared_norm': str(norm), 'momentum_impulse': GAMMA,
                        'support': '[0,9/20] x [0,1]^2',
                        'detector_support': '[11/20,1] x [0,1]^2',
                        'detector_effect': '(I+sin(phi(g)))/2', 'baseline_probability': '1/2'},
        'continuum': {'smooth_comparison_tail_mass': tail.encode(),
                      'coherent_mode_norm': bnorm.encode(),
                      'localisation_probability_error_upper': cut_error.encode(),
                      'window': ['19/20', '1'], 'closed_time_cells': 100,
                      'band_induced_response': response.encode(),
                      'compact_induced_response': I(response.lo-cut_error.hi,
                                                    response.hi+cut_error.hi).encode()},
        'levels': [build_level(n, context) for n in LEVELS], 'source_pins': pins,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    data = canonical(build())
    if args.check:
        if OUTPUT.read_bytes() != data:
            raise SystemExit('source scalar receipt drift')
    else:
        OUTPUT.write_bytes(data)
    print('source common scalar: exact outward generation PASS')


if __name__ == '__main__':
    main()
