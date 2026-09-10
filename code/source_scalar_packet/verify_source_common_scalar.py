"""Independent edge-flux and reverse-elimination scalar certificate replay.

No producer is imported. Only the separately tested rational interval arithmetic
is shared. Every output bound, tail, time cell and current source pin is checked.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb, isqrt
from pathlib import Path
import json

from source_scalar_intervals import I, PI, SCALE, sin, exp, nonnegative

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE/'source_common_scalar_receipt.json'
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


def require(value, message):
    if not value:
        raise ValueError(message)


def exact(got, expected, path='receipt'):
    require(type(got) is type(expected), path+': wrong type')
    if isinstance(expected, dict):
        require(got.keys() == expected.keys(), path+': wrong keys')
        for k in expected:
            exact(got[k], expected[k], path+'/'+k)
    elif isinstance(expected, list):
        require(len(got) == len(expected), path+': wrong length')
        for i, (a, b) in enumerate(zip(got, expected)):
            exact(a, b, path+'/'+str(i))
    else:
        require(got == expected, path+': independently recomputed value differs')


def unique(items):
    out = {}
    for key, value in items:
        require(key not in out, 'duplicate JSON key')
        out[key] = value
    return out


def forbidden(token):
    raise ValueError('floating/nonfinite JSON tokens are forbidden')


def load(path=OUTPUT):
    raw = Path(path).read_bytes()
    require(len(raw) <= 500_000, 'receipt size limit')
    return json.loads(raw.decode('utf-8'), object_pairs_hook=unique,
                      parse_float=forbidden, parse_constant=forbidden)


def summands(values):
    return sum(values, I.of(0))


def coeffs():
    # Independent Laurent/binomial expansion of sin(x)(1+cos(x))^8.
    def choose(j):
        return comb(16, j) if 0 <= j <= 16 else 0
    raw = [F(choose(8+k-1)-choose(8+k+1), 256) for k in range(1, 10)]
    norm = sum(v*v for v in raw)
    return raw, norm, [I.of(a)/I.of(norm).sqrt() for a in raw]


def characteristic_bound(b, a, dv, dz):
    # Vacuum characteristic factor plus the bilinear phase difference.
    return dv*(2*b+dv)/4 + dv*a + dv*dz + b*dz


def reference():
    raw, norm, aa = coeffs()
    omega = [(1+(k*k+2)*PI**2).sqrt() for k in range(1, 10)]
    b = nonnegative(summands(a*a/(2*w) for a, w in zip(aa, omega))).sqrt()
    # Integrate the cosine polynomial for |f|^2; this is separate from the
    # producer's pairwise sine-product antiderivatives.
    cc = [I.of(0)]*19
    for k in range(1, 10):
        for l in range(1, 10):
            cc[abs(k-l)] += aa[k-1]*aa[l-1]
            cc[k+l] -= aa[k-1]*aa[l-1]
    tail = nonnegative(cc[0]*F(11, 20)-summands(
        cc[k]*sin(k*PI*F(9, 20))/(k*PI) for k in range(1, 19)))
    v = (tail/2).sqrt().upper()
    cut = characteristic_bound(b.upper(), 4*b.upper(), v, 4*v).upper()
    outcomes = []
    for j in range(100):
        t = I(I.of(F(1900+j, 2000)).lo, I.of(F(1901+j, 2000)).hi)
        phase = summands(4*((-1)**k)*a*a*sin(w*t)/w
                         for k, (a, w) in enumerate(zip(aa, omega)))
        outcomes.append(exp(-b*b/2)*sin(phase)/2)
    response = I(min(x.lo for x in outcomes), max(x.hi for x in outcomes))
    return raw, norm, aa, omega, b, tail, cut, response


def reverse_solve(diag, edge, rhs):
    """Eliminate from the right, opposite to the producer's Thomas solve."""
    d, r = list(diag), list(rhs)
    for j in range(len(d)-2, -1, -1):
        require(d[j+1].lo > 0, 'nonpositive reverse pivot')
        d[j] -= edge[j]*edge[j]/d[j+1]
        r[j] -= edge[j]*r[j+1]/d[j+1]
    require(d[0].lo > 0, 'nonpositive first pivot')
    y = [r[0]/d[0]]
    for j in range(1, len(d)):
        y.append((r[j]-edge[j-1]*y[-1])/d[j])
    return y


def verify_source_coordinates(q, pairs):
    # Verify the source six-axis control symbolically in Z[phi]. For each
    # coordinate control basis, sum z_i times the six primitive icosahedral axes.
    axes = (((-1, 0), (0, 1), (0, 0)), ((1, 0), (0, 1), (0, 0)),
            ((0, 0), (-1, 0), (0, 1)), ((0, 0), (1, 0), (0, 1)),
            ((0, 1), (0, 0), (-1, 0)), ((0, 1), (0, 0), (1, 0)))
    for v in range(6):
        a, b = [0]*3, [0]*3
        (a if v < 3 else b)[v % 3] = 1
        z = [b[1]-a[0], b[1]+a[0], b[2]-a[1], b[2]+a[1], b[0]-a[2], b[0]+a[2]]
        require(sum(z) % 2 == 0, 'source control is not D6')
        result = [[sum(z[i]*axes[i][j][t] for i in range(6)) for t in range(2)]
                  for j in range(3)]
        require(result == [[2*a[j], 2*b[j]] for j in range(3)], 'source Gram control mismatch')
    require(sorted(b for a, b in pairs) == list(range(q)), 'golden census mismatch')
    for a, b in pairs:
        # floor(b phi) independently bounded by squaring integer sqrt(5)b.
        require(a == -((b+isqrt(5*b*b))//2),
                'wrong golden integer control')


def level(n, context):
    fib = [0, 1]
    for _ in range(n):
        fib.append(fib[-1]+fib[-2])
    q = fib[n]
    root5 = I.of(5).sqrt()
    pairs = [[-((b+isqrt(5*b*b))//2), b] for b in range(q)]
    def value(pair):
        a, b = pair
        return (2*a+b+b*root5)/2
    pairs.sort(key=lambda pair: value(pair).lo)
    verify_source_coordinates(q, pairs)
    x = [value(v) for v in pairs]+[I.of(1)]
    lengths = [x[i+1]-x[i] for i in range(q)]
    require(all(v.lo > 0 and (q*v*v).hi < SCALE for v in lengths), 'edge/radius failure')
    mass = [(lengths[i]+lengths[i+1])/2 for i in range(q-1)]
    modes = [[I.of(2).sqrt()*sin(k*PI*z) for z in x[1:-1]] for k in range(1, 10)]
    def inner(u, v):
        return summands(m*a*b for m, a, b in zip(mass, u, v))
    rr = []
    for k, row in enumerate(modes, 1):
        extended = [I.of(0)]+row+[I.of(0)]
        flux = [(extended[i+1]-extended[i])/lengths[i] for i in range(q)]
        rr.append([(flux[i]-flux[i+1])/mass[i]-(k*PI)**2*row[i] for i in range(q-1)])
    gram = [[inner(a, b) for b in modes] for a in modes]
    g1 = gram[0][0]
    G = [[g1*g1*v for v in row] for row in gram]
    rowbounds = [summands((v-int(i == j)).abs() for j, v in enumerate(row))
                 for i, row in enumerate(G)]
    d = I(max(v.hi for v in rowbounds), max(v.hi for v in rowbounds))
    require(d.hi < SCALE, 'noninvertible sampled Gram bound')
    norms = [nonnegative(gram[k][k]).sqrt() for k in range(9)]
    _, _, aa, omega, b, tail, cut, response = context
    base_diag = [1/lengths[i]+1/lengths[i+1] for i in range(q-1)]
    edge = [-1/v for v in lengths[1:-1]]
    energies, col = [], []
    for k, om in enumerate(omega):
        sigma = 1+om*om
        diagonal = [v+sigma*m for v, m in zip(base_diag, mass)]
        qs = []
        for j in (k, 0):
            y = reverse_solve(diagonal, edge, [m*r for m, r in zip(mass, rr[j])])
            # Independent energy identity y^T(K+sigma W)y = rhs^T y.
            yy = [I.of(0)]+y+[I.of(0)]
            energy = summands((yy[i+1]-yy[i])**2/lengths[i] for i in range(q))
            energy += sigma*inner(y, y)
            qs.append(nonnegative(energy))
        energies.append([a.encode() for a in qs])
        col.append(qs[0].sqrt()*norms[0]**2+2*qs[1].sqrt()*norms[k]*norms[0])
    c = 1/(1-d).sqrt()-1
    sd = ((1+d).sqrt()*c).upper()
    r = summands(v*v for v in col).sqrt()
    op = (r/(1-d).sqrt()+2*omega[-1]*sd).upper()
    correction = r*c+2*omega[-1]*sd
    kf = (summands(v*a.abs() for v, a in zip(col, aa))+correction).upper()
    kv = (summands(v*a.abs()/w.sqrt() for v, a, w in zip(col, aa, omega))
          +correction*I.of(2).sqrt()*b).upper()
    db = ((kf/2+sd)/I.of(2).sqrt()).upper()
    alpha = 4*b.upper()
    da = (4*db+4*kv/I.of(2).sqrt()).upper()
    band = characteristic_bound(b.upper(), alpha, db, da).upper()
    fv = [summands(a*row[i] for a, row in zip(aa, modes)) for i in range(q-1)]
    gv = [summands(((-1)**k)*a*row[i] for k, (a, row) in enumerate(zip(aa, modes)))
          for i in range(q-1)]
    left, right = [], []
    for z in x[1:-1]:
        c, dcut = I.of(F(9, 20)), I.of(F(11, 20))
        require(z.hi < c.lo or z.lo > c.hi, 'ambiguous preparation cutoff')
        require(z.hi < dcut.lo or z.lo > dcut.hi, 'ambiguous detector cutoff')
        left.append(z.hi < c.lo)
        right.append(z.lo > dcut.hi)
    ft = nonnegative(g1*g1*summands(m*v*v for keep, m, v in zip(left, mass, fv) if not keep))
    gt = nonnegative(g1*g1*summands(m*v*v for keep, m, v in zip(right, mass, gv) if not keep))
    energy = nonnegative(8*g1*g1*summands(m*v*v for keep, m, v in zip(left, mass, fv) if keep))
    dc = (gt/2).sqrt().upper()
    ac = 4*(ft/2).sqrt().upper()
    graphcut = characteristic_bound((b+db).upper(), (alpha+4*db).upper(), dc, ac).upper()
    total = (band+graphcut+cut).upper()
    grapherror = (band+graphcut).upper()
    signal = max(0, -response.hi-grapherror.hi)
    return {
        'fibonacci_level': n, 'q': q, 'dynamic_oscillators': (q-1)**3,
        'tensor_algorithm_largest_matrix': q-1,
        'ordered_source_coordinates_Qphi': pairs,
        'dimensionless_max_gap': I(max(v.lo for v in lengths), max(v.hi for v in lengths)).encode(),
        'all_field_edges_within_old_radius': True,
        'gram_row_sum_error_upper': d.encode(), 'full_resolvent_column_energies': energies,
        'full_frequency_residual_upper': op.encode(),
        'full_preparation_frequency_residual_upper': kf.encode(),
        'full_coherent_frequency_residual_upper': kv.encode(),
        'sampling_isometry_error_upper': sd.encode(),
        'sampled_preparation_tail_mass': ft.encode(), 'sampled_detector_tail_mass': gt.encode(),
        'compact_preparation_energy': energy.encode(),
        'band_probability_error_upper': band.encode(), 'graph_localisation_error_upper': graphcut.encode(),
        'graph_vs_compact_continuum_error_upper': I(min(SCALE, total.lo), min(SCALE, total.hi)).encode(),
        'graph_compact_induced_response': I(max(-SCALE//2, response.lo-grapherror.hi),
                                            min(SCALE//2, response.hi+grapherror.hi)).encode(),
        'resolved_graph_signal_lower': I(signal, signal).encode()[0],
        'error_smaller_than_resolved_signal': total.hi < -response.hi-grapherror.hi,
    }


def verify(data):
    require(type(data) is dict, 'receipt must be an object')
    pins = {p: {'sha256': sha256((ROOT/p).read_bytes()).hexdigest(),
                'bytes': (ROOT/p).stat().st_size} for p in PINS}
    raw, norm, aa, omega, b, tail, cut, response = context = reference()
    expected = {
        'schema': 'oph.source-common-massive-scalar.v1',
        'scope': {
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
        },
        'units': {'source_length_L': '2/sqrt(phi+2)', 'coordinates': 'x/L', 'time': 'c*t/L',
                  'hbar': '1', 'dimensionless_mass': '1', 'causal_radius': 'L/sqrt(q)',
                  'causal_layer_duration': 'L/(c*sqrt(q))'},
        'preparation': {'sine_coefficients': [str(v) for v in raw],
                        'coefficient_squared_norm': str(norm), 'momentum_impulse': 4,
                        'support': '[0,9/20] x [0,1]^2',
                        'detector_support': '[11/20,1] x [0,1]^2',
                        'detector_effect': '(I+sin(phi(g)))/2', 'baseline_probability': '1/2'},
        'continuum': {'smooth_comparison_tail_mass': tail.encode(), 'coherent_mode_norm': b.encode(),
                      'localisation_probability_error_upper': cut.encode(),
                      'window': ['19/20', '1'], 'closed_time_cells': 100,
                      'band_induced_response': response.encode(),
                      'compact_induced_response': I(response.lo-cut.hi, response.hi+cut.hi).encode()},
        'source_pins': pins,
    }
    require(data.keys() == expected.keys() | {'levels'}, 'wrong receipt fields')
    exact({k: v for k, v in data.items() if k != 'levels'}, expected)
    require(type(data['levels']) is list and len(data['levels']) == 4, 'wrong refinement census')
    for observed, n in zip(data['levels'], (7, 9, 11, 13)):
        last = level(n, context)
        exact(observed, last, 'level/'+str(n))
    require(last['error_smaller_than_resolved_signal'], 'fine packet does not resolve the signal')
    require(F(last['graph_vs_compact_continuum_error_upper'][1]) < F('0.050055'), 'error theorem mismatch')
    require(F(last['resolved_graph_signal_lower']) > F('0.137754'), 'signal theorem mismatch')
    return {'verified': True, 'levels': 4, 'resolved_q': 233,
            'full_tensor_oscillators': 232**3, 'largest_tridiagonal': 232,
            'uniform_window': ['19/20', '1'],
            'error_upper': last['graph_vs_compact_continuum_error_upper'][1],
            'resolved_signal_lower': last['resolved_graph_signal_lower']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('receipt', nargs='?', type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(verify(load(args.receipt)), sort_keys=True))


if __name__ == '__main__':
    main()
