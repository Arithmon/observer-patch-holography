"""Rigorous continuous-action pressure enclosing all recorded scalar times.

The existing observer records remain split-integrator records. This package
computes a separate exact-Taylor enclosure for continuous dynamics of the
same supplied 64-mode action. It does not infer a native or equilibrium EoS.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
RER = HERE.parents[2]
PARENT = 'code/source_scalar_execution/source_scalar_execution_receipt.json'
sys.path.insert(0, str(RER / 'code/source_scalar_execution'))
import scalar_execution_algebra as algebra
import source_scalar_execution as source
import verify_source_scalar_execution as parent_verifier

Q, parse, rb = algebra.Q, algebra.parse, algebra.rational_bounds
SCHEMA = 'oph.continuous_scalar_pressure.v1'
SCOPE = {'same_supplied_finite_action': True, 'all_22_model_times_retained': True,
         'continuous_pressure_rigorously_enclosed': True, 'native_repair_eos': False,
         'equilibrium_thermodynamics': False, 'new_observer_history': False,
         'recorded_split_values_identified_with_exact_continuous_values': False,
         'physical_clock_calibration': False, 'continuum_limit': False}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def tails():
    dq = 4 * F(3, 2)**163 * 614**81 / factorial(163)
    dv = 4 * F(3, 2)**162 * 614**81 / factorial(162)
    dp = 5*dv + dv*dv/2 + F(616, 3)*(7*dq+dq*dq/2)
    assert 0 < dq < dv < F(1, 10**28)
    return dq, dv, dp


def bound(value, radius=F()):
    return source.widen(rb(value), radius)


def operator_checks(masses, rows, velocity):
    dense = [dict(row) for row in rows]
    for i, mass in enumerate(masses):
        assert mass.sign() > 0
        total = Q()
        for j, coefficient in rows[i]:
            assert masses[i]*coefficient == masses[j]*dense[j].get(i, Q())
            if i != j:
                assert coefficient.sign() <= 0
            total += coefficient if coefficient.sign() >= 0 else -coefficient
        assert (614-total).sign() >= 0
        # The weighted gradient form is an edge-square sum plus nonnegative
        # Dirichlet boundary diagonals; this certifies A>=I.
        assert (sum((c for j, c in rows[i]), Q())-1).sign() >= 0
    assert (16-source.inner(masses, velocity, velocity)).sign() >= 0
    assert (Q(F(3, 2))-21*source.TAU).sign() > 0


def observable(masses, rows, q, v):
    aq = source.action(rows, q)
    K = source.inner(masses, v, v)/2
    G = source.inner(masses, q, [a-x for a, x in zip(aq, q)])/2
    U = source.inner(masses, q, q)/2
    return K, G, U, K-G/3-U


def build():
    packet = parent_verifier.load(RER/PARENT)
    replay = parent_verifier.verify(packet, root=RER)
    masses = [parse(x) for x in packet['mass_Qphi']]
    rows = [[(j, parse(c)) for j, c in row] for row in packet['action_rows_Qphi']]
    velocity = [parse(x) for x in packet['initial_velocity_Qphi']]
    layers = [[parse(x) for x in layer] for layer in packet['traces']['ascending_intervention']['layers']]
    operator_checks(masses, rows, velocity)
    powers = [velocity]
    for _ in range(80):
        powers.append(source.action(rows, powers[-1]))
    dq, dv, dp = tails()
    H0 = source.inner(masses, velocity, velocity)/2
    h0lo = F(rb(H0)[0])
    assert h0lo > 0
    samples = []
    for step in range(22):
        q, v = [Q() for _ in range(64)], [Q() for _ in range(64)]
        for k, vector in enumerate(powers):
            cq = source.TAU*F((-1)**k*step**(2*k+1), 245**k*factorial(2*k+1))
            cv = F((-1)**k*step**(2*k), 245**k*factorial(2*k))
            for i, x in enumerate(vector):
                q[i] += cq*x
                v[i] += cv*x
        assert (49-source.inner(masses, q, q)).sign() >= 0
        assert (25-source.inner(masses, v, v)).sign() >= 0
        K, G, U, p = observable(masses, rows, q, v)
        old, split_q = layers[step:step+2]
        split_aq = source.action(rows, split_q)
        split_v = [(x-y)/source.TAU-source.TAU*a/2 for x, y, a in zip(split_q, old, split_aq)]
        _, _, _, split_p = observable(masses, rows, split_q, split_v)
        # At t=0 both initial Taylor polynomials are exact, so no tail is needed.
        radius = F() if step == 0 else dp
        signed_error = split_p-p
        samples.append({
            'step': step, 'model_time_Qphi': (step*source.TAU).encode(),
            'polynomial_q_Qphi': [x.encode() for x in q],
            'polynomial_v_Qphi': [x.encode() for x in v],
            'polynomial_K_Qphi': K.encode(), 'polynomial_G_Qphi': G.encode(),
            'polynomial_U_Qphi': U.encode(), 'polynomial_pressure_Qphi': p.encode(),
            'continuous_pressure_interval': bound(p, radius),
            'continuous_w_interval': bound(p/H0, radius/h0lo),
            'split_pressure_Qphi': split_p.encode(),
            'split_minus_continuous_pressure_interval': bound(signed_error, radius),
            'pressure_tail_bound': str(radius),
        })
    return {
        'schema': SCHEMA, 'scope': SCOPE,
        'assumptions': 'supplied q5 Dirichlet scalar action, preparation, model time, V=1 and fixed-original-q,p metric work',
        'polynomial': {'q_degree': 161, 'v_degree': 160, 'terms': 81,
                       'all_operator_modes_retained': 64},
        'bounds': {'operator_lower': '1', 'operator_upper': '614', 'gradient_operator_upper': '613',
                   'initial_velocity_norm_upper': '4', 'time_upper': '3/2',
                   'q_polynomial_norm_upper': '7', 'v_polynomial_norm_upper': '5',
                   'q_tail': str(dq), 'v_tail': str(dv), 'pressure_tail': str(dp),
                   'pressure_tail_formula': '5*dv+dv^2/2+(616/3)*(7*dq+dq^2/2)'},
        'continuous_energy_Qphi': H0.encode(), 'continuous_energy_interval': rb(H0),
        'baseline_control': {'continuous_q_v_and_pressure': 'identically zero',
                             'w': None, 'all_22_times_retained_in_parent': True},
        'schedule_control': 'both independently replayed intervention schedules give identical layer states',
        'samples': samples, 'parent_replay': replay,
        'source_sha256': {p: digest(RER/p) for p in (PARENT,
            'code/source_scalar_execution/source_scalar_execution.py',
            'code/source_scalar_execution/scalar_execution_algebra.py',
            'code/source_scalar_execution/verify_source_scalar_execution.py')},
        'producer_sha256': digest(Path(__file__)),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE/'continuous_pressure_receipt.json')
    args = parser.parse_args()
    receipt = build()
    args.out.write_text(json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(args.out), 'samples': len(receipt['samples']),
                      'final_pressure_interval': receipt['samples'][-1]['continuous_pressure_interval'],
                      'final_w_interval': receipt['samples'][-1]['continuous_w_interval']}, indent=2))
