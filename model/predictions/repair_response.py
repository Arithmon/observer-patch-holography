"""Native feedback response from passive clock/field records, without a fit.

A finite amplitude, occupation-controlled repair has an exact stationary
likelihood ratio. The identity covers every subsequent fixed record event,
including stopped clocks and timeouts. Parameters count native events.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as sp

from ..laboratory_clock import event_geometry as inherited
from ..observer import native
from ..observer.scalar_response.experiment import compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*inherited.SOURCE_PATHS,
    'model/predictions/repair_response.py', 'tools/strict_json.py'))))
PULSE = (0, 1)
# The impossible negative-current and coin-one positive-current labels are
# present in the readout alphabet and receive exact zero probability.
MARKS = ((0, -1, 0), *((1, coin, current) for coin in (0, 1) for current in (-1, 0, 1)))
READ_PORTS = (0, 1, 2)


def _integer(value, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('an original integer inside the declared range is required')
    return value


@lru_cache(maxsize=1)
def source():
    kernel = inherited.transfer.native_counts()
    if PULSE not in kernel.seams:
        raise ArithmeticError('the pulse seam is absent from the native carrier')
    return kernel


def pulse_weights(state):
    """Postpulse stationary density for each retained trigger/coin/current mark.

    These formulas use only the disclosed endpoint bits. No counterfactual
    future, hidden state search or fitted response parameter is consulted.
    """
    if (type(state) is not tuple or len(state) != 12 or
            any(type(x) is not int or x not in (0, 1) for x in state)):
        raise ValueError('twelve original binary occupations are required')
    a, b = (state[j] for j in PULSE)
    result = {mark: F(0) for mark in MARKS}
    result[0, -1, 0] = F(1-a)
    result[1, 0, 1] = F((1-a)*b, 2)
    result[1, 1, 0] = F(a, 2)
    result[1, 0, 0] = F(a*b, 2)
    if sum(result.values()) != 1+F(b-a, 2):
        raise ArithmeticError('the complete pulse density lost its exact identity')
    return result


def pulse_distribution(*, reverse_forecast=False, prior=None):
    """Full actual pulse instrument, or its passive-record density forecast."""
    states = source().states
    if prior is None:
        prior = (F(1, len(states)),)*len(states)
    if (type(prior) is not tuple or len(prior) != len(states)
            or any(type(p) is not F or p < 0 for p in prior) or sum(prior) != 1):
        raise ValueError('a complete normalized exact native prior is required')
    if type(reverse_forecast) is not bool:
        raise ValueError('a Boolean calculation mode is required')
    prior = tuple(inherited.exact_core._fraction(p) for p in prior)
    result = defaultdict(F)
    for state, weight in zip(states, prior, strict=True):
        if reverse_forecast:
            for mark, density in pulse_weights(state).items():
                result[mark, state] += weight*density
        elif state[PULSE[0]] == 0:
            result[(0, -1, 0), state] += weight
        else:
            for coin in (0, 1):
                after, (_, current) = native.native(state, PULSE, coin)
                mark = (1, coin, current)
                before_bits = tuple(state[j] for j in READ_PORTS)
                if pulse_readback(before_bits, mark) != tuple(after[j] for j in READ_PORTS):
                    raise ArithmeticError('the pulse record lost the disclosed initial fields')
                result[mark, after] += weight/2
    return {key: value for key, value in result.items() if value}


def pulse_readback(before_bits, mark):
    """Postpulse bits from the pre-pulse disclosure and retained pulse current."""
    if (type(before_bits) is not tuple or len(before_bits) != 3
            or any(type(x) is not int or x not in (0, 1) for x in before_bits)
            or type(mark) is not tuple or len(mark) != 3
            or any(type(x) is not int for x in mark) or mark not in MARKS):
        raise ValueError('three disclosed bits and a complete pulse mark are required')
    trigger, coin, current = mark
    if trigger != before_bits[0]:
        raise ValueError('the pulse trigger must read the disclosed port-zero bit')
    if not trigger:
        return before_bits
    # Checking the deterministic native endpoint prevents a fabricated mark
    # from acquiring observer-access status by a formal current subtraction.
    pair = before_bits[:2]+(0,)*10
    after, expected = native.native(pair, PULSE, coin)
    if expected[1] != current:
        raise ValueError('the retained pulse mark is incompatible with the disclosed bits')
    return after[0], after[1], before_bits[2]


def _reading(mark, initial, tape, terminal, status, clock_b):
    """Every field is reconstructed from disclosed bits and retained currents."""
    fields = list(initial)
    for edge, current in tape:
        for i, port in enumerate(READ_PORTS):
            fields[i] += current*(int(port == edge[1])-int(port == edge[0]))
    if tuple(fields) != tuple(terminal[j] for j in READ_PORTS):
        raise ArithmeticError('the visible current tape did not reconstruct the field')
    a_ticks = mark[0]+int(status == 'tick')
    b_ticks = mark[0]*int(clock_b in PULSE)+sum(clock_b in edge for edge, _ in tape)
    return mark, initial, tape, status, (a_ticks, b_ticks), tuple(fields)


def stopped_law(cap=2, *, clock_b=1, reverse_forecast=False, intervention=True):
    """Complete union tape until the first postpulse port-zero tick or timeout.

    The pulse counts as a tick iff the trigger is one. The independent
    attempt watchdog limits execution but its counter is not a clock readout.
    """
    cap = _integer(cap, 1, 2)
    clock_b = _integer(clock_b, 1, 2)
    if type(intervention) is not bool or type(reverse_forecast) is not bool:
        raise ValueError('Boolean intervention and forecast choices are required')
    if reverse_forecast and not intervention:
        raise ValueError('the forecast concerns the declared intervention')
    initial_law = (pulse_distribution(reverse_forecast=reverse_forecast) if intervention else
                   {((0, -1, 0), state): F(1, 66) for state in source().states})
    populations = {}
    for (mark, state), mass in initial_law.items():
        endpoint = tuple(state[j] for j in READ_PORTS)
        before = (endpoint[0]+mark[2], endpoint[1]-mark[2], endpoint[2])
        # In the intervention arm before is the originally disclosed field
        # record, recoverable from endpoint+mark. No extra field read occurs.
        initial = pulse_readback(before, mark) if intervention else endpoint
        key = mark, initial, (), state
        populations[key] = mass
    result = defaultdict(F)
    for step in range(cap):
        future = defaultdict(F)
        for (mark, initial, tape, state), mass in populations.items():
            for edge in source().seams:
                for coin in (0, 1):
                    terminal, (_, current) = native.native(state, edge, coin)
                    read = tape+((edge, current),) if any(j in edge for j in READ_PORTS) else tape
                    if 0 in edge:
                        result[_reading(mark, initial, read, terminal, 'tick', clock_b)] += mass/60
                    elif step+1 == cap:
                        result[_reading(mark, initial, read, terminal, 'timeout', clock_b)] += mass/60
                    else:
                        future[mark, initial, read, terminal] += mass/60
        populations = future
    if sum(result.values()) != 1:
        raise ArithmeticError('the complete stopped law is not normalized')
    return dict(result)


def density_response(attempts):
    """All twelve mean-field responses, from endpoint sums and passive covariance."""
    attempts = _integer(attempts, 0, 8)
    kernel = source()
    distribution = defaultdict(F)
    for (_, state), mass in pulse_distribution().items():
        distribution[state] += mass
    row = tuple(distribution[state] for state in kernel.states)
    for _ in range(attempts):
        row = inherited.transfer._row_times(row, kernel.counts)
    actual = tuple(sum((mass*(state[c]-F(1, 6)) for state, mass in zip(kernel.states, row)), F())
                   for c in range(12))
    laplacian = sp.zeros(12)
    for a, b in kernel.seams:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    transfer = (sp.eye(12)-laplacian/60)**attempts
    predicted = tuple(F(5, 66)*F(transfer[1, c]-transfer[0, c]) for c in range(12))
    if actual != predicted:
        raise ArithmeticError('the finite-amplitude response differs from passive covariance')
    return actual


def discounted_response(z=1):
    """Geometrically discounted response is the same native spatial Green dipole."""
    if type(z) not in (int, F) or z <= 0:
        raise ValueError('a positive exact discount parameter is required')
    z = inherited.exact_core._fraction(z)
    kernel = source()
    c = sp.Matrix(kernel.counts)
    mass = defaultdict(F)
    for (_, state), value in pulse_distribution().items():
        mass[state] += value
    delta = sp.Matrix(1, 66, [sp.Rational(mass[s]-F(1, 66)) for s in kernel.states])
    field = sp.Matrix(66, 12, lambda i, j: sp.Rational(kernel.states[i][j])-sp.Rational(1, 6))
    oracle = delta*((60+sp.Rational(z))*sp.eye(66)-c).inv(method='DM')*field
    laplacian = sp.zeros(12)
    for a, b in kernel.seams:
        laplacian[a, a] += 1; laplacian[b, b] += 1
        laplacian[a, b] -= 1; laplacian[b, a] -= 1
    dipole = sp.zeros(1, 12); dipole[0, 1], dipole[0, 0] = 1, -1
    predicted = sp.Rational(5, 66)*dipole*(laplacian+sp.Rational(z)*sp.eye(12)).inv()
    if oracle != predicted:
        raise ArithmeticError('the full configuration resolvent lost the spatial response')
    return tuple(F(value) for value in predicted)


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(item) for item in value]
    return value


def law_digest(law):
    rows = [[_json(key), str(value)] for key, value in sorted(law.items())]
    return hashlib.sha256(json.dumps(rows, separators=(',', ':')).encode()).hexdigest()


def run():
    actual, predicted = pulse_distribution(), pulse_distribution(reverse_forecast=True)
    if actual != predicted:
        raise ArithmeticError('native pulse and passive-record forecast disagree')
    rows = []
    for cap in (1, 2):
        for clock_b in (1, 2):
            law = stopped_law(cap, clock_b=clock_b)
            if law != stopped_law(cap, clock_b=clock_b, reverse_forecast=True):
                raise ArithmeticError('the complete clock/field response forecast failed')
            mark_mass = tuple(sum((v for key, v in law.items() if key[0] == mark), F()) for mark in MARKS)
            timeout = sum((v for key, v in law.items() if key[3] == 'timeout'), F())
            if timeout != F(5, 6)**cap:
                raise ArithmeticError('the watchdog tail differs from native incidence')
            signal = sum((v for key, v in law.items() if key[3] == 'tick'
                          and key[-1][1] == 1 and key[-2][1] >= 1), F())
            baseline = stopped_law(cap, clock_b=clock_b, intervention=False)
            def future_signal(key):
                return (key[3] == 'tick' and key[-1][1] == 1
                        and key[-2][1]-key[0][0]*int(clock_b in PULSE) >= 1)
            response = sum((v for key, v in law.items() if future_signal(key)), F())
            passive = sum((v for key, v in baseline.items() if future_signal(key)), F())
            covariance = sum((v*F(key[1][1]-key[1][0], 2) for key, v in baseline.items()
                              if future_signal(key)), F())
            if response-passive != covariance or covariance < 0:
                raise ArithmeticError('the held-out joint event violates its passive response forecast')
            rows.append({'attempt_cap': cap, 'second_clock_port': clock_b,
                'nonzero_complete_labels': len(law), 'law_sha256': law_digest(law),
                'pulse_mark_probabilities': mark_mass, 'timeout_probability': timeout,
                'joint_signal': 'completed port-zero tick, final port-one field one, second clock at least one',
                'joint_signal_probability': signal,
                'matched_continuation_signal': 'completed tick, final port-one field one, at least one subsequent B tick',
                'passive_signal_probability': passive, 'feedback_signal_probability': response,
                'signal_increase_predicted_from_passive_covariance': covariance,
                'all_current_zero_events_retained': True,
                'pulse_is_charged_to_each_incident_clock': True})
    chosen = next(i for i, state in enumerate(source().states) if state[0] == 1 and state[1] == 0)
    prior = tuple(F(int(i == chosen)) for i in range(66))
    nonstationary = pulse_distribution(prior=prior)
    invalid_forecast = pulse_distribution(prior=prior, reverse_forecast=True)
    if nonstationary == invalid_forecast:
        raise ArithmeticError('the nonstationary premise control failed')
    return _json({'schema': 'oph.native_repair_response.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'pulse': {'seam': PULSE, 'policy': 'repair iff disclosed initial port-zero occupation equals one',
            'complete_marks_trigger_coin_current': MARKS,
            'no_pulse_coin_sentinel': -1,
            'sum_density': '1+(n_1-n_0)/2',
            'passive_future_event_response': 'Delta Pr(F)=E[(n_1-n_0) 1_F]/2'},
        'stopped_joint_forecasts': rows,
        'mean_field_responses': {str(n): density_response(n) for n in (0, 1, 2, 4)},
        'discounted_spatial_response': {str(z): discounted_response(z) for z in (1, 2)},
        'controls': {'nonstationary_prior_rejects_stationary_weight_formula': True,
            'nonstationary_actual_total_mass': sum(nonstationary.values()),
            'nonstationary_invalid_forecast_total_mass': sum(invalid_forecast.values()),
            'nonstationary_invalid_forecast_is_not_renormalized': True,
            'future_law_includes_trigger_coin_current': True},
        'scope': {'stationary_uniform_fixed_number_preparation': True,
            'source_readouts_and_future_continuation_fixed_before_comparison': True,
            'response_parameter_fitted': False, 'pulse_has_measured_physical_duration': False,
            'unmarked_response_requires_continuation_independent_of_pulse_mark': True,
            'marked_response_may_use_a_separately_declared_continuation_for_each_mark': True,
            'arbitrary_continuation_theorem_requires_identical_continuation_in_both_arms': True,
            'matched_public_intervention_dataset_available': False,
            'new_discriminator_against_SM_GR_established': False,
            'Lorentzian_or_CFT_limit_established': False}, 'all_pass': True})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(result, strict_json.load(args.check))
        print('Native passive-record response and complete clock/field forecasts replayed.')
    else:
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
