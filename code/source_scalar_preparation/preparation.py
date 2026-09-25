"""Exact finite local-force preparation bounds on the pinned scalar action.

This certifies conditional quantum channels and errors, never sampled outcomes.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / 'source_scalar_execution'))
from scalar_execution_algebra import Q, parse, rational_bounds

PARENTS = {
    'code/source_scalar_execution/source_scalar_execution_receipt.json':
        '6981fcd4e13fbc33f371224228341a995f34b4733aa70c2ec0b3f6ee5975c7af',
    'code/source_scalar_finite_instrument/finite_instrument_receipt.json':
        'a9e922b62896419bb96ac24d9023c24f1c7f4f0cc40ca84c601ece12577a09ad',
}
OUTPUT = HERE / 'preparation_receipt.json'
SCALE = 10**15
DELTA_MIN, DELTA_MAX = F(1, 10000), F(1, 500)
VACUUM_ERROR = FORCE_L1_ERROR = F(1, 10**6)
FILES = ('preparation.py', 'verify_preparation.py', 'test_preparation.py', 'README.md')


def raw(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode('ascii')


def sqrt_upper(value):
    if value.sign() < 0:
        raise ValueError('negative squared norm')
    lo, hi = 0, SCALE
    while (Q(F(hi, SCALE)**2) - value).sign() < 0:
        hi *= 2
    while hi-lo > 1:
        mid = (hi+lo)//2
        if (Q(F(mid, SCALE)**2)-value).sign() < 0:
            lo = mid
        else:
            hi = mid
    return F(hi, SCALE)


def produce():
    parents = {}
    for path, digest in PARENTS.items():
        data = (ROOT/path).read_bytes()
        if sha256(data).hexdigest() != digest:
            raise ValueError('immutable parent ' + path)
        parents[path] = json.loads(data)
    scalar, instrument = parents.values()
    m = list(map(parse, scalar['mass_Qphi']))
    v = list(map(parse, scalar['initial_velocity_Qphi']))
    a = [[(j, parse(x)) for j, x in row] for row in scalar['action_rows_Qphi']]
    av = [sum((x*v[j] for j, x in row), Q()) for row in a]
    v0 = sum((mass*x*x for mass, x in zip(m,v)), Q())
    v2 = sum((mass*x*x for mass, x in zip(m,av)), Q())
    norm_av = sqrt_upper(v2)
    pulse = DELTA_MAX**2*norm_av/24
    baseline = VACUUM_ERROR+FORCE_L1_ERROR
    intervention = baseline+pulse
    support = [i for i,x in enumerate(v) if x.sign()]
    inherited = F(instrument['parameters']['preparation_trace_distance_max'])
    first = instrument['schedule'][0]
    start = first['operations'][0]['start_offset']
    first_start = F(first['center_time_interval'][0])+F(start[0])*F(instrument['parameters']['duration_max'])+F(start[1])
    operations = [
        {'site': i, 'velocity_Qphi': v[i].encode(),
         'force_coefficient_m_v_Qphi': (m[i]*v[i]).encode(),
         'start_over_delta': '-1/2', 'stop_over_delta': '1/2',
         'control_reads': ['mass', 'velocity', 'duration'],
         'on_write': 'force=-m*v/delta multiplying q in H',
         'off_write': 'force=0'} for i in support]
    return {
        'schema': 'oph.scalar-finite-local-preparation.v1',
        'parents': PARENTS,
        'scope': {
            'same_original64_mode_action': True,
            'free_field_active_throughout_preparation': True,
            'local_force_support_equals_original16_site_intervention': True,
            'original_vacuum_at_pulse_start_supplied': True,
            'force_control_law_and_error_budget_supplied': True,
            'virtual_center_state_is_actual_midpulse_state': False,
            'trace_comparison_valid_after_pulse_end_only': True,
            'all_bounded_later_instruments_share_trace_error_bound': True,
            'source_selected_quantization_or_vacuum': False,
            'physical_clock_or_quantum_outcome_custody': False,
            'actual_quantum_hardware_execution': False,
            'noisy_energy_bound_from_trace_distance': False,
        },
        'law': {
            'equation': 'q_ddot + A q = v/delta + e(t) during [-delta/2,delta/2]',
            'Hamiltonian': 'H0 - sum_i m_i*v_i*q_i/delta',
            'reference': 'free evolution from D(momentum=M*v)|vacuum> at virtual t=0',
            'effective_velocity': 'sinc(delta*sqrt(A)/2)*v',
            'ideal_trace_bound': 'delta^2*||A v||_M/(24*sqrt(2)) <= delta^2*||A v||_M/24',
            'force_error_contract': 'integral ||e(t)||_M dt <= force_L1_error',
            'force_trace_bound': 'force_L1_error/sqrt(2) <= force_L1_error',
            'initial_error_contract': 'trace distance from original vacuum at pulse start <= vacuum_error',
            'baseline_force': '0, with same declared residual force error contract',
        },
        'parameters': {'delta_min':str(DELTA_MIN),'delta_max':str(DELTA_MAX),
                       'vacuum_error':str(VACUUM_ERROR),'force_L1_error':str(FORCE_L1_ERROR)},
        'moments': {'v_mass_squared_Qphi':v0.encode(),'Av_mass_squared_Qphi':v2.encode(),
                    'Av_norm_upper':str(norm_av),
                    'ideal_added_energy_upper':rational_bounds(v0/2,SCALE)[1]},
        'bounds': {'pulse_trace_error_upper':str(pulse),
                   'baseline_trace_error_upper':str(baseline),
                   'intervention_trace_error_upper':str(intervention),
                   'inherited_per_branch_preparation_budget':str(inherited),
                   'both_branches_fit_inherited_budget':intervention<=inherited,
                   'first_pointer_preparation_start_lower':str(first_start),
                   'preparation_end_to_first_pointer_gap_lower':str(first_start-DELTA_MAX/2)},
        'instrument_composition': {
            'all21_unconditional_marginals_retained': True,
            'independent_shots_assumed': False,
            'resolved_steps':instrument['summary']['resolved_steps'],
            'parent_error_budgets_remain_valid':intervention<=inherited and baseline<=inherited,
            'step16_response_lower':instrument['rows'][15]['noisy_paired_response_abs_lower'],
            'step16_total_error_upper':instrument['rows'][15]['total_paired_response_error_upper'],
        },
        'control_schedule':operations,
        'resource_model': {'parallel_force_ports':len(support),
                           'coefficient_reads':3*len(support),
                           'on_off_writes':2*len(support),
                           'elapsed_preparation_time':'delta, not 16*delta',
                           'coefficient_bit_precision':'exact Q(phi) law plus declared integrated force error; no hardware precision claim',
                           'phase_and_outcome_records':'none synthesized'},
        'source_pins':{str((HERE/f).relative_to(ROOT)).replace('\\','/'):sha256((HERE/f).read_bytes()).hexdigest() for f in FILES},
    }


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');ap.add_argument('--check',action='store_true');args=ap.parse_args()
    value=produce();data=raw(value)
    if args.write: OUTPUT.write_bytes(data)
    if args.check and OUTPUT.read_bytes()!=data: raise SystemExit('preparation producer parity failure')
    print(json.dumps(value['bounds'],indent=2))
