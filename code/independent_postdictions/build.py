"""Build compact historical arithmetic; no natural held-out comparison."""
import itertools
import json
from pathlib import Path

from .replay import context, alpha_replay, measured_pixel_diagnostic, tau_replay
from .comparison import compare, reference_log_likelihood_ratios


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def calculations():
    c = context()
    fmt = lambda x: c.nstr(x, 40)
    inputs = json.loads((HERE/'inputs.json').read_text())
    reference = c.mpf(inputs['codata']['alpha_inverse']['value'])
    sigma_alpha = c.mpf(inputs['codata']['alpha_inverse']['sigma'])
    alpha = alpha_replay()
    diagnostics = []
    for row in alpha['rows']:
        d = reference-c.mpf(row['alpha_inverse'])
        diagnostics.append(dict(mode=row['mode'], signed_reference_minus_output=fmt(d),
                                absolute_ppm=fmt(abs(d)/reference*1000000),
                                absolute_in_reference_sigma_only=fmt(abs(d)/sigma_alpha),
                                minimum_symmetric_remainder_for_two_sigma_overlap=fmt(max(abs(d)-2*sigma_alpha, 0))))
    e, m = [c.mpf(inputs['codata'][k]['value']) for k in ('electron', 'muon')]
    se, sm = [c.mpf(inputs['codata'][k]['sigma']) for k in ('electron', 'muon')]
    tau = tau_replay(str(e), str(m))
    corners = [tau_replay(str(e+i*se), str(m+j*sm)) for i, j in itertools.product((-1, 1), repeat=2)]
    # tau_replay intentionally uses its own fixed context; use a macroscopic
    # symmetric difference then Richardson, well above its precision floor.
    def derivative(which):
        def central(h):
            args = [e, m]
            args[which] += h
            high = tau_replay(*map(str, args))
            args[which] -= 2*h
            low = tau_replay(*map(str, args))
            return (high-low)/(2*h)
        h = c.mpf('1e-12')
        return (4*central(h/2)-central(h))/3
    derivative_e, derivative_m = derivative(0), derivative(1)
    sigma_max = abs(derivative_e)*se+abs(derivative_m)*sm
    measured, measured_sigma = c.mpf(inputs['tau']['value']), c.mpf(inputs['tau']['sigma'])
    distance = abs(tau-measured)
    Q = lambda es, ms, ts: (es+ms+ts)/(c.sqrt(es)+c.sqrt(ms)+c.sqrt(ts))**2
    q_corners = [Q(e+i*se, m+j*sm, measured+k*measured_sigma)
                 for i, j, k in itertools.product((-1, 1), repeat=3)]
    tau_data = dict(central_mev=fmt(tau), corners_mev=list(map(fmt, corners)),
                    excluded_small_root_mev=fmt(tau_replay(str(e), str(m), ordered=False)),
                    measured_Q=fmt(Q(e, m, measured)), measured_Q_corner_range=list(map(fmt, (min(q_corners), max(q_corners)))),
                    outward_mev=[fmt(c.floor(min(corners)*10**6)/10**6), fmt(c.ceil(max(corners)*10**6)/10**6)],
                    jacobian=list(map(fmt, (derivative_e, derivative_m))),
                    first_order_sigma_max_mev=fmt(sigma_max),
                    residual_in_sigma_with_arbitrary_correlations=[fmt(distance/(measured_sigma+sigma_max)),
                                                                  fmt(distance/(measured_sigma-sigma_max))],
                    covariance_convention='first-order Cauchy bounds; no independence, Gaussian confidence level or theory-error distribution assumed')
    historical = compare({k: inputs['tau'][k] for k in ('value', 'sigma', 'unit')})
    ratios = reference_log_likelihood_ratios(inputs['tau']['value'], inputs['tau']['sigma'], fmt(tau))
    trunk = json.loads((ROOT/'code/P_derivation/runtime/p_closure_trunk_current.json').read_text())['fixed_point_candidate']
    defect = (c.mpf(trunk['P'])-(1+c.sqrt(5))/2)*c.mpf(trunk['alpha_inv'])/c.sqrt(c.pi)-1
    asym = c.mpf(alpha['rows'][2]['alpha_inverse'])
    return dict(alpha=alpha, comparison_diagnostics=diagnostics,
                measured_pixel_diagnostic=measured_pixel_diagnostic(str(reference)),
                approximate_trunk=dict(printed_inverse=trunk['alpha_inv'], relative_pair_defect=fmt(defect),
                                       printed_minus_converged_asymptotic=fmt(c.mpf(trunk['alpha_inv'])-asym)),
                tau=tau_data, historical_tau=historical, normal_error_reference_diagnostic=ratios,
                conclusions=dict(reproduction='AGREEMENT_UNDER_DECLARED_CONVENTIONS', natural_comparison='NOT_READY',
                                 selected_target='FZ-10', new_natural_outcomes=0, global_chance_probability=None,
                                 corners_are_joint_68_percent_interval=False,
                                 tau_independently_tests_P=False, oph_versus_koide_log_likelihood_ratio='0'))


if __name__ == '__main__':
    from .verify import pins
    packet = dict(schema='oph-independent-postdictions-v1', sources=pins(), calculations=calculations())
    (HERE/'receipt.json').write_text(json.dumps(packet, indent=2, sort_keys=True)+'\n', encoding='ascii', newline='\n')
    print('Built independent historical replay and NOT_READY decision')
