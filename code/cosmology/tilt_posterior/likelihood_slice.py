"""Conditional ACT DR6 likelihood slice; no full cosmological fit or p-value.

Online reproduction uses pinned measured data, official ACT window/covariance
code and CAMB. Offline verification recomputes both retained quadratic scores
from the complete observed covariance and predicted bandpowers, without CAMB.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import inspect
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
RECEIPT = HERE / 'likelihood_slice.json'
EDGE_NS = 463690317905960407 / 480000000000000000
ARCHIVE_SHA = '3f057c2569211ada03759530b74848b322edc3d68d66b6b8c2db0679547dbbd8'
FITS_SHA = 'c887dc9178d81f5e5e0ce76eca0cd3c9f089056630ab78a1cc4bb28ff8751c29'
ACT_COMMIT = '627aeafb88ae5ad1aa66b406bea2d65cfa66a27d'
BACKGROUND = dict(ombh2=.0226, omch2=.118, cosmomc_theta=.010417, tau=.0544)
ACCURACY = dict(lmax=9000, lens_potential_accuracy=8, min_l_logl_sampling=6000)
BOUNDS = dict(ns=[.9, 1.1], logA=[2.6, 3.5], P_act=[.9, 1.1])
ELL_CUTS = dict(TT=[600, 6500], TE=[600, 6500], EE=[600, 6500])
PHYSICS_PARAMETERS = dict(mnu=.06, nnu=3.044, num_massive_neutrinos=1,
    neutrino_hierarchy='degenerate', omk=0., TCMB=2.7255, Alens=1.,
    w=-1., wa=0., cs2=1., dark_energy_model='fluid')
FIXED_PHYSICS = dict(parameters=PHYSICS_PARAMETERS, helium='CAMB BBN consistency',
    recombination='Recfast', reionization='TanhReionization',
    nonlinear='NonLinear_lens', nonlinear_model='Halofit:mead2020',
    tensor_to_scalar_ratio=0., scalar_running=0., scalar_running_running=0.)
RUNTIME = dict(camb='1.6.6', sacc='0.16', numpy='2.2.6', scipy='1.15.3')
# Independently projected from the pinned measured FITS using SACC0.16.
OBSERVATION_PROJECTION_SHA256 = 'f683ebd198cd4ec0ec5b64c681cc0d824ff8ef9688b9ee939d438c9cd73b9d29'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def observation_projection(like):
    """Retain the complete data/covariance and bind their ordered windows."""
    menu = []
    for row in like.spec_meta:
        window = row['window']
        menu.append(dict(spectrum=row['pol'], indices=row['idx'].tolist(),
            ell=row['ell'].tolist(), tracers=[row['tracer1'], row['tracer2']],
            window_shape=list(window.weight.shape),
            window_sha256=canonical_digest(dict(ell=window.values.tolist(), weights=window.weight.tolist()))))
    return dict(archive_sha256=ARCHIVE_SHA, fits_sha256=FITS_SHA, official_commit=ACT_COMMIT,
        official_source_sha256=ACT_SOURCE_SHA, spectra=[row['pol'] for row in like.spec_meta],
        bins_per_spectrum=[len(row['idx']) for row in like.spec_meta], ell_cuts=like.ell_cuts,
        bandpower_metadata=menu, data=like.data_vec.tolist(),
        covariance_lower=like.covmat[np.tril_indices(len(like.data_vec))].tolist())


def finite_number(value):
    return type(value) in (int, float) and math.isfinite(value)


def validate_point(ns, logA, polarization):
    values = dict(ns=ns, logA=logA, P_act=polarization)
    if any(not finite_number(value) or not BOUNDS[key][0] <= value <= BOUNDS[key][1]
           for key, value in values.items()):
        raise ValueError('finite cosmological and calibration parameters inside the declared bounds required')


def chi_square(data, covariance, theory):
    y, c, t = np.asarray(data, dtype=float), np.asarray(covariance, dtype=float), np.asarray(theory, dtype=float)
    if y.ndim != 1 or t.shape != y.shape or c.shape != (len(y), len(y)):
        raise ValueError('complete aligned data, theory and covariance required')
    if not all(np.isfinite(x).all() for x in (y, c, t)) or not np.allclose(c, c.T, atol=0, rtol=1e-13):
        raise ValueError('finite symmetric covariance required')
    # Full Cholesky whitening preserves all TT/TE/EE correlations.
    residual = np.linalg.solve(np.linalg.cholesky(c), y-t)
    return float(residual @ residual)


def unpack_covariance(lower, n):
    if len(lower) != n*(n+1)//2:
        raise ValueError('complete covariance triangle required')
    c = np.zeros((n, n)); idx = np.tril_indices(n); c[idx] = lower
    c += np.tril(c, -1).T
    return c


def projected_theory(like, cl, polarization):
    result = np.zeros_like(like.data_vec)
    for m in like.spec_meta:
        # ACT convention: D_l in microkelvin squared, TE/P and EE/P^2.
        result[m['idx']] = m['window'].weight.T @ (
            cl[m['pol']][m['window'].values] / polarization**m['pol'].count('e'))
    return result


def check(packet):
    """Verify custody, scope and quadratic scores, not the CAMB transfer map.

    The optional online producer independently recalculates that map at
    both retained points. Neither route certifies a global optimum.
    """
    expected = {'schema', 'scope', 'source_sha256', 'background', 'accuracy', 'calibration_A_act',
        'pivot_Mpc_inverse', 'parameter_bounds', 'runtime', 'fixed_physics', 'observations',
        'fits', 'delta_chi2', 'p_value', 'FZ15_verdict', 'full_cosmological_profile',
        'native_primordial_generator_identified'}
    if type(packet) is not dict or set(packet) != expected:
        raise ValueError('complete declared conditional-likelihood packet required')
    if packet['schema'] != 'oph.act_dr6.conditional_likelihood_slice.v1':
        raise ValueError('conditional likelihood schema differs')
    if packet['source_sha256'] != digest(__file__):
        raise ValueError('likelihood analysis source differs')
    if packet['observations']['fits_sha256'] != FITS_SHA or packet['observations']['archive_sha256'] != ARCHIVE_SHA:
        raise ValueError('measured ACT provenance differs')
    if (canonical_digest(packet['background']) != canonical_digest(BACKGROUND)
            or canonical_digest(packet['accuracy']) != canonical_digest(ACCURACY)):
        raise ValueError('conditional background or accuracy differs')
    if (not finite_number(packet['calibration_A_act']) or packet['calibration_A_act'] != 1.
            or not finite_number(packet['pivot_Mpc_inverse']) or packet['pivot_Mpc_inverse'] != .05
            or canonical_digest(packet['parameter_bounds']) != canonical_digest(BOUNDS)
            or canonical_digest(packet['fixed_physics']) != canonical_digest(FIXED_PHYSICS)):
        raise ValueError('fixed physical, pivot or calibration conditions differ')
    runtime = packet['runtime']
    if runtime != RUNTIME:
        raise ValueError('declared runtime versions differ')
    if packet['scope'] != 'fixed_background_ACT_only_likelihood_slice' or packet['p_value'] is not None or packet['FZ15_verdict'] is not None:
        raise ValueError('this diagnostic is neither a full fit nor a significance verdict')
    if packet['full_cosmological_profile'] is not False or packet['native_primordial_generator_identified'] is not False:
        raise ValueError('the conditional diagnostic cannot promote physical or profile claims')
    obs = packet['observations']; n = len(obs['data'])
    if n != 135 or obs['spectra'] != ['tt', 'te', 'ee'] or obs['bins_per_spectrum'] != [45, 45, 45]:
        raise ValueError('complete measured TT/TE/EE menu required')
    if canonical_digest(obs) != OBSERVATION_PROJECTION_SHA256:
        raise ValueError('observed data, covariance, ordered windows or provenance projection differs')
    cov = unpack_covariance(obs['covariance_lower'], n)
    if set(packet['fits']) != {'fixed_tilt', 'free_tilt'}:
        raise ValueError('exactly the fixed and free conditional tilt fits are required')
    scores = {}
    for key, fit in packet['fits'].items():
        if set(fit) != {'ns', 'logA', 'P_act', 'A_act', 'chi2', 'official_chi2',
                'predicted_bandpowers', 'optimizer_success', 'optimizer_reported_chi2',
                'optimizer_origin', 'full_recalculation_delta_chi2'}:
            raise ValueError('complete declared conditional fit required')
        validate_point(fit['ns'], fit['logA'], fit['P_act'])
        if not finite_number(fit['A_act']) or fit['A_act'] != 1.:
            raise ValueError('ACT absolute calibration must stay fixed')
        numeric = ['chi2', 'official_chi2', 'full_recalculation_delta_chi2', 'optimizer_reported_chi2']
        if any(not finite_number(fit[name]) or fit[name] < 0 for name in numeric):
            raise ValueError('finite nonnegative likelihood and recalculation scores required')
        if (fit['optimizer_success'] is not True or fit['full_recalculation_delta_chi2'] > .002
                or fit['optimizer_origin'] not in ('recomputed', 'supplied_points_recalculated')):
            raise ValueError('optimizer or independent CAMB recalculation failed')
        if key == 'fixed_tilt' and abs(fit['ns']-EDGE_NS) > 1e-15:
            raise ValueError('fixed tilt was refitted')
        score = chi_square(obs['data'], cov, fit['predicted_bandpowers'])
        if abs(score-fit['chi2']) > 1e-8 or abs(score-fit['official_chi2']) > 1e-8:
            raise ValueError('complete-covariance score differs')
        if abs(score-fit['optimizer_reported_chi2']) > .002:
            raise ValueError('reported optimized point differs from its full recalculation')
        scores[key] = score
    delta = scores['fixed_tilt']-scores['free_tilt']
    if not finite_number(packet['delta_chi2']) or delta < 0 or abs(delta-packet['delta_chi2']) > 1e-8:
        raise ValueError('nested conditional slice comparison differs')
    return scores, delta


def reproduce(data_root, official_code, fit_file=None):
    """External libraries are optional for the ordinary offline check."""
    import camb
    import scipy
    import sacc
    from scipy.optimize import minimize
    runtime = dict(camb=camb.__version__, sacc=sacc.__version__, numpy=np.__version__, scipy=scipy.__version__)
    if runtime != RUNTIME:
        raise ValueError('reproduction requires the declared CAMB/SACC/NumPy/SciPy versions')
    fits = data_root / 'data/ACTDR6CMBonly/v1.0/dr6_data_cmbonly.fits'
    if digest(fits) != FITS_SHA:
        raise ValueError('the measured ACT FITS hash differs; bundled mock data are inadmissible')
    code = official_code / 'act_dr6_cmbonly/act_dr6_cmbonly.py'
    if digest(code) != ACT_SOURCE_SHA:
        raise ValueError('official ACT implementation differs from the pinned revision')
    sys.path.insert(0, str(official_code))
    from act_dr6_cmbonly import ACTDR6CMBonly
    if Path(inspect.getfile(ACTDR6CMBonly)).resolve() != code.resolve():
        raise ValueError('the imported ACT implementation is not the pinned supplied source')
    like = object.__new__(ACTDR6CMBonly)
    like.packages_path = str(data_root); like.log = logging.getLogger('ACT')
    like.initialize()
    observations = observation_projection(like)
    if canonical_digest(observations) != OBSERVATION_PROJECTION_SHA256:
        raise ValueError('online observed SACC projection differs from the independent pin')
    params = camb.set_params(**BACKGROUND, **ACCURACY, **PHYSICS_PARAMETERS, As=2.1e-9, ns=EDGE_NS)
    actual = dict(parameters=PHYSICS_PARAMETERS, helium='CAMB BBN consistency',
        recombination=type(params.Recomb).__name__, reionization=type(params.Reion).__name__,
        nonlinear=str(params.NonLinear),
        nonlinear_model=type(params.NonLinearModel).__name__+':'+params.NonLinearModel.halofit_version,
        tensor_to_scalar_ratio=float(params.InitPower.r), scalar_running=float(params.InitPower.nrun),
        scalar_running_running=float(params.InitPower.nrunrun))
    if actual != FIXED_PHYSICS:
        raise ValueError('resolved CAMB physical defaults differ from the declared model')
    # Required: changing ns/As must recompute nonlinear-lensing corrections.
    transfers = camb.get_transfer_functions(params, only_time_sources=True)
    def spectra(ns, logA, full=False):
        power = camb.initialpower.InitialPowerLaw(As=1e-10*np.exp(logA), ns=ns, pivot_scalar=.05)
        if full:
            fresh = camb.set_params(**BACKGROUND, **ACCURACY, **PHYSICS_PARAMETERS, As=power.As, ns=ns)
            result = camb.get_results(fresh)
        else:
            transfers.power_spectra_from_transfer(power); result = transfers
        raw = result.get_cmb_power_spectra(CMB_unit='muK', lmax=9000)['total']
        return {key: raw[:, index] for key, index in [('tt', 0), ('ee', 1), ('te', 3)]}
    def objective(x, fixed):
        ns, logA, pol = (fixed, *x) if fixed is not None else x
        if not (.9 <= ns <= 1.1 and 2.6 <= logA <= 3.5 and .9 <= pol <= 1.1):
            return 1e10
        return like.chi_square(spectra(ns, logA), 1., pol)
    fits_out = {}
    supplied = json.loads(fit_file.read_text()) if fit_file else None
    for name, fixed, start in [('fixed_tilt', EDGE_NS, [np.log(21), 1.]),
                               ('free_tilt', None, [.974, np.log(21), 1.])]:
        if supplied is None:
            step = [.01, .005] if fixed is not None else [.006, .01, .005]
            simplex = np.array([start]+[np.array(start)+np.eye(len(start))[i]*step[i] for i in range(len(start))])
            fit = minimize(objective, start, args=(fixed,), method='Nelder-Mead',
                           options=dict(maxiter=180, xatol=2e-6, fatol=1e-5, initial_simplex=simplex))
            x, success, reported = [float(v) for v in fit.x], bool(fit.success), float(fit.fun)
        else:
            row = supplied['edge' if fixed is not None else 'free']
            x, success, reported = row['x'], row['success'], row['chi2']
        if success is not True: raise ArithmeticError('conditional optimizer did not converge')
        if type(x) is not list or len(x) != (2 if fixed is not None else 3) or any(not finite_number(v) for v in x):
            raise ValueError('the declared optimized parameter vector is required')
        ns, logA, pol = map(float, (fixed, *x) if fixed is not None else x)
        validate_point(ns, logA, pol)
        if not finite_number(reported) or reported < 0:
            raise ValueError('a finite nonnegative optimizer score is required')
        cached, full = spectra(ns, logA), spectra(ns, logA, full=True)
        theory = projected_theory(like, full, pol)
        score = chi_square(like.data_vec, like.covmat, theory)
        fits_out[name] = dict(ns=float(ns), logA=float(logA), P_act=float(pol), A_act=1.,
            chi2=score, official_chi2=float(like.chi_square(full, 1., pol)),
            predicted_bandpowers=theory.tolist(), optimizer_success=success,
            optimizer_reported_chi2=reported,
            optimizer_origin='supplied_points_recalculated' if supplied is not None else 'recomputed',
            full_recalculation_delta_chi2=float(abs(score-like.chi_square(cached, 1., pol))))
    packet = dict(schema='oph.act_dr6.conditional_likelihood_slice.v1',
        scope='fixed_background_ACT_only_likelihood_slice', source_sha256=digest(__file__),
        background=BACKGROUND, accuracy=ACCURACY, calibration_A_act=1.,
        pivot_Mpc_inverse=.05, parameter_bounds=BOUNDS,
        runtime=runtime,
        fixed_physics=actual, observations=observations,
        fits=fits_out, delta_chi2=fits_out['fixed_tilt']['chi2']-fits_out['free_tilt']['chi2'],
        p_value=None, FZ15_verdict=None, full_cosmological_profile=False,
        native_primordial_generator_identified=False)
    check(packet)
    return packet


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--official-code', type=Path)
    parser.add_argument('--fit-file', type=Path, help='optional already optimized points, independently recalculated')
    parser.add_argument('--output', type=Path, default=RECEIPT)
    args = parser.parse_args()
    if args.check:
        scores, delta = check(json.loads(args.output.read_text()))
        print('ACT conditional slice:', scores, 'delta_chi2=', delta)
    else:
        if args.data_root is None or args.official_code is None:
            parser.error('online reproduction needs --data-root and --official-code')
        packet = reproduce(args.data_root, args.official_code, args.fit_file)
        args.output.write_text(json.dumps(packet, indent=2)+'\n', encoding='utf-8', newline='\n')
        print('wrote', args.output)


ACT_SOURCE_SHA = 'b1614ff2f69697376212eec857bcde7a94a91d18bd8ab476d118c66fa6d6c964'

if __name__ == '__main__':
    main()
