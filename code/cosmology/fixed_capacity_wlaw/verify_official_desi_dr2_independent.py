#!/usr/bin/env python3
"""Independent centered-Decimal replay of the retained official DESI summaries.

No producer import. Input chain identities are fixed here independently of the
receipt, so changing both a receipt pin and its input cannot manufacture a pass.
This checks numerical accounting, not posterior convergence or model validity.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RECEIPT = HERE / 'runtime/official_desi_dr2_fz13_retrospective.json'
# Published v1.0 chain pins, unchanged by the accounting repair.
PINS = {'DESI_DR2_BAO+CMB': [('cmb_chain.1.txt', 'c228de7bbaec19ddb22eec25c3dd7c40ef218976c020a7b55dd1c78dc3a638c5'),
                      ('cmb_chain.2.txt', '01bb30f43b3207d8575cc16354159fceff2ad4deffacc87964b0aef7f8e8ee44'),
                      ('cmb_chain.3.txt', 'db12623c8c69c03ad219b28bbc517416fc941cbb1e717054344f30b4e43a4adc'),
                      ('cmb_chain.4.txt',
                       '339312d9b7d3027147c38433fb4335cb9a12585776c4dfffc592c7c43af9a5ff')],
 'DESI_DR2_BAO+CMB+DESY5': [('desy5_chain.1.txt',
                             '8c783ebf283a205b7f569ce36a6694a32646ced5e28fd3b4683742733f6165e0'),
                            ('desy5_chain.2.txt',
                             'cd4f2ff3a66aa92aceecd47c8452520c2de09d887f231fc0df033cd68597a886'),
                            ('desy5_chain.3.txt',
                             'f7f8dbf28ff23d371e0b987b938d321adb920f1a09e97f746b3c0bb2d6247a01'),
                            ('desy5_chain.4.txt',
                             '4a3607867d34890832f431c4d61947e5572a72ff1a407141f4b8cb8d3d7ec769')],
 'DESI_DR2_BAO+CMB+PantheonPlus': [('pantheon_chain.1.txt',
                                    'db81d299d59051ae5e1e8f67952320e08a1644a4a1f8d19267705aee32798877'),
                                   ('pantheon_chain.2.txt',
                                    '442390266f6a3bfe88e9c7cd8e29d1e7cff47fc8582f3a3af80a6b40b9f83939'),
                                   ('pantheon_chain.3.txt',
                                    '1c92edf523df34784633f6852f7564f3d953203ae939d5d2e8cd4c3fd6f580ae'),
                                   ('pantheon_chain.4.txt',
                                    'b17dc02689c3a30dd34a96d91ac35180006f17d10b82331396ef55ce999594af')],
 'DESI_DR2_BAO+CMB+Union3': [('union3_chain.1.txt',
                              'f95d091203f615e1654ea9f4fea332db08fa66c59a5fe188b4915e67e1d73b11'),
                             ('union3_chain.2.txt',
                              '23aac58326257d207b8ee3699d21ed908e46238e42041d8c2ccf12fd85c8b839'),
                             ('union3_chain.3.txt',
                              '9e0998c5685e7552b94fd45e985dd137604f34692bebbc6fb922a1588669ab6a'),
                             ('union3_chain.4.txt',
                              '7710ababc2ab0c7fb5f3f67c8b790409703bce18707a2f99a1cd3fea531b4c20')],
 'base_lcdm_capacity_display': [('lcdm_chain.1.txt',
                                 '00f3766f7a7b6370d21323886cd72869087b2b1346a04d729c8f3bc9e65ef698'),
                                ('lcdm_chain.2.txt',
                                 '33b154eebdf4e9dca3b8f02ed2680120879d35c10b32fef42261a490104e1dc1'),
                                ('lcdm_chain.3.txt',
                                 'd4717e7e5a13de851c86f24c87213faccef2b5f8747900274ab509d9dfa40aa2'),
                                ('lcdm_chain.4.txt',
                                 'c827cd767a4864ca28aa15c902bda32004e803050d4be330e25aefddd78b5c36')]}


# Canonical scientific context reviewed at a14a3dac. This binds provenance,
# domains, formulas and interpretation; numerical outputs are replayed below.
# The digest excludes all computed statistics and per-file records, whose
# independent checks must still pass. It is not a replacement for replay.
CONTRACT_SHA256 = "ba85bfdad9e987b30c1cb0c02dfb8a2d9f1003ff3da6326e1a4e5777c133a726"


def _keys(value, keys, path):
    if not isinstance(value, dict) or value.keys() != set(keys):
        raise ValueError(f'{path}: field set mismatch')


def scientific_context(receipt):
    """Extract every declaration, refusing unrecognized receipt fields."""
    metadata = ('schema', 'producer', 'source', 'arithmetic',
                'epistemic_status', 'subset_definition')
    _keys(receipt, (*metadata, 'datasets', 'base_lcdm_capacity_display'), 'receipt')
    _keys(receipt['datasets'], set(PINS)-{'base_lcdm_capacity_display'}, 'datasets')
    context = {key: receipt[key] for key in metadata}
    context['datasets'] = {}
    for name, dataset in receipt['datasets'].items():
        _keys(dataset, ('source_directory', 'chains', 'combined',
                       'chain_range_for_monotone_subset_mass', 'rare_tail_resolution',
                       'fixed_capacity_point_gaussian_diagnostic'), name)
        tail = dataset['rare_tail_resolution']
        _keys(tail, ('classification', 'combined_raw_tail_rows', 'warning'), name+'/tail')
        diagnostic = dataset['fixed_capacity_point_gaussian_diagnostic']
        _keys(diagnostic, ('classification', 'status', 'mahalanobis_squared',
                          'chi2_2dof_survival', 'log_chi2_2dof_survival',
                          'two_sided_normal_sigma_equivalent'), name+'/diagnostic')
        context['datasets'][name] = {
            'source_directory': dataset['source_directory'],
            'tail_warning': tail['warning'],
            'diagnostic_classification': diagnostic['classification'],
        }
    base = receipt['base_lcdm_capacity_display']
    _keys(base, ('source_directory', 'model_scope', 'sample_level_formula',
                 'constants', 'chains', 'combined', 'classification'), 'base LCDM')
    context['base_lcdm_capacity_display'] = {key: value for key, value in base.items()
                                          if key not in ('chains', 'combined')}
    return context


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'duplicate JSON key: {key}')
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError(f'non-finite JSON constant: {value}')


def compare(actual, expected, path='summary'):
    """No absolute floor; exact counts/zeros/shape, relative 1e-12 otherwise."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or actual.keys() != expected.keys():
            raise ValueError(f'{path}: field set mismatch')
        for key, value in expected.items():
            compare(actual[key], value, f'{path}/{key}')
    elif expected is None or isinstance(expected, (int, str)):
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(f'{path}: exact value mismatch')
    else:
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise ValueError(f'{path}: expected finite number')
        value = Decimal(actual)
        if not value.is_finite() or abs(value-expected) > abs(expected)*Decimal('1e-12'):
            raise ValueError(f'{path}: numerical mismatch ({value} != {expected})')


def control(rows, base=False):
    """Two-pass centered moments, independent of the producer's raw sums."""
    weights = [row[0] for row in rows]
    values = [list(row[1:]) for row in rows]
    if base:
        for row in values:
            h0, omega = row
            row.append(3*omega*(h0*1000/Decimal('3.0856775814913673e22') /
                               Decimal(299792458))**2*Decimal('1.616255e-35')**2)
    total = sum(weights)
    means = [sum(w*x[i] for w,x in zip(weights,values))/total for i in range(len(values[0]))]
    centered = [[v-m for v,m in zip(row,means)] for row in values]
    cov = [[sum(w*x[i]*x[j] for w,x in zip(weights,centered))/total
            for j in range(len(means))] for i in range(len(means))]
    corr = cov[0][1]/(cov[0][0]*cov[1][1]).sqrt() if cov[0][0] and cov[1][1] else None
    result = dict(raw_rows=len(rows), expanded_posterior_weight=total,
                  weight_concentration_ess_not_autocorrelation_corrected=total**2/sum(w*w for w in weights))
    if base:
        for i,name in enumerate(('H0_km_s_Mpc','OmegaLambda','Lambda_lP2')):
            result[name] = dict(weighted_mean=means[i], weighted_std=cov[i][i].sqrt())
        quantiles = {}
        ordered = sorted((row[2],w) for row,w in zip(values,weights))
        for p in map(Decimal, ('.025','.16','.5','.84','.975')):
            cumulative = Decimal(0)
            for value,w in ordered:
                cumulative += w
                if cumulative >= p*total:
                    quantiles[f'{p:.3f}'] = value
                    break
        result['H0_OmegaLambda_weighted_correlation'] = corr
        result['Lambda_lP2'].update(fractional_std_about_weighted_mean=cov[2][2].sqrt()/means[2],
                                  weighted_step_cdf_quantiles=quantiles)
    else:
        # Integer-denominator endpoint comparison, separately counted complements.
        admitted = [x >= -1 and 3*x+2*y >= -3 for _,x,y in rows]
        result.update(w0_mean=means[0], wa_mean=means[1], w0_std=cov[0][0].sqrt(),
                      wa_std=cov[1][1].sqrt(), w0_wa_covariance=cov[0][1], w0_wa_correlation=corr,
                      raw_rows_in_monotone_subset=sum(admitted))
        for name, keep in (
            ('w_ge_minus_one_for_0_le_z_le_2', admitted),
            ('capacity_loss_somewhere_for_0_le_z_le_2', [not v for v in admitted]),
            ('w0_gt_minus_one', [x > -1 for _,x,_ in rows]),
            ('wa_nonnegative', [y >= 0 for _,_,y in rows])):
            result['posterior_mass_'+name] = sum(w for w,k in zip(weights,keep) if k)/total
    return result, means, cov


def replay(receipt, data_dir):
    context = json.dumps(scientific_context(receipt), sort_keys=True,
                         separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode('utf-8')
    if hashlib.sha256(context).hexdigest() != CONTRACT_SHA256:
        raise ValueError('scientific context differs from the reviewed v3 contract')
    fields = 0
    with localcontext() as ctx:
        ctx.prec = 160
        for name, pins in PINS.items():
            base = name == 'base_lcdm_capacity_display'
            dataset = receipt[name] if base else receipt['datasets'][name]
            if len(dataset['chains']) != 4:
                raise ValueError(f'{name}: expected four chains')
            combined = []
            for index, (chain, (filename, expected_hash)) in enumerate(zip(dataset['chains'], pins), 1):
                if chain['file'] != filename or chain['sha256'] != expected_hash or type(chain['chain']) is not int or chain['chain'] != index:
                    raise ValueError(f'{name}: source pin mismatch')
                payload = (data_dir/filename).read_bytes()
                if hashlib.sha256(payload).hexdigest() != expected_hash or type(chain['bytes']) is not int or chain['bytes'] != len(payload):
                    raise ValueError(f'{filename}: source bytes mismatch')
                lines = payload.decode('utf-8').splitlines()
                header = lines[0].lstrip('#').split()
                columns = [header.index(c) for c in (('weight','H0','omegal') if base else ('weight','w','wa'))]
                rows = [tuple(Decimal(parts[i]) for i in columns)
                        for line in lines[1:] if line.strip() and not line.startswith('#')
                        for parts in [line.split()]]
                expected, _, _ = control(rows, base)
                observed = {k:v for k,v in chain.items() if k not in ('chain','file','bytes','sha256')}
                compare(observed, expected, filename)
                fields += len(expected)
                combined.extend(rows)
            expected, mean, cov = control(combined, base)
            compare(dataset['combined'], expected, name+'/combined')
            fields += len(expected)
            if not base:
                d0, da = mean[0]+1, mean[1]
                # Solve by elimination; do not reuse the producer's inverse numerator.
                slope = cov[0][1]/cov[0][0]
                q = d0*d0/cov[0][0] + (da-slope*d0)**2/(cov[1][1]-slope*cov[0][1])
                diagnostic = dataset['fixed_capacity_point_gaussian_diagnostic']
                if diagnostic['status'] != 'available':
                    raise ValueError(name+': expected available diagnostic on these official chains')
                for key, value in [('mahalanobis_squared',q), ('log_chi2_2dof_survival',-q/2),
                                   ('chi2_2dof_survival',(-q/2).exp())]:
                    compare(diagnostic[key], value, name+'/'+key)
                # Independent high-precision erfc, not the producer's inverse normal.
                import mpmath
                mp = mpmath.mp.clone(); mp.dps = 80
                sigma = diagnostic['two_sided_normal_sigma_equivalent']
                if not isinstance(sigma, (int,float)) or isinstance(sigma,bool) or sigma < 0:
                    raise ValueError(name+': invalid normal sigma')
                log_tail = Decimal(str(mp.log(mp.erfc(mp.mpf(sigma)/mp.sqrt(2)))))
                if not log_tail.is_finite() or abs(log_tail+q/2) > abs(q)*Decimal('1e-12'):
                    raise ValueError(name+': normal sigma mismatch')
                masses = [c['posterior_mass_w_ge_minus_one_for_0_le_z_le_2'] for c in dataset['chains']]
                if dataset['chain_range_for_monotone_subset_mass'] != [min(masses),max(masses)]:
                    raise ValueError(name+': chain range mismatch')
                tail = dataset['rare_tail_resolution']
                classification = ('resolved_in_all_four_chains' if all(c['raw_rows_in_monotone_subset'] > 0
                                  for c in dataset['chains']) else 'at_least_one_chain_has_no_raw_tail_row')
                if tail['classification'] != classification or type(tail['combined_raw_tail_rows']) is not int or tail['combined_raw_tail_rows'] != expected['raw_rows_in_monotone_subset']:
                    raise ValueError(name+': rare-tail classification mismatch')
    return fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, default=RECEIPT)
    args = parser.parse_args()
    receipt = json.loads(args.receipt.read_text(encoding='utf-8'),
                         object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    fields = replay(receipt, args.data_dir)
    print(f'PASS: 20 pinned chains, 25 summaries ({fields} top-level fields), and four Gaussian diagnostics')


if __name__ == '__main__':
    main()
