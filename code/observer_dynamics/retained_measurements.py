"""Retained scale statistics: independent moments, not a raw-array replay."""
import json
import math
from fractions import Fraction as Q

try:
    from .archive_adapter import SIM, sha, load
except ImportError:
    from archive_adapter import SIM, sha, load

RTOL = Q(1, 10**9)
DATA = 'data/codex_audit_20260925/'


def need(condition, message):
    if not condition:
        raise ValueError(message)


def number(value):
    need(type(value) is int or (type(value) is float and math.isfinite(value)),
         'Retained statistic must be finite numeric data, not a Boolean')
    return Q(value)


def integer(value, minimum=1):
    need(type(value) is int and value >= minimum, 'Retained positive integer required')
    return value


def close(actual, expected, *, scale=None):
    actual = number(actual) if not isinstance(actual, Q) else actual
    expected = number(expected) if not isinstance(expected, Q) else expected
    magnitude = max(abs(actual), abs(expected)) if scale is None else abs(scale)
    need(abs(actual-expected) <= RTOL*magnitude, 'Retained arithmetic mismatch')


def response_statistics(responses, average, deviation, error=None):
    count = len(responses)
    mean = sum(responses)/count
    variance = sum((a-mean)**2 for a in responses)/(count-1)
    close(average, mean, scale=max(map(abs, responses)))
    sd = number(deviation)
    need(sd >= 0, 'Negative schedule standard deviation')
    close(sd**2, variance)
    if error is not None:
        se = number(error)
        need(se >= 0, 'Negative schedule standard error')
        close(se**2, variance/count)


def sqrt_bounds(value):
    """Outward rational bounds, absolute width <= 2^-60 on [0, 1]."""
    units = 1 << 60
    floor = math.isqrt((value.numerator << 120)//value.denominator)
    lower = Q(floor, units)
    return lower, lower if lower**2 == value else Q(floor+1, units)


def ensemble_bounds(responses, powers, noise):
    """Check a common vector population, with at least two residual dimensions.

    a_i=<x,y_i>/||x||^2, p_i=||y_i||^2/||x||^2 and n is the unbiased
    schedule variance in the same units. The residual vectors have squared
    lengths p_i-a_i^2. Their sum must have a length in the polygon interval.
    This checks realizability of these moments, not native bounded loads or
    the unretained arrays. Comparisons allow RTOL on the arithmetic scale.
    """
    need(isinstance(responses, (list, tuple)) and isinstance(powers, (list, tuple))
         and len(responses) == len(powers) and len(powers) >= 2,
         'Matched schedule population required')
    a, p = ([number(x) if not isinstance(x, Q) else x for x in values]
            for values in (responses, powers))
    n = number(noise) if not isinstance(noise, Q) else noise
    need(all(x >= 0 for x in p) and n >= -RTOL*max(p),
         'Negative schedule power or noise exceeds roundoff')
    residuals = []
    for response, power in zip(a, p):
        residual = power-response**2
        error = RTOL*max(power, response**2)
        need(residual >= -error,
             'Impossible schedule covariance: response exceeds power')
        residuals.append((max(Q(0), residual-error), max(Q(0), residual+error)))
    count = len(p)
    total = count*sum(p)
    subtraction = count*(count-1)*n
    response_sum = sum(a)**2
    resultant = total-subtraction-response_sum
    scale = max(upper for _, upper in residuals)
    if scale:
        # Propagate moment uncertainty BEFORE the ill-conditioned square root
        # at a zero residual. The retained power may have rounded away a real
        # residual. These outward rational roots avoid new floating losses.
        lengths = [(sqrt_bounds(lo/scale)[0], sqrt_bounds(hi/scale)[1])
                   for lo, hi in residuals]
        perimeter = sum(hi for _, hi in lengths)
        minimum = max(Q(0), max(lo-(perimeter-hi) for lo, hi in lengths))
        lower, upper = minimum**2*scale, perimeter**2*scale
    else:
        lower = upper = Q(0)
    tolerance = RTOL*max(total, subtraction, response_sum, upper)
    need(lower-tolerance <= resultant <= upper+tolerance,
         'Impossible schedule covariance: residual population cannot realize noise')


def geometry_factor(scale, carriers):
    groups = integer(scale['groups'], 4)
    low, high = integer(scale['cells_min']), integer(scale['cells_max'])
    arithmetic, harmonic = (number(scale[k]) for k in
                            ('cells_arithmetic_mean', 'cells_harmonic_mean'))
    need(low <= high and harmonic > 0 and arithmetic > 0
         and low-harmonic <= RTOL*max(low, harmonic)
         and harmonic-arithmetic <= RTOL*max(harmonic, arithmetic)
         and arithmetic-high <= RTOL*max(arithmetic, high),
         'Invalid retained group geometry')
    need(groups*arithmetic <= carriers+RTOL*carriers, 'Groups exceed carrier population')
    name = scale['name']
    if name.startswith('ancestor_'):
        size = int(name.removeprefix('ancestor_'))
        need(groups*size == carriers and low == high == arithmetic == harmonic == size,
             'Ancestral geometry differs from parent population')
    else:
        resolution = int(name.split('_')[1])
        need(low >= 4 and groups <= resolution**3, 'Invalid retained cube geometry')
    cells = arithmetic if name.endswith('_weighted') else harmonic
    return Q(groups-1, 12*groups)/cells


def check_scale(scale, seeds, carriers, coefficient):
    factor = geometry_factor(scale, carriers)
    equilibrium = coefficient*factor
    close(scale['equilibrium_variance'], equilibrium)
    rows = scale['per_schedule']
    need([integer(r['seed'], 0) for r in rows] == seeds, 'Scale lost or reordered schedule')
    initial = number(rows[0]['initial_variance'])
    need(initial > 0, 'Positive initial variance required')
    close(scale['initial_variance_over_iid_null'], initial/(Q(35, 12)*factor))
    responses, powers, variances = [], [], []
    for row in rows:
        close(row['initial_variance'], initial)
        variance = number(row['terminal_variance'])
        need(variance > 0, 'Positive terminal variance required')
        response = number(row['response'])
        power = variance/initial
        residual = power-response**2
        close(row['power_ratio'], power)
        close(row['coherence_squared'], response**2/power)
        close(row['residual_power_ratio'], residual, scale=max(power, response**2))
        close(row['variance_over_exchangeable_equilibrium'], variance/equilibrium)
        number(row['legacy_unweighted_global_centered_response'])
        responses.append(response); powers.append(power); variances.append(variance)
    response_statistics(responses, scale['response_mean'], scale['response_schedule_sd'],
                        scale['response_schedule_se'])
    total = sum(variances)/len(variances)
    noise = number(scale['conditional_schedule_variance_unbiased'])
    common = number(scale['cross_schedule_common_power_unbiased'])
    close(scale['total_terminal_variance'], total)
    close(noise+common, total, scale=max(total, abs(noise), abs(common)))
    close(scale['schedule_noise_fraction'], noise/total)
    ensemble_bounds(responses, powers, noise/initial)


def check_angular(row, schedules):
    initial = number(row['C_initial'])
    need(initial > 0, 'Positive initial angular power required')
    responses = [number(x) for x in row['per_schedule_response']]
    powers = [number(x) for x in row['per_schedule_power_ratio']]
    need(len(responses) == len(powers) == schedules and all(p > 0 for p in powers),
         'Complete positive angular schedule population required')
    mean_power = sum(powers)/schedules
    response_statistics(responses, row['cross_response_mean'], row['cross_response_schedule_sd'])
    close(row['total_power_ratio_mean'], mean_power)
    close(row['C_terminal_mean'], initial*mean_power)
    close(row['coherence_squared_mean'], sum(a*a/p for a, p in zip(responses, powers))/schedules)
    close(row['orthogonal_residual_power_ratio_mean'],
          sum(p-a*a for a, p in zip(responses, powers))/schedules,
          scale=max(mean_power, sum(a*a for a in responses)/schedules))
    noise = number(row['conditional_schedule_noise_power'])
    common = number(row['cross_schedule_common_power'])
    close(noise+common, initial*mean_power,
          scale=max(initial*mean_power, abs(noise), abs(common)))
    ensemble_bounds(responses, powers, noise/initial)


def scale_ensemble(value=None):
    value = json.loads((SIM/(DATA+'scale_ensemble.json')).read_text()) if value is None else value
    spec_path = SIM/(DATA+'scale_ensemble_spec.json')
    spec = json.loads(spec_path.read_text())
    need(value['schema'] == 'oph.sim.codex-scale-ensemble.v1'
         and json.dumps(value['spec'], sort_keys=True) == json.dumps(spec, sort_keys=True)
         and value['spec_sha256'] == sha(spec_path)
         and value['script_sha256'] == sha(SIM/'scripts/codex_scale_ensemble.py'),
         'Retained scale specification or producer identity differs')
    need([integer(r['level']) for r in value['levels']] == spec['levels'], 'Unexpected tower levels')
    crosspins = rows = angular_rows = 0
    for level in value['levels']:
        path = SIM/(DATA+f"inputs/data/receipts/L{level['level']}/receipt.json")
        need(sha(path) == level['source_receipt_sha256'], 'Parent settlement receipt hash differs')
        parent = json.loads(path.read_text()); entries = parent['integer_law']['entries']
        carriers, ports = parent['carriers'], parent['ports']
        need(integer(level['carriers']) == carriers and integer(level['ports']) == ports == 12*carriers,
             'Retained parent population differs')
        need(integer(level['initial_load_seed']) == 20260909+level['level']
             and level['initial_load_dtype'] == 'int64'
             and level['geometry_file_sha256'] == parent['geometry']['arrays']['cell_points']['sha256'],
             'Retained preparation or geometry differs')
        remainder = parent['expected']['component_total'][0] % ports
        close(level['balanced_fraction'], Q(remainder, ports))
        coefficient = Q(remainder*(ports-remainder), ports*(ports-1))
        need(len(entries) == integer(level['schedules']) == len(level['inputs']) == 16, 'Missing schedules')
        seeds = [r['seed'] for r in entries]
        for source, retained in zip(entries, level['inputs']):
            need(integer(retained['seed'], 0) == source['seed']
                 and source['terminal_state']['sha256'] == retained['array_sha256'],
                 'Terminal-state digest crosspin differs')
            crosspins += 1
        names = [f'ancestor_{size}' for size in spec['primary_ancestral_group_cells']
                 if carriers % size == 0 and carriers//size >= 20]
        names += [f'cube_{r}_{mode}' for r in spec['cube_control_resolutions']
                  for mode in ('weighted', 'unweighted')]
        need([s['name'] for s in level['scales']] == names, 'Missing or reordered scale population')
        for scale in level['scales']:
            check_scale(scale, seeds, carriers, coefficient)
            rows += 1
        by_name = {s['name']: s for s in level['scales']}
        for resolution in spec['cube_control_resolutions']:
            weighted, unweighted = (by_name[f'cube_{resolution}_{mode}']
                                    for mode in ('weighted', 'unweighted'))
            need(all(weighted[k] == unweighted[k] for k in
                     ('groups', 'cells_min', 'cells_max', 'cells_arithmetic_mean', 'cells_harmonic_mean')),
                 'Weighted/unweighted cube geometry differs')
        grids = {f'z{nz}_phi{nphi}': nz*nphi for nz, nphi in spec['angular_grids']}
        need(set(level['angular']) == set(grids), 'Missing angular grids')
        for key, pixels in grids.items():
            grid = level['angular'][key]
            low, high = integer(grid['counts_min']), integer(grid['counts_max'])
            need(low*pixels <= carriers <= high*pixels, 'Invalid angular pixel population')
            need([integer(r['ell']) for r in grid['rows']] == list(range(1, spec['angular_lmax']+1)),
                 'Missing or reordered angular multipoles')
            for row in grid['rows']:
                check_angular(row, len(seeds))
                angular_rows += 1
    return {'terminal_digest_crosspins': crosspins, 'scale_rows': rows,
            'angular_rows': angular_rows, 'moment_relative_tolerance': float(RTOL),
            'raw_terminal_arrays_rehashed': False}


def orbit_census():
    # Burnside coefficients from the actual pinned carrier rotation action.
    from oph_exact import carrier
    rotations=carrier.rotations(); fixed=[0]*13
    if len(rotations)!=60: raise ValueError('Expected the sixty carrier rotations')
    for rotation in rotations:
        seen=set();polynomial=[1]
        for start in range(12):
            if start in seen:continue
            current=start;length=0
            while current not in seen:
                seen.add(current);length+=1;current=int(rotation[current])
            result=polynomial+[0]*length
            for i,c in enumerate(polynomial): result[i+length]+=c
            polynomial=result
        fixed=[a+b for a,b in zip(fixed,polynomial,strict=True)]
    if any(n%60 for n in fixed): raise ValueError('Burnside quotient not integral')
    coefficients=[n//60 for n in fixed]
    total_states=0
    for level in (6,8,9,10):
        record=json.loads((SIM/f'data/census_L{level}.json').read_text())
        if record['orbits_total']!=sum(coefficients) or [record['orbits_by_k'][str(k)] for k in range(13)]!=coefficients:
            raise ValueError('Orbit census disagrees with carrier action')
        for state in record['states']:
            if sum(state['raised_count_histogram'])!=record['carriers']: raise ValueError('Census does not cover all carriers')
            total_states+=1
    return {'orbits':sum(coefficients),'by_raised_count':coefficients,'retained_states':total_states,'raw_terminal_arrays_recounted':False}


def conditional_measured_comparison():
    module=load('codex/observer_cmb/measured_comparison.py')
    columns,tables,cal,residuals,fresh=module.prepare()
    retained=json.loads((SIM/'codex/observer_cmb/measured_comparison.json').read_text())
    outputs=retained.pop('outputs_sha256')
    module.observer._compare(retained,fresh)
    for channel,rows in residuals.items():
        name=f'measured_residuals_{channel}.txt';path=SIM/'codex/observer_cmb'/name
        if path.read_bytes()!=module.table_bytes(rows) or sha(path)!=outputs[name]:
            raise ValueError('Changed conditional measured residual table')
    return {'channels':list(residuals),'standard_transfer_recomputed':False,'retained_spectrum_and_residuals_recomputed':True,
            'amplitude_tilt_background':'supplied inputs','plot_hashes_rechecked':False}
