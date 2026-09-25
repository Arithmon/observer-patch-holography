"""Cross-check retained large-run readouts without loading raw terminal states."""
import json
import math
from pathlib import Path
from statistics import mean, stdev
try:
    from .archive_adapter import SIM, sha, load
except ImportError:
    from archive_adapter import SIM, sha, load


def close(a,b):
    if not math.isclose(a,b,rel_tol=1e-9,abs_tol=1e-12):
        raise ValueError(f"Retained arithmetic mismatch: {a} != {b}")


def scale_ensemble(value=None):
    value = json.loads((SIM/'data/codex_audit_20260925/scale_ensemble.json').read_text()) if value is None else value
    if [r['level'] for r in value['levels']] != [8,9,10]: raise ValueError('Unexpected tower levels')
    crosspins=rows=0
    for level in value['levels']:
        path=SIM/f"data/codex_audit_20260925/inputs/data/receipts/L{level['level']}/receipt.json"
        if sha(path)!=level['source_receipt_sha256']: raise ValueError('Parent settlement receipt hash differs')
        parent=json.loads(path.read_text());entries=parent['integer_law']['entries']
        if len(entries)!=16 or level['schedules']!=16 or len(level['inputs'])!=16: raise ValueError('Missing schedules')
        seeds=[r['seed'] for r in entries]
        for source, retained in zip(entries,level['inputs'],strict=True):
            if source['seed']!=retained['seed'] or source['terminal_state']['sha256']!=retained['array_sha256']:
                raise ValueError('Terminal-state digest crosspin differs')
            crosspins+=1
        for scale in level['scales']:
            rows+=1
            records=scale['per_schedule']
            if [r['seed'] for r in records]!=seeds: raise ValueError('Scale lost or reordered schedule')
            responses=[r['response'] for r in records]
            close(mean(responses),scale['response_mean']);close(stdev(responses),scale['response_schedule_sd'])
            close(scale['response_schedule_sd']/4,scale['response_schedule_se'])
            close(scale['total_terminal_variance'],mean(record['terminal_variance'] for record in records))
            close(scale['total_terminal_variance'],scale['conditional_schedule_variance_unbiased']+scale['cross_schedule_common_power_unbiased'])
            close(scale['schedule_noise_fraction'],scale['conditional_schedule_variance_unbiased']/scale['total_terminal_variance'])
            for record in records:
                close(record['variance_over_exchangeable_equilibrium'],record['terminal_variance']/scale['equilibrium_variance'])
                close(record['power_ratio'],record['terminal_variance']/record['initial_variance'])
                close(record['residual_power_ratio'],record['power_ratio']-record['response']**2)
    return {'terminal_digest_crosspins':crosspins,'scale_rows':rows,'raw_terminal_arrays_rehashed':False}


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
