#!/usr/bin/env python3
"""Run the fixed offline audit controls; no search, network or model discovery.

Requires an ordinary, unoptimized Python (assertions enabled). Verifies package
hashes before executing seven prescribed controls. Each subprocess has a
120-second ceiling; failure/timeout stops the suite. Receipts are written to a
separate directory, never replacing retained evidence. Finite checks do not
prove all-level model membership, physical interpretation, or exclusion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parent

def main():
    if not __debug__:
        raise SystemExit('Assertions must be enabled; do not use python -O.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,
                        help='empty directory outside this package; default: fresh temporary directory')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'provenance/package_sha256.json').read_text())
    for rel,expected in manifest['files'].items():
        actual=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
        if actual!=expected: raise SystemExit(f'Package hash mismatch: {rel}')
    fixture=ROOT/'data/response.json'
    raw=fixture.read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='226ab001894f357d3732d8a898fbae80e28b9e16ca41104a6f293190b441e053'
    assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()=='e56b96a9429ac618aea9de7cd5ffa53a13354ff6'
    provenance=json.loads((ROOT/'provenance/fixed_sources.json').read_text())
    assert provenance['target_revision']=='0510578e317de5da0d485d09792dc8bcad2e2f3c'
    matched=[x for x in provenance['verified_blobs'] if x['path']=='code/source_selection_model/response.json']
    assert len(matched)==1 and matched[0]['sha256']==hashlib.sha256(raw).hexdigest()
    jobs=[('channel_author',[]),('channel_independent',[]),('joint_author',['data/response.json']),
          ('joint_independent',[]),('kinetic_control',[]),('attempt_independent',[]),('compact_classifier',[])]
    if args.output_dir is None:
        temporary_root = Path(tempfile.gettempdir()).resolve()
        if temporary_root == ROOT or ROOT in temporary_root.parents:
            raise SystemExit('Output directory must be outside the input package.')
        output_dir = Path(tempfile.mkdtemp(prefix='oph-programme-obstruction-', dir=temporary_root))
    else:
        output_dir = args.output_dir.resolve()
    if output_dir == ROOT or ROOT in output_dir.parents:
        raise SystemExit('Output directory must be outside the input package.')
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise SystemExit('Output directory must be empty; retained results are never overwritten.')
    output_dir.mkdir(parents=True,exist_ok=True)
    summary={'target_revision':'0510578e317de5da0d485d09792dc8bcad2e2f3c','python':sys.version,
             'hash_manifest_verified':True,'scope':'fixed finite controls; analytic membership is reviewed separately; no physical observations',
             'stopping_rule':'seven prescribed controls, stop on any failure or per-script120s ceiling; no adaptive expansion','runs':[]}
    for name,extra in jobs:
        start=time.monotonic()
        result=subprocess.run([sys.executable,'-B',str(ROOT/(name+'.py')),*extra],capture_output=True,text=True,timeout=120,cwd=ROOT)
        if result.returncode: raise SystemExit(f'{name} failed ({result.returncode}): {result.stderr}')
        packet=json.loads(result.stdout)
        destination=output_dir/(name+'.json');destination.write_text(result.stdout)
        summary['runs'].append({'script':name+'.py','exit_code':result.returncode,
                 'seconds':round(time.monotonic()-start,6),'receipt':destination.name,
                 'receipt_sha256':hashlib.sha256(destination.read_bytes()).hexdigest()})
    (output_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'all_seven_controls_passed':True,'immutable_fixture_and_package_hashes_verified':True,'receipt_file_count':len(jobs)+1,'full_summary':str(output_dir/'summary.json'),'scope':summary['scope']},indent=2))

if __name__=='__main__':main()
