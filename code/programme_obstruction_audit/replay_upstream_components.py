#!/usr/bin/env python3
"""Optional replay on an exact pinned upstream source fixture with SymPy1.14.0.

Supply --source-root containing the source paths retained in fixed_sources.json
and additional_source_pins.json. Does not fetch, alter or generate source bytes.
Checks fixed-source custody before seven finite components;
not top-level upstream verifier/hostile suite/producer regeneration/Lean build.
"""
import argparse,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 if not __debug__:raise SystemExit('Do not use -O; assertions must be enabled.')
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-root',type=Path,required=True);a=p.parse_args();source=a.source_root.resolve()
 pins=json.loads((ROOT/'provenance/fixed_sources.json').read_text())['verified_blobs']
 extra=json.loads((ROOT/'provenance/additional_source_pins.json').read_text())
 relevant=[x for x in pins+extra if x['path'].startswith('code/source_selection_model/') or x['path'].startswith('code/a5_closure/')]
 for entry in relevant:
  raw=(source/entry['path']).read_bytes();assert hashlib.sha256(raw).hexdigest()==entry['sha256'],entry['path']
 sys.path.insert(0,str(source/'code'))
 from source_selection_model import channels,finite_model,geometry,grammar,verify_response
 import sympy,mpmath
 assert sympy.__version__=='1.14.0' and mpmath.__version__=='1.3.0'
 support=geometry.tower();actual={
 'response':verify_response.verify(verify_response.strict_load(source/'code/source_selection_model/response.json')),
 'algebra':finite_model.algebra_controls(),'constraint_grammar':grammar.controls(),'central_channels':channels.controls(),
 'support':geometry.check_tower(support),'records':finite_model.record_controls(verify_response.strict_load(source/'code/source_selection_model/records.json')),
 'refinement':finite_model.refinement_controls(support)}
 retained=json.loads((ROOT/'receipts/upstream_components.json').read_text())
 assert json.dumps(actual,sort_keys=True)==json.dumps(retained['components'],sort_keys=True)
 print(json.dumps({'target_revision':retained['target_revision'],'python':sys.version,'sympy':sympy.__version__,'mpmath':mpmath.__version__,
  'all_seven_components_equal_retained_receipt':True,'components':actual,'source_pins_verified':len(relevant),
  'scope':'bounded finite components only; all-level model membership is analytic, not inferred from finite support levels'},indent=2))

if __name__=='__main__':main()
