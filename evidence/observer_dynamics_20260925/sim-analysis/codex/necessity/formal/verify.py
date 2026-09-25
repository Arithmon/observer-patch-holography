"""Compile the existing cocycle and local witness without modifying RER build state."""
from __future__ import annotations
import hashlib,json,os,re,subprocess,tempfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[3]
LEAN=WORKSPACE/'reverse-engineering-reality'/'Lean'
MODULE=Path('ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltCocycle')
DEPENDENCY=LEAN/(str(MODULE)+'.lean')
SOURCE=HERE/'NecessityCore.lean'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    version=subprocess.check_output(['lake','env','lean','--version'],cwd=LEAN,text=True).strip()
    leanpath=subprocess.check_output(['lake','env','printenv','LEAN_PATH'],cwd=LEAN,text=True).strip()
    logs=[]
    source_before,dependency_before=sha(SOURCE),sha(DEPENDENCY)
    with tempfile.TemporaryDirectory(prefix='oph-necessity-proof-') as temp:
        obj=Path(temp)/(str(MODULE)+'.olean');obj.parent.mkdir(parents=True)
        commands=[['lake','env','lean','-o',str(obj),str(DEPENDENCY)],['lean',str(SOURCE)]]
        env=os.environ.copy();env['LEAN_PATH']=temp+os.pathsep+leanpath
        for command in commands:
            p=subprocess.run(command,cwd=LEAN,env=env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            logs.append(p.stdout)
            if p.returncode:
                (HERE/'failed_build.log').write_text('\n'.join(logs))
                raise SystemExit(p.stdout)
    log='\n'.join(logs)
    if 'sorryAx' in log or 'error:' in log: raise SystemExit('Unexpected unproved/error output')
    allowed={'propext','Classical.choice','Quot.sound'}
    for group in re.findall(r'depends on axioms: \[(.*?)\]',log,re.S):
        actual={name.strip() for name in group.split(',') if name.strip()}
        if not actual <= allowed: raise SystemExit('Unexpected custom proof axioms: '+str(actual-allowed))
    names=re.findall(r'^theorem (\w+)',SOURCE.read_text(),re.M)
    for name in names:
        if f"'Codex.Necessity.{name}' depends on axioms:" not in log:
            raise SystemExit('Missing theorem dependency report: '+name)
    if (sha(SOURCE),sha(DEPENDENCY))!=(source_before,dependency_before):
        raise SystemExit('Proof source changed during verification')
    (HERE/'build.log').write_text(log)
    receipt={'schema':'oph.kernel-necessity-cores.v1','lean_version':version,'exit_code':0,
             'dependency':str(DEPENDENCY.relative_to(WORKSPACE)),'dependency_sha256':sha(DEPENDENCY),
             'source_sha256':sha(SOURCE),'verifier_sha256':sha(Path(__file__)),
             'log_sha256':sha(HERE/'build.log'),
             'theorem_count':len(re.findall(r'^theorem ',SOURCE.read_text(),re.M)),
             'scope':'Kernel-checked finite algebra, cocycle freedom and pointwise filter inequality; no formalization of stochastic CLT or cosmology identification',
             'allowed_axioms':['propext','Classical.choice','Quot.sound'],
             'no_sorryAx_in_printed_axioms':True}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(version);print('Both Lean files compiled; no sorryAx in printed theorem dependencies')

if __name__=='__main__':main()
