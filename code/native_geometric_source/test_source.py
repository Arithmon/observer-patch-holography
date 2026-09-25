"""Adversarial checks of the scientific interpretation and exact covariance data."""
import copy
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
import pytest
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('native_geometry_verifier',HERE/'verify.py')
verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)

@pytest.fixture(scope='module')
def receipt():
    return json.loads((verifier.EVIDENCE/'receipt.json').read_text())

def test_independent_exact_enumeration_and_native_replay():
    assert verifier.verify()['status']=='PASS'

def apply_mutation(r,kind):
    finite=r['finite_controls'][0]; readout=finite['readouts']['load']
    if kind=='hidden_gain': finite['kappa']='1/6'
    elif kind=='source_promotion': r['decision']='GEOMETRIC_CURVATURE_DERIVED'
    elif kind=='observed_target': r['observational_data_read']=True
    elif kind=='second_preparation_missing': r['finite_controls'].pop()
    elif kind=='equal_diagonal_false_shape': finite['unit_total_variance_covariances']['local_drive']=copy.deepcopy(finite['unit_total_variance_covariances']['load'])
    elif kind=='separated_covariance': readout['instant_covariance'][0][3]='0'
    elif kind=='cross_scale': readout['cross_scale_covariance'][0][1]='0'
    elif kind=='clock': readout['two_attempts_per_record_covariance']['4']=copy.deepcopy(readout['record_average_covariance']['4'])
    elif kind=='projection': readout['density_projection_B'][0][0]='1'
    elif kind=='nonlinear_residual': r['finite_controls'][1]['readouts']['local_mismatch']['residual_covariance'][0][0]='0'
    elif kind=='invented_volume': r['geometry_traces'][0]['physical_collar_volume_measured']=True
    elif kind=='changed_geometry': r['geometry_traces'][0]['geometry_max_log_ratio']=0.1
    elif kind=='geometry_hash': r['geometry_traces'][0]['geometry_hash']='0'*64
    elif kind=='omitted_refinement': r['geometry_traces']=r['geometry_traces'][:2]
    elif kind=='native_final_load': r['geometry_traces'][0]['final_loads'][0]+=1
    elif kind=='source_pin': r['source_sha256']['oph_exact/carrier.py']='0'*64
    else: raise ValueError(kind)

@pytest.mark.parametrize('kind',[
    'hidden_gain','source_promotion','observed_target','second_preparation_missing',
    'equal_diagonal_false_shape','separated_covariance','cross_scale','clock',
    'projection','nonlinear_residual','invented_volume','changed_geometry',
    'geometry_hash','omitted_refinement','native_final_load','source_pin'])
def test_false_green_mutations_fail(receipt,tmp_path,kind):
    r=copy.deepcopy(receipt);apply_mutation(r,kind)
    path=tmp_path/'mutated.json';path.write_text(json.dumps(r))
    with pytest.raises(ValueError): verifier.verify(path)


def test_native_custody_checks_survive_python_optimization(receipt,tmp_path):
    r=copy.deepcopy(receipt);r['geometry_traces'][0]['geometry_hash']='0'*64
    path=tmp_path/'mutated.json';path.write_text(json.dumps(r))
    result=subprocess.run([sys.executable,'-O',str(HERE/'check_native_geometry.py'),str(path),str(verifier.EVIDENCE/'spec.json')],capture_output=True,text=True)
    assert result.returncode!=0 and 'native geometry hash' in result.stderr
