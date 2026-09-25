"""False-green mutations and an independent driven-oscillator control."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.linalg import expm

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_scalar_preparation',HERE/'verify_preparation.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)

def receipt():return json.loads((HERE/'preparation_receipt.json').read_text())

def test_independent_certificate():
    assert checker.check()['verified']

@pytest.mark.parametrize('mutation',[
    lambda d:d['moments'].__setitem__('Av_mass_squared_Qphi',['0','0']),
    lambda d:d['bounds'].__setitem__('pulse_trace_error_upper','0'),
    lambda d:d['bounds'].__setitem__('baseline_trace_error_upper','0'),
    lambda d:d['bounds'].__setitem__('preparation_end_to_first_pointer_gap_lower','1'),
    lambda d:d['parameters'].__setitem__('delta_max','1/100'),
    lambda d:d['law'].__setitem__('Hamiltonian','H0 + force'),
    lambda d:d['control_schedule'].pop(),
    lambda d:d['control_schedule'][0].__setitem__('site',32),
    lambda d:d['control_schedule'][0].__setitem__('force_coefficient_m_v_Qphi',['0','0']),
    lambda d:d['control_schedule'][0].__setitem__('start_over_delta','0'),
    lambda d:d['control_schedule'][0].__setitem__('off_write','retain previous force'),
    lambda d:d['resource_model'].__setitem__('on_off_writes',16),
    lambda d:d['resource_model'].__setitem__('coefficient_bit_precision','hardware precision certified'),
    lambda d:d['scope'].__setitem__('original_vacuum_at_pulse_start_supplied',False),
    lambda d:d['scope'].__setitem__('virtual_center_state_is_actual_midpulse_state',True),
    lambda d:d['scope'].__setitem__('actual_quantum_hardware_execution',True),
    lambda d:d['instrument_composition'].__setitem__('resolved_steps',[16]),
    lambda d:d['instrument_composition'].__setitem__('independent_shots_assumed',True),
])
def test_misbound_or_promoted_receipt_rejected(mutation):
    d=deepcopy(receipt());mutation(d)
    with pytest.raises(ValueError):checker.check(d)

@pytest.mark.parametrize('delta',[0.0001,0.002,0.03])
def test_force_solution_has_centered_sinc_preparation(delta):
    # Independent augmented exponential solves a coupled classical force ODE.
    # This checks the exact mean law, not quantum outcome sampling.
    a=np.array([[4.,-1.],[-1.,2.]])
    v=np.array([0.7,0.0])
    generator=np.zeros((5,5));generator[:2,2:4]=np.eye(2)
    generator[2:4,:2]=-a;generator[2:4,4]=v/delta
    initial=np.array([0.,0.,0.,0.,1.])
    actual=(expm(delta*generator)@initial)[:4]
    lam,u=np.linalg.eigh(a);omega=np.sqrt(lam)
    effective=u@(np.sinc(delta*omega/(2*np.pi))*(u.T@v))
    flow=generator[:4,:4]
    expected=expm(delta*flow/2)@np.r_[np.zeros(2),effective]
    assert np.max(np.abs(actual-expected))<2e-14
    # At pulse END compare to the ideal kick at virtual center time zero.
    reference=expm(delta*flow/2)@np.r_[np.zeros(2),v]
    diff=actual-reference
    distance_squared=(np.sum(omega*(u.T@diff[:2])**2)+np.sum((u.T@diff[2:])**2/omega))/2
    trace_distance=np.sqrt(-np.expm1(-distance_squared))
    sharp_bound=delta**2*np.linalg.norm(a@v)/(24*np.sqrt(2))
    assert trace_distance<=sharp_bound+1e-14
    # An impulse placed at the start has the wrong free-evolution clock.
    wrong=expm(delta*flow)@np.r_[np.zeros(2),v]
    assert np.linalg.norm(actual-wrong)>delta*np.linalg.norm(v)/3
