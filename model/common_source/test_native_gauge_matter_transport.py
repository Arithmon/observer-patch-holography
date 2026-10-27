"""Tests of literal native-current linked CAR transport and its boundaries."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as s

from . import native_gauge_matter_transport as lane


@pytest.fixture(scope='module')
def packet():
    return lane.build()


def test_every_local_source_factor_retains_vacuum_and_double_occupation():
    for z in (s.Integer(1),-s.Integer(1),s.I,-s.I):
        factor=lane.local_factor('hopping',F(s.re(z)),F(s.im(z)),F(1,7))
        assert factor['actual4'].shape==(4,4)
        assert np.linalg.norm(factor['actual4']-factor['target4'],2)<float(lane.LOCAL_ADMISSION)
        assert abs(factor['actual4'][0,0]-1)<1e-11
        assert abs(factor['actual4'][3,3]-1)<1e-11
    swap=lane.local_factor('fswap')
    assert abs(swap['actual4'][3,3]+1)<1e-11
    bosonic=swap['target4'].copy();bosonic[3,3]=1
    assert np.linalg.norm(swap['actual4']-bosonic,2)>1.9


def test_routed_nonadjacent_bond_matches_all_sixteen_fock_columns():
    word=lane.FockWord(4);t=F(2,7);J=s.Rational(lane.field.native_parameters()['exchange'])
    word.term('hopping',0,3,-s.I*J,t)
    h=np.zeros((4,4),complex);h[0,3]=-1j*float(J);h[3,0]=h[0,3].conjugate()
    target=expm(-1j*float(t)*lane.field._fock_hamiltonian(h))
    assert np.linalg.norm(word.complete_matrix()-target,2)<float(word.error_bound())
    assert word.complete_matrix().shape==(16,16)
    assert any(word.catalogue[i]['kind']=='fswap' for i,_ in word.order)


def test_actual_current_paths_retain_all_auxiliary_levels(packet):
    report,arrays=packet
    rows=report['same_current_carrier']['canonical_current_requests']
    assert len(rows)>=12
    for row in rows:
        value=arrays[row['array']]
        assert value.shape==(32,32)
        assert np.linalg.norm(value.conj().T@value-np.eye(32),2)<1e-10
        assert row['selected_to_auxiliary_complement_error_diagnostic']<1e-10


def test_same_source_current_is_physically_second_quantized_without_species_table():
    path=((1,F(1,11)),(9,F(1,13)),(11,F(1,17)))
    word=lane.FockWord(15);word.current(path,families=1)
    expected=lane.current_path(path)['ideal_matter15']
    assert np.linalg.norm(word.ideal_one_body()-expected,2)<1e-10
    assert lane.current_path(path)['matter15'].shape==(15,15)
    assert word.record()['complete_Fock_dimension']=='32768'


def test_all_twelve_physical_car_current_words_match_the_same_full_exterior_action():
    for index in range(12):
        path=((index,F(1,8)),);word=lane.FockWord(15)
        word.current(path,families=1)
        assert np.linalg.norm(word.ideal_one_body()-lane.current_path(path)['ideal_matter15'],2)<1e-10
        assert word.record()['native_nine_bit_carrier_occupation_wiring_checked']


def test_linked_bond_joins_native_echo_and_complete_matter_current(packet):
    report,arrays=packet
    h=arrays['bond_one_body_generator']
    assert h.shape==(90,90)
    assert np.linalg.norm(arrays['bond_ideal_one_body_word']-expm(-.125j*h),2)<1e-10
    word=report['native_bond_word']
    assert word['physical_CAR_modes']==90
    assert word['complete_Fock_dimension']==str(1<<90)
    assert F(word['unrestricted_Fock_operator_error_upper_conditional'])>0
    assert word['all_occupation_sectors_and_JW_routing_retained']


def test_nonflat_background_has_independent_local_gauge_covariance(packet):
    report,arrays=packet
    h,g=arrays['lattice_hamiltonian'],arrays['whole_site_gauge']
    assert h.shape==(360,360)
    assert np.linalg.norm(h-h.conj().T,2)<1e-12
    assert np.linalg.norm(arrays['gauged_hamiltonian']-g@h@g.conj().T,2)<1e-10
    assert np.linalg.norm(arrays['plaquette']-np.eye(45),2)>1e-3
    assert np.linalg.norm(g@h@g.conj().T-h,2)>1e-3
    assert report['local_covariance']['independent_site_transformations']==8


def test_complete_taste_count_is_not_a_chiral_sm_claim(packet):
    report,_=packet;species=report['complete_species']
    assert species['free_reduced_momentum_spinor_dimension']==8*15*3
    assert species['Dirac_tastes_per_internal_component']==2
    assert species['Dirac_fields_counting_internal_and_family']==90
    assert not species['chiral_SM_projection_performed']
    assert not report['scope']['gauge_field_dynamics_Gauss_constraint_and_gauge_measure_derived']
    assert not report['scope']['chiral_interacting_continuum_or_measured_spectrum_derived']


def test_full_fock_budget_uses_complete_local_errors_and_factor_count():
    word=lane.FockWord(4)
    for _ in range(7):word.append(lane.local_factor('fswap'),(0,1))
    total=7*lane.LOCAL_ADMISSION
    assert word.error_bound()==total/(1-total)
    assert np.linalg.norm(word.complete_matrix()-lane.field.attachment._little_apply(
        lane.local_factor('fswap')['target4'],(0,1),4,np.eye(16)),2)<float(word.error_bound())


@pytest.mark.parametrize('path',[[(-1,F(1))],[(12,F(1))],[(0,.1)],[(0,True)],[(True,F(1))]])
def test_current_path_rejects_ambiguous_or_nonexact_requests(path):
    with pytest.raises(ValueError):lane.current_path(path)


def test_nonadjacent_ungraded_tensor_hopping_is_rejected():
    with pytest.raises(ValueError,match='adjacent CAR'):
        lane.FockWord(4).append(lane.local_factor('hopping',F(1),F(0),F(1,9)),(0,3))


def test_complete_retained_arrays_revalidate_without_replacing_diagnostics(packet):
    report,arrays=packet
    lane.verify_arrays(lane._json_copy(report),arrays)


@pytest.mark.parametrize('tamper',['missing','extra','nonfinite','real_dtype','wrong_double_fill'])
def test_retained_array_custody_rejects_missing_or_changed_complete_factors(packet,tamper):
    report,original=packet;arrays=dict(original)
    key=next(k for k in arrays if k.startswith('bond_factor_') and k.endswith('_actual4'))
    if tamper=='missing':del arrays[key]
    elif tamper=='extra':arrays['extra']=np.zeros((1,1),complex)
    elif tamper=='nonfinite':
        arrays[key]=arrays[key].copy();arrays[key][0,0]=np.nan
    elif tamper=='real_dtype':arrays[key]=arrays[key].real.copy()
    else:
        arrays[key]=arrays[key].copy();arrays[key][3,3]+=1e-4
    with pytest.raises(ValueError):lane.verify_arrays(lane._json_copy(report),arrays)


def test_retained_fock_budget_cannot_be_reduced_to_one_particle(packet):
    report,arrays=packet;bad=lane._json_copy(report)
    bad['native_bond_word']['complete_Fock_dimension']='90'
    with pytest.raises(ValueError,match='Fock domain'):lane.verify_arrays(bad,arrays)
    bad=lane._json_copy(report);bad['native_bond_word']['unrestricted_Fock_operator_error_upper_conditional']='0'
    with pytest.raises(ValueError,match='budget'):lane.verify_arrays(bad,arrays)


def test_typed_report_compare_rejects_boolean_counts_and_nonfinite_diagnostics():
    with pytest.raises(ValueError,match='type'):lane._compare({'count':True},{'count':1})
    with pytest.raises(ValueError):lane._compare({'error':float('nan')},{'error':0.})
    with pytest.raises(ValueError):lane._compare({'count':1,'extra':0},{'count':1})


def test_source_drift_is_rejected_before_array_loading(tmp_path,monkeypatch):
    path=tmp_path/'receipt.json';path.write_text('{"source_bindings": {}}')
    monkeypatch.setattr(lane,'source_bindings',lambda:{'source.py':'new'})
    with pytest.raises(ValueError,match='source bindings'):lane.check(path)
