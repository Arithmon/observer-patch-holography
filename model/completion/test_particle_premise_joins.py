"""Cross-component identities and scope boundaries of the same-action push."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from model.bound_pair.fiber import PairFiber,clock_subspace
from model.bound_pair.mass import curvature_ratio_interval
from model.interacting_field.lattice import symbol
from model.particle_clock import clock

ROOT=Path(__file__).resolve().parents[2]


def test_selected_issues_and_pending_prs_do_not_accept_or_merge_physics():
    c=json.loads((ROOT/'model/completion/particle_premise_campaign.json').read_text())
    assert c['selected_issues']==[1139,1140,1142,1137]
    assert not c['integrated_gate_accepted'] and c['physical_comparison_count_added']==0
    assert not c['whole_program_frozen']
    assert {r['number'] for r in c['pending_complementary_prs']}=={1153,1128}
    assert all(r['status']=='DRAFT_NOT_MERGED' for r in c['pending_complementary_prs'])
    assert not c['common_action']['infinite_volume_is_spatial_continuum']


def test_particle_fiber_uses_original_field_symbol_not_a_fitted_surrogate():
    f=PairFiber([0,np.pi,np.pi],side=4)
    for k,block in zip(f.k,f.blocks,strict=True):
        independent=np.kron(symbol(list(k)),np.eye(2))+np.kron(np.eye(2),symbol(list(f.momentum-k)))
        assert np.max(abs(block-independent))<1e-15


def test_same_number_energy_reference_and_native_corner_gibbs_map():
    s=clock_subspace()
    h=s['restriction'];beta=np.log(25/8)
    thermal=expm(-beta*h);thermal/=np.trace(thermal)
    assert np.max(abs(thermal-np.diag([25/33,8/33])))<1e-15
    assert max(s['residuals'].values())<1e-12
    receipt=json.loads((ROOT/'model/particle_clock/results.json').read_text())
    assert receipt['field_dictionary']['same_field_source_module']=='model/bound_pair/fiber.py'
    assert receipt['source_sha256']['model/bound_pair/fiber.py']==hashlib.sha256(
        (ROOT/'model/bound_pair/fiber.py').read_bytes()).hexdigest()
    assert not receipt['error_bound_conditioned_on_herald_or_training_acceptance']
    assert F(receipt['complete_average_half_trace_upper'])<F(19,1000)
    assert receipt['field_dictionary']['supplied_source_to_field_isometry']
    assert not receipt['field_dictionary']['faithful_on_full_field_or_local_regional_algebras']


def test_common_unit_change_cannot_repair_curvature_rest_energy_mismatch():
    lo,hi=curvature_ratio_interval()
    assert F(17,16)<lo<hi<F(14,13)
    # Independent arbitrary clock/rod conversion on the certified endpoints.
    for tau,ell in ((F(7,3),F(5,2)),(F(1,100),F(11))):
        energy=F(1)/tau; speed=F(1,8)*ell/tau
        for ratio in (lo,hi):
            curvature=ratio*F(1,64)*ell**2/tau
            assert energy*curvature/speed**2==ratio


def test_receipt_does_not_assign_species_or_topological_defect():
    d=json.loads((ROOT/'model/bound_pair/results.json').read_text())
    for key in ('physical_species_identified','nontrivial_topological_defect_assigned',
                'full_native_field_encoding_executed','integrated_gate_accepted'):
        assert d[key] is False
    assert d['physical_comparison_count_added']==0
    for path,sha in d['source_bindings'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
