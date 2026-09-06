"""Independent replay and false-green controls for smooth spatial consistency."""
from copy import deepcopy
import json
from pathlib import Path
import sys

import numpy as np
import pytest
import sympy as s

sys.path.insert(0,str(Path(__file__).resolve().parent))
import whitney_spatial_consistency as producer
import verify_whitney_spatial_consistency as verifier


def test_recorded_refinements_independently_replay():
    summary=verifier.verify(verifier.load())
    assert summary["tetrahedra"]==[20,160,1280]
    assert summary["vertices"]==[13,55,309]
    assert summary["finite_algebra_lean_declarations"]==5
    assert summary["continuum_trajectory_claimed"] is False
    assert summary["uniform_bound_certified_by_numerics"] is False


@pytest.mark.parametrize("mutation,message",[
    (lambda r:r.update(scope="CONTINUUM_PHYSICS_COMPLETE"),"scope"),
    (lambda r:r.update(continuum_trajectory_claimed=True),"promotion"),
    (lambda r:r.update(observer_history_claimed=True),"promotion"),
    (lambda r:r.update(uniform_bound_certified_by_numerics=True),"promotion"),
    (lambda r:r["source_pins"].pop(next(iter(r["source_pins"]))),"pin set"),
    (lambda r:r["source_pins"].update({next(iter(r["source_pins"])):"0"*64}),"source pin"),
    (lambda r:r["refinements"].pop(),"scheduled refinements"),
    (lambda r:r["parameters"]["levels"].__setitem__(0,False),"parameters"),
    (lambda r:r["quadrature"].update(time_legendre_order=4.0),"declared quadrature"),
    (lambda r:r["polynomials"]["Ax"][0].__setitem__(1,0.0),"manufactured fields"),
    (lambda r:r["polynomials"]["Ax"][0].__setitem__(1,False),"manufactured fields"),
    (lambda r:r["polynomials"]["Ax"][0].__setitem__(0,0.3),"manufactured fields"),
    (lambda r:r["refinements"][0].update(level=False),"exact refinement count"),
    (lambda r:r["refinements"][0].update(level=1e-16),"exact refinement count"),
    (lambda r:r["refinements"][0].update(vertices=13.0),"exact refinement count"),
    (lambda r:r["refinements"][0].update(action=float("nan")),"finite numeric"),
    (lambda r:r["refinements"][0].update(action=-100.),"independent quadrature"),
    (lambda r:r["refinements"][0].update(first_variation=100.),"independent quadrature"),
    (lambda r:r["refinements"][0].update(volume=1.),"independent geometry"),
])
def test_receipt_mutations_fail(mutation,message):
    receipt=deepcopy(verifier.load())
    mutation(receipt)
    with pytest.raises(ValueError,match=message):
        verifier.verify(receipt)


@pytest.mark.parametrize("old,new",[
    ('"continuum_trajectory_claimed":false',
     '"continuum_trajectory_claimed":true,"continuum_trajectory_claimed":false'),
    ('"vertices":13','"vertices":1,"vertices":13'),
])
def test_actual_receipt_duplicate_keys_rejected_before_last_value_wins(tmp_path,old,new):
    original=verifier.OUTPUT.read_text(encoding="utf-8")
    assert old in original
    altered=original.replace(old,new,1)
    # The default parser conceals the contradictory first occurrence.
    assert json.loads(altered)==json.loads(original)
    path=tmp_path/"duplicate-receipt.json"
    path.write_text(altered,encoding="utf-8")
    with pytest.raises(ValueError,match="duplicate JSON key"):
        verifier.load(path)


@pytest.mark.parametrize("constant",["NaN","Infinity","-Infinity"])
def test_nonfinite_json_constants_rejected_by_loader(tmp_path,constant):
    path=tmp_path/"nonfinite.json"
    path.write_text('{"value":'+constant+'}',encoding="utf-8")
    with pytest.raises(ValueError,match="nonfinite JSON constant"):
        verifier.load(path)


def test_quadratic_gradient_distinguishes_original_and_whitney_paths():
    h,x,y,z=s.symbols("h x y z",positive=True)
    barycentric=s.Matrix([1-(x+y+z)/h,x/h,y/h,z/h])
    nodal=[0,h*h,0,0]
    edges=s.Matrix(4,4,lambda i,j:nodal[j]-nodal[i])
    interpolated=edges*barycentric
    original=s.Matrix([x*x-value for value in nodal])
    assert s.simplify((barycentric.T*interpolated)[0])==0
    assert s.simplify((barycentric.T*original)[0]-(x*x-h*x))==0
    assert original!=interpolated
    assert (x*x-h*x).subs(x,h/2)==-h*h/4


def test_omitted_dressing_derivative_fails_independent_action_difference():
    actual=producer.evaluate(0)
    omitted=producer.evaluate(0,omit_dressing_variation=True)
    step=.002
    minus2,minus1,plus1,plus2=[producer.evaluate(0,amplitude=shift)["action"]
                             for shift in (-2*step,-step,step,2*step)]
    derivative=(minus2-8*minus1+8*plus1-plus2)/(12*step)
    assert abs(actual["first_variation"]-derivative)<1e-9
    assert abs(omitted["first_variation"]-derivative)>1e-5


def test_full_dressing_derivative_at_fields_and_covariant_derivatives():
    vertices,cells=producer.original_mesh()
    xyz=vertices[cells[:2]]
    lam=np.array([[.1,.2,.3,.4],[.4,.3,.2,.1]])
    center=producer.reconstructed(.17,xyz,lam)
    step=1e-4
    minus=producer.reconstructed(.17,xyz,lam,-step)
    plus=producer.reconstructed(.17,xyz,lam,step)
    # Explicitly exercise E, B, scalar, temporal and spatial covariant fields.
    for position,variation in zip(range(5),range(5,10),strict=True):
        finite_difference=(plus[position]-minus[position])/(2*step)
        assert np.allclose(center[variation],finite_difference,atol=1e-9,rtol=1e-9)


def test_affine_potential_curl_is_reproduced_independently(monkeypatch):
    vertices,cells=producer.original_mesh()
    xyz=vertices[cells[:3]]
    matrix=np.array([[.2,-.3,.4],[.1,.2,-.2],[-.1,.3,.5]])
    shift=np.array([.2,-.1,.3])
    def affine(prefix,time,points,derivative=None):
        if derivative is not None:
            return np.zeros_like(points)
        return points@matrix.T+shift
    monkeypatch.setattr(producer,"vector",affine)
    coefficients=producer.sample_edges("A",0.,xyz)
    gradients,_=producer.element_geometry(xyz)
    reconstructed=sum(2*coefficients[:,edge,None]*np.cross(gradients[:,i],gradients[:,j])
                      for edge,(i,j) in enumerate(producer.EDGES))
    expected=np.array([matrix[2,1]-matrix[1,2],matrix[0,2]-matrix[2,0],matrix[1,0]-matrix[0,1]])
    assert np.allclose(reconstructed,expected,atol=1e-13,rtol=1e-13)
    assert np.linalg.norm(expected)>.1


def test_nodal_point_sampling_is_unbounded_in_three_dimensional_h1():
    r=s.symbols("r",real=True)
    bump=(1-r*r)**3
    l2=4*s.pi*s.integrate(bump**2*r*r,(r,0,1))
    gradient=4*s.pi*s.integrate(s.diff(bump,r)**2*r*r,(r,0,1))
    assert l2>0 and gradient>0 and bump.subs(r,0)==1
    epsilon=s.symbols("epsilon",positive=True)
    norm_squared=epsilon**3*l2+epsilon*gradient
    assert s.limit(norm_squared,epsilon,0,dir="+")==0
    # The nodal interpolant on any fixed surrounding mesh retains value1,
    # while the input H1 norm tends to zero. No bounded point-value operator.
    assert s.limit(1/s.sqrt(norm_squared),epsilon,0,dir="+")==s.oo


def test_source_reads_are_explicit_utf8_under_windows_default(monkeypatch):
    original=Path.read_text
    def windows_default(path,*args,**kwargs):
        if not args and "encoding" not in kwargs:
            kwargs["encoding"]="cp1252"
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,"read_text",windows_default)
    assert producer.original_mesh()[0].shape==(13,3)
    assert len(list(verifier.independent_meshes()))==3
    assert verifier.polynomial_data()==json.loads(json.dumps(producer.POLYNOMIALS))
