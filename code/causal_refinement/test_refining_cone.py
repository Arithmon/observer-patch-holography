import copy
import importlib.util
from pathlib import Path
import pytest

HERE=Path(__file__).resolve().parent
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    out=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(out)
    return out
v=module("independent_refining_causal_cone",HERE/"verify_refining_cone.py")
p=module("proposed_refining_causal_cone",HERE/"refining_cone.py")

def test_frozen_control_and_reproduction():
    packet=v.load()
    assert packet==p.produce()
    result=v.verify(packet)
    assert result["exact_history_width"]==27
    assert result["native_physical_spacetime_selected"] is False

@pytest.mark.parametrize("field",["writer","value","register"])
def test_provenance_corruption(field):
    packet=p.produce()
    read=packet["history"][27]["reads"][0]
    read[field]=0 if field=="value" else [0,9,9,9]
    with pytest.raises(ValueError):v.verify(packet)

def test_missing_parent():
    packet=p.produce()
    packet["history"][27]["reads"].pop()
    with pytest.raises(ValueError):v.verify(packet)

def test_read_before_write():
    packet=p.produce()
    packet["history"][0],packet["history"][27]=packet["history"][27],packet["history"][0]
    with pytest.raises(ValueError,match="committed prefix"):v.verify(packet)

@pytest.mark.parametrize("index",range(16))
def test_every_path_boundary_is_checked(index):
    packet=p.produce()
    packet["cone_paths"][index]["path"][1][0]+=2
    with pytest.raises(ValueError):v.verify(packet)

def test_wrong_feedback():
    packet=p.produce()
    packet["history"][28]["write_value"]+=1
    with pytest.raises(ValueError,match="consumed records"):v.verify(packet)

def test_no_physical_promotion():
    packet=p.produce()
    packet["scope"]["native_OPH_manifold"]=True
    with pytest.raises(ValueError):v.verify(packet)

@pytest.mark.parametrize("payload",['{"a":1,"a":2}','{"a":NaN}','{"a":1.1}'])
def test_strict_json(tmp_path,payload):
    path=tmp_path/"bad.json"
    path.write_text(payload,encoding="utf-8")
    with pytest.raises(ValueError):v.load(path)

def test_boolean_radius_is_not_an_integer():
    packet=p.produce()
    packet["cone_paths"][0]["radius"]=True
    with pytest.raises(ValueError):v.verify(packet)
