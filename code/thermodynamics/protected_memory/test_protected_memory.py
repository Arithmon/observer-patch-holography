"""Adversarial exact replay of conserved records with finite memory tails."""
import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path

import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("protected_memory_independent_tests",HERE/"verify_protected_memory.py")
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)


@pytest.fixture(scope="module")
def packet():return v.load()


def test_whole_exact_replay(packet):
    result=v.verify(packet)
    assert result["accepted"] is True
    assert result["protected_fibres"]==2
    assert F(result["two_step_minorization"])>0
    assert result["native_source_attachment"] is False
    assert result["historical_obstruction_overturned"] is False


@pytest.mark.parametrize("key",["native_source_attachment","historical_obstruction_overturned",
    "global_ergodicity","physical_clock_or_conductivity","continuum_limit","empirical_evidence","new_Lean_formalization"])
def test_no_physical_or_historical_promotion(packet,key):
    p=copy.deepcopy(packet);p["scope"][key]=True
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("kind",["reference","local_kernel","transition","coupling","equilibrium","record_leak"])
def test_actual_law_controls(packet,kind):
    p=copy.deepcopy(packet)
    if kind=="reference":p["reference"][0]="1/3"
    elif kind=="local_kernel":p["coordinate_kernels"]["u"][0][1]="1/9"
    elif kind=="transition":p["transition"]=copy.deepcopy(p["equilibrium_projection"])
    elif kind=="coupling":p["stationary_pair_coupling"][0][2]="0"
    elif kind=="equilibrium":p["equilibrium_projection"]=copy.deepcopy(p["transition"])
    else:p["transition"][0][4]="1/8"
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("kind",["epsilon","gap","minorization","inverse","solution","current","global_centering"])
def test_nontrivial_algebra_controls(packet,kind):
    p=copy.deepcopy(packet)
    if kind=="epsilon":p["two_step_minorization"]["epsilon"]="1"
    elif kind=="gap":p["two_step_minorization"]["gap_lower"]="1/2"
    elif kind=="minorization":p["two_step_minorization"]["remainder"][0][0]="-1"
    elif kind=="inverse":p["poisson_resolvent"][0][0]="0"
    elif kind=="solution":p["poisson_solutions"][0]=["0"]*8
    elif kind=="current":p["currents"][0]=["0"]*8
    else:p["currents"][0]=copy.deepcopy(p["global_centering_counterexample"]["current"])
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("kind",["matrix","lag","partial","remainder","bound","history","persistent"])
def test_memory_and_tail_controls(packet,kind):
    p=copy.deepcopy(packet);g=p["green_kubo"]
    if kind=="matrix":g["matrix"][0][1]="0"
    elif kind=="lag":g["correlations"][1][0][0]="0"
    elif kind=="partial":g["partial_sum"][0][0]="0"
    elif kind=="remainder":g["exact_remainder"][0][0]="0"
    elif kind=="bound":g["tail_upper"][0][0]="0"
    elif kind=="history":p["distribution_history"][-1]["distribution"][0]="0"
    else:p["global_centering_counterexample"]["persistent_correlation"]="0"
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("value",[False,0.0,"0",-1])
def test_state_categories_are_exact(packet,value):
    p=copy.deepcopy(packet);p["model"]["states"][0][0]=value
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("value",[True,1,"2/60","nan","1/0","0.0333333333333333"])
def test_canonical_rational_values(packet,value):
    p=copy.deepcopy(packet);p["reference"][0]=value
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("kind",["source","parent","missing_lag","count_bool","read_footprint"])
def test_custody_and_complete_inventory(packet,kind):
    p=copy.deepcopy(packet)
    if kind=="source":p["source_pins"]["protected_memory.py"]["sha256"]="0"*64
    elif kind=="parent":p["parent_pins"]["Lean/Thermodynamics/GreenKubo.lean"]["bytes"]+=1
    elif kind=="missing_lag":p["green_kubo"]["correlations"].pop()
    elif kind=="count_bool":p["green_kubo"]["cutoff"]=True
    else:p["model"]["operations"][1]["reads"]=["c"]
    with pytest.raises(ValueError):v.verify(p)


@pytest.mark.parametrize("text",['{"a":1,"a":2}','{"x":{"a":1,"a":2}}','{"a":NaN}',
    '{"a":Infinity}','{"a":1e999}','{"a":0.25}'])
def test_loader_rejects_ambiguous_or_nonexact_numbers(tmp_path,text):
    p=tmp_path/"invalid.json";p.write_text(text,encoding="utf-8")
    with pytest.raises(ValueError):v.load(p)


def test_global_centering_is_insufficient_independently(packet):
    pi=[F(x) for x in packet["reference"]]
    T=[[F(x) for x in row] for row in packet["transition"]]
    f=[F(-2,3)]*4+[F(1,3)]*4
    assert sum(p*x for p,x in zip(pi,f))==0
    assert [sum(t*x for t,x in zip(row,f)) for row in T]==f
    assert sum(p*x*x for p,x in zip(pi,f))==F(2,9)


def test_conditionals_are_distinct_nonproduct_laws(packet):
    weights=packet["model"]["conditional_reference_weights"]
    assert weights[0]!=weights[1]
    assert all(w[0]*w[3]!=w[1]*w[2] for w in weights)
