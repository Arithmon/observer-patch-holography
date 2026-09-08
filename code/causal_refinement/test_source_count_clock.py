"""Exact bounds and forged-record controls for retrospective count clocks."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import random
import types

import pytest

HERE = Path(__file__).resolve().parent


def module(name):
    path = HERE / (name+".py")
    obj = types.ModuleType(name)
    obj.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), obj.__dict__)
    return obj


producer = module("source_count_clock")
verifier = module("verify_source_count_clock")


@pytest.fixture(scope="module")
def source():
    return producer.source_fixture()


def digest(trace):
    chain = bytes(32)
    for item in trace:
        raw = (json.dumps(item, sort_keys=True, separators=(",", ":"),
                          ensure_ascii=True, allow_nan=False)+"\n").encode("ascii")
        chain = sha256(chain+raw).digest()
    return chain.hex()


def test_fresh_receipt():
    stored = verifier.load()
    assert verifier.verify(stored) == verifier.verify(producer.report())


def test_exact_root_random_rationals():
    rng = random.Random(805729)
    for _ in range(100):
        value = F(rng.randrange(10**15), rng.randrange(1,10**10))
        lo, hi = producer.fourth_root_bounds(value)
        assert lo**4 <= value <= hi**4
        assert hi-lo <= F(1, 1 << 80)
        assert (lo*(1 << 80)).denominator == (hi*(1 << 80)).denominator == 1


@pytest.mark.parametrize("root", [F(0),F(1),F(3,16),F(11),F(1,1<<40)])
def test_exact_fourth_powers(root):
    assert producer.fourth_root_bounds(root**4) == (root,root)


@pytest.mark.parametrize("value", [True,False,1.0,-1,"-1",None,[],{},float("inf")])
def test_invalid_root_inputs(value):
    with pytest.raises((ValueError, TypeError)):
        producer.fourth_root_bounds(value)


@pytest.mark.parametrize("bits", [True,0,-1,257,80.0,"80"])
def test_invalid_precision(bits):
    with pytest.raises(ValueError):
        producer.fourth_root_bounds(2,bits)


@pytest.mark.parametrize("n,m", [(0,1),(1,0),(-1,1),(1,True),(2.0,3)])
def test_invalid_counts(n,m):
    with pytest.raises(ValueError):
        producer.count_clock(n,m)


@pytest.mark.parametrize("en,em", [(-1,0),(0,-1),(0,81),(0,82),(True,0),(0,1.0)])
def test_error_certificate_fails_closed(en,em):
    with pytest.raises(ValueError):
        producer.certified_ratio_bounds(16,81,en,em)


def test_finite_error_enclosure_extremes():
    for n,m,en,em in [(16,81,1,2),(1,16,2,1),(16,81,0,0)]:
        lo,hi = producer.certified_ratio_bounds(n,m,en,em)
        for v in (max(0,n-en),n,n+en):
            for w in (m-em,m,m+em):
                assert lo**4 <= F(v,w) <= hi**4


def test_ancestry_counts_match_independent_parent(source):
    trace,order,row,counts = source
    assert len(order) == 500 and counts == [2,41,80]
    assert digest(trace) == row["forward_execution"]["audit_trace_sha256"]
    assert producer.interval_count(order,(0,0),(0,0)) == 1
    with pytest.raises(ValueError):
        producer.interval_count(order,(1,0),(0,0))


@pytest.mark.parametrize("anchor", [(False,0),(0.0,0),(0,True),(0,-1),(0,),"00",None])
def test_invalid_query_anchors(source,anchor):
    _,order,_,_ = source
    with pytest.raises(ValueError):
        producer.interval_count(order,anchor,(1,0))
    with pytest.raises(ValueError):
        producer.interval_count(order,(0,0),anchor)


def test_reversed_within_layer_schedule(source):
    trace,order,row,_ = source
    reverse = sorted(trace, key=lambda e: (e[0][0],-e[0][1]))
    expected = row["within_layer_reversed_execution"]["audit_trace_sha256"]
    assert producer.authenticated_order(reverse,expected) == order


def test_opaque_labels_need_no_numeric_clock(source):
    trace,order,_,_ = source
    other = deepcopy(trace)
    rename = lambda event: (10000-17*event[0],20000-13*event[1])
    for item in other:
        item[0] = list(rename(item[0]))
        for read in item[1]:
            read[2] = list(rename(read[2]))
        item[2][2] = list(rename(item[2][2]))
    result = producer.authenticated_order(other,digest(other))
    assert result == {rename(k): frozenset(rename(p) for p in v) for k,v in order.items()}


@pytest.mark.parametrize("mutation", ["value","writer","version","boolean_writer","duplicate_read","write_version","write_writer","duplicate_event"])
def test_resealed_invalid_records_fail(source,mutation):
    trace,_,_,_ = source
    changed = deepcopy(trace)
    event = changed[125]
    if mutation == "value": event[1][0][3] += 1
    elif mutation == "writer": event[1][0][2] = [999,999]
    elif mutation == "version": event[1][0][1] += 1
    elif mutation == "boolean_writer": event[1][0][2][0] = False
    elif mutation == "duplicate_read": event[1].append(deepcopy(event[1][0]))
    elif mutation == "write_version": event[2][1] += 1
    elif mutation == "write_writer": event[2][2][0] = True
    elif mutation == "duplicate_event": changed.append(deepcopy(changed[0]))
    with pytest.raises(ValueError):
        producer.authenticated_order(changed,digest(changed))


def test_missing_real_parent_cannot_use_source_commitment(source):
    trace,_,row,_ = source
    changed = deepcopy(trace)
    changed[125][1].pop()
    with pytest.raises(ValueError,match="commitment"):
        producer.authenticated_order(changed,row["forward_execution"]["audit_trace_sha256"])
    # A different valid history can have a different commitment. Hashing alone
    # does not establish that the declared read law was executed.
    producer.authenticated_order(changed,digest(changed))


def test_add_audit_chain_changes_semantic_order(source):
    trace,_,row,_ = source
    changed = deepcopy(trace)
    changed[1][1].append([0,1,[0,0],1])
    with pytest.raises(ValueError,match="commitment"):
        producer.authenticated_order(changed,row["forward_execution"]["audit_trace_sha256"])
    order = producer.authenticated_order(changed,digest(changed))
    assert (0,0) in order[(0,1)]


def test_raw_count_is_not_inertially_additive():
    whole = producer.count_clock(16,1)[0]
    half = producer.count_clock(1,1)[0]
    assert whole == 2*half
    assert 16 != 2*1


@pytest.mark.parametrize("mutation", ["count","floatcount","boolcount","root","extra","scope","signature","parenthash","tracehash","error","interior","precision"])
def test_receipt_forgeries_rejected(mutation):
    item = verifier.load()
    if mutation == "count": item["interval_counts"][1] += 1
    elif mutation == "floatcount": item["authenticated_events"] = 500.0
    elif mutation == "boolcount": item["interval_counts"][0] = True
    elif mutation == "root": item["ratio_to_first_interval_bounds"][1] = ["2","2"]
    elif mutation == "extra": item["undocumented"] = True
    elif mutation == "scope": item["scope"]["finite_clock_accuracy_certified"] = True
    elif mutation == "signature": item["scope"]["external_signature_verified"] = True
    elif mutation == "parenthash": item["source_receipt_sha256"] = "0"*64
    elif mutation == "tracehash": item["source_trace_sha256"] = "0"*64
    elif mutation == "error": item["conditional_error_controls"][0]["count_error_bounds"][0] = "0"
    elif mutation == "interior": item["continuum_diamond_inside_source_window"][2] = True
    elif mutation == "precision": item["dyadic_bits"] = True
    with pytest.raises(ValueError):
        verifier.verify(item)


@pytest.mark.parametrize("raw", ['{"a":1,"a":2}', '{"a":1.0}', '{"a":1e9999}', '{"a":NaN}', '{"a":Infinity}'])
def test_strict_json_parser(tmp_path,raw):
    path = tmp_path / "forged.json"
    path.write_bytes(raw.encode("utf-8"))
    with pytest.raises(ValueError):
        verifier.load(path)


def test_independent_verifier_never_imports_producer(monkeypatch):
    original = verifier.fresh
    called = []
    def check(path):
        called.append(path.name)
        assert path.name == "verify_source_net_causet.py"
        return original(path)
    monkeypatch.setattr(verifier,"fresh",check)
    verifier.verify(verifier.load())
    assert called == ["verify_source_net_causet.py"]


def test_fresh_import_bypasses_same_size_source_cache(tmp_path):
    path = tmp_path / "parent.py"
    path.write_bytes(b"VALUE = 1\n")
    assert verifier.fresh(path).VALUE == 1
    path.write_bytes(b"VALUE = 2\n")
    assert verifier.fresh(path).VALUE == 2
