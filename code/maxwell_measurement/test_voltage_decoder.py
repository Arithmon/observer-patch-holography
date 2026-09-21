"""Finite interface checks only; no device readings or measured outcome fixture."""
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("maxwell_voltage_decoder",HERE/"voltage_decoder.py")
decoder=importlib.util.module_from_spec(spec)
spec.loader.exec_module(decoder)


def test_complete_immutable_alphabet_and_physical_margin():
    book=decoder.codebook()
    values=list(map(Q,book["symbols"]))
    assert len(values)==123 and len(set(values))==123
    assert min(b-a for a,b in zip(values,values[1:]))==Q(23,42896)
    # 0.5 V per model unit: nearest-neighbour midpoint is 23/171584 V away.
    assert Q(book["maximum_endpoint_distance_from_symbol"])/2 < Q(23,171584)
    for value in values:
        result=decoder.decode_interval(str(value-Q(1,20000)),str(value+Q(1,20000)),book)
        assert Q(result["symbol"])==value and result["measurement_verdict"] is False


def test_no_nearest_reference_snapping():
    book=decoder.codebook()
    value=Q(book["symbols"][10])
    with pytest.raises(ValueError,match="no alphabet symbol"):
        decoder.decode_interval(str(value+Q(1,100000)),str(value+Q(1,50000)),book)


def test_measured_wrong_symbol_cannot_be_forced_to_prediction():
    book=decoder.codebook()
    with pytest.raises(ValueError,match="disagrees"):
        decoder.decode_interval("1/2","1/2",book,"-1/2")


def test_ambiguous_interval_fails_even_with_expected_value():
    book=decoder.codebook()
    values=list(map(Q,book["symbols"]))
    left,right=min(zip(values,values[1:]),key=lambda pair:pair[1]-pair[0])
    with pytest.raises(ValueError,match="ambiguous"):
        decoder.decode_interval(str(left),str(right),book,str(left))


@pytest.mark.parametrize("interval,message",[
    (("1/100","-1/100"),"reversed"),
    (("34135/10724",str(Q(34135,10724)+Q(1,1000))),"too wide"),
    (("-3/20000","0"),"code margin"),
    (("999","1000"),"no alphabet symbol"),
    (("0.0","0"),"canonical"),
    ((False,"0"),"rational string"),
])
def test_bad_intervals(interval,message):
    with pytest.raises(ValueError,match=message):
        decoder.decode_interval(*interval,decoder.codebook())


@pytest.mark.parametrize("change",["extra","missing","reordered","radius","parent","source","boolean"])
def test_rehashed_codebook_changes_fail(change):
    book=copy.deepcopy(decoder.codebook())
    if change=="extra": book["symbols"].append("123")
    if change=="missing": book["symbols"].pop()
    if change=="reordered": book["symbols"].reverse()
    if change=="radius": book["maximum_endpoint_distance_from_symbol"]="1"
    if change=="parent": book["parent"]["sha256"]="0"*64
    if change=="source": book["source_pins"]["code/maxwell_measurement/voltage_decoder.py"]="0"*64
    if change=="boolean": book["uses_measured_validation_data"]=0
    with pytest.raises(ValueError,match="changed frozen"):
        decoder.decode_interval("0","0",book)


@pytest.mark.parametrize("raw",[b'{"x":0,"x":1}',b'{"x":0.5}',b'{"x":NaN}'])
def test_ambiguous_json(raw):
    with pytest.raises(ValueError): decoder.strict_json(raw)


def test_cli_preparation_is_not_experiment_freeze(tmp_path):
    path=tmp_path/"codebook.json"
    command=[sys.executable,str(HERE/"voltage_decoder.py"),"prepare-codebook","--output",str(path)]
    result=subprocess.run(command,capture_output=True,text=True)
    assert result.returncode==0 and '"experiment_preregistered": false' in result.stdout
    assert subprocess.run(command,capture_output=True).returncode==2


def test_cli_external_pin_and_output_bind_same_codebook_bytes(tmp_path):
    path=tmp_path/"codebook.json"
    raw=decoder.canonical(decoder.codebook())
    path.write_bytes(raw)
    pin=hashlib.sha256(raw).hexdigest()
    command=[sys.executable,str(HERE/"voltage_decoder.py"),"decode","--codebook",str(path),
             "--codebook-sha256",pin,"--lower","0","--upper","0"]
    result=subprocess.run(command,capture_output=True,text=True)
    assert result.returncode==0
    assert json.loads(result.stdout)["codebook_sha256"]==pin
    command[6]="0"*64
    result=subprocess.run(command,capture_output=True,text=True)
    assert result.returncode==2 and "external codebook pin mismatch" in result.stdout
    raw=json.dumps(decoder.codebook(),indent=2).encode()
    path.write_bytes(raw)
    command[6]=hashlib.sha256(raw).hexdigest()
    result=subprocess.run(command,capture_output=True,text=True)
    assert result.returncode==2 and "canonical codebook bytes" in result.stdout
