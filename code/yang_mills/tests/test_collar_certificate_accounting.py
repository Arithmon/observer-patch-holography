"""Original-input controls for the two finite collar certificate interfaces."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import finite_collar_gap_certificate as calibration
import verify_collar_gap_certificate as witness


def entry():
    return {
        "id": "a",
        "rooted_neighborhood": "one declared binary root",
        "boundary_sector_flags": [],
        "local_alphabets": ["binary"],
        "potential_templates": [],
        "complete_repaired_readback": "all other declared registers",
        "conditional_kernel": "supplied binary row pairs",
        "rate_lower": "3/4",
        "influences": [],
        "refinement_targets": ["a"],
    }


def influence(bound="1/4", multiplicity=1):
    # Independent binary TV: distance((1,0),(1-b,b)) = b exactly.
    b = Q(bound)
    return {"target_type": "a", "multiplicity": multiplicity, "upper": bound,
            "conditional_rows": [["1", "0"], [str(1-b), str(b)]]}


def payload(kind, entries=None):
    entries = copy.deepcopy(entries if entries is not None else [entry()])
    if kind == "witness":
        return {"schema": "oph.yang_mills.collar_gap.v1",
                "scope": "theorem_contract_witness", "types": entries}
    if kind == "explicit":
        return {"schema": calibration.SCHEMA, "scope": "calibration", "type_table": entries}
    assert kind == "family" and len(entries) == 1
    template = entries[0]
    del template["id"], template["refinement_targets"]
    for item in template["influences"]:
        del item["target_type"]
    return {"schema": calibration.SCHEMA, "scope": "calibration",
            "calibration_type_family": {"count": 3, "id_prefix": "a", "template": template}}


def check(kind, data):
    return witness.verify(data) if kind == "witness" else calibration.validate(data)


KINDS = ("witness", "explicit", "family")


@pytest.mark.parametrize("bounds", [("1/8", "1/4"), ("1/4", "1/8"), ("1/8",)*3])
@pytest.mark.parametrize("kind", KINDS)
def test_every_declared_influence_contributes(kind, bounds):
    e = entry()
    e["influences"] = [influence(b) for b in bounds]
    data = payload(kind, [e])
    before = copy.deepcopy(data)
    result = check(kind, data)
    eta = sum(map(Q, bounds))
    assert Q(result["eta_upper"]) == eta
    assert Q(result["gap_lower"]) == Q(3, 4) * (1-eta)
    assert all(Q(v) == eta for v in result["row_sums"].values())
    assert data == before


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("bounds", [("3/4", "1/4"), ("3/4", "1/2"), ("1/2", "3/4")])
def test_all_influences_decide_contractivity(kind, bounds):
    e = entry()
    e["influences"] = [influence(b) for b in bounds]
    with pytest.raises(ValueError):
        check(kind, payload(kind, [e]))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("bad", [True, False, 0.9, -0.9, 1.5, "2"])
def test_multiplicity_is_an_exact_nonnegative_integer(kind, bad):
    e = entry()
    e["influences"] = [influence("1/4", bad)]
    with pytest.raises(ValueError):
        check(kind, payload(kind, [e]))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("field", ["rate_lower", "upper", "probability"])
def test_boolean_is_not_exact_rational_evidence(kind, field):
    e = entry()
    e["influences"] = [influence("0")]
    if field == "rate_lower":
        e[field] = True
    elif field == "upper":
        e["influences"][0][field] = False
    else:
        e["influences"][0]["conditional_rows"] = [[True, False], [True, False]]
    with pytest.raises(ValueError):
        check(kind, payload(kind, [e]))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("extra", [{}, {"physical_source_provenance": None, "continuum_receipts": None},
                                  {"physical_source_provenance": {"verified": True},
                                   "continuum_receipts": ["asserted"]}])
def test_metadata_cannot_promote_a_fixture_to_a_physical_receipt(kind, extra):
    data = payload(kind)
    data.update(extra, scope="physical_source_receipt")
    with pytest.raises(ValueError):
        check(kind, data)


@pytest.mark.parametrize("kind,field,value", [
    (kind, field, value) for kind in KINDS
    for field, value in [("rooted_neighborhood", None), ("conditional_kernel", ""),
                         ("local_alphabets", []), ("influences", {}),
                         ("refinement_targets", "a")]
    if not (kind == "family" and field == "refinement_targets")
])
def test_malformed_source_fields_do_not_certify(kind, field, value):
    e = entry()
    e[field] = value
    with pytest.raises(ValueError):
        check(kind, payload(kind, [e]))


@pytest.mark.parametrize("bad", [True, False, 1.5, "3"])
def test_family_count_is_an_exact_positive_integer(bad):
    data = payload("family")
    data["calibration_type_family"]["count"] = bad
    with pytest.raises(ValueError):
        calibration.validate(data)


def test_family_expansion_does_not_alias_independent_types():
    data = payload("family")
    expanded = calibration.expand_types(data)
    expanded[0]["local_alphabets"].append("unexpected")
    assert expanded[1]["local_alphabets"] == ["binary"]
    assert data["calibration_type_family"]["template"]["local_alphabets"] == ["binary"]


def test_explicit_and_compact_tables_cannot_shadow_each_other():
    data = payload("family")
    data["type_table"] = [entry()]
    with pytest.raises(ValueError):
        calibration.validate(data)


@pytest.mark.parametrize("kind", ["witness", "explicit"])
def test_cli_rejects_duplicate_json_evidence_before_writing(tmp_path, kind):
    data = payload(kind)
    text = json.dumps(data).replace('"influences": []',
                                  '"influences": [null], "influences": []')
    path = tmp_path / "ambiguous.json"
    path.write_text(text, encoding="utf-8")
    output = tmp_path / "result.json"
    args = ([str(HERE / "verify_collar_gap_certificate.py"), str(path)] if kind == "witness"
            else [str(HERE / "finite_collar_gap_certificate.py"), "certify",
                  "--manifest", str(path), "--output", str(output)])
    run = subprocess.run([sys.executable, *args], capture_output=True, text=True)
    assert run.returncode != 0
    assert "duplicate JSON key" in run.stderr
    assert not output.exists()


@pytest.mark.parametrize("kind", KINDS)
def test_zero_influence_and_nonunit_rate_are_supported(kind):
    result = check(kind, payload(kind))
    assert result["physical_clay_receipt"] is False
    assert Q(result["gap_lower"]) == Q(3, 4)
    assert Q(result["approximate_tensorization_upper"]) == 1
