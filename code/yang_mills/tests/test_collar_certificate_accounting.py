"""Original-input controls for the two finite collar certificate interfaces."""

from __future__ import annotations

import copy
import itertools
import json
import math
import random
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


@pytest.mark.parametrize("kind", ["witness", "explicit"])
def test_original_multistate_laws_against_event_distance(kind):
    """Independent TV oracle: maximize probability differences over all events."""
    rng = random.Random(306)
    accepted = refused = 0
    for trial in range(240):
        entries = [entry() for _ in range(1 + trial % 3)]
        expected_rows, rates = {}, []
        for index, e in enumerate(entries):
            e["id"] = str(index)
            e["refinement_targets"] = [str(j) for j in range(len(entries))]
            rate = Q(1 + rng.randrange(20), 1 + rng.randrange(10))
            e["rate_lower"] = str(rate)
            rates.append(rate)
            total = Q(0)
            for number in range(trial % 4):
                size = 2 + (trial + number) % 3
                a = [rng.randrange(1, 8) for _ in range(size)]
                b = [rng.randrange(1, 8) for _ in range(size)]
                left, right = [[Q(v, sum(row)) for v in row] for row in (a, b)]
                distance = max(abs(sum((left[j]-right[j] for j in range(size) if mask >> j & 1), Q(0)))
                               for mask in range(1 << size))
                bound = distance + Q(trial % 2, 20)
                count = rng.randrange(4)
                e["influences"].append({"target_type": str((index+1) % len(entries)),
                                        "multiplicity": count, "upper": str(bound),
                                        "conditional_rows": [[str(v) for v in row] for row in (left, right)]})
                total += count * bound
            expected_rows[e["id"]] = total
        eta = max(expected_rows.values())
        if eta >= 1:
            with pytest.raises(ValueError, match="< 1"):
                check(kind, payload(kind, entries))
            refused += 1
        else:
            result = check(kind, payload(kind, entries))
            assert {key: Q(value) for key, value in result["row_sums"].items()} == expected_rows
            assert Q(result["c_floor"]) == min(rates)
            assert Q(result["eta_upper"]) == eta
            assert Q(result["gap_lower"]) == min(rates)*(1-eta)
            assert Q(result["approximate_tensorization_upper"]) == 1/(1-eta)
            accepted += 1
    assert accepted > 50 and refused > 50


@pytest.mark.parametrize("kind", ["witness", "explicit"])
@pytest.mark.parametrize("correlation", [Q(0), Q(1, 4), Q(-3, 4), 1-Q(1, 2**100), Q(1, 2**100)-1])
@pytest.mark.parametrize("rate", [Q(3, 7), Q(2**1100), Q(1, 2**1100)])
def test_bound_attained_by_an_explicit_two_spin_heat_bath(kind, correlation, rate):
    """Construct the full operator, then check a complete exact eigenbasis."""
    states = list(itertools.product((-1, 1), repeat=2))
    law = {(s, t): (1+correlation*s*t)/4 for s, t in states}
    generator = [[Q(0) for _ in states] for _ in states]
    entries = []
    for site in range(2):
        e = entry()
        e.update(id=str(site), refinement_targets=[str(site)], rate_lower=str(rate))
        conditional = []
        for other in (-1, 1):
            pair = [(v, other) if site == 0 else (other, v) for v in (-1, 1)]
            mass = sum(law[state] for state in pair)
            conditional.append([law[state]/mass for state in pair])
        distance = abs(conditional[0][0]-conditional[1][0])
        e["influences"] = [{"target_type": str(1-site), "upper": str(distance),
                             "conditional_rows": [[str(v) for v in row] for row in conditional]}]
        entries.append(e)
        for i, state in enumerate(states):
            generator[i][i] += rate
            compatible = [j for j, target in enumerate(states) if target[1-site] == state[1-site]]
            mass = sum(law[states[j]] for j in compatible)
            for j in compatible:
                generator[i][j] -= rate*law[states[j]]/mass
    vectors = [[Q(1) for _ in states], [s+t for s, t in states],
               [s-t for s, t in states], [s*t-correlation for s, t in states]]
    eigenvalues = [Q(0), rate*(1-correlation), rate*(1+correlation), 2*rate]
    # Four independent vectors exclude an untested smaller eigenvalue.
    determinant = sum((-1)**sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4))
                      * math.prod(vectors[i][p[i]] for i in range(4))
                      for p in itertools.permutations(range(4)))
    assert determinant != 0
    for vector, eigenvalue in zip(vectors, eigenvalues, strict=True):
        assert [sum(a*b for a, b in zip(row, vector, strict=True)) for row in generator] == [
            eigenvalue*v for v in vector]
    result = check(kind, payload(kind, entries))
    assert Q(result["gap_lower"]) == min(eigenvalues[1:])


@pytest.mark.parametrize("kind", KINDS)
def test_large_multiplicity_and_small_positive_gap_remain_exact(kind):
    e = entry()
    e["influences"] = [influence(str(Q(1, 2**100)), 2**100-1)]
    result = check(kind, payload(kind, [e]))
    assert Q(result["gap_lower"]) == Q(3, 2**102)
    assert Q(result["approximate_tensorization_upper"]) == 2**100
    e["influences"][0]["multiplicity"] += 1
    with pytest.raises(ValueError, match="< 1"):
        check(kind, payload(kind, [e]))


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("bad", ["understated_tv", "negative", "unnormalized", "empty_rows", "bad_target"])
def test_invalid_rows_are_rejected_even_at_zero_multiplicity(kind, bad):
    e = entry()
    item = influence("1/4", 0)
    if bad == "understated_tv":
        item["upper"] = "1/8"
    elif bad == "negative":
        item["conditional_rows"] = [["-1", "2"], ["-1", "2"]]
    elif bad == "unnormalized":
        item["conditional_rows"] = [["1", "1"], ["1", "1"]]
    elif bad == "empty_rows":
        item["conditional_rows"] = [[], []]
    else:
        item["target_type"] = "missing"
    e["influences"] = [item]
    data = payload(kind, [e])
    if kind == "family" and bad == "bad_target":
        # A compact self-type must not silently overwrite an explicit target.
        data["calibration_type_family"]["template"]["influences"][0]["target_type"] = "missing"
    with pytest.raises(ValueError):
        check(kind, data)


@pytest.mark.parametrize("key,bad", [("physical_clay_receipt", 0), ("active_type_count", True),
                                    ("eta_upper", "1/4"), ("type_hashes", {}),
                                    ("manifest_sha256", "0"*64)])
def test_replay_checks_types_and_all_bindings(tmp_path, key, bad):
    data = payload("explicit")
    source, receipt = tmp_path / "source.json", tmp_path / "receipt.json"
    source.write_text(json.dumps(data), encoding="utf-8")
    calibration.certify(source, receipt)
    calibration.verify(source, receipt)
    retained = json.loads(receipt.read_text(encoding="utf-8"))
    retained[key] = bad
    receipt.write_text(json.dumps(retained), encoding="utf-8")
    with pytest.raises(ValueError, match="exactly recompute"):
        calibration.verify(source, receipt)


def test_original_calibration_receipt_reproduces_byte_for_byte(tmp_path):
    source = HERE / "manifests" / "atomic_4d_ising_calibration.json"
    retained = HERE / "receipts" / "atomic_4d_ising_calibration.receipt.json"
    rebuilt = tmp_path / "rebuilt.json"
    calibration.certify(source, rebuilt)
    assert rebuilt.read_bytes() == retained.read_bytes()
    calibration.verify(source, retained)


@pytest.mark.parametrize("bad", [{"gap_lower": "1"}, {"c_floor": True},
                                {"eta_upper": False}, {"gap_lwoer": "1"}, [], None])
def test_witness_checks_every_supplied_expectation(bad):
    data = payload("witness")
    data["expected"] = bad
    with pytest.raises(ValueError):
        witness.verify(data)


@pytest.mark.parametrize("kind", ["witness", "explicit"])
@pytest.mark.parametrize("case", ["duplicate_type", "missing_id", "unknown_refinement", "missing_rate"])
def test_type_table_membership_and_completeness(kind, case):
    e = entry()
    entries = [e]
    if case == "duplicate_type":
        entries.append(copy.deepcopy(e))
    elif case == "missing_id":
        del e["id"]
    elif case == "unknown_refinement":
        e["refinement_targets"] = ["absent"]
    else:
        del e["rate_lower"]
    with pytest.raises(ValueError):
        check(kind, payload(kind, entries))


@pytest.mark.parametrize("kind", ["witness", "explicit"])
@pytest.mark.parametrize("bad_json", ["duplicate", "NaN", "Infinity"])
def test_file_entry_points_reject_ambiguous_or_nonfinite_evidence(tmp_path, kind, bad_json):
    source, output = tmp_path / "source.json", tmp_path / "receipt.json"
    data = payload(kind)
    text = json.dumps(data)
    if bad_json == "duplicate":
        text = text.replace('"scope":', '"scope": "physical_source_receipt", "scope":')
    else:
        text = text[:-1] + ', "extra_metadata": ' + bad_json + '}'
    source.write_text(text, encoding="utf-8")
    output.write_text("retained evidence", encoding="utf-8")
    command = ([str(HERE / "verify_collar_gap_certificate.py"), str(source)] if kind == "witness"
               else [str(HERE / "finite_collar_gap_certificate.py"), "certify",
                     "--manifest", str(source), "--output", str(output)])
    run = subprocess.run([sys.executable, *command], capture_output=True, text=True, cwd=tmp_path)
    assert run.returncode != 0
    assert "ValueError" in run.stderr or "DuplicateKeyError" in run.stderr
    assert output.read_text(encoding="utf-8") == "retained evidence"


def test_replay_rejects_duplicate_receipt_keys(tmp_path):
    source, output = tmp_path / "source.json", tmp_path / "receipt.json"
    source.write_text(json.dumps(payload("explicit")), encoding="utf-8")
    calibration.certify(source, output)
    text = output.read_text(encoding="utf-8").replace('"gap_lower":', '"gap_lower": "10", "gap_lower":')
    output.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key"):
        calibration.verify(source, output)


@pytest.mark.parametrize("kind", KINDS)
def test_valid_cli_runs_outside_the_repository(tmp_path, kind):
    source, output = tmp_path / "source.json", tmp_path / "receipt.json"
    e = entry()
    e["influences"] = [influence("1/8"), influence("1/4")]
    source.write_text(json.dumps(payload(kind, [e])), encoding="utf-8")
    command = ([str(HERE / "verify_collar_gap_certificate.py"), str(source)] if kind == "witness"
               else [str(HERE / "finite_collar_gap_certificate.py"), "certify",
                     "--manifest", str(source), "--output", str(output)])
    run = subprocess.run([sys.executable, *command], capture_output=True, text=True, cwd=tmp_path)
    assert run.returncode == 0, run.stderr
    result = json.loads(run.stdout if kind == "witness" else output.read_text(encoding="utf-8"))
    assert Q(result["gap_lower"]) == Q(15, 32)
    assert result["physical_clay_receipt"] is False
    if kind != "witness":
        run = subprocess.run([sys.executable, str(HERE / "finite_collar_gap_certificate.py"),
                              "verify", "--manifest", str(source), "--receipt", str(output)],
                             capture_output=True, text=True, cwd=tmp_path)
        assert run.returncode == 0, run.stderr
