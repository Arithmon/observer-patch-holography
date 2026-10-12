"""Endpoint comparison checks independent arithmetic and evidence boundaries."""
from copy import deepcopy
from decimal import localcontext
from fractions import Fraction
import hashlib
import json

import pytest

from model.predictions import endpoint_audit as endpoint


def exact_interval(box):
    return Fraction(box["lo"]), Fraction(box["hi"])


def test_known_endpoint_disagrees_without_promoting_incomplete_transport():
    result = endpoint.audit()
    assert result["verdicts"] == {
        "inherited_incomplete_transport_contract": "physical_completion_pending",
        "literal_complete_zero_correction_hypothesis": "inconsistent_with_declared_comparison_band",
        "OPH_as_a_whole": "not_determined",
    }
    assert result["comparison_class"] == "retrospective_known_data"
    assert result["arithmetic"]["source_interval_disjoint_from_comparison_band"]
    scope = result["interpretation"]
    for field in ("fit_applied", "P_recalibrated_from_target", "prospective_prediction", "empirical_validation_of_OPH"):
        assert scope[field] is False
    assert "no sigma significance" in scope["theory_uncertainty"]


def test_exact_residual_and_outward_ppm_against_certificate():
    result = endpoint.audit()
    certificate = json.loads((endpoint.ROOT / endpoint.response.CERT).read_text())
    source = certificate["modes"][endpoint.MODE]["certified_enclosure"]["alpha_inv"]
    slo, shi = exact_interval(source)
    target = Fraction(137035999177, 10**9)
    half_width = Fraction(105, 10**9)
    arithmetic = result["arithmetic"]
    assert exact_interval(arithmetic["target_comparison_band"]) == (target-half_width, target+half_width)
    assert exact_interval(arithmetic["target_centre_minus_source_interval"]) == (target-shi, target-slo)
    required = exact_interval(arithmetic["correction_needed_to_enter_comparison_band"])
    assert required == (target-half_width-shi, target+half_width-slo)
    assert required[0] > Fraction(41, 1000)
    assert Fraction(arithmetic["minimum_absolute_correction_to_enter_band"]) == required[0]
    plo, phi = exact_interval(arithmetic["relative_inverse_alpha_deficit_ppm"])
    assert plo <= 10**6*(target-shi)/target <= 10**6*(target-slo)/target <= phi
    assert Fraction(300) < plo < phi < Fraction(301)
    # The printed P-detector point is lower precision than its parent enclosure.
    receipt = json.loads((endpoint.ROOT / endpoint.RECEIPT).read_text())
    assert exact_interval(result["source"]["forward_P_interval"]) == exact_interval(receipt["source_branch"]["P_imported_interval"])
    assert result["source"]["inverse_alpha_interval"] != {"lo": receipt["source_branch"]["inverse_alpha"], "hi": receipt["source_branch"]["inverse_alpha"]}


def test_global_decimal_context_cannot_change_audit_arithmetic():
    reference = endpoint.audit()
    with localcontext() as context:
        context.prec = 6
        context.rounding = "ROUND_UP"
        assert endpoint.audit() == reference


def test_target_fixture_is_cited_bound_and_returned_without_aliasing():
    target = endpoint.target_fixture()
    original = deepcopy(target)
    output = endpoint.audit(target)
    assert target == original
    assert target["citation"] == "https://physics.nist.gov/cuu/pdf/wall_2022.pdf"
    assert target["value"] == "137.035999177"
    assert target["standard_uncertainty"] == "0.000000021"
    raw = json.dumps(target, sort_keys=True, separators=(",", ":")).encode()
    assert hashlib.sha256(raw).hexdigest() == output["target_sha256"]
    output["target"]["value"] = "1"
    assert endpoint.target_fixture() == original


@pytest.mark.parametrize("field,value", [
    ("value", "NaN"), ("value", "Infinity"), ("value", "0"),
    ("value", -1), ("value", 137.035999177), ("value", "nonsense"),
    ("standard_uncertainty", "-0.1"), ("standard_uncertainty", "0"),
    ("comparison_half_width_in_standard_uncertainties", "0"),
    ("value", "136.994835177412937"),
    ("standard_uncertainty", "1"),
    ("comparison_half_width_in_standard_uncertainties", "5000000"),
    ("citation", "https://example.invalid/target"),
    ("exposure", "prospective"),
])
def test_invalid_or_mutated_target_cannot_rescue_source(field, value):
    target = endpoint.target_fixture()
    target[field] = value
    with pytest.raises(ValueError):
        endpoint.audit(target)


@pytest.mark.parametrize("target", [[], {}, None, "137.035999177"])
def test_target_parser_refuses_missing_schema(target):
    # audit(None) deliberately selects its pinned default; parser itself is strict.
    with pytest.raises(ValueError):
        endpoint._target(target)


def test_source_mutation_fails_before_comparison(tmp_path):
    for name in endpoint.SOURCE_PINS:
        copied = tmp_path / name
        copied.parent.mkdir(parents=True, exist_ok=True)
        copied.write_bytes((endpoint.ROOT / name).read_bytes())
    assert endpoint.audit(root=tmp_path) == endpoint.audit()
    certificate_path = tmp_path / endpoint.response.CERT
    certificate = json.loads(certificate_path.read_text())
    box = certificate["modes"][endpoint.MODE]["certified_enclosure"]["alpha_inv"]
    box["lo"], box["hi"] = "137.035999176", "137.035999178"
    certificate_path.write_text(json.dumps(certificate))
    with pytest.raises(ValueError, match="Source pin mismatch"):
        endpoint.audit(root=tmp_path)


def test_complete_provenance_and_parent_certificate_are_not_rebranded():
    result = endpoint.audit()
    assert result["source"]["pins"] == endpoint.SOURCE_PINS
    assert result["source"]["mode"] == "thomson_structured_running"
    assert result["source"]["artifact_hashes_verified"] is True
    assert result["source"]["certificate_replayed_here"] is False
    assert "not physical theory uncertainty" in result["source"]["scope"]
    for path, digest in result["source"]["pins"].items():
        assert hashlib.sha256((endpoint.ROOT / path).read_bytes()).hexdigest() == digest
