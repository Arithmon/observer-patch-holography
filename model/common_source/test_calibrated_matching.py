"""Input ancestry, convention and adversarial retained-output controls."""
from copy import deepcopy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import PurePosixPath, PureWindowsPath

import pytest

from tools import strict_json
from . import calibrated_matching as matching


@pytest.mark.parametrize("flavor,root", [(PurePosixPath, "/repo"), (PureWindowsPath, "C:/repo")])
def test_custody_paths_are_identical_across_host_path_flavors(flavor, root):
    base = flavor(root)
    assert matching.repository_key(base/"model"/"common_source"/"calibrated_matching.py", base) == (
        "model/common_source/calibrated_matching.py")


def test_complete_retained_external_comparisons_and_unfitted_miss():
    result = matching.check(strict_json.load(matching.OUTPUT))
    assert result["summary"]["all_fitted_inputs_reproduced_within_quoted_one_sigma"]
    assert not result["summary"]["W_within_one_quoted_sigma"]
    row = result["comparisons"]["calibrated"]["SMDR_MW_PDG"]
    assert not row["used_for_calibration"]
    assert -1.51 < row["residual_over_quoted_uncertainty"] < -1.50
    assert result["summary"]["rounded_muon_and_electron_miss_their_precision"]


def test_one_sigma_is_comparison_not_hidden_acceptance_tolerance():
    outputs = deepcopy(strict_json.load(matching.OUTPUT)["external_outputs"]["calibrated"])
    target, sigma, _ = matching.TARGETS["SMDR_Melectron_pole"]
    with localcontext() as ctx:
        ctx.prec = 60
        outputs["SMDR_Melectron_pole"] = str(Decimal(target)+Decimal(sigma))
        assert matching.comparison(outputs)["SMDR_Melectron_pole"]["absolute_residual_le_one_quoted_uncertainty"]
        outputs["SMDR_Melectron_pole"] = str(Decimal(target)+Decimal(sigma)*(1+Decimal("1e-20")))
        assert not matching.comparison(outputs)["SMDR_Melectron_pole"]["absolute_residual_le_one_quoted_uncertainty"]


def test_gut_and_quartic_normalizations_are_explicit():
    packet = strict_json.load(matching.PACKET)
    data = matching.smdr_inputs(packet)
    assert F(data["lambda"])*4 == F(packet["matched_SM"]["lambda_paper"])
    squared = F(data["gp"])**2/F(packet["matched_SM"]["g1_GUT"])**2
    assert abs(squared-F(3, 5)) < F(1, 10**50)
    assert data["Q"] == "91.1876"
    assert data["Delta_alpha_had_5_MZ"] == "0.02783"


def test_changed_calibration_is_not_silently_detached_from_forward_output(monkeypatch):
    packet = strict_json.load(matching.PACKET)
    original = matching.smdr_inputs
    packet["matched_SM"]["yukawa"]["charged_lepton"][0] = "0.01"
    monkeypatch.setattr(matching, "smdr_inputs", lambda **kw: original(packet, **kw))
    with pytest.raises(ValueError, match="receipt, inputs"):
        matching.check(strict_json.load(matching.OUTPUT))


@pytest.mark.parametrize("mutation", ("omit", "bool", "nan", "double", "duplicate"))
def test_complete_external_output_contract(mutation):
    outputs = deepcopy(strict_json.load(matching.OUTPUT)["external_outputs"]["calibrated"])
    if mutation == "omit":
        outputs.pop("SMDR_MW_PDG")
    elif mutation == "bool":
        outputs["SMDR_Melectron_pole"] = True
    elif mutation == "nan":
        outputs["SMDR_Melectron_pole"] = "NaN"
    raw = "OPH_MATCHING_JSON="+json.dumps(outputs)
    if mutation == "double":
        raw += "\n"+raw
    elif mutation == "duplicate":
        raw = raw.replace('{', '{"SMDR_MW_PDG":"80",', 1)
    with pytest.raises(ValueError):
        matching.parse_output(raw)


def test_tampered_scope_or_summary_is_rejected():
    result = strict_json.load(matching.OUTPUT)
    result["role"] = "derived_particle_spectrum"
    with pytest.raises(ValueError):
        matching.check(result)
    result = strict_json.load(matching.OUTPUT)
    result["summary"]["W_within_one_quoted_sigma"] = True
    with pytest.raises(ValueError):
        matching.check(result)
