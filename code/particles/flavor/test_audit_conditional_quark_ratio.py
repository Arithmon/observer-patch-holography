"""Guards for a downstream diagnostic of the frozen conditional quark chart."""
import copy
import importlib.util
import json
from decimal import Decimal
from pathlib import Path

import pytest


PATH = Path(__file__).with_name("audit_conditional_quark_ratio.py")
SPEC = importlib.util.spec_from_file_location("conditional_quark_ratio_diagnostic", PATH)
diagnostic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(diagnostic)
D = Decimal


def load_inputs():
    return json.loads(diagnostic.PACKET.read_text()), json.loads(diagnostic.REFERENCE.read_text())


def test_receipt_replays_and_all_fixed_ablations_miss_the_lower_edge():
    result = diagnostic.build()
    assert result == json.loads(diagnostic.OUTPUT.read_text())
    comparison = result["comparison"]
    native = result["native"]
    assert D(comparison["full_below_lower_by"]) > 0
    assert D(comparison["lower_order_below_lower_by"]) > 0
    gap = D(comparison["full_below_lower_by"])
    assert all(abs(D(row["change_from_full"])) < gap
               for row in native["single_term_ablations"].values())
    assert native["single_term_ablations"]["common_delta_g"]["change_from_full"] == "0.000000000000"
    assert result["status"] == "retrospective_conditional_sensitivity_not_model_selection"


def test_reconstruction_rejects_mutated_mass_or_promotion():
    packet, _ = load_inputs()
    broken = copy.deepcopy(packet)
    broken["full_rscc"]["native_mass_chart_GeV"]["s"] = "0.094"
    with pytest.raises(ValueError, match="exponents and masses disagree"):
        diagnostic.native_diagnostic(broken)
    broken = copy.deepcopy(packet)
    broken["guards"]["source_selected_flavor_law"] = True
    with pytest.raises(ValueError, match="conditional interpretation"):
        diagnostic.native_diagnostic(broken)


def test_reference_changes_only_downstream_comparison():
    packet, reference = load_inputs()
    native = diagnostic.native_diagnostic(packet)
    changed = copy.deepcopy(reference)
    changed["light_quantities"][1]["central"] = "28"
    assert diagnostic.native_diagnostic(packet) == native
    assert diagnostic.reference_diagnostic(native, changed)["full_below_lower_by"] != diagnostic.reference_diagnostic(native, reference)["full_below_lower_by"]


def test_separate_pdg_marginals_are_not_substituted_for_direct_ratio():
    packet, reference = load_inputs()
    comparison = diagnostic.reference_diagnostic(diagnostic.native_diagnostic(packet), reference)
    direct = D(comparison["pdg_direct_quoted_90_percent"]["central"])
    from_separate = D(comparison["displayed_marginal_central_arithmetic_only"]["from_u_d_s"])
    assert abs(direct - from_separate) > D("0.2")
    assert comparison["boundary"].startswith("PDG mass and ratio rows are separate marginals")
