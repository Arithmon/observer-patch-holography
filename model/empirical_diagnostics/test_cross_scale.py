"""Arithmetic and claim controls for the fixed spectral cross-scale audit."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path

import mpmath
import pytest

from . import cross_scale as audit


def test_exact_sqrt_and_pi_enclosures():
    for value in (F(0), F(4, 9), F(2), F(5), F(17, 13)):
        low, high = audit.sqrt_bounds(value)
        assert low**2 <= value <= high**2
        assert high-low <= F(1, 10**60)
    assert audit.sqrt_bounds(F(4, 9), digits=1) == (F(6, 10), F(7, 10))
    assert audit.sqrt_bounds(F(4)) == (F(2), F(2))
    lo, hi = audit.pi_bounds()
    known_lo = F("3.14159265358979323846264338327950288419716939937510")
    assert known_lo < lo < hi < known_lo + F(1, 10**50)
    assert hi-lo < F(1, 10**65)


@pytest.mark.parametrize("value,digits", [(F(-1), 60), (2, 60), (F(2), 0), (F(2), True)])
def test_invalid_enclosure_inputs_rejected(value, digits):
    with pytest.raises(ValueError):
        audit.sqrt_bounds(value, digits)


def test_prediction_has_no_file_or_empirical_input(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("forward prediction attempted file access")
    monkeypatch.setattr(Path, "read_text", forbidden)
    monkeypatch.setattr(Path, "read_bytes", forbidden)
    prediction = audit.spectral_prediction()
    assert prediction["observation_inputs"] == []
    mp = mpmath.mp.clone()
    mp.dps = 100
    alpha = (79+48*mp.sqrt(2)+42*mp.sqrt(3)+3*mp.sqrt(6))/31104
    p = (1+mp.sqrt(5))/2 + mp.sqrt(mp.pi)*alpha
    for name, value in (("alpha_interval", alpha), ("inverse_alpha_interval", 1/alpha),
                        ("P_interval", p), ("scalar_tilt_interval", 1-p/48)):
        low, high = map(F, prediction[name])
        assert low < F(str(value)) < high


@pytest.fixture(scope="module")
def result():
    return audit.build()


def test_actual_CODATA_and_ACT_comparison(result):
    assert result["coupling"]["within_one_reported_standard_uncertainty"]
    standardized = result["coupling"]["signed_difference_over_reported_standard_uncertainty"]
    assert F("0.5983") < F(standardized["lo"]) <= F(standardized["hi"]) < F("0.5984")
    assert result["cosmology"]["forward_tail_weight_interval"] == [3114, 3114]
    assert result["cosmology"]["multiplicity"] == 508985
    assert result["cosmology"]["forward_and_inherited_have_same_tail_weight"]
    assert not result["cosmology"]["posterior_mass_is_p_value"]
    contract = result["cosmology"]["physical_comparison_contract"]
    assert contract["pivot_Mpc_inverse"] == "0.05"
    assert contract["ns_prior"]["bounds"] == ["0.85", "1.1"]
    assert contract["documentation"] in result["custody"]["input_sha256"]


def test_agreeing_scalar_does_not_promote_missing_physical_join(result):
    status = result["joint_status"]
    assert status["declared_chiral_prescriptions_passing_both_lepton_ratios"] == []
    assert len(status["declared_chiral_prescriptions_failing_both_ratio_requirement"]) == 4
    assert "spectral_scalar_to_Thomson_photon_residue" in status["missing_contracts"]
    assert status["joint_likelihood_or_significance"] is None
    assert not status["all_public_datasets_reproduced"]
    assert not status["source_only_measured_universe_established"]


def test_target_cannot_be_moved_to_improve_agreement():
    target = audit.endpoint_audit.target_fixture()
    target["value"] = "137.035999189566155"
    with pytest.raises(ValueError, match="pinned"):
        audit.comparison(audit.spectral_prediction(), target)


def test_joint_contract_rejects_vacuous_or_malformed_mass_success(result):
    chiral = audit.strict_json.loads((audit.ROOT/audit.CHIRAL).read_text())
    quantum = audit.strict_json.loads((audit.ROOT/audit.QUANTUM).read_text())
    for comparisons in ([], [{"prescription": "x", "both_ratios_pass": "true"}]):
        changed = deepcopy(chiral)
        changed["separate_empirical_comparison"]["comparisons"] = comparisons
        with pytest.raises(ValueError):
            audit.joint_status(alpha_identified=True, same_cell_proved=True,
                               chiral=changed, quantum=quantum)


def test_retained_receipt_matches_current_sources(result):
    assert audit.strict_json.loads(audit.OUTPUT.read_text(encoding="utf-8")) == result


def test_canonical_model_verifier_replays_cross_scale_comparison():
    from model.verify import REPLAYS
    assert ("model.empirical_diagnostics.cross_scale", None,
            "model/empirical_diagnostics/cross_scale_results.json") in REPLAYS
