import math

import numpy as np
import pytest

from . import action, history as h


@pytest.fixture(scope="module")
def executed():
    return h.execute()


def test_full_chiral_rest_history_matches_actual_source_fock_and_clock(executed):
    report,arrays = executed
    assert report["normal_records"] == 64 and report["gaussian_final_field_records"] == 16
    assert abs(report["total_record_probability"]-1) < 2e-11
    assert report["complete_probe_CQ_distance_to_independent_continuous"] < 2e-11
    assert report["full90_mode_and_probe_final_residual"] < 2e-11
    assert report["old_record_marginal_residual"] < 2e-11
    assert report["work_preparations"] == 1 and report["work_resets"] == 0
    assert arrays["probe_outputs"].shape == (64,256,4)
    assert arrays["gaussian_covariances"].shape == (16,180,180)
    assert report["layout"]["modes_discarded"] == 0
    assert report["layout"]["neutral_zero_modes_retained"] == 6


def test_model_electric_charge_is_conserved_at_completed_record_checkpoints(executed):
    report,arrays = executed
    q = np.array(action.spin_electric_charge().diagonal(),float).ravel()
    q = q[report["layout"]["initial_fermionic_wire_order"]]
    for covariance in arrays["gaussian_covariances"]:
        occupations = np.array([(1+covariance[2*j,2*j+1])/2 for j in range(90)])
        assert abs(q@occupations) < 2e-11
        assert np.linalg.eigvalsh(1j*covariance).max() < 1+2e-10


def test_gap_inverse_uses_only_record_means_and_declared_clock(executed,monkeypatch):
    report,_ = executed
    def reject(*args,**kwargs):
        raise AssertionError("decoder consulted source Hamiltonian")
    monkeypatch.setattr(action,"pairing_kernel",reject)
    monkeypatch.setattr(h.model,"parameters",reject)
    probabilities = {tuple(map(int,key)):value for key,value in report["source_record_probabilities"].items()}
    result = h.infer_ratios(probabilities,report["clock_calibration"],report["phase_window_contract"])
    np.testing.assert_allclose(result["native_family_order_gap_over_clock"],
                               report["independent_analytic_rest_gaps_over_clock"],atol=2e-12,rtol=0)
    ratios = sorted(result["native_family_order_gap_over_clock"])
    assert ratios[-1]/ratios[0] < 1.5


def test_clock_and_old_readouts_remain_in_the_same_reference_bearing_history(executed):
    report,arrays = executed
    reference = sum(state.conj().T@state for state in arrays["probe_outputs"])
    np.testing.assert_allclose(reference,np.diag([1,5,5,25])/36,atol=2e-11,rtol=0)
    assert abs(math.fsum(arrays["gaussian_weights"])-1) < 2e-11
    assert report["observer_inverse"]["finite_shot_or_physical_service_errors_certified"] is False


def test_broken_clock_is_inconclusive(executed):
    report,_ = executed
    probabilities = {tuple(map(int,key)):.8*value for key,value in report["source_record_probabilities"].items()}
    with pytest.raises(ValueError,match="calibrated circle"):
        h.infer_ratios(probabilities,report["clock_calibration"],report["phase_window_contract"])


def test_frequency_aliases_need_a_prior_contract(executed):
    report,_ = executed
    probabilities = {tuple(map(int,key)):value for key,value in report["source_record_probabilities"].items()}
    with pytest.raises(ValueError,match="prior admitted"):
        h.infer_ratios(probabilities,report["clock_calibration"])
