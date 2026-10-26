"""Native transition-current response, complete outcomes and clock dictionary."""
from fractions import Fraction as F

import pytest
import sympy as s

from . import current_response as current
from ...observer import native


def test_native_branches_fix_continuity_stationarity_and_noise_scale():
    data = current.source()
    assert data["native_branches_checked"] == 66*20*2
    assert data["current_scale"] == s.Rational(19, 150)
    assert data["selector_survival"] == s.Rational(27, 32)
    assert data["selector_timeout"] == s.Rational(5, 32)
    assert data["transition"]*data["reference"] == data["reference"]
    assert s.ones(1, 12)*data["incidence"] == s.zeros(1, 20)


@pytest.mark.parametrize("lag", [1, 2, 3])
def test_full_native_states_give_the_reduced_transition_current_response(lag):
    assert current.full_state_covariance(lag) == current.current_covariance(lag)
    assert current.current_covariance(lag) == current.current_covariance(lag).T


def test_collective_local_charge_attaches_to_the_same_transition_current():
    data = current.source()
    assert data["collective_current_scale"] == s.Rational(1, 36)
    assert s.ones(1, 12)*data["collective_occupations"] == s.Matrix([[state[0] for state in data["states"]]])
    for lag in (0, 1, 2, 3):
        expected = s.Rational(25, 114)*current.current_covariance(lag)
        assert current.current_covariance(lag, collective=True) == expected
        if lag:
            assert current.full_state_covariance(lag, collective=True) == expected


def test_cycle_noise_is_a_positive_form_and_retains_the_finite_sum_remainder():
    data, result = current.source(), current.cycle_noise()
    c, d, b = result["weighted_projector"], data["weights"], data["incidence"]
    assert result["incidence_rank"] == 9 and result["cycle_rank"] == 11
    assert len(result["active_ports"]) == 10
    assert c == c.T and c*d.inv()*c == c
    assert b*c == s.zeros(12, 20)
    # x^T C x = (Cx)^T D^-1(Cx) >= 0, using exact positive D.
    assert all(value > 0 for value in d.diagonal())
    assert all(F(value) > 0 for value in result["finite_remainder_traces"])
    # A gradient readout has no zero-frequency current circulation.
    potential = s.Matrix(range(12))
    assert (b.T*potential).T*c*(b.T*potential) == s.zeros(1)


def test_current_is_a_transition_instrument_not_its_mean_operator_squared():
    data, report = current.source(), current.modular_dictionary()
    assert report["operator_variance"] == "19/300"
    assert report["record_variance"] == "19/150"
    assert not report["repair_proposal_current_is_modular_charge_derivative"]
    first, second = data["first_effects"][0], data["second_effects"][0]
    assert any(b != a*a for a, b in zip(first, second, strict=True))
    assert all(b == 2*a*a for a, b in zip(first, second, strict=True))
    # The actual reference and local charge algebra are diagonal, so their
    # modular commutator vanishes exactly despite nonzero native record noise.
    reference = s.diag(*data["reference"])
    q = s.diag(*first)
    assert reference*q == q*reference


def test_source_flag_protocol_retains_zero_currents_and_both_outcomes():
    result = current.fixed_seam_protocol()
    law = result["complete_joint_current_flag_law"]
    assert result["native_histories"] == 1056
    assert sum((F(row["probability"]) for row in law), F()) == 1
    assert sum((F(row["probability"]) for row in law if row["first_current"] == 0), F()) == F(131, 150)
    assert result["acceptance_probability"] == "1/2"
    assert result["accepted_next_current_mean"] == "-19/600"
    assert result["rejected_next_current_mean"] == "19/600"
    assert result["current_lag_one_covariance"] == "-19/300"
    assert result["unconditional_signed_flag_current_mean"] == "-19/600"
    assert not result["reset_of_work_between_events"]
    assert not result["conditioning_is_a_force_perturbation_or_conductivity_measurement"]


def test_unbiased_reversed_and_nonreset_auxiliary_controls():
    assert current.fixed_seam_protocol(0)["unconditional_signed_flag_current_mean"] == "0"
    assert current.fixed_seam_protocol(-1)["unconditional_signed_flag_current_mean"] == "19/600"
    expected = current.fixed_seam_protocol()["complete_joint_current_flag_law"]
    for auxiliary in (1, 20, len(native.states())-1):
        assert current.fixed_seam_protocol(auxiliary_index=auxiliary)["complete_joint_current_flag_law"] == expected


def test_actual_stopped_current_covariance_keeps_survival_probability():
    data = current.source()
    # Recompute the two-event moment directly from accepted native branches.
    seam = data["seams"][0]
    weight = data["selector_survival"]*data["probabilities"][0]
    exact = F()
    for probability, state in zip(current.local_preparation.reference(), native.states(binary=True), strict=True):
        for coin0 in (0, 1):
            middle, (_, j0) = native.native(state, seam, coin0)
            for coin1 in (0, 1):
                _, (_, j1) = native.native(middle, seam, coin1)
                exact += probability*F(weight*weight)*j0*j1/4
    predicted = data["selector_survival"]**2*current.current_covariance(1)[0, 0]
    assert s.Rational(exact) == predicted
    assert predicted != current.current_covariance(1)[0, 0]


def test_invalid_protocol_or_lag_is_rejected():
    for lag in (-1, True, 1.5):
        with pytest.raises(ValueError):
            current.current_covariance(lag)
    for sign in (True, 2, "1"):
        with pytest.raises(ValueError):
            current.fixed_seam_protocol(sign)


def test_finite_reference_budget_supports_a_complete_attempted_trial_readout():
    report = current.report()
    preparation = report["finite_preparation_error"]
    assert F(preparation["bounded_signed_current_mean_error_upper"]) < F("0.0000015")
    assert preparation["finite_preparation_bound_below_signal"]
    contract = report["finite_attempted_trial_contract"]
    assert contract["ideal_biased_mean"] == "-19/3960"
    assert contract["attempted_preparations"] == 2_000_000
    assert F(contract["Hoeffding_exponent_lower"]) > 5
    assert contract["preparation_failures_count_as_trials"]
    assert not contract["exact_reset_or_actual_trial_independence_assumed"]
    assert not contract["complete_sampled_campaign_executed"]
    assert not report["physical_scope"]["Chain_scalar_identified_with_current_noise_or_photon_residue"]
    assert not report["collective_local_charge_current"]["protected_bit_is_part_of_raw_current_record"]


def test_nonzero_field_frequencies_follow_complete_records_without_initial_loads():
    result = current.finite_field_readout()
    assert result["complete_native_histories"] == 66*16
    assert result["open_endpoint_histories"] > 0
    assert result["total_probability_exact"] == "1"
    assert not result["periodic_endpoint_assumed"]
    assert not result["hidden_initial_occupation_used_by_readout"]
    left, right = result["same_zero_current_word_distinct_DC_witness"]
    assert left != right
    assert not result["DC_occupation_recovered"]


def test_fourier_boundary_correction_works_on_an_open_native_path():
    seam = current.source()["seams"][0]
    # All four admissible records represent different addressed events; the
    # final event follows x_(N-1), so it cannot affect its Fourier coordinate.
    records = ((seam, 1), (seam, 0), (seam, 0), (seam, 0))
    altered_final = records[:-1]+((seam, -1),)
    b = current.source()["incidence"][:, 0]
    assert current.density_fourier(records, 1) == b
    assert current.density_fourier(altered_final, 1) == b
    # Explicit open path, without relying on the reconstruction formula.
    initial = s.zeros(12, 1)
    initial[seam[0]] = 1
    initial[0] = 1
    snapshots = [initial, initial-b, initial-b, initial-b]
    assert sum((s.I**n*x for n, x in enumerate(snapshots)), s.zeros(12, 1)) == current.density_fourier(records, 1)


def test_loop_read_is_a_martingale_difference_not_a_photon_mode():
    data, result = current.source(), current.frequency_dictionary()
    h = s.Matrix(list(map(s.Rational, result["bounded_example_readout_weights"])))
    assert max(map(abs, h)) == 1
    assert h.T*data["current_mean_map"] == s.zeros(1, 66)
    assert result["bounded_example_variance"] == "19/1620"
    for lag in (1, 2, 3):
        assert current.full_state_covariance(lag)*h == s.zeros(20, 1)
    assert result["cycle_record_is_martingale_difference_for_fixed_IID_selector"]
    assert not result["white_two_point_law_implies_independent_or_Gaussian_records"]


@pytest.mark.parametrize("cosine", [-1, F(0), F(1, 2), 1])
def test_frequency_Ward_and_collective_normalization(cosine):
    data = current.source()
    spectrum = current.frequency_spectrum(cosine)
    assert spectrum == spectrum.T
    assert current.frequency_spectrum(cosine, collective=True) == s.Rational(25, 114)*spectrum
    if cosine == 1:
        assert spectrum == current.cycle_noise()["matrix"]
    else:
        b, m = data["incidence"], data["mean_map"]
        density = data["current_scale"]*(s.eye(12)-m*m)*(s.eye(12)-2*cosine*m+m*m).inv()
        assert b*spectrum*b.T == 2*(1-cosine)*density


def test_modal_spectrum_identity_has_only_positive_relaxation_poles():
    m, c = s.symbols("m c", real=True)
    matrix_factor = 1-(1-m)*(c-m)/(1-2*m*c+m*m)
    displayed = (1+m)*(1-c)/((1-m)**2+2*m*(1-c))
    assert s.factor(matrix_factor-displayed) == 0
    data = current.source()
    assert min(2*data["mean_map"][i, i]-1 for i in range(12)) == s.Rational(19, 27)


def test_invalid_frequency_and_incomplete_record_contracts_fail_closed():
    seam = current.source()["seams"][0]
    records = ((seam, 0),)*4
    for frequency in (0, 4, True, 1.5):
        with pytest.raises(ValueError):
            current.density_fourier(records, frequency)
    for cosine in (True, 0.5, s.nan, 2):
        with pytest.raises(ValueError):
            current.frequency_spectrum(cosine)
    with pytest.raises(ValueError):
        current.density_fourier((((0, 1), 0),)*4, 1)
