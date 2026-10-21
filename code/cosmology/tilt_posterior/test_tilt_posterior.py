"""Independent expanded-sample and original-decimal posterior controls."""
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import statistics

import pytest

PATH = Path(__file__).with_name("audit.py")
SPEC = importlib.util.spec_from_file_location("oph_tilt_audit", PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_weighted_results_match_independently_expanded_sample():
    rows = [(F(3), F(91, 100)), (F(1), F(99, 100)), (F(4), F(101, 100)),
            (F(2), F(97, 100))]
    expanded = sorted(x for w, x in rows for _ in range(int(w)))
    target = F(99, 100)
    result = mod.summary(rows, {"control": target})
    assert result["ns_mean"] == float(statistics.mean(expanded))
    assert result["ns_posterior_sd"] == math.sqrt(float(statistics.pvariance(expanded)))
    for q in mod.QUANTILES:
        expected = expanded[max(0, math.ceil(q * len(expanded)) - 1)]
        assert result["step_cdf_quantiles"][f"{float(q):.3f}"] == float(expected)
    observed = result["positions"]["control"]
    assert F(observed["posterior_mass_le_target_exact"]) == F(
        sum(x <= target for x in expanded), len(expanded))
    assert observed["weight_equal"] == "1"


def test_burn_in_is_per_chain_stored_rows_not_multiplicity():
    rows = mod.parse_rows([["1000", "0.9"], ["1", "0.95"],
                           ["1", "1.0"], ["1", "1.05"], ["1", "1.1"]])
    kept, dropped = mod.discard_burn_in(rows)
    assert dropped == 2
    assert kept == rows[2:]
    assert mod.summary(kept, {})["total_multiplicity"] == "3"


def test_original_decimal_tail_boundary_survives_float_rounding():
    rows = mod.parse_rows([["2", "0.99999999999999999999"],
                           ["3", "1"], ["5", "1.00000000000000000001"]])
    assert len({float(x) for _, x in rows}) == 1
    result = mod.summary(rows, {"boundary": F(1)})["positions"]["boundary"]
    assert result["weight_strictly_below"] == "2"
    assert result["weight_equal"] == "3"
    assert result["weight_strictly_above"] == "5"
    assert result["posterior_mass_le_target_exact"] == "1/2"


@pytest.mark.parametrize("bad", [[], [["0", "0.96"]], [["-1", "0.96"]],
                                 [["0.5", "0.96"]], [["1", "NaN"]],
                                 [["1", "inf"]], [["1", "3"]],
                                 [["1/2", "0.96"]]])
def test_invalid_or_non_cobaya_samples_fail(bad):
    with pytest.raises(ValueError):
        mod.parse_rows(bad)


def test_prediction_is_frozen_rational_not_refitted():
    assert mod.EDGE_TILT == F(463690317905960407, 480000000000000000)
    assert mod.predictions()["edge_center"] == mod.EDGE_TILT


def test_ordering_keeps_excursions_separate_from_multiplicity():
    rows = [(F(100), F(9, 10)), (F(2), F(9, 10)), (F(1), F(11, 10)),
            (F(4), F(9, 10)), (F(1), F(11, 10))]
    result = mod.ordering_diagnostics(rows, {"boundary": F(1)})
    assert result["ns_changes_between_stored_rows"] == 3
    assert result["contiguous_tail_excursions"]["boundary"] == 2
    assert result["effective_sample_size_estimator"] is False


def test_unequal_chain_pooling_and_omission_match_expanded_observations():
    inputs = [
        [["30", "0.7"], ["2", "0.95"], ["1", "1.02"]],
        [["10", "0.8"], ["1", "1.01"], ["7", "0.94"], ["2", "1.02"]],
        [["9", "0.8"], ["1", "0.9"], ["1", "1.01"], ["2", "0.96"], ["4", "1.02"]],
        [["20", "0.7"], ["3", "0.93"], ["1", "1.03"]],
    ]
    payload = {"schema": "oph.act_dr6.ns_projection.v1",
               "archive_sha256": mod.SOURCE_ARCHIVE_SHA256,
               "burn_in_fraction": "0.5",
               "chains": [{"member": f"synthetic-{i}", "source_sha256": "0" * 64,
                           "rows": rows} for i, rows in enumerate(inputs)]}
    result = mod.analyze(payload)
    expanded = [[F(x) for w, x in rows[len(rows)//2:] for _ in range(int(w))]
                for rows in inputs]
    for omitted, observed in enumerate(result["leave_one_chain_out"]):
        samples = [x for i, values in enumerate(expanded) if i != omitted for x in values]
        for key, target in mod.predictions().items():
            expected = F(sum(x <= target for x in samples), len(samples))
            assert observed["posterior_mass_le_target"][key] == float(expected)
    all_samples = [x for samples in expanded for x in samples]
    assert result["pooled"]["ns_mean"] == float(statistics.mean(all_samples))


def test_retained_actual_data_replay():
    result = mod.build()
    assert result == json.loads(mod.RECEIPT.read_text())
    pooled = result["pooled"]
    assert pooled["stored_rows"] == 379175
    assert pooled["total_multiplicity"] == "508985"
    # Separate archive-based float parser produced these integer event totals.
    edge = pooled["positions"]["edge_center"]
    assert edge["weight_strictly_below"] == "3114"
    assert edge["weight_equal"] == "0"
    assert F(edge["posterior_mass_le_target_exact"]) == F(3114, 508985)
    assert result["total_tail_excursions_across_chains"]["edge_center"] == 85
    assert result["analysis"]["tail_mass_is_p_value"] is False
    assert result["analysis"]["FZ15_prospective_verdict"] is None


def test_rewritten_projection_is_rejected_before_analysis(tmp_path):
    altered = tmp_path / "altered.gz"
    altered.write_bytes(b"arbitrary replacement archive")
    with pytest.raises(ValueError, match="projection hash differs"):
        mod.build(altered)


SLICE_PATH = Path(__file__).with_name('likelihood_slice.py')
SLICE_SPEC = importlib.util.spec_from_file_location('oph_act_likelihood_slice', SLICE_PATH)
slice_mod = importlib.util.module_from_spec(SLICE_SPEC)
SLICE_SPEC.loader.exec_module(slice_mod)


@pytest.fixture
def slice_packet(monkeypatch):
    """Synthetic validation fixture; its observation pin is local to this test.

    Production retains the independently measured-SACC projection pin.
    No synthetic observation is admitted by that production pin.
    """
    import copy
    n = 135
    observations = dict(archive_sha256=slice_mod.ARCHIVE_SHA, fits_sha256=slice_mod.FITS_SHA,
        official_commit=slice_mod.ACT_COMMIT, official_source_sha256=slice_mod.ACT_SOURCE_SHA,
        spectra=['tt', 'te', 'ee'], bins_per_spectrum=[45]*3, ell_cuts=slice_mod.ELL_CUTS,
        bandpower_metadata=[dict(spectrum='synthetic_unit_control', indices=list(range(n)), ell=list(range(n)))],
        data=[0.]*n, covariance_lower=[float(i == j) for i in range(n) for j in range(i+1)])
    fits = {}
    for name, theory, ns in [('fixed_tilt', [1.]*n, slice_mod.EDGE_NS), ('free_tilt', [0.]*n, .97)]:
        score = float(n) if name == 'fixed_tilt' else 0.
        fits[name] = dict(ns=ns, logA=3., P_act=1., A_act=1., chi2=score, official_chi2=score,
            predicted_bandpowers=theory, optimizer_success=True, optimizer_reported_chi2=score,
            optimizer_origin='recomputed', full_recalculation_delta_chi2=0.)
    packet = copy.deepcopy(dict(schema='oph.act_dr6.conditional_likelihood_slice.v1',
        scope='fixed_background_ACT_only_likelihood_slice', source_sha256=slice_mod.digest(SLICE_PATH),
        background=slice_mod.BACKGROUND, accuracy=slice_mod.ACCURACY, calibration_A_act=1.,
        pivot_Mpc_inverse=.05, parameter_bounds=slice_mod.BOUNDS, runtime=slice_mod.RUNTIME,
        fixed_physics=slice_mod.FIXED_PHYSICS, observations=observations, fits=fits, delta_chi2=135.,
        p_value=None, FZ15_verdict=None, full_cosmological_profile=False,
        native_primordial_generator_identified=False))
    monkeypatch.setattr(slice_mod, 'OBSERVATION_PROJECTION_SHA256', slice_mod.canonical_digest(observations))
    return packet


def test_conditional_slice_whitens_the_complete_covariance():
    score = slice_mod.chi_square([3., 0.], [[4., 1.], [1., 2.]], [1., 1.])
    assert score == pytest.approx(16/7, abs=1e-14)
    diagonal = slice_mod.chi_square([3., 0.], [[4., 0.], [0., 2.]], [1., 1.])
    assert abs(score-diagonal) > .5


def test_conditional_slice_uses_ACT_TE_and_EE_polarization_calibration():
    from types import SimpleNamespace
    import numpy as np
    window = SimpleNamespace(weight=np.ones((1, 1)), values=np.array([0]))
    like = SimpleNamespace(data_vec=np.zeros(3), spec_meta=[
        dict(pol=pol, idx=np.array([i]), window=window) for i, pol in enumerate(('tt', 'te', 'ee'))])
    predicted = slice_mod.projected_theory(like, {pol: np.array([4.]) for pol in ('tt', 'te', 'ee')}, 2.)
    assert predicted.tolist() == [4., 2., 1.]


def test_conditional_slice_validation_recomputes_both_scores(slice_packet):
    scores, delta = slice_mod.check(slice_packet)
    assert scores == {'fixed_tilt': 135., 'free_tilt': 0.}
    assert delta == 135.


@pytest.mark.parametrize('key', ['ns', 'logA', 'P_act', 'A_act', 'chi2', 'official_chi2',
                                  'optimizer_reported_chi2', 'full_recalculation_delta_chi2'])
@pytest.mark.parametrize('bad', [float('nan'), float('inf'), True])
def test_conditional_slice_rejects_nonfinite_or_boolean_fit_numbers(slice_packet, key, bad):
    slice_packet['fits']['fixed_tilt'][key] = bad
    with pytest.raises(ValueError):
        slice_mod.check(slice_packet)


@pytest.mark.parametrize('key,bad', [('schema', 'wrong'), ('p_value', .04), ('FZ15_verdict', 'confirmed'),
    ('full_cosmological_profile', True), ('native_primordial_generator_identified', True),
    ('calibration_A_act', True), ('pivot_Mpc_inverse', .1), ('delta_chi2', float('nan'))])
def test_conditional_slice_rejects_promotions_and_changed_units(slice_packet, key, bad):
    slice_packet[key] = bad
    with pytest.raises(ValueError):
        slice_mod.check(slice_packet)


@pytest.mark.parametrize('key', ['data', 'covariance_lower', 'bandpower_metadata', 'official_commit',
                                  'official_source_sha256', 'ell_cuts'])
def test_conditional_slice_binds_observation_arrays_and_order_to_the_pin(slice_packet, key):
    if key in ('data', 'covariance_lower'):
        slice_packet['observations'][key][0] += .1
    elif key == 'bandpower_metadata':
        slice_packet['observations'][key][0]['indices'].reverse()
    elif key == 'ell_cuts':
        slice_packet['observations'][key]['TT'][0] = 601
    else:
        slice_packet['observations'][key] = '0'*len(slice_packet['observations'][key])
    with pytest.raises(ValueError, match='projection differs'):
        slice_mod.check(slice_packet)


@pytest.mark.parametrize('change', ['theory', 'optimization_score', 'recalculation', 'runtime',
                                    'physics', 'bounds', 'fit_set', 'refitted_tilt'])
def test_conditional_slice_rejects_failed_reproduction_contracts(slice_packet, change):
    fit = slice_packet['fits']['fixed_tilt']
    if change == 'theory': fit['predicted_bandpowers'][0] += .1
    elif change == 'optimization_score': fit['optimizer_reported_chi2'] += .01
    elif change == 'recalculation': fit['full_recalculation_delta_chi2'] = .01
    elif change == 'runtime': slice_packet['runtime']['camb'] = '0.0'
    elif change == 'physics': slice_packet['fixed_physics']['parameters']['Alens'] = True
    elif change == 'bounds': fit['P_act'] = 1.2
    elif change == 'fit_set': slice_packet['fits']['unrequested_alternative'] = fit
    elif change == 'refitted_tilt': fit['ns'] += .001
    with pytest.raises(ValueError):
        slice_mod.check(slice_packet)


def test_conditional_slice_keeps_the_exact_existing_tilt_commitment():
    assert slice_mod.EDGE_NS == float(mod.EDGE_TILT)


def test_conditional_slice_rejects_invalid_covariance_and_nonfinite_theory():
    import numpy as np
    with pytest.raises(ValueError):
        slice_mod.chi_square([1., 2.], [[1., .3], [.2, 1.]], [0., 0.])
    with pytest.raises(np.linalg.LinAlgError):
        slice_mod.chi_square([1., 2.], [[1., 2.], [2., 1.]], [0., 0.])
    with pytest.raises(ValueError):
        slice_mod.chi_square([1., 2.], [[1., 0.], [0., 1.]], [float('nan'), 0.])
    with pytest.raises(ValueError, match='triangle'):
        slice_mod.unpack_covariance([1., 2.], 2)


def test_retained_measured_cmb_likelihood_slice():
    packet = json.loads(slice_mod.RECEIPT.read_text())
    scores, delta = slice_mod.check(packet)
    assert 171.9 < scores['fixed_tilt'] < 172.0
    assert 167.8 < scores['free_tilt'] < 167.9
    assert 4.0 < delta < 4.2
    assert packet['fits']['fixed_tilt']['ns'] == float(mod.EDGE_TILT)
    assert packet['full_cosmological_profile'] is False
