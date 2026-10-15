"""Receipt-level scientific scope and actual offline replay regressions."""
from copy import deepcopy
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from . import experiment
from .data import FOLDER, load_observations
from ...observer.scalar_response.experiment import compare_receipt


FROZEN_PREVIOUS_ARTIFACTS = {
    'protocol.py': 'b07cd72b87aa73f464ec90b79d53077f8c5da322acced1d58058bcd49a19dc7d',
    'data/occupation_source_receipt.json.gz': '1c3163dfffa5f19c9b6cb18c96d3d3e47dd378a1bcb0b531369f59c1fa2e65f7',
    'occupation_results.json': 'd846be5d7f6cdd79014ebc501048664f5e761d3affcdbd1cb358fd3738cb43ed',
}


def _previous_hashes():
    result = {}
    for name in FROZEN_PREVIOUS_ARTIFACTS:
        raw = (FOLDER.parent/name).read_bytes()
        if name.endswith('.gz'):
            raw = gzip.decompress(raw)
        result[name] = hashlib.sha256(raw).hexdigest()
    return result


@pytest.fixture(scope='module')
def receipt():
    before = _previous_hashes()
    result = experiment.run()
    assert before == _previous_hashes() == FROZEN_PREVIOUS_ARTIFACTS
    return result


def test_main_receipt_preserves_original_experiment_and_occupation_evidence(receipt):
    assert _previous_hashes() == FROZEN_PREVIOUS_ARTIFACTS
    assert receipt['checks']['original_noisy_origin_preserved'] is True
    # The source's observed central value must not be normalized to 1+0j.
    origin = next(row['value'] for row in receipt['calibration']['axis_means']
                  if row['index'] == [0, 0])
    assert origin['real'] == pytest.approx(.984, abs=1e-15)
    assert origin['imag'] == pytest.approx(.112, abs=1e-15)


def test_real_data_result_matches_independently_computed_numeric_regressions(receipt):
    # Reference values obtained directly from the retained 10x10 complex
    # subarray and its two outer-product axes, independently of run()/score().
    result = receipt['heldout_comparison']
    assert result['complex_rmse'] == pytest.approx(.10803244966731063, abs=1e-14)
    assert result['spin_mean_rmse'] == pytest.approx(.07639047774794971, abs=1e-14)
    assert result['spin_probability_rmse'] == pytest.approx(.038195238873974856, abs=1e-14)
    assert receipt['zero_response_baseline']['complex_rmse'] == pytest.approx(.18116246851928244, abs=1e-14)
    assert result['withheld_points'] == len(result['rows']) == 100
    assert result['real_statistics'] == 200
    assert receipt['calibration']['real_statistics'] == 42
    axis_keys = {tuple(row['index']) for row in receipt['calibration']['axis_means']}
    withheld_keys = {tuple(row['index']) for row in result['rows']}
    assert len(axis_keys) == 21 and len(withheld_keys) == 100
    assert axis_keys.isdisjoint(withheld_keys)
    assert result['box_violation_count'] == 0
    assert result['verdict'] == 'not_rejected_under_assumed_sampling'


def test_all_published_coordinate_conventions_remain_visible(receipt):
    conventions = receipt['unfitted_nominal_conventions']
    assert [(row['extent'], row['row_direction']) for row in conventions] == [
        (4, 1), (4, -1), (5, 1), (5, -1)]
    # Scientific failures remain in an otherwise passing software receipt.
    assert [row['box_violation_count'] for row in conventions] == [11, 0, 9, 2]
    assert receipt['all_pass'] is True
    for row in conventions:
        assert row['parameters_fitted'] is False
        assert row['physical_coordinate_identification_established'] is False
        assert row['empirical_confirmation'] is False
    assert receipt['orientation']['coordinate_direction_known'] is False
    assert receipt['orientation']['coordinate_extent_resolved'] is False


def test_software_success_never_depends_on_passing_the_scientific_comparison(monkeypatch, receipt):
    changed = deepcopy(load_observations())
    changed['validation'] = {key: 1+1j for key in changed['validation']}
    monkeypatch.setattr(experiment, 'load_observations', lambda: changed)
    result = experiment.run()
    assert result['calibration'] == receipt['calibration']
    assert result['all_pass'] is True
    assert result['heldout_comparison']['box_violation_count'] > 0
    assert result['heldout_comparison']['verdict'] == 'conditional_box_violation'
    assert result['heldout_comparison']['empirical_confirmation'] is False
    assert result['scope']['independent_validation_of_oph'] is False


def test_published_file_provenance_is_separate_from_trial_authentication(receipt):
    assert receipt['checks']['source_files_match_published_checksums'] is True
    provenance = receipt['provenance']
    assert provenance['doi'] == '10.5281/zenodo.18752287'
    assert provenance['raw_shot_records_included'] is False
    assert provenance['accepted_counts_included'] is False
    for name, metadata in provenance['files'].items():
        content = (FOLDER/'data'/name).read_bytes()
        assert hashlib.md5(content).hexdigest() == metadata['published_md5']
        assert hashlib.sha256(content).hexdigest() == metadata['sha256']
        assert len(content) == metadata['bytes']
    scope = receipt['scope']
    assert scope['published_experimental_data'] is True
    assert scope['data_kind'] == 'published_experimental_spin_means'
    assert scope['initialization_postselection'] is True
    assert scope['retrospective_withheld_points_not_prospective_acquisition'] is True
    for name in ('acquisition_independence_or_accepted_counts_authenticated',
                 'all_acquisition_attempts_available', 'physical_clock_calibrated',
                 'seven_mode_scalar_dynamics_tested', 'complete_native_history_tested',
                 'native_action_selected', 'oph_specific_prediction',
                 'source_execution_of_this_laboratory_preparation', 'independent_validation_of_oph'):
        assert scope[name] is False
    assert receipt['heldout_comparison']['provenance_verified'] is False
    assert receipt['heldout_comparison']['processed_means_establish_sampling_premise'] is False


def test_interface_and_protocol_bindings_reference_actual_bytes(receipt):
    assert set(receipt['source_interface_hashes']) == {'source.py', 'protocol.py', 'occupation_readout.py'}
    for name, digest in receipt['source_interface_hashes'].items():
        assert hashlib.sha256((FOLDER.parent/name).read_bytes()).hexdigest() == digest
    assert receipt['protocol_sha256'] == hashlib.sha256((FOLDER/'protocol.json').read_bytes()).hexdigest()


def test_two_conditional_sampling_certificates_are_combined_without_claiming_service_calibration(receipt):
    certificates = receipt['concentration_certificates']
    assert certificates['calibration']['simultaneous_real_statistics'] == 42
    assert certificates['validation']['simultaneous_real_statistics'] == 200
    joint = sum((F(value['familywise_error_upper']) for value in certificates.values()), F(0))
    assert F(receipt['conditional_joint_failure_upper_exact']) == joint <= F(1, 20)
    for value in certificates.values():
        assert F(value['familywise_error_upper']) <= F(1, 40)
        assert value['samples_per_real_mean_assumed'] == 250
        assert value['coverage'] == 'conditional_on_independent_bounded_raw_observations'
    result = receipt['heldout_comparison']
    assert result['calibration_failure_budget_included'] is True
    assert result['physical_service_error_budget_included'] is False
    assert result['empirical_confirmation'] is False


def test_receipt_is_json_safe_and_replays_the_retained_result(receipt):
    serialized = json.dumps(receipt, allow_nan=False)
    assert json.loads(serialized) == receipt
    retained = json.loads((FOLDER/'results.json').read_text())
    compare_receipt(retained, receipt)


def test_actual_cli_checks_retained_result_offline():
    root = Path(__file__).resolve().parents[3]
    completed = subprocess.run([sys.executable, '-m',
        'model.scalar_observations.experimental_readout', '--check',
        str(FOLDER/'results.json')], cwd=root, capture_output=True, text=True, check=True)
    result = json.loads(completed.stdout)
    assert result['all_pass'] is True
    assert result['calibration_points'] == 21 and result['withheld_points'] == 100
    assert result['conditional_box_violations'] == 0
    assert result['seven_mode_dynamics_tested'] is False
    assert result['coordinate_metadata_discrepancy_resolved'] is False
