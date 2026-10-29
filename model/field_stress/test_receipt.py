"""Receipt source custody and explicit conditional-premise boundaries."""
from pathlib import Path
from tools import strict_json
from .experiment import source_bindings


def test_retained_receipt_binds_actual_sources_and_theorem_consumers():
    receipt = strict_json.load(Path(__file__).with_name('results.json'))
    assert receipt['source_bindings'] == source_bindings()
    assert receipt['all_pass']
    assert not receipt['EinsteinRegisterData_inhabited']
    assert not receipt['integrated_gate_accepted']
    assert receipt['physical_comparison_count_added'] == 0
    assert set(receipt['premise_attachment']) == {f'PR-{i}' for i in range(23,29)}
    assert receipt['energy_continuity']['operator_identities_exact']
    assert not receipt['same_action_entropy']['density_preparation_source_executed']
    assert not receipt['same_action_entropy']['physical_clock_or_temperature_calibrated']
