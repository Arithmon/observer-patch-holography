"""Hostile receipt checks retain observables while allowing counted roundoff."""
from copy import deepcopy
import math

import pytest

from .experiment import ROOT, source_bindings
from .replay import compare_band_receipt
from ..interacting_field.replay import LOCAL_RESIDUAL_REPLAY_ATOL
from tools import strict_json


@pytest.fixture
def packet():
    return strict_json.load(ROOT/'model/band_clock/results.json')


def test_count_scaled_roundoff_is_checked_without_widening_observables(packet):
    changed = deepcopy(packet)
    resources = changed['source_execution']['requested_resources']
    calls = sum(resources[key] for key in ('completed_hopping_blocks',
                'completed_contact_blocks', 'routing_FSWAP_calls', 'fresh_pointer_blanks'))
    resources['floating_local_residual_sum_with_routing'] = calls*LOCAL_RESIDUAL_REPLAY_ATOL/2
    resources['floating_residual_telescoping_estimate'] = math.expm1(
        resources['floating_local_residual_sum_with_routing'])
    compare_band_receipt(packet, changed)
    changed['source_execution']['complete_pointer_probabilities'][0] += 1e-7
    with pytest.raises(ValueError):
        compare_band_receipt(packet, changed)


@pytest.mark.parametrize('mutation', ['count', 'bool_count', 'negative_residual',
        'cap', 'inconsistent_expm1', 'nonfinite', 'scope', 'custody'])
def test_mutated_receipts_fail_closed(packet, mutation):
    changed = deepcopy(packet)
    resources = changed['source_execution']['requested_resources']
    if mutation == 'count':
        resources['routing_FSWAP_calls'] += 1
    elif mutation == 'bool_count':
        resources['fresh_pointer_blanks'] = True
    elif mutation == 'negative_residual':
        resources['floating_local_residual_sum_with_routing'] = -1e-12
    elif mutation == 'cap':
        resources['floating_local_residual_sum_with_routing'] = 1e-5
        resources['floating_residual_telescoping_estimate'] = math.expm1(1e-5)
    elif mutation == 'inconsistent_expm1':
        resources['floating_residual_telescoping_estimate'] = 0.
    elif mutation == 'nonfinite':
        resources['floating_local_residual_sum_with_routing'] = float('nan')
    elif mutation == 'scope':
        changed['unfulfilled_joins']['physical_time_energy_mass_calibration'] = True
    else:
        key = next(iter(changed['source_bindings']))
        changed['source_bindings'][key] = '0'*64
    with pytest.raises(ValueError):
        compare_band_receipt(packet, changed)


def test_check_cli_uses_count_scaled_checker_and_rejects_changed_observable(packet, tmp_path, monkeypatch):
    import json
    import sys
    from . import __main__ as cli
    monkeypatch.setattr(cli, 'run', lambda: packet)
    changed = deepcopy(packet)
    changed['source_execution']['complete_pointer_probabilities'][0] += .01
    receipt = tmp_path/'hostile.json'
    receipt.write_text(json.dumps(changed))
    monkeypatch.setattr(sys, 'argv', ['band_clock', '--check', str(receipt)])
    with pytest.raises(ValueError):
        cli.main()


def test_custody_covers_transitive_local_imports_and_eager_initializers():
    bound = source_bindings()
    for path in ('model/bound_pair/__init__.py', 'model/interacting_field/__init__.py',
                 'model/particle_clock/__init__.py', 'model/local_repair_family/learned_flow.py',
                 'model/local_repair_family/modular_clock.py', 'model/band_clock/replay.py',
                 'model/interacting_field/replay.py', 'tools/strict_json.py'):
        assert path in bound
