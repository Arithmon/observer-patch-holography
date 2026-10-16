"""Receipt replay must recompute science, not trust success flags or digests."""
import copy
from fractions import Fraction as F
import json
import subprocess
import sys

import pytest

from ..observer.scalar_response.experiment import compare_receipt
from . import receipt


@pytest.fixture(scope='module')
def produced():
    return receipt.run()


def test_actual_producer_roundtrips_strict_replay(produced):
    restored = json.loads(json.dumps(produced, allow_nan=False))
    compare_receipt(restored, produced)
    assert restored['elementary_exact_instance']['complete_record_count'] == 64
    assert F(restored['elementary_exact_instance']['baseline_feedback_total_variation']) == F(5, 99)
    assert F(restored['effective_logarithm']['federation_operator_remainder_upper']) < F(1, 10**11)


@pytest.mark.parametrize('mutate', [
    lambda r: r['elementary_exact_instance'].__setitem__('states', 1),
    lambda r: r['elementary_exact_instance'].__setitem__('complete_record_count', True),
    lambda r: r['elementary_exact_instance'].__setitem__('baseline_feedback_total_variation', '0'),
    lambda r: r['effective_logarithm'].__setitem__('federation_operator_remainder_upper', '0'),
    lambda r: r.__setitem__('source_sha256', {}),
    lambda r: r.pop('federation_complete_instance'),
])
def test_scientific_receipt_tampering_is_rejected(produced, mutate):
    hostile = copy.deepcopy(produced)
    mutate(hostile)
    with pytest.raises((ValueError, AssertionError)):
        compare_receipt(hostile, produced)


def test_cli_writes_real_production_types_then_checks_them(tmp_path):
    target = tmp_path/'receipt.json'
    argv = [sys.executable, '-m', 'model.observer_time_fields']
    subprocess.run([*argv, '--write', str(target)], check=True, capture_output=True, text=True)
    subprocess.run([*argv, '--check', str(target)], check=True, capture_output=True, text=True)
    data = json.loads(target.read_text(encoding='utf-8'))
    data['effective_logarithm']['degree'] = 0
    target.write_text(json.dumps(data), encoding='utf-8')
    rejected = subprocess.run([*argv, '--check', str(target)], capture_output=True, text=True)
    assert rejected.returncode != 0
    assert 'degree' in rejected.stderr


@pytest.mark.parametrize('value', [float('nan'), 0.0, {1: 'wrong key'}, object()])
def test_exact_receipt_does_not_silently_coerce_unknown_values(value):
    with pytest.raises(TypeError):
        receipt.json_value(value)
