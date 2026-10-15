import json

import pytest

from .__main__ import compare, read_receipt


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -float('inf')])
@pytest.mark.parametrize('side', ['actual', 'expected', 'both'])
def test_nonfinite_values_cannot_pass_numeric_replay(value, side):
    actual = value if side in ('actual', 'both') else .2
    expected = value if side in ('expected', 'both') else .2
    with pytest.raises(ValueError, match='nonfinite'):
        compare({'runs': [{'bound': actual}]}, {'runs': [{'bound': expected}]})


@pytest.mark.parametrize('constant', ['NaN', 'Infinity', '-Infinity', '1e999', '-1e999'])
def test_loader_rejects_nonfinite_json_numbers(tmp_path, constant):
    path = tmp_path/'receipt.json'
    path.write_text('{"bound": '+constant+'}')
    # Python's default decoder accepts these. The receipt decoder must not.
    assert 'bound' in json.loads(path.read_text())
    with pytest.raises(ValueError, match='[Nn]onfinite'):
        read_receipt(path)


def test_finite_replay_still_checks_changes():
    compare({'bound': .125}, {'bound': .125})
    with pytest.raises(ValueError, match='numerical result changed'):
        compare({'bound': .126}, {'bound': .125})


@pytest.mark.parametrize('payload', [
    '{"bound": 0.1, "bound": 0.2}',
    '{"runs": [{"bound": 0.1, "bound": 0.2}]}',
])
def test_loader_rejects_duplicate_keys_at_every_depth(tmp_path, payload):
    path = tmp_path/'receipt.json'
    path.write_text(payload)
    with pytest.raises(ValueError, match='duplicate JSON key'):
        read_receipt(path)


def test_finite_producer_output_roundtrips_through_receipt_loader(tmp_path):
    from .experiment import run
    result = run(dimension=8, steps=0)
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(result, allow_nan=False))
    compare(result, read_receipt(path))
