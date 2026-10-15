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


@pytest.mark.parametrize('constant', ['NaN', 'Infinity', '-Infinity'])
def test_loader_rejects_nonstandard_json_nonfinite_constants(tmp_path, constant):
    path = tmp_path/'receipt.json'
    path.write_text('{"bound": '+constant+'}')
    # Python's default decoder accepts these. The receipt decoder must not.
    assert 'bound' in json.loads(path.read_text())
    with pytest.raises(ValueError, match='Nonfinite JSON'):
        read_receipt(path)


def test_finite_replay_still_checks_changes():
    compare({'bound': .125}, {'bound': .125})
    with pytest.raises(ValueError, match='numerical result changed'):
        compare({'bound': .126}, {'bound': .125})
