"""Receipt checking never trusts an uploaded bound or drops unexpected fields."""
import copy
import json
from pathlib import Path
import sys

import pytest

from . import thermal_certificate as certificate


@pytest.fixture
def reference(monkeypatch):
    # This tests the CLI acceptance boundary. Mathematical tests independently
    # exercise the producer; model.verify replays the full expensive receipt.
    result=json.loads(Path(__file__).with_name('thermal_results.json').read_text())
    monkeypatch.setattr(certificate,'run',lambda: copy.deepcopy(result))
    return result


@pytest.mark.parametrize('mutation',[
    lambda r:r['cold_protocol']['histories'].pop(),
    lambda r:r['cold_protocol']['histories'][0]['record_alphabet'].pop(),
    lambda r:r['cold_protocol']['histories'][0].update(retained_reference_dimension=True),
    lambda r:r.update(complete_source_to_infinite_output_error_upper_under_numerical_guards='0'),
    lambda r:r['boundary'].update(physical_seconds_or_universal_empirical_agreement_established=True),
    lambda r:r['source_sha256'].clear(),
    lambda r:r.update(extra_claim='all premises proved'),
    lambda r:r['source_execution']['histories'][0]['probabilities'].__setitem__(0,float('nan')),
])
def test_semantic_mutations_fail_closed(reference,mutation,tmp_path,monkeypatch):
    bad=copy.deepcopy(reference);mutation(bad)
    path=tmp_path/'bad.json';path.write_text(json.dumps(bad))
    monkeypatch.setattr(sys,'argv',['thermal','--check',str(path)])
    with pytest.raises((ValueError,ArithmeticError)):
        certificate.main()


def test_duplicate_fields_and_overwrite_fail(reference,tmp_path,monkeypatch):
    path=tmp_path/'bad.json';path.write_text('{"schema":"one","schema":"two"}')
    monkeypatch.setattr(sys,'argv',['thermal','--check',str(path)])
    with pytest.raises(ValueError):certificate.main()
    monkeypatch.setattr(sys,'argv',['thermal','--write',str(path)])
    with pytest.raises(FileExistsError):certificate.main()


def test_outward_display_does_not_round_down():
    from fractions import Fraction as F
    for value in (F(1,3),F(-1,3),F(1,10**30),F(17,8)):
        assert F(certificate.upper(value))>=value
    with pytest.raises(ValueError):certificate.upper(True)
