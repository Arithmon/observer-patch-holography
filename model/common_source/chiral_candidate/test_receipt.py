"""Custody and state-weight checks for the retained chiral source experiment."""
from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from . import run


HERE = Path(__file__).resolve().parent


def test_retained_candidate_verifies():
    run.verify(HERE/"results.json")


@pytest.mark.parametrize("damage",["state_bytes","weight","nonfinite","source","mass_claim"])
def test_damaged_chiral_candidate_is_rejected(tmp_path,damage):
    document = deepcopy(json.loads((HERE/"results.json").read_text()))
    artifact = tmp_path/document["artifact"]["filename"]
    shutil.copyfile(HERE/document["artifact"]["filename"],artifact)
    if damage == "state_bytes":
        with artifact.open("ab") as handle:
            handle.write(b"changed")
    elif damage in ("weight","nonfinite"):
        document["continuing_source_history"]["source_record_probabilities"]["000000"] = (
            .5 if damage == "weight" else float("nan"))
    elif damage == "source":
        del document["source_bindings"]["model/common_source/chiral_candidate/history.py"]
    else:
        document["verdict"]["observed_standard_model_mass_spectrum_derived"] = True
    path = tmp_path/"damaged.json"
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError):
        run.verify(path)
