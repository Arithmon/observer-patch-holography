"""Checks of the retained full reference-bearing field output artifact."""
from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from . import run


HERE = Path(__file__).resolve().parent


def test_retained_particle_audit_and_quantum_outputs_verify():
    run.verify(HERE/"results.json")


@pytest.mark.parametrize("damage", ["state_bytes", "probability", "nonfinite", "missing_source", "physical_promotion", "kinetic_promotion"])
def test_damaged_particle_output_is_rejected(tmp_path, damage):
    document = deepcopy(json.loads((HERE/"results.json").read_text()))
    artifact = tmp_path/document["artifact"]["filename"]
    shutil.copyfile(HERE/document["artifact"]["filename"], artifact)
    if damage == "state_bytes":
        with artifact.open("ab") as handle:
            handle.write(b"changed")
    elif damage in ("probability", "nonfinite"):
        document["continuing_field_read"]["source_record_probabilities"]["000000"] = (
            .5 if damage == "probability" else float("nan"))
    elif damage == "missing_source":
        del document["source_bindings"]["model/common_source/particle_identification/response.py"]
    elif damage == "physical_promotion":
        document["verdict"]["observed_standard_model_species_identified"] = ["electron"]
    elif damage == "kinetic_promotion":
        document["finite_inverse_response_kinetics"]["response_inferred_from_observer_records"] = True
    path = tmp_path/"receipt.json"
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError):
        run.verify(path)
