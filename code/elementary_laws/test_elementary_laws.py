"""Mutation controls for the elementary-law receipt."""
import copy
import importlib.util
from fractions import Fraction as Q
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


producer = module("elementary_laws_producer", HERE / "build.py")
verifier = module("elementary_laws_verifier", HERE / "verify.py")


def test_rebuild_matches_committed_receipt_and_verifies():
    receipt = producer.build()
    assert receipt == verifier.load(HERE / "receipt.json")
    assert verifier.verify(receipt)["ok"]


def test_wrong_aberration_sign_is_rejected():
    receipt = producer.build()
    row = receipt["celestial"]["boosts"][0]["directions"][1]
    ct, v = Q(row["cos_theta"]), Q(receipt["celestial"]["boosts"][0]["velocity"])
    row["cos_theta_boosted"] = str((ct + v) / (1 + v * ct))
    with pytest.raises(ValueError, match="aberration"):
        verifier.verify(receipt)


def test_swapped_statistics_are_rejected():
    receipt = producer.build()
    row = receipt["occupation"][0]
    row["fermi_dirac"], row["bose_einstein"] = row["bose_einstein"], row["fermi_dirac"]
    with pytest.raises(ValueError, match="occupation"):
        verifier.verify(receipt)


def test_zeno_survival_below_bound_is_rejected():
    receipt = producer.build()
    row = receipt["lueders_zeno"]["zeno_sequence"][0]
    row["one_minus_survival"] = str(Q(4, row["n"]) + 1)
    with pytest.raises(ValueError, match="zeno bound"):
        verifier.verify(receipt)


def test_entropy_increasing_repair_is_rejected():
    receipt = copy.deepcopy(producer.build())
    block = receipt["coarse_graining"]
    block["collision_sum_after"] = str(Q(block["collision_sum_before"]) - Q(1, 64))
    with pytest.raises(ValueError, match="collision"):
        verifier.verify(receipt)


def test_fragment_pin_is_enforced():
    receipt = producer.build()
    receipt["fragment_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="fragment pin"):
        verifier.verify(receipt)
