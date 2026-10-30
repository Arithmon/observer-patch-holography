"""Hostile finite/algebraic controls; these do not test all continuous maps."""

import copy
import json
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import equivariant_degree_crosscheck as producer
import verify_equivariant_degree_crosscheck as independent


@pytest.fixture(scope="module")
def packet():
    return json.loads(producer.RECEIPT.read_text())


def test_independent_e8_character_and_basepoint_controls(packet):
    result = independent.verify(packet)
    assert result["minimum_homogeneous_degrees"] == {"Q": 1, "Q_prime": 7}
    assert result["checked_through"] == 94


def test_common_orientation_reversal_and_marked_generator_square(packet):
    for residues in packet["necessary_residues_mod_30"].values():
        assert sorted((-x) % 30 for x in residues) == residues
    for prime in (2, 3, 5):
        for a in range(1, prime):
            for b in range(1, prime):
                for k in range(1, prime):
                    assert (b*pow(a, -1, prime)) % prime == ((k*b)*pow(k*a, -1, prime)) % prime
    assert 2*pow(1, -1, 5) % 5 == 2
    assert 4*pow(2, -1, 5) % 5 == 2  # coherent generator square
    assert 4*pow(1, -1, 5) % 5 != 2  # target-only square changes the typing


@pytest.mark.parametrize("mutation", ["wrong_Qprime", "wrong_class", "omit_lift", "multiplicity",
    "wrong_numerator", "wrong_graph", "common_zero", "ambient_sign", "universal_promotion",
    "physical_promotion", "false_provenance", "implementation_hash"])
def test_rejects_algebraic_and_scope_forgery(packet, mutation):
    bad = copy.deepcopy(packet)
    if mutation == "wrong_Qprime":
        bad["characters"][7] = bad["characters"][1]
    elif mutation == "wrong_class":
        bad["class_Q_traces"][5] = [-1, 1]
    elif mutation == "omit_lift":
        bad["class_sizes"][1] = 0
    elif mutation == "multiplicity":
        bad["multiplicities"]["Q_prime"][3] = 1
    elif mutation == "wrong_numerator":
        bad["hilbert_numerator_exponents"]["Q_prime"] = [1, 11, 19, 29]
    elif mutation == "wrong_graph":
        bad["mckay_adjacency"][6][7] = 0
    elif mutation == "common_zero":
        bad["hms"]["denominator"] = bad["hms"]["numerator"]
    elif mutation == "ambient_sign":
        bad["minimum_absolute_degree_from_congruences"]["twin"] = 1
    elif mutation == "universal_promotion":
        bad["continuous_theorem"] = "CONTINUOUS_TWIN_DEGREE_FLOOR_PROVED"
    elif mutation == "physical_promotion":
        bad["physical_selection"] = True
    elif mutation == "false_provenance":
        bad["bernhard_original_proof_and_scripts"] = "CHECKED"
    else:
        bad["implementation_sha256"]["equivariant_degree_crosscheck.py"] = "0"*64
    with pytest.raises(ValueError):
        independent.verify(bad)
