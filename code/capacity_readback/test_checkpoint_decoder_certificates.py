"""Independent integer enumeration and hostile replay of decoder evidence."""
from copy import deepcopy
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import subprocess
import sys

import pytest

from correctable_public_record_capacity import (
    approximate_public_capacity, compound_confusability_graph,
    maximum_independent_set, support_relation_semigroup, tv_robustness_bound,
)
from verify_checkpoint_decoder import verify_decoder_witness


def channel(counts):
    return {"rows": {str(x): {str(y): Fraction(n, sum(row)) for y, n in enumerate(row)}
                     for x, row in enumerate(counts)}}


def integer_oracle(family, epsilon):
    """Flat enumeration using integer counts, including unassigned outputs."""
    inputs, outputs = len(family[0]), len(family[0][0])
    for size in range(inputs, 0, -1):
        for code in combinations(range(inputs), size):
            for counts in family:
                good = False
                for decoder in product((-1, *code), repeat=outputs):
                    if all(Fraction(sum(n for y, n in enumerate(counts[x]) if decoder[y] != x),
                                    sum(counts[x])) <= epsilon for x in code):
                        good = True
                        break
                if not good:
                    break
            else:
                return size
    raise AssertionError("valid kernels must admit a singleton")


BINARY = tuple(((a, 4-a), (b, 4-b)) for a, b in product(range(5), repeat=2))
TERNARY_ROWS = tuple(row for row in product(range(3), repeat=3) if sum(row) == 2)


@pytest.mark.parametrize("first,second", tuple(combinations_with_replacement(range(len(BINARY)), 2)))
def test_every_pair_of_binary_quarter_count_channels(first, second):
    family = [BINARY[first], BINARY[second]]
    kernels = [channel(counts) for counts in family]
    for epsilon in [Fraction(k, 4) for k in range(5)]:
        result = approximate_public_capacity(["0", "1"], kernels, epsilon)
        assert result["capacity"] == integer_oracle(family, epsilon)
        assert verify_decoder_witness(["0", "1"], kernels, epsilon,
                                      json.loads(json.dumps(result)))


@pytest.mark.parametrize("rows", tuple(product(TERNARY_ROWS, repeat=3)))
def test_every_three_input_three_output_half_count_channel(rows):
    kernels = [channel(rows)]
    previous = 0
    for epsilon in [Fraction(0), Fraction(1, 2), Fraction(1)]:
        result = approximate_public_capacity(["0", "1", "2"], kernels, epsilon)
        assert result["capacity"] == integer_oracle([rows], epsilon)
        assert result["capacity"] >= previous
        previous = result["capacity"]
        assert verify_decoder_witness(["0", "1", "2"], kernels, epsilon, result)
        if epsilon == 0:
            graph = compound_confusability_graph(["0", "1", "2"], kernels)
            assert result["capacity"] == len(maximum_independent_set(graph))


@pytest.mark.parametrize("size", [1, 3, 12, 24])
def test_reversible_codes_need_no_exponential_decoder_enumeration(size):
    records = [str(i) for i in range(size)]
    kernels = [{"rows": {str(i): {str((i+shift) % size): 1} for i in range(size)}}
               for shift in [0, 1]]
    result = approximate_public_capacity(records, kernels, Fraction(1, 10),
                                        max_vertices=size, max_decoder_nodes=2)
    assert result["capacity"] == size
    assert verify_decoder_witness(records, kernels, Fraction(1, 10), result)
    # Identity and the known cyclic shift have different inverse decoders.
    if size > 1:
        assert result["decoder_certificates"][0]["decoder"] != result["decoder_certificates"][1]["decoder"]


def test_complete_retained_source_family_has_replayable_24_record_decoders():
    from source_derived_public_checkpoint_packet import build_source_derived_packet
    packet = build_source_derived_packet()
    records = sorted(packet["reachability_witnesses"])
    kernels = packet["global_checkpoint_kernels"]
    assert len(kernels) == 40
    result = approximate_public_capacity(records, kernels, 0,
                                        max_vertices=24, max_decoder_nodes=40)
    assert result["capacity"] == 24
    assert verify_decoder_witness(records, kernels, 0, json.loads(json.dumps(result)))


def test_replay_scope_is_attainability_not_maximality():
    # A singleton is achievable in a two-record identity channel. Its decoder
    # is valid evidence of that fact, but is not evidence that capacity is one.
    records = ["0", "1"]
    kernels = [channel(((1, 0), (0, 1)))]
    singleton = [{"rows": {"0": kernels[0]["rows"]["0"]}}]
    result = approximate_public_capacity(["0"], singleton, 0)
    assert result["capacity"] == 1
    assert verify_decoder_witness(records, kernels, 0, result)
    assert approximate_public_capacity(records, kernels, 0)["capacity"] == 2


def test_exhausted_search_makes_no_capacity_claim():
    with pytest.raises(ValueError, match="budget exhausted"):
        approximate_public_capacity(["0", "1"], [channel(((1, 1), (1, 1)))],
                                    Fraction(1, 2), max_decoder_nodes=1)


def read_ratio(item):
    return Fraction(int(item["numerator"]), int(item["denominator"]))


def test_near_unit_float_row_interpretation_is_disclosed_exactly():
    kernels = [{"rows": {"0": {"0": .9, "1": .1}, "1": {"0": .1, "1": .9}}}]
    result = approximate_public_capacity(["0", "1"], kernels, .1)
    certificate = result["decoder_certificates"][0]
    mass = Fraction(.9)+Fraction(.1)
    assert mass != 1
    assert read_ratio(certificate["input_row_masses"]["0"]) == mass
    assert read_ratio(certificate["worst_input_error"]) == Fraction(.1)/mass
    assert verify_decoder_witness(["0", "1"], kernels, .1, result)


def test_unrepresentable_positive_error_is_not_displayed_as_perfect_success():
    e = Fraction(1, 10**400)
    kernels = [{"rows": {"0": {"0": 1-e, "1": e}, "1": {"0": e, "1": 1-e}}}]
    result = approximate_public_capacity(["0", "1"], kernels, e)
    assert result["capacity"] == 2
    assert 0 < result["worst_input_error_by_channel"][0]
    assert result["worst_input_success_by_channel"][0] < 1
    assert read_ratio(result["decoder_certificates"][0]["worst_input_error"]) == e
    assert verify_decoder_witness(["0", "1"], kernels, e, result)
    assert approximate_public_capacity(["0", "1"], kernels, tv_robustness_bound(0, e))["capacity"] == 2


def test_relation_closure_retains_delayed_erasure():
    records = ["0", "1", "2"]
    kernel = {"rows": {"0": {"1": 1}, "1": {"2": 1}, "2": {"2": 1}}}
    single = compound_confusability_graph(records, [kernel])
    assert len(maximum_independent_set(single)) == 2
    relations = support_relation_semigroup(records, [kernel])
    assert relations == {frozenset({("0", "1"), ("1", "2"), ("2", "2")}),
                         frozenset((x, "2") for x in records)}
    # After two checkpoints every record has the same output: capacity one.


@pytest.mark.parametrize("mutation", ["constant_pass", "empty_code", "duplicate_code", "wrong_capacity",
    "wrong_decoder", "missing_output", "wrong_error", "wrong_mass", "wrong_epsilon",
    "missing_channel", "infinite_display", "false_perfect_display", "wrong_scope"])
def test_independent_replay_rejects_forged_certificates(mutation):
    kernels = [channel(((3, 1), (1, 3)))]
    epsilon = Fraction(1, 4)
    result = approximate_public_capacity(["0", "1"], kernels, epsilon)
    changed = deepcopy(result)
    certificate = changed["decoder_certificates"][0]
    zero = {"numerator": "0", "denominator": "1"}
    if mutation == "constant_pass": changed = {"status": "PASS"}
    elif mutation == "empty_code": changed["code_witness"] = []
    elif mutation == "duplicate_code": changed["code_witness"] = ["0", "0"]
    elif mutation == "wrong_capacity": changed["capacity"] = 99
    elif mutation == "wrong_decoder": certificate["decoder"] = {"0": "1", "1": "0"}
    elif mutation == "missing_output": certificate["decoder"].pop("0")
    elif mutation == "wrong_error": certificate["worst_input_error"] = zero
    elif mutation == "wrong_mass": certificate["input_row_masses"]["0"] = zero
    elif mutation == "wrong_epsilon": changed["epsilon_exact"] = zero
    elif mutation == "missing_channel": changed["decoder_certificates"] = []
    elif mutation == "infinite_display": changed["worst_input_error_by_channel"] = [float("inf")]
    elif mutation == "false_perfect_display": changed["worst_input_success_by_channel"] = [1.]
    elif mutation == "wrong_scope": changed["certificate_scope"] = "optimal_capacity_proved"
    assert not verify_decoder_witness(["0", "1"], kernels, epsilon, changed)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1, 2, True, "0.1", [], None])
def test_invalid_error_tolerances_are_rejected(bad):
    kernels = [channel(((1, 0), (0, 1)))]
    with pytest.raises(ValueError):
        approximate_public_capacity(["0", "1"], kernels, bad)
    with pytest.raises(ValueError):
        tv_robustness_bound(0, bad)
    assert not verify_decoder_witness(["0", "1"], kernels, bad, {})


def test_verifier_replays_with_producer_absent_under_optimized_python(tmp_path):
    import shutil
    shutil.copy2(Path(__file__).with_name("verify_checkpoint_decoder.py"), tmp_path/'verify_checkpoint_decoder.py')
    channels = [{"rows": {"0": {"0": .75, "1": .25}, "1": {"0": .25, "1": .75}}}]
    result = approximate_public_capacity(["0", "1"], channels, .25)
    (tmp_path/'input.json').write_text(json.dumps([channels, result]))
    script = """
import json, sys
from verify_checkpoint_decoder import verify_decoder_witness
channels, result = json.load(open('input.json'))
if not verify_decoder_witness(['0','1'], channels, .25, result):
    raise SystemExit('valid witness rejected')
result['decoder_certificates'][0]['decoder'] = {'0':'0','1':'0'}
if verify_decoder_witness(['0','1'], channels, .25, result):
    raise SystemExit('forged decoder accepted')
if 'checkpoint_channels' in sys.modules or 'correctable_public_record_capacity' in sys.modules:
    raise SystemExit('producer imported')
"""
    completed = subprocess.run([sys.executable, '-O', '-c', script], cwd=tmp_path,
                               capture_output=True, text=True, timeout=20)
    assert completed.returncode == 0, completed.stdout+completed.stderr
