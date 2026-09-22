import copy
import json
import subprocess
import sys
from itertools import product

import pytest

from .certificates import algebra_certificate, clock_certificate, resource_certificate
from .checker import check_boundary_lowering, check_program, decode_is_lossless, execute, validate
from .compiler import cases, compile_table
from .receipt import EVIDENCE, ROOT, compute, verify


def test_all_finite_reference_functions():
    assert all(len(gate) in (1, 3) for _, program in cases() for gate in program["gates"])
    results = [check_program(p, name in ("retained_copy", "clean_conjunction")) for name, p in cases()]
    assert len(results) == 263
    assert sum(r["interventions"] for r in results) > 1000
    assert all(r["executed_gates"] == r["forward_gates"] + r["inverse_gates"] + r["intervention_gates"] for r in results)


@pytest.mark.parametrize("n", range(6))
def test_xor_into_nonblank_outputs_and_clean_work(n):
    table = [[i % 2, int(i == 0)] for i in range(2**n)]
    p = compile_table(n, 2, table)
    for i in range(2**n):
        source = [(i >> j) & 1 for j in range(n)]
        for output in product((0, 1), repeat=2):
            start = source + list(output) + [0] * p["work"]
            final = execute(start, p["gates"])[-1]
            assert final == source + [a ^ b for a, b in zip(output, table[i])] + [0] * p["work"]


def test_cached_version_substitution_fails_intermediate_interventions():
    logical = dict(cases())["retained_copy"]
    forged = copy.deepcopy(logical)
    forged["gates"][2] = [3, 0, 2]
    # Both functions agree for every initial preparation; that check is insufficient.
    check_program(logical)
    check_program(forged)
    assert check_boundary_lowering(logical, logical)
    with pytest.raises(ValueError, match="intermediate intervention"):
        check_boundary_lowering(logical, forged)


@pytest.mark.parametrize("bad", [None, {}, [], {"n": 0}, True])
def test_junk_circuit_rejected(bad):
    with pytest.raises(ValueError):
        check_program(bad)


@pytest.mark.parametrize("gate", [[True], [-1], [999], [0, 0], [0, 1, 2, 3], [], "NOT", [0.0]])
def test_invalid_operands_rejected(gate):
    p = compile_table(3, 1, [[1]] * 8)
    p["gates"].append(gate)
    with pytest.raises(ValueError):
        check_program(p)


@pytest.mark.parametrize("mutation", ["empty", "drop", "retarget", "dirty", "table", "bool", "extra"])
def test_plausible_but_false_programs_rejected(mutation):
    p = compile_table(4, 1, [[int(i == 15)] for i in range(16)])
    if mutation == "empty": p["gates"] = []
    elif mutation == "drop": p["gates"].pop()
    elif mutation == "retarget": p["gates"][0] = [4]
    elif mutation == "dirty": p["gates"].append([5])
    elif mutation == "table": p["table"][15] = [0]
    elif mutation == "bool": p["n"] = True
    else: p["trusted"] = True
    with pytest.raises(ValueError):
        check_program(p)


def test_compiler_rejects_incomplete_nonbinary_tables():
    for table in ([], [[1]], [[True], [0]], [[2], [0]], [[0, 1], [1, 0]]):
        with pytest.raises(ValueError):
            compile_table(1, 1, table)


def test_exact_clock_algebra_and_resource_controls():
    assert clock_certificate()["invariant_quadratic_dimension"] == 1
    assert algebra_certificate()["matrix_units_checked"] == 1728
    assert len(resource_certificate()) == 9


def test_lossless_tuple_and_additive_collision():
    words = list(product((0, 1), repeat=5))
    assert decode_is_lossless(words, tuple, tuple)
    assert not decode_is_lossless(words, sum, lambda total: (1,) * total + (0,) * (5-total))
    assert sum((0, 1)) == sum((1, 0))


def test_incoming_cut_excludes_the_locally_available_record():
    odd = next(row for row in resource_certificate() if row["n"] == 3)
    assert odd["reads_per_receiver"] == 27
    assert odd["remote_bits_per_receiver"] == 26
    assert odd["remote_bits_per_receiver"] < odd["n"]**3


@pytest.fixture(scope="module")
def replay():
    return compute()


def test_committed_receipt(replay):
    assert verify(EVIDENCE) == replay


@pytest.mark.parametrize("mutation", ["empty", "source", "case", "count", "bool", "phase", "clock", "resource", "extra"])
def test_forged_receipts_rejected(tmp_path, replay, mutation):
    forged = copy.deepcopy(replay)
    if mutation == "empty": forged = {}
    elif mutation == "source": forged["sources"]["code/rg_principle/checker.py"] = "0" * 64
    elif mutation == "case": forged["circuits"].pop("boolean3_000")
    elif mutation == "count": forged["circuits"]["retained_copy"]["executed_gates"] = 0
    elif mutation == "bool": forged["clock"]["invariant_quadratic_dimension"] = True
    elif mutation == "phase": forged["algebra"]["phase_effect_probabilities"] = [0, 0]
    elif mutation == "clock": forged["clock"]["polyhedral_equal_radius_squared_gap"] = "0"
    elif mutation == "resource": forged["resources"][0]["reads_per_receiver"] = 1
    else: forged["physical_source_certified"] = True
    path = tmp_path / "forged.json"
    path.write_text(json.dumps(forged), encoding="utf-8")
    with pytest.raises(ValueError):
        verify(path)


@pytest.mark.parametrize("content", ['{"schema":0,"schema":1}', '{"x":NaN}', '[]', '{'])
def test_invalid_receipt_json(tmp_path, content):
    path = tmp_path / "invalid.json"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValueError):
        verify(path)


def test_cli_rejects_junk_even_with_python_optimization(tmp_path):
    path = tmp_path / "empty.json"
    path.write_text("{}", encoding="utf-8")
    result = subprocess.run([sys.executable, "-O", "code/rg_principle/receipt.py", "verify", "--path", str(path)],
                            cwd=ROOT, capture_output=True, text=True)
    assert result.returncode != 0
    assert "receipt differs" in result.stderr
