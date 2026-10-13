"""Independent integral/dilation oracles and adversarial complete-output checks."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys
from types import ModuleType

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import eval_hermite

from . import halfline, infinite_history, sharp_output as s


def _hermite(n, x):
    return eval_hermite(n, x)*math.exp(-x*x/2)/math.sqrt(2**n*math.factorial(n)*math.sqrt(math.pi))


def _integral_projector(size):
    return np.array([[quad(lambda x: _hermite(m, x)*_hermite(n, x), 0, np.inf,
                           epsabs=3e-13)[0] for n in range(size)] for m in range(size)])


def _numeric(value):
    return complex(*map(float, value))


def _as_fock(vector, size, references):
    result = np.zeros(size*len(references), complex)
    for (n0, n1, gauge), coefficient in vector.items():
        index = size*references.index((n1, gauge))+n0
        result[index] = _numeric(coefficient)*math.sqrt(math.factorial(n0)*math.factorial(n1))
    return result


def _contains(interval, value, tolerance=5e-12):
    return float(interval[0])-tolerance <= value <= float(interval[1])+tolerance


@pytest.mark.parametrize("outcome", (0, 1))
def test_complete_overlap_and_projected_masses_match_independent_coordinate_integrals(outcome):
    source = {(0, 0, 0): (F(2, 3), F(1, 7)), (1, 0, 0): (F(-2, 5), F(3, 4)),
              (3, 0, 0): (F(1, 11), F(-1, 9)), (2, 2, 3): (F(2, 7), F(1, 8)),
              (0, 3, 1): (F(-1, 4), F(2, 9))}
    other = {(2, 0, 0): (F(1, 2), F(-1, 3)), (1, 2, 3): (F(2, 9), F(3, 5)),
             (0, 3, 1): (F(1, 7), F(-4, 9))}
    references = [(0, 0), (2, 3), (3, 1)]
    source_fock = _as_fock(source, 4, references)
    other_fock = _as_fock(other, 4, references)
    positive = _integral_projector(4)
    projector = positive if outcome else np.eye(4)-positive
    operator = np.kron(np.eye(3), projector)
    overlap = np.vdot(source_fock, operator@other_fock)
    exact = s.halfline_overlap_coefficients(source, other, outcome)
    for component, target in ((0, overlap.real), (1, overlap.imag)):
        assert _contains(s.polynomial_interval(tuple(z[component] for z in exact)), target)
    normalization = F(7, 5)
    masses = s.projected_masses(source, 3, outcome, normalization)
    projected = (operator@source_fock).reshape(3, 4)
    # The last spectator has n1=3, outside the rectangular d=3 code.
    inside = np.linalg.norm(projected[:2, :3])**2/float(normalization)
    total = np.vdot(source_fock, operator@source_fock).real/float(normalization)
    assert _contains(s.polynomial_interval(masses["total"]), total)
    assert _contains(s.polynomial_interval(masses["inside"]), inside)
    assert _contains(s.polynomial_interval(masses["outside"]), total-inside)
    assert s._add(masses["inside"], masses["outside"]) == masses["total"]


@pytest.mark.parametrize("outcome", (0, 1))
def test_complete_distance_matches_independent_projector_dilation_and_trace_norm(outcome):
    # The integral matrix A is the compression of the infinite projector.
    # [[A,sqrt(A-A²)],[sqrt(A-A²),I-A]] is an actual finite projector with
    # identical overlaps on this polynomial support. Its extra rows retain
    # the output tail; taking A alone would test the wrong quantum output.
    a = _integral_projector(4)
    eigenvalues, eigenvectors = np.linalg.eigh(a-a@a)
    cross = (eigenvectors*np.sqrt(np.maximum(0., eigenvalues)))@eigenvectors.T
    projector = np.block([[a, cross], [cross, np.eye(4)-a]])
    assert np.linalg.norm(projector@projector-projector) < 1e-12
    if outcome == 0:
        projector = np.eye(8)-projector
    left = {(0, 0, 0): (F(1), F(2, 5)), (3, 0, 0): (F(-1, 7), F(1, 4)),
            (1, 2, 3): (F(3, 8), F(-1, 9))}
    right = {(0, 0, 0): (F(1, 2), F(-1, 3)), (2, 0, 0): (F(-1, 5), F(4, 7)),
             (0, 2, 3): (F(2, 5), F(1, 6))}
    references = [(0, 0), (2, 3)]
    left_fock = _as_fock(left, 4, references).reshape(2, 4)
    right_fock = _as_fock(right, 4, references).reshape(2, 4)
    # Two reference labels and two distinct initial normalizations.
    first = np.hstack((left_fock, np.zeros((2, 4))))/math.sqrt(7/3)
    second = np.hstack((right_fock, np.zeros((2, 4))))/math.sqrt(5/2)
    first = (first@projector.T).reshape(-1)
    second = second.reshape(-1)
    difference = np.outer(first, first.conj())-np.outer(second, second.conj())
    expected = np.sum(np.abs(np.linalg.eigvalsh(difference)))/2
    interval = s.pure_output_distance_interval(left, right, outcome, F(7, 3), F(5, 2))
    assert _contains(interval, expected)
    assert interval[1]-interval[0] < F(1, 10**45)


def test_vacuum_probabilities_hide_large_complete_output_error():
    vacuum = {(0, 0, 0): (F(1), F(0))}
    actual = F(0)
    for outcome in (0, 1):
        finite = {(0, 0, 0): (F(1, 2), F(0)),
                  (1, 0, 0): (F(2*outcome-1, 2), F(0))}
        assert halfline.norm_squared(finite) == F(1, 2)
        assert halfline.positive_weight_coefficients(vacuum) == (F(1, 2), F(0))
        actual += s.pure_output_distance_interval(vacuum, finite, outcome)[0]
    # Equal outcome probabilities give no certificate of output-state fidelity.
    assert F(438, 1000) < actual < F(439, 1000)
    control = s.vacuum_certificate((2, 16, 256))
    assert F(control[0]["complete_output_half_trace_distance_lower"]) > F(385, 1000)
    assert F(control[-1]["complete_output_half_trace_distance_lower"]) > F(124, 1000)


def test_phase_and_reference_coherence_are_checked_without_changing_probabilities():
    left = {(0, 0, 0): (F(1), F(0)), (1, 1, 3): (F(0), F(1, 3))}
    right = {(0, 0, 0): (F(1, 2), F(0)), (1, 0, 0): (F(1, 2), F(0)),
             (1, 1, 3): (F(0), F(1, 6))}
    wrong = dict(right)
    wrong[1, 0, 0] = (F(-1, 2), F(0))
    assert halfline.norm_squared(right) == halfline.norm_squared(wrong)
    good = s.pure_output_distance_interval(left, right, 1)
    bad = s.pure_output_distance_interval(left, wrong, 1)
    assert bad[0] > good[1]+F(1, 10)
    # A common global phase leaves the output density unchanged.
    phase = {key: (-imag, real) for key, (real, imag) in right.items()}
    assert s.pure_output_distance_interval(left, phase, 1) == good


def test_unnormalized_zero_tiny_and_supernormalized_branches_keep_their_mass():
    source = {(0, 0, 0): (F(1), F(0))}
    target = {(0, 0, 0): (F(1, 2), F(0)), (1, 0, 0): (F(1, 2), F(0))}
    base = s.pure_output_distance_interval(source, target, 1)
    for factor in (F(0), F(1, 10**20), F(7)):
        first = {key: tuple(factor*x for x in z) for key, z in source.items()}
        second = {key: tuple(factor*x for x in z) for key, z in target.items()}
        result = s.pure_output_distance_interval(first, second, 1)
        assert result[0] <= factor*factor*base[1]
        assert result[1] >= factor*factor*base[0]
        assert result[1]-result[0] < F(1, 10**50)
    assert s.pure_output_distance_interval({}, {}, 0) == (F(0), F(0))
    with pytest.raises(ValueError):
        s.pure_output_distance_interval(source, target, 1, F(0))


def test_signed_interval_evaluation_and_roots_are_outward_even_near_cancellation():
    low, high = halfline.inverse_sqrt_two_pi_interval()
    center = (low+high)/2
    coefficients = (center*center, -2*center, F(1))
    interval = s.polynomial_interval(coefficients)
    assert interval[0] <= 0 <= interval[1]
    assert interval[1] < F(1, 10**55)
    for squared in (F(0), F(2), F(1, 10**130), F(3, 7)):
        a, b = s.root_interval((squared, squared))
        assert a*a <= squared <= b*b
        assert b-a <= F(1, s.ROOT_SCALE)
    with pytest.raises(ArithmeticError, match="nonnegative"):
        s._nonnegative_interval((F(-1),))


@pytest.mark.parametrize("value", (True, 0, -1, 514, 2., "2", None))
def test_invalid_code_dimensions_fail(value):
    with pytest.raises(ValueError):
        s.projected_masses({}, value, 0)


@pytest.mark.parametrize("value", (True, -1, 2, 0., "0", None))
def test_invalid_signs_fail(value):
    with pytest.raises(ValueError):
        s.halfline_overlap_coefficients({}, {}, value)


@pytest.fixture(scope="module")
def actual_history():
    return infinite_history.certificate(retain_polynomials=True)


def test_actual_history_keeps_complete_records_and_subtracts_the_imported_error(actual_history):
    rows = s.history_support_bounds(actual_history, (2, 16))
    assert [row["feedback_enabled"] for row in rows] == [True, False]
    for row in rows:
        error = F(row["inherited_complete_output_error_upper"])
        assert F(0) < error < F(1, 10**6)
        for entry in row["dimensions"]:
            assert [branch["record"] for branch in entry["retained_records"]] == [[0, 0], [0, 1], [1, 0], [1, 1]]
            polynomial = F(entry["polynomial_complete_output_lower"])
            final = F(entry["infinite_history_complete_output_lower"])
            assert final <= polynomial-error+F(1, 10**15)
            assert F(23, 100) < final < F(4, 10)


def test_actual_two_level_complete_comparison_is_not_an_outcome_only_comparison(actual_history):
    result = s.finite_two_level_comparison(actual_history)
    assert result["onsite_code_dimension"] == 2
    for row in result["configurations"]:
        assert [branch["record"] for branch in row["retained_records"]] == [[0, 0], [0, 1], [1, 0], [1, 1]]
        continuous = tuple(map(F, row["continuous_complete_output_distance_interval"]))
        word = tuple(map(F, row["ideal_source_word_complete_output_distance_interval"]))
        assert F(45, 100) < word[0] <= continuous[0] <= continuous[1] <= word[1] < F(46, 100)
        assert continuous[1]-continuous[0] < F(1, 10**6)
    assert not result["numerical_execution_or_physical_service_error_included"]


@pytest.mark.parametrize("mutation", ("missing_first", "duplicate_first", "missing_configuration",
                                      "hash", "mass", "error", "boolean_first", "negative_error"))
def test_corrupted_retained_history_cannot_certify_a_complete_output(actual_history, mutation):
    modified = deepcopy(actual_history)
    first = modified["configurations"][0]
    branch = first["retained_first_records"][0]
    if mutation == "missing_first":
        first["retained_first_records"].pop()
    elif mutation == "duplicate_first":
        first["retained_first_records"][1]["first_record"] = 0
    elif mutation == "missing_configuration":
        modified["configurations"].pop()
    elif mutation == "hash":
        branch["pre_readout_polynomial"][0][-1] = "17"
    elif mutation == "mass":
        branch["approximate_pre_readout_norm_squared"] = "1"
    elif mutation == "error":
        first["complete_output_half_trace_norm_error_upper"] = "0"
    elif mutation == "boolean_first":
        branch["first_record"] = False
    else:
        branch["second_segment_vector_errors"]["Taylor"] = "-1"
    with pytest.raises(ValueError):
        s.retained_history(modified)


def test_duplicate_polynomial_coordinates_are_rejected_even_with_a_matching_hash(actual_history):
    modified = deepcopy(actual_history)
    branch = modified["configurations"][0]["retained_first_records"][0]
    rows = branch["pre_readout_polynomial"]
    rows.append(rows[0])
    branch["pre_readout_polynomial_sha256"] = hashlib.sha256(
        json.dumps(rows, separators=(",", ":")).encode()).hexdigest()
    with pytest.raises(ValueError, match="duplicate"):
        s.retained_history(modified)


def test_no_certificate_is_created_by_renormalizing_a_bad_approximation():
    coarse = infinite_history.certificate(cutoff=2, preparation_cutoff=0, degree=1,
                                          steps_per_unit=2, bits=32, retain_polynomials=True)
    bounds = s.history_support_bounds(coarse, (2,))
    for configuration in bounds:
        assert F(configuration["inherited_complete_output_error_upper"]) > 1
        assert F(configuration["dimensions"][0]["infinite_history_complete_output_lower"]) == 0


@pytest.fixture
def retained_receipt():
    return json.loads(Path(s.__file__).with_name("sharp_output_results.json").read_text(encoding="utf-8"))


def test_cli_compares_exact_receipt_and_refuses_an_existing_write_path(monkeypatch, tmp_path, capsys, retained_receipt):
    receipt = tmp_path/"receipt.json"
    receipt.write_text(infinite_history.canonical(retained_receipt), encoding="utf-8")
    monkeypatch.setattr(s, "run", lambda: deepcopy(retained_receipt))
    monkeypatch.setattr(sys, "argv", ["sharp_output", "--check", str(receipt)])
    s.main()
    assert json.loads(capsys.readouterr().out)["schema"] == "oph.common_source.sharp_output.v1"
    monkeypatch.setattr(sys, "argv", ["sharp_output", "--write", str(receipt)])
    with pytest.raises(FileExistsError):
        s.main()


@pytest.mark.parametrize("mutation", ("missing_higher_record", "probabilities_only", "source_hash",
                                      "support_lower", "pointer_eta", "evolution_error",
                                      "bool_for_integer", "floating_dimension", "duplicate_key"))
def test_cli_rejects_binding_scientific_receipt_mutations(monkeypatch, tmp_path, retained_receipt, mutation):
    corrupted = deepcopy(retained_receipt)
    larger = corrupted["larger_code_complete_output_comparisons"][0]
    configuration = larger["configurations"][0]
    if mutation == "missing_higher_record":
        configuration["retained_record_comparisons"].pop()
    elif mutation == "probabilities_only":
        configuration.pop("complete_finite_to_infinite_trace_distance_interval")
        configuration["sign_probability"] = "1/2"
    elif mutation == "source_hash":
        corrupted["source_sha256"][next(iter(corrupted["source_sha256"]))] = "0"*64
    elif mutation == "support_lower":
        corrupted["actual_history_code_support_bounds"][0]["dimensions"][0]["infinite_history_complete_output_lower"] = "0"
    elif mutation == "pointer_eta":
        larger["pointer_certificate"]["positive_projector_operator_error_upper"] = "0"
    elif mutation == "evolution_error":
        larger["pre_readout_trace_distance_upper"] = "0"
    elif mutation == "bool_for_integer":
        corrupted["scope"]["native_preparation_current"] = False
    elif mutation == "floating_dimension":
        larger["dimension"] = float(larger["dimension"])
    text = json.dumps(corrupted)
    if mutation == "duplicate_key":
        text = '{"schema":"duplicate",'+text[1:]
    receipt = tmp_path/"corrupted.json"
    receipt.write_text(text, encoding="utf-8")
    monkeypatch.setattr(s, "run", lambda: deepcopy(retained_receipt))
    monkeypatch.setattr(sys, "argv", ["sharp_output", "--check", str(receipt)])
    with pytest.raises(ValueError):
        s.main()


@pytest.mark.parametrize("mutation", ("omit_factorial_metric", "discard_infinite_output_tail",
                                      "omit_initial_normalization"))
def test_independent_scientific_oracles_kill_production_code_mutations(monkeypatch, mutation):
    # Execute a separate, in-memory copy of the production module. None of
    # the tracked source or receipt files is edited by the mutation audit.
    code = Path(s.__file__).read_text(encoding="utf-8")
    if mutation == "omit_factorial_metric":
        old = "F(math.factorial(n1), math.factorial(m))/normalization"
        assert code.count(old) == 1
        code = code.replace(old, "F(1)/normalization")
        check = lambda: test_complete_overlap_and_projected_masses_match_independent_coordinate_integrals(1)
    elif mutation == "discard_infinite_output_tail":
        old = "p = (r, (2*outcome-1)*s)"
        assert code.count(old) == 1
        code = code.replace(old, 'p = projected_masses(infinite_vector, 2, outcome, normalization)["inside"]')
        check = test_vacuum_probabilities_hide_large_complete_output_error
    else:
        old = "h.positive_weight_coefficients(infinite_vector, normalization)"
        assert code.count(old) == 1
        code = code.replace(old, "h.positive_weight_coefficients(infinite_vector, F(1))")
        old = "-F(4)/(normalization*finite_normalization)"
        assert code.count(old) == 1
        code = code.replace(old, "-F(4)/finite_normalization")
        check = lambda: test_complete_distance_matches_independent_projector_dilation_and_trace_norm(1)
    mutant = ModuleType("model.common_source._mutated_sharp_output")
    mutant.__package__ = "model.common_source"
    mutant.__file__ = s.__file__
    exec(compile(code, s.__file__, "exec"), mutant.__dict__)
    monkeypatch.setattr(sys.modules[__name__], "s", mutant)
    with pytest.raises(AssertionError):
        check()
