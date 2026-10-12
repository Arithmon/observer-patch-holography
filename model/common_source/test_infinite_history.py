"""Independent Hamiltonian and adversarial checks of the infinite history."""
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import halfline, infinite_history as history, scalar_coupling


def _rational(vector):
    return {key: (F(real, vector.denominator), F(imag, vector.denominator))
            for key, (real, imag) in vector.coefficients.items()}


def _norm_squared(vector):
    return halfline.norm_squared(_rational(vector))


@pytest.mark.parametrize("dimension", (2, 4))
@pytest.mark.parametrize("coupling", (F(0), F(1, 4)))
def test_full_hamiltonian_compresses_to_the_existing_common_source_action(dimension, coupling):
    # Normalize the Bargmann monomials before comparing with the separately
    # constructed quadrature/CAR/electric matrix in the retained model.
    size = 4 * dimension**2
    actual = np.zeros((size, size), complex)
    for n0 in range(dimension):
        for n1 in range(dimension):
            for gauge in range(4):
                column = gauge + 4 * n0 + 4 * dimension * n1
                source_weight = math.factorial(n0) * math.factorial(n1)
                images = history._hamiltonian_numerator({(n0, n1, gauge): (1, 0)}, coupling)
                for (m0, m1, output), (real, imag) in images.items():
                    if m0 >= dimension or m1 >= dimension:
                        continue
                    target_weight = math.factorial(m0) * math.factorial(m1)
                    row = output + 4 * m0 + 4 * dimension * m1
                    actual[row, column] = complex(real, imag) / 128 * math.sqrt(
                        target_weight / source_weight)
    assert np.linalg.norm(actual - scalar_coupling.ideal_hamiltonian(dimension, float(coupling)), 2) < 1e-13
    assert np.linalg.norm(actual - actual.conj().T, 2) < 1e-14


@pytest.mark.parametrize("cutoff", (0, 2, 4, 8))
def test_projected_squeezed_preparation_has_the_exact_common_normalization(cutoff):
    polynomial = history._initial(cutoff)
    vector = _rational(polynomial)
    for (n0, n1, gauge), coefficient in vector.items():
        pair = (n0 + n1) // 2
        assert gauge == 0 and (n0 + n1) % 2 == 0 and n0 + n1 <= cutoff
        assert coefficient == (F((-1)**(pair + n1) * math.comb(2 * pair, n0),
                                 36**pair * math.factorial(pair)), F(0))
    expected_norm = sum((F(math.comb(2 * pair, pair), 324**pair)
                         for pair in range(cutoff // 2 + 1)), F(0))
    assert _norm_squared(polynomial) == expected_norm
    assert F(history._norm_numerator(polynomial.coefficients), polynomial.denominator**2) == expected_norm
    # Finite polynomials are NOT renormalized to one in their coefficients.
    assert vector[0, 0, 0] == (F(1), F(0))
    if cutoff:
        assert expected_norm > 1


def test_taylor_power_keeps_excursions_out_of_and_back_into_the_retained_space():
    vacuum = history._Polynomial({(0, 0, 0): (1, 0)})
    first = history._hamiltonian_numerator(vacuum.coefficients)
    second = history._hamiltonian_numerator(first)
    # <0|(128H)^2|0> = ||128H|0>||^2, including creation of one/two
    # oscillator quanta before returning. Projecting each power to N=0
    # would silently delete the last three contributions below.
    expected = 146**2 + 2 * 128**2 + 4 * 9**2 + 18**2 + 32**2
    assert expected == 55756
    assert second[0, 0, 0] == (expected, 0)
    assert history._norm_numerator(first) == expected
    projected_first = {key: value for key, value in first.items() if sum(key[:2]) == 0}
    incorrect = history._hamiltonian_numerator(projected_first)[0, 0, 0][0]
    assert expected - incorrect == 1672
    polynomial, remainder = history._taylor(vacuum, F(1, 2), 1)
    assert _rational(polynomial)[0, 0, 0] == (F(1), -F(73, 128))
    assert _rational(remainder)[0, 0, 0] == (-F(55756, 131072), F(0))
    assert max(sum(key[:2]) for key in polynomial.coefficients) == 2
    assert max(sum(key[:2]) for key in remainder.coefficients) == 4


def test_support_discard_and_rounding_charge_the_actual_factorial_weighted_losses():
    polynomial = history._Polynomial({(0, 0, 0): (1, -2), (2, 0, 1): (3, 4),
                                     (0, 3, 2): (-5, 6)}, 7)
    normalization = F(7, 3)
    rounded, discard, rounding = history._cut_and_round(polynomial, 2, 1, normalization)
    original, output = _rational(polynomial), _rational(rounded)
    omitted = {key: value for key, value in original.items() if sum(key[:2]) > 2}
    difference = {key: tuple(original[key][part] - output.get(key, (F(0), F(0)))[part]
                             for part in (0, 1))
                  for key in original if sum(key[:2]) <= 2}
    exact_discard_squared = halfline.norm_squared(omitted) / normalization
    exact_rounding_squared = halfline.norm_squared(difference) / normalization
    for bound, squared in ((discard, exact_discard_squared), (rounding, exact_rounding_squared)):
        assert bound > 0 and bound * bound >= squared
        assert max(F(0), bound - F(1, history.ROOT_SCALE))**2 <= squared
    assert set(output) <= {key for key in original if sum(key[:2]) <= 2}
    assert _norm_squared(rounded) != normalization


def test_first_record_and_feedback_preserve_total_mass_without_branch_normalization():
    initial = history._Polynomial({(0, 0, 0): (2, 1), (1, 0, 1): (1, -3),
                                   (0, 2, 2): (2, -1), (1, 1, 3): (-4, 2)}, 13)
    for enabled in (False, True):
        branches = [history._first_record(initial, bit, enabled) for bit in (0, 1)]
        assert sum(map(_norm_squared, branches), F(0)) == _norm_squared(initial)
        for bit, branch in enumerate(branches):
            original_mass = sum((math.factorial(n0) * math.factorial(n1) *
                                 F(real * real + imag * imag, 13**2)
                                 for (n0, n1, gauge), (real, imag) in initial.coefficients.items()
                                 if gauge & 1 == bit), F(0))
            assert _norm_squared(branch) == original_mass
    branch = history._first_record(initial, 1, True)
    assert branch.denominator == 65
    assert branch.coefficients[1, 0, 1] == (3, -9)
    assert branch.coefficients[1, 0, 0] == (-12, -4)


def test_inadequate_history_does_not_certify_and_keeps_excess_raw_mass():
    result = history.certificate(cutoff=2, preparation_cutoff=0, degree=1,
                                 steps_per_unit=2, bits=32)
    assert result["nonzero_response_certified"] is False
    lower, upper = map(F, result["infinite_feedback_contrast_interval"])
    assert lower <= 0 <= upper
    assert result["maximum_expanded_total_occupation"] > result["parameters"]["retained_total_occupation"]
    preparation = F(result["initial_trace_distance_upper"])
    for row in result["configurations"]:
        error = F(row["normalized_preparation_to_polynomial_vector_error_upper"])
        # Unnormalized phi=(1+epsilon)psi already has event-weight error
        # 2epsilon+epsilon^2; the half trace norm alone cannot bound it.
        assert F(row["event_probability_error_upper"]) >= preparation + 2 * error + error**2
        assert F(row["complete_output_half_trace_norm_error_upper"]) >= preparation + error + error**2 / 2
        assert tuple(map(F, row["infinite_positive_probability_interval"])) == (F(0), F(1))
        mass = sum((F(branch["approximate_pre_readout_norm_squared"])
                    for branch in row["retained_first_records"]), F(0))
        assert mass > 1  # The weak Taylor polynomial is never normalized away.
        assert [branch["first_record"] for branch in row["retained_first_records"]] == [0, 1]


def test_disconnected_control_preserves_scalar_parity_and_cannot_certify_a_signal():
    result = history.certificate(cutoff=4, preparation_cutoff=2, degree=6,
                                 steps_per_unit=4, bits=64, coupling=F(0),
                                 retain_polynomials=True)
    assert result["nonzero_response_certified"] is False
    normalization = F(result["initial_polynomial_norm_squared"])
    for row in result["configurations"]:
        lower, upper = map(F, row["infinite_positive_probability_interval"])
        assert lower <= F(1, 2) <= upper
        for branch in row["retained_first_records"]:
            vector = {(n0, n1, gauge): (F(real), F(imag))
                      for n0, n1, gauge, real, imag in branch["pre_readout_polynomial"]}
            assert all((n0 + n1) % 2 == 0 for n0, n1, _ in vector)
            constant, cross = halfline.positive_weight_coefficients(vector, normalization)
            assert cross == 0
            assert 2 * constant == F(branch["approximate_pre_readout_norm_squared"])


@pytest.mark.parametrize("name,value", (
    ("cutoff", True), ("preparation_cutoff", True), ("degree", True),
    ("steps_per_unit", True), ("bits", True), ("cutoff", 1), ("cutoff", 65),
    ("preparation_cutoff", 3), ("preparation_cutoff", 18), ("degree", 0),
    ("steps_per_unit", 3), ("bits", 31), ("bits", 257),
    ("coupling", True), ("coupling", .25), ("coupling", F(1, 8)),
    ("retain_polynomials", 1),
))
def test_invalid_history_parameters_fail_before_computation(name, value):
    with pytest.raises(ValueError):
        history.certificate(**{name: value})
