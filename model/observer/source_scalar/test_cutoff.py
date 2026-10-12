"""Independent boundary, Gaussian-history and rational cutoff checks."""
from decimal import Decimal, localcontext
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import expm

from . import cutoff, field


@pytest.mark.parametrize("value", (F(0), F(1, 10), F(1), F(4), F(10), F(40)))
def test_directed_exponential_against_independent_high_precision(value):
    bound = cutoff.exponential_upper(value)
    with localcontext() as context:
        context.prec = 100
        expected = (Decimal(value.numerator)/Decimal(value.denominator)).exp()
        actual = Decimal(bound.numerator)/Decimal(bound.denominator)
        assert 0 <= actual-expected <= Decimal("2e-30")


@pytest.mark.parametrize("q", (F(5, 2), F(4), F(6)))
@pytest.mark.parametrize("components", ((F(0), F(1, 2), F(1, 4)),
                                        (F(1, 2), F(0), F(1, 4)),
                                        (F(1, 4), F(1, 2), F(0))))
def test_number_moment_dominates_independent_bargmann_series(q, components):
    # Expand the displaced squeezed normal mode by its creation recurrence;
    # this does not call the closed generating-function expression.
    x, y, common = components
    alpha = complex(float(x), float(y))
    tau = 1/9
    linear = alpha+tau*alpha.conjugate()
    amplitude = (1-tau*tau)**.25*np.exp(-abs(alpha)**2/2-tau*alpha.conjugate()**2/2)
    previous = 0j
    moment = abs(amplitude)**2
    for occupied in range(120):
        following = (linear*amplitude-tau*math.sqrt(occupied)*previous)/math.sqrt(occupied+1)
        previous, amplitude = amplitude, following
        moment += abs(amplitude)**2*float(q)**(occupied+1)
    moment *= math.exp((float(q)-1)*float(common*common))
    radius = x*x+y*y+common*common
    assert moment <= float(cutoff.number_moment_bound(q, radius))*(1+2e-14)


@pytest.mark.parametrize("d", (2, 4, 8))
def test_cross_boundary_generator_constants_on_an_independent_larger_space(d):
    larger = d+2
    a = np.zeros((larger, larger))
    for n in range(1, larger):
        a[n-1, n] = math.sqrt(n)
    q = (a+a.T)/math.sqrt(2)
    p = (a-a.T)/(1j*math.sqrt(2))
    q0, q1 = np.kron(np.eye(larger), q), np.kron(q, np.eye(larger))
    p0, p1 = np.kron(np.eye(larger), p), np.kron(p, np.eye(larger))
    h = (p0@p0+p1@p1+q0@q0+q1@q1)/2+(9/64)*(q0-q1)@(q0-q1)
    inside = [n0+larger*n1 for n1 in range(d) for n0 in range(d)]
    outside = [index for index in range(larger*larger) if index not in inside]
    numbers = np.array([index % larger+index//larger for index in outside])
    h_residual = h[np.ix_(inside, outside)]/numbers
    q_residual = q0[np.ix_(inside, outside)]/np.sqrt(numbers)
    pointer_residual = q1[np.ix_(inside, outside)]/(2*np.sqrt(numbers))
    assert np.linalg.norm(h_residual, 2) <= 9/32+2e-15
    assert np.linalg.norm(q_residual, 2) <= 1/math.sqrt(2)+2e-15
    assert np.linalg.norm(pointer_residual, 2) <= 1/(2*math.sqrt(2))+2e-15


def _gaussian_coefficients(d, displacement):
    """Site-basis Bargmann recurrence for the exact displaced coupled vacuum."""
    b = np.array([[-1., 1.], [1., -1.]])/18
    displacement = np.asarray(displacement, dtype=complex)
    linear = displacement-b@displacement.conj()
    amplitudes = np.zeros((d, d), dtype=complex)
    amplitudes[0, 0] = (80/81)**.25*np.exp(
        -np.vdot(displacement, displacement)/2
        +displacement.conj()@b@displacement.conj()/2)
    for n1 in range(1, d):
        amplitudes[n1, 0] = (linear[1]*amplitudes[n1-1, 0]
                            +(b[1, 1]*math.sqrt(n1-1)*amplitudes[n1-2, 0] if n1 > 1 else 0))/math.sqrt(n1)
    for n1 in range(d):
        for n0 in range(1, d):
            amplitudes[n1, n0] = (linear[0]*amplitudes[n1, n0-1]
                +(b[0, 0]*math.sqrt(n0-1)*amplitudes[n1, n0-2] if n0 > 1 else 0)
                +(b[0, 1]*math.sqrt(n1)*amplitudes[n1-1, n0-1] if n1 > 0 else 0))/math.sqrt(n0)
    return amplitudes.reshape(-1)


@pytest.mark.parametrize("d", (4, 8, 16))
def test_normalized_preparation_versus_projected_infinite_vacuum(d):
    target, _ = field.vacuum_vector(d)
    projected_reference = _gaussian_coefficients(d, (0., 0.))
    certificate = cutoff.occupation_cutoff_bound(d)
    assert np.linalg.norm(target-projected_reference) <= float(F(certificate["preparation_vector_error_upper"]))+1e-15


@pytest.mark.parametrize("d", (4, 8))
def test_complete_quantum_instrument_against_independent_gaussian_history(d):
    # An independent analytic Gaussian Fock expansion checks each complete
    # terminal branch. The numerical comparison is a regression, not the proof.
    q, _ = field.quadratures(d)
    initial, _ = field.vacuum_vector(d)
    finite = expm(-1j*field.physical_hamiltonian(d))@np.kron(np.eye(d), expm(1j*q))@initial
    finite_arms = [np.kron(expm(-.5j*sign*q), np.eye(d))@finite for sign in (1, -1)]

    position = np.array([(math.sin(1)+.8*math.sin(1.25))/2,
                         (math.sin(1)-.8*math.sin(1.25))/2])
    momentum = np.array([(math.cos(1)+math.cos(1.25))/2,
                         (math.cos(1)-math.cos(1.25))/2])
    alpha = (position+1j*momentum)/math.sqrt(2)
    larger = 32
    exact_arms = []
    for sign in (1, -1):
        beta = np.array([0., -sign*1j/(2*math.sqrt(2))])
        weyl_phase = np.exp(1j*np.imag(np.sum(beta*alpha.conj())))
        exact_arms.append(np.exp(-9j/8)*weyl_phase*_gaussian_coefficients(larger, alpha+beta))
    finite_embedded = []
    for arm in finite_arms:
        embedded = np.zeros((larger, larger), dtype=complex)
        embedded[:d, :d] = arm.reshape(d, d)
        finite_embedded.append(embedded.reshape(-1))
    vector_error = math.sqrt(sum(np.linalg.norm(a-b)**2 for a, b in zip(exact_arms, finite_embedded))/2)
    certificate = cutoff.occupation_cutoff_bound(d)
    assert vector_error <= float(F(certificate["full_stinespring_vector_error_upper"]))

    instrument_distance = 0.
    for outcome in (1, -1):
        exact_branch = (exact_arms[0]-1j*outcome*exact_arms[1])/2
        finite_branch = (finite_embedded[0]-1j*outcome*finite_embedded[1])/2
        p, r = np.vdot(exact_branch, exact_branch).real, np.vdot(finite_branch, finite_branch).real
        overlap = abs(np.vdot(exact_branch, finite_branch))**2
        # Trace norm of the difference of two *unnormalized* rank-one states.
        instrument_distance += .5*math.sqrt(max(0., (p+r)**2-4*overlap))
        if outcome == 1:
            mean = (math.sin(1)-.8*math.sin(1.25))/2
            expected_probability = (1+math.exp(-9/40)*math.sin(mean))/2
            assert abs(p-expected_probability) < 3e-15
    assert instrument_distance <= float(F(certificate["instrument_trace_distance_upper"]))


def test_nontrivial_certificates_and_converging_occupation_refinement():
    bounds = {d: F(cutoff.occupation_cutoff_bound(d)["instrument_trace_distance_upper"])
              for d in (4, 8, 16, 32, 64)}
    assert bounds[4] < F(4, 5)
    assert bounds[8] < F(78, 1000)
    assert bounds[16] < F(196, 10**6)
    assert bounds[32] < F(21, 10**11)
    assert bounds[64] < F(1, 10**22)
    assert all(bounds[a] > bounds[b] for a, b in ((4, 8), (8, 16), (16, 32), (32, 64)))


def test_parity_uniform_history_and_no_conditional_probability_division():
    positive = cutoff.occupation_cutoff_bound(8)
    negative = cutoff.occupation_cutoff_bound(8, current=-1)
    zero = cutoff.occupation_cutoff_bound(8, current=0)
    for key in ("instrument_trace_distance_upper", "full_stinespring_vector_error_upper",
                "final_outside_probability_upper"):
        assert positive[key] == negative[key]
        assert F(zero[key]) < F(positive[key])
    assert zero["impulse_projected_residual_integral_upper"] == "0"
    at_zero = cutoff.occupation_cutoff_bound(8, time=F(0))
    assert at_zero["evolution_projected_residual_integral_upper"] == "0"
    assert F(at_zero["instrument_trace_distance_upper"]) < F(positive["instrument_trace_distance_upper"])


def test_reference_probability_interval_against_independent_high_precision():
    def sine(value):
        term, result, power = value, value, 1
        while abs(term) > Decimal("1e-95"):
            term *= -value*value/Decimal((power+1)*(power+2))
            result += term
            power += 2
        return result

    def decimal(value):
        fraction = F(value)
        return Decimal(fraction.numerator)/Decimal(fraction.denominator)

    interval = cutoff.reference_probability_interval()
    with localcontext() as context:
        context.prec = 100
        mean = (sine(Decimal(1))-Decimal(4)*sine(Decimal(5)/4)/5)/2
        attenuation = (-Decimal(9)/40).exp()
        expected = (1+attenuation*sine(mean))/2
        lower, upper = decimal(interval["plus_probability_lower"]), decimal(interval["plus_probability_upper"])
        assert lower <= expected <= upper
        assert upper-lower < Decimal("1e-29")
        assert decimal(interval["mean_lower"]) <= mean <= decimal(interval["mean_upper"])
        assert decimal(interval["attenuation_lower"]) <= attenuation <= decimal(interval["attenuation_upper"])
        assert lower > Decimal("0.5164215017481")


def test_reference_interval_parity_and_certified_implemented_signal_margin():
    from .bounds import split_error_bound

    positive = cutoff.reference_probability_interval(1)
    negative = cutoff.reference_probability_interval(-1)
    zero = cutoff.reference_probability_interval(0)
    assert F(negative["plus_probability_lower"]) == 1-F(positive["plus_probability_upper"])
    assert F(negative["plus_probability_upper"]) == 1-F(positive["plus_probability_lower"])
    assert zero["plus_probability_lower"] == zero["plus_probability_upper"] == "1/2"
    error = F(cutoff.occupation_cutoff_bound(16)["instrument_trace_distance_upper"])
    error += F(split_error_bound(16, 128)["detector_probability_error_upper"])
    assert F(positive["plus_probability_lower"])-error > F(5149, 10000)
    assert F(negative["plus_probability_upper"])+error < F(4851, 10000)


def test_executed_receipt_composes_cutoff_without_claiming_unexecuted_programs():
    from . import experiment

    # Run the actual source executor, including the finest implemented cutoff;
    # its certified interval must distinguish the fair-detector alternative.
    result = experiment.run(trials=1000)
    assert result["all_pass"]
    reference = result["infinite_two_oscillator_reference"]["plus_probability"]
    for level in result["levels"]:
        occupation = F(level["continuous_instrument_occupation_certificate"]["instrument_trace_distance_upper"])
        for row in level["step_refinement"]:
            finite = F(row["finite_split_certificate"]["detector_probability_error_upper"])
            combined = row["ideal_infinite_instrument_certificate"]
            error = F(combined["instrument_trace_distance_upper"])
            assert error == min(F(1), occupation+finite)
            # Preserve the existing finite-only fields and distinguish the
            # complete instrument from its parity-symmetric classical records.
            assert F(row["ideal_classical_history_tv_error_upper"]) == F(5, 33)*finite
            assert F(combined["classical_history_tv_error_upper_under_parity"]) == F(5, 33)*error
            assert F(combined["native_current_error_upper_under_parity"]) == F(5, 132)*error
            lower = F(combined["ideal_finite_plus_probability_lower"])
            upper = F(combined["ideal_finite_plus_probability_upper"])
            measured = row["histories"][-1]["branches"][0]["probability"]
            assert abs(measured-reference) < float(error)
            assert float(lower) < measured < float(upper)
            assert combined["ideal_nonzero_response_certified"] == (lower > F(1, 2))
    refinements = result["occupation_refinement"]
    assert [row["dimension"] for row in refinements] == [4, 8, 16, 32, 64]
    assert [row["executed"] for row in refinements] == [True, True, True, False, False]
    assert all("no source program emitted or executed" in row["execution_scope"] for row in refinements[3:])
    assert "certified numerical execution error" in result["scope"]["open"]


@pytest.mark.parametrize("kwargs", ({"d": True}, {"d": 3}, {"d": 0}, {"d": 4, "time": 1.},
    {"d": 4, "time": F(2)}, {"d": 4, "time": F(-1)}, {"d": 4, "current": True}, {"d": 4, "current": 2}))
def test_invalid_scope_is_rejected(kwargs):
    with pytest.raises(ValueError):
        cutoff.occupation_cutoff_bound(**kwargs)
