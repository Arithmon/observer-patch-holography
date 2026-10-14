"""Complete projector histories and independently exponentiated finite checks."""
from fractions import Fraction as F
from math import factorial

import numpy as np
from scipy.linalg import expm

from model.relativistic_scalar import exact, exact_history as history


def test_projectors_are_exact_complete_and_orthogonal():
    state = history.Polynomial({exact.VACUUM: (3, -2),
                                (4,)+(0,)*7: (1, 7),
                                (2,)+(0,)*7: (2, 1)}, 13)
    yes, no = history.project(state, 1), history.project(state, 0)
    assert history.norm_squared(history.difference(history.project(yes, 1), yes)) == 0
    assert history.norm_squared(history.project(no, 1)) == 0
    assert history.norm_squared(yes)+history.norm_squared(no) == history.norm_squared(state)
    assert history.norm_squared(history.apply_feedback(yes, 1)) == history.norm_squared(yes)
    assert history.norm_squared(history.project(history.vacuum(), 1)) == F(1, 25)


def _independent_history(time, coupling, feedback, instrument_sign):
    # Five levels at the origin and only the vacuum elsewhere.  Construct
    # P H P from a larger oscillator, including q^4 paths outside P.
    d, large = 5, 10
    a = np.diag(np.sqrt(np.arange(1, large)), 1)
    b4 = np.linalg.matrix_power(a+a.T, 4)[:d, :d]
    h = np.diag(2*np.arange(d, dtype=float)) + float(coupling)*(b4-3*np.eye(d))/384
    u = expm(-1j*float(time)*h)
    e = np.zeros(d, dtype=complex)
    e[0], e[4] = 1/5, instrument_sign*1j*np.sqrt(factorial(4))/5
    p = np.outer(e, e.conj())
    phase = np.eye(d, dtype=complex)
    phase[4, 4] = 1j
    vacuum = np.eye(d)[0]
    out = {}
    for bit in (0, 1):
        state = (p if bit else np.eye(d)-p) @ u @ vacuum
        if bit and feedback:
            state = phase @ state
        state = u @ state
        for final in (0, 1):
            branch = (p if final else np.eye(d)-p) @ state
            out[bit, final] = float(np.vdot(branch, branch).real)
    return out


def test_all_finite_probabilities_enclose_independent_exact_evolution():
    time, cutoff = F(1, 128), (5,)+(1,)*7
    for coupling, feedback, sign in ((F(1), True, -1), (F(1), False, -1),
                                    (F(0), True, -1), (F(1), True, 1)):
        result = history.certificate(time, local_cutoff=cutoff,
                                     coupling=coupling, feedback=feedback, instrument_sign=sign)
        independent = _independent_history(time, coupling, feedback, sign)
        assert len(result["branches"]) == 4
        for row in result["branches"]:
            lo, hi = map(float, row["probability_interval"])
            assert lo <= independent[tuple(row["records"])] <= hi
        lo, hi = result["normalization_interval"]
        assert lo <= 1 <= hi
        assert result["raw_polynomial_mass"] != 1


def test_zero_time_and_infinite_oscillator_history():
    zero = history.certificate(F(0))
    assert zero["complete_retained_vector_error_upper"] == 0
    assert zero["raw_polynomial_mass"] == 1
    result = history.certificate(F(1, 256))
    assert result["scope"]["oscillator_space"] == "infinite"
    assert result["normalization_interval"][0] <= 1 <= result["normalization_interval"][1]
    assert 0 < result["complete_retained_vector_error_upper"] < F(1, 1000)


def test_complete_cutoff_comparison_is_nonvacuous():
    result = history.comparison_certificate((8,)+(2,)*7, time=F(1, 256))
    assert 0 < result["complete_history_tv_upper"] < F(1, 100)
    assert result["polynomial_vector_difference_norm_squared"] > 0
