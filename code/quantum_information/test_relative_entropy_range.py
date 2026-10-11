"""Original-entry controls for small quantum relative entropy contractions."""

from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from quantum_information.states import relative_entropy


def _qubit_oracle(source, reference):
    """Independent two-root formula; no producer spectrum or remainder."""
    with mp.workdps(450):
        a, d = (mp.mpf(float(source[i, i].real)) for i in range(2))
        z = source[0, 1]
        radius = mp.mpf(float(z.real))**2 + mp.mpf(float(z.imag))**2
        large = (a+d+mp.sqrt((a-d)**2+4*radius))/2
        small = (a*d-radius)/large
        assert small > 0
        entropy = small*mp.log(small)+large*mp.log(large)
        pairing = mp.fsum(mp.mpf(float(source[i, i].real))
                         * mp.log(mp.mpf(float(reference[i, i].real)))
                         for i in range(2))
        return float(entropy-pairing)


@pytest.mark.parametrize("phase", (1., 1j, 1+1j))
@pytest.mark.parametrize("exponent", (-537, -538, -539))
def test_normalized_states_retain_log_amplified_subnormal_coherence(phase, exponent):
    p = 2.**-40
    reference = np.diag([p, 1-p]).astype(complex)
    source = reference.copy()
    source[0, 1] = phase*2.**exponent
    source[1, 0] = source[0, 1].conjugate()
    assert sum(Fraction(float(x.real)) for x in source.diagonal()) == 1
    expected = _qubit_oracle(source, reference)
    assert expected > 0
    assert relative_entropy(source, reference) == expected


@pytest.mark.parametrize("phase", (1., 1j))
@pytest.mark.parametrize("exponent", (-538, -539))
def test_fallback_uses_original_subnormal_population_not_rounded_eigenvalue(phase, exponent):
    p = np.nextafter(0., 1.)
    source = np.array([[p, phase*2.**exponent],
                       [np.conj(phase)*2.**exponent, 1.]])
    reference = np.diag([p, 1.])
    expected = _qubit_oracle(source, reference)
    assert expected > 0
    assert relative_entropy(source, reference) == expected


@pytest.mark.parametrize("phase", (1., 1j))
def test_independent_blocks_add_before_the_final_subnormal_rounding(phase):
    p, t = 2.**-40, 2.**-539
    block = np.array([[p, phase*t], [np.conj(phase)*t, .5-p]])
    reference_block = np.diag([p, .5-p])
    source = np.kron(np.eye(2), block)
    reference = np.kron(np.eye(2), reference_block)
    # Scale the exact block before applying the oracle. D(cA||cB)=cD(A||B),
    # so this compares the complete real sum, without rounding each block.
    expected = _qubit_oracle(2*block, 2*reference_block)
    assert expected > 0
    assert relative_entropy(source, reference) == expected


def test_exact_hermitian_averaging_precedes_zero_comparison():
    tiny = np.nextafter(0., 1.)
    reference = np.diag([1., *([tiny]*6)])
    source = reference.copy()
    for i in (1, 3, 5):
        source[i, i+1] = tiny
    # Each exact Hermitian block has eigenvalues 3t/2,t/2, versus t,t.
    with mp.workdps(450):
        expected = float(3*mp.mpf(tiny)
                         *(mp.mpf('1.5')*mp.log(mp.mpf('1.5'))
                           + mp.mpf('.5')*mp.log(mp.mpf('.5'))))
    assert expected == tiny
    assert relative_entropy(source, reference) == expected
    antisymmetric = reference.copy()
    antisymmetric[1, 2], antisymmetric[2, 1] = tiny, -tiny
    assert relative_entropy(antisymmetric, reference) == 0


def test_complex_singular_reference_uses_exact_support_containment():
    p, t = 2.**-40, 2.**-538
    q = (1-p)/2
    reference = np.array([[p, 0, 0], [0, q, -1j*q], [0, 1j*q, q]])
    source = reference.copy()
    source[0, 1] = source[1, 0] = t
    source[0, 2], source[2, 0] = -1j*t, 1j*t
    # The last two coordinates have support vector (1,i)/sqrt(2).
    # Its coupling has squared modulus2t², with one exact kernel spectator.
    reduced = np.array([[p, (1+1j)*t], [(1-1j)*t, 1-p]])
    expected = _qubit_oracle(reduced, np.diag([p, 1-p]))
    with mp.workdps(17):
        assert relative_entropy(source, reference) == expected
        assert mp.mp.dps == 17
    escaped = reference.copy()
    escaped[1, 2], escaped[2, 1] = -1j*(q-p), 1j*(q-p)
    assert relative_entropy(escaped, reference) == float("inf")


@pytest.mark.parametrize("phase", (1., 1j))
def test_nonzero_divergence_below_output_range_is_refused(phase):
    reference = np.diag([.25, .75]).astype(complex)
    source = reference.copy()
    source[0, 1] = phase*2.**-600
    source[1, 0] = source[0, 1].conjugate()
    # Strict Klein inequality proves D>0 for these unequal trace-one states;
    # the original-entry formula confirms that it rounds below binary64.
    assert _qubit_oracle(source, reference) == 0
    with pytest.raises(ValueError, match="relative entropy.*unresolved"):
        relative_entropy(source, reference)


def test_normal_result_exact_zero_and_negative_trace_offset_controls():
    reference = np.diag([2.**-40, 1-2.**-40])
    source = reference.astype(complex)
    source[0, 1], source[1, 0] = 2.**-30*(1+1j), 2.**-30*(1-1j)
    assert relative_entropy(source, reference) == pytest.approx(
        _qubit_oracle(source, reference), rel=3e-14, abs=0)
    assert relative_entropy(source, source) == 0
    tiny = np.nextafter(0., 1.)
    a, b = np.diag([1., 2*tiny]), np.diag([1., 4*tiny])
    with mp.workdps(450):
        expected = float(2*mp.mpf(tiny)*mp.log(mp.mpf('.5')))
    assert expected < 0
    assert relative_entropy(a, b) == expected


@pytest.mark.parametrize("bad", (
    [[True, 0.], [0., 0.]],
    np.ma.array(np.eye(2)/2, mask=[[True, False], [False, False]]),
))
@pytest.mark.parametrize("source", (True, False))
def test_source_level_evaluation_rejects_hidden_boolean_or_missing_evidence(bad, source):
    pair = (bad, np.eye(2)/2) if source else (np.eye(2)/2, bad)
    with pytest.raises(ValueError):
        relative_entropy(*pair)


@pytest.mark.skipif(np.finfo(np.longdouble).nmant <= 52,
                    reason="longdouble carries no extra bits on this platform")
@pytest.mark.parametrize("source", (True, False))
def test_extended_positive_input_cannot_be_erased_before_original_value_fallback(source):
    tiny = np.longdouble(np.nextafter(0., 1.))/2
    state = np.diag(np.array([1., tiny], dtype=np.longdouble))
    pair = (state, np.diag([1., 0.])) if source else (np.diag([1., 0.]), state)
    with pytest.raises(ValueError, match="loses information"):
        relative_entropy(*pair)
    ordinary = np.diag(np.array([.25, .75], dtype=np.longdouble))
    assert relative_entropy(ordinary, ordinary) == 0
