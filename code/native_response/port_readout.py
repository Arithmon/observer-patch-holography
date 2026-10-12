"""Occupied-port records on the binary two-load shell.

The twelve effects E_p=N_p/2 are a POVM, not central projections on the
native matrix algebra. The instrument's output is the genuine direct sum
of twelve matrix algebras B(N_p H), each of dimension 11*d as a Hilbert
space when H=C^66 tensor C^d. Its sectors record a randomly selected
occupied port. These finite operations do not establish A1 compliance.

Forgetting the record multiplies |x><y| tensor A by |x intersect y|/2.
It preserves the tracial state and every load-diagonal private operator,
but disturbs native load coherences. Matrix routines are floating-point
evaluators; their input validation is not an interval certificate.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from numbers import Integral, Rational, Real

import numpy as np

from quantum_information.algebras import operator
from quantum_information.states import density_matrix


BINARY_LOADS = tuple(combinations(range(12), 2))


def _matrix(value):
    if np.ma.isMaskedArray(value):
        raise ValueError("occupied-port operators must not be masked arrays")
    return operator(value)


def _positive_integer(value, name):
    if (isinstance(value, (bool, np.bool_)) or
            not isinstance(value, Integral) or value <= 0):
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def sector_indices(port, private_dimension=1):
    """Indices of N_p H in lexicographic load, then private, ordering."""
    if (isinstance(port, (bool, np.bool_)) or not isinstance(port, Integral)
            or not 0 <= port < 12):
        raise ValueError("port must be an integer from 0 to 11")
    d = _positive_integer(private_dimension, "private dimension")
    return np.array([d * i + j for i, load in enumerate(BINARY_LOADS)
                     if port in load for j in range(d)], dtype=int)


def occupied_port_effects(private_dimension=1):
    """Return E_p=N_p/2; every state has exactly two occupied ports."""
    d = _positive_integer(private_dimension, "private dimension")
    result = []
    for port in range(12):
        diagonal = np.zeros(66 * d)
        diagonal[sector_indices(port, d)] = 0.5
        result.append(np.diag(diagonal))
    return tuple(result)


def occupied_port_kraus(private_dimension=1):
    """Rectangular Kraus maps H -> N_p H, with the outcome kept classical."""
    d = _positive_integer(private_dimension, "private dimension")
    identity = np.eye(66 * d)
    return tuple(identity[sector_indices(port, d)] / np.sqrt(2)
                 for port in range(12))


def luders_branches(value):
    """Linear instrument on operators; positive inputs give positive blocks.

    For a normalized density matrix, use measure_occupied_port to validate
    positivity and normalization as well as the native-space dimensions.
    """
    value = _matrix(value)
    if len(value) % 66:
        raise ValueError("native operator dimension must be 66*d")
    d = len(value) // 66
    return tuple(value[np.ix_(indices, indices)] / 2
                 for indices in (sector_indices(port, d) for port in range(12)))


def measure_occupied_port(state):
    """Validate a state and return twelve subnormalized record blocks.

    No branch is postselected or normalized separately. The blocks retain
    the full private factor and eleven compatible load configurations.
    """
    value = _matrix(state)
    if len(value) % 66:
        raise ValueError("native state dimension must be 66*d")
    return luders_branches(density_matrix(value))


def _port_blocks(blocks):
    if not isinstance(blocks, (list, tuple)) or len(blocks) != 12:
        raise ValueError("twelve occupied-port operator blocks required")
    values = tuple(_matrix(block) for block in blocks)
    n = len(values[0])
    if n % 11 or any(block.shape != (n, n) for block in values):
        raise ValueError("occupied-port blocks must share dimension 11*d")
    return values, n // 11


def expand_port_record(blocks):
    """Embed each compressed branch into its labelled full native block.

    The result has the address carrier's shape (12, 66*d, 66*d). Composing
    expand_port_record(measure_occupied_port(decode(sigma))) implements
    J^A(sigma)_p=N_p decode(sigma) N_p/2 in that same carrier. It preserves
    the new occupied-port label, while the preceding decode forgets the old
    address. It does not store an endpoint history or reverse measurement.
    """
    values, d = _port_blocks(blocks)
    result = np.zeros((12, 66 * d, 66 * d), complex)
    for port, block in enumerate(values):
        indices = sector_indices(port, d)
        result[port][np.ix_(indices, indices)] = block
    return result


def forget_port_record(blocks):
    """Embed and sum twelve operator blocks, without undoing their readout."""
    values, d = _port_blocks(blocks)
    result = np.zeros((66 * d, 66 * d), complex)
    for port, block in enumerate(values):
        indices = sector_indices(port, d)
        result[np.ix_(indices, indices)] += block
    return result


def unread_port_readout(value):
    """The forgotten-record channel, not the identity on native states."""
    return forget_port_record(luders_branches(value))


def _joint_probability(probability, upper):
    if (np.ma.isMaskedArray(probability) or
            isinstance(probability, (bool, np.bool_))):
        raise ValueError("joint probability must be a finite real scalar")
    if isinstance(probability, Rational):
        q = Fraction(probability)
    elif isinstance(probability, Real) and np.isfinite(probability):
        q = float(probability)
        if q != probability:
            raise ValueError("joint probability loses precision in binary64 conversion")
        upper = float(upper)
    else:
        raise ValueError("joint probability must be a finite real scalar")
    if not 0 <= q <= upper:
        raise ValueError(f"joint probability outside the tracial endpoint range [0,{upper}]")
    return q


def correlation_from_same_port(probability):
    """Decode c_n=(144*q_pp(n)-1)/5 from fresh endpoint port records.

    q_pp is a JOINT probability, not a conditional repeat probability.
    The identity assumes initial normalized trace on the binary two-load
    shell, evolution preserving that trace, and two fresh executions of
    this instrument. Importing the N4
    moments requires the unmodified declared evolution between them, with
    no intermediate readout. This scalar conversion cannot verify those
    experimental premises. Rational inputs retain exact arithmetic.
    """
    q = _joint_probability(probability, Fraction(1, 24))
    return (144 * q - 1) / 5


def correlation_from_equal_ports(probability):
    """Decode the mean port correlation from P(tag_initial=tag_final).

    Under the same endpoint protocol, c_mean=(12*P_equal-1)/5. Proper
    carrier covariance makes every local c_n equal to this mean; without
    that symmetry, substituting this mean into the local memory statistic
    has a different null. Record storage for both endpoints is supplied.
    """
    q = _joint_probability(probability, Fraction(1, 2))
    return (12 * q - 1) / 5


def corrected_memory_from_same_port(one_attempt, two_attempts):
    """Convert endpoint probabilities to the corrected local memory test.

    Uniform IID native/N4 seam attempts, their tracial stationarity and the
    no-intermediate-readout protocol are required. A positive floating or
    empirical result is not a positivity or statistical certificate.
    """
    c1 = correlation_from_same_port(one_attempt)
    c2 = correlation_from_same_port(two_attempts)
    return c2 - c1**2 - (1 - c1)**2 / 10


def corrected_memory_from_equal_ports(one_attempt, two_attempts):
    """The same corrected statistic using all equal-tag outcomes.

    This requires the preceding endpoint premises AND proper-carrier
    covariance, as supplied by the declared uniform native/N4 process.
    """
    c1 = correlation_from_equal_ports(one_attempt)
    c2 = correlation_from_equal_ports(two_attempts)
    return c2 - c1**2 - (1 - c1)**2 / 10
