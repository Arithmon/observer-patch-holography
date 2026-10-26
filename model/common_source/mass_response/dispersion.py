"""Independent analytic controls for the source-defined dimer continuation.

The four-mode interacting ring and the free three-axis extension are different
contracts.  The former has an exact finite charge-transfer gap.  The latter
has an all-volume structural dispersion gap and a controlled continuum limit.
Neither formula identifies a Standard Model particle or the pole of the
homogeneous interacting infinite-volume theory.

All returned floats are numerical evaluations of the displayed exact formulas,
not interval-arithmetic or physical-error certificates.  No measured mass is
an input.  Mapping the native reference log odds to hopping coefficients is
an explicitly chosen coherent architecture.
"""
from __future__ import annotations

import math
from numbers import Real

import numpy as np


SOURCE_A = math.log(25 / 8)
SOURCE_B = math.log(8 / 5)
SOURCE_J = SOURCE_A - SOURCE_B
SOURCE_OMEGA = SOURCE_A + SOURCE_B


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        value = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must have a finite execution representation") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite real number")
    return value


def _level(value):
    if type(value) is not int or value < 1:
        raise ValueError("a positive integer refinement level is required")
    return value


def _momentum(value):
    try:
        values = tuple(value)
    except TypeError as exc:
        raise ValueError("three finite real momentum components are required") from exc
    if len(values) != 3:
        raise ValueError("three finite real momentum components are required")
    return tuple(_real(x, "momentum") for x in values)


def _hoppings(intra, inter):
    intra, inter = _real(intra, "intracell hopping"), _real(inter, "intercell hopping")
    if intra <= 0 or inter <= 0:
        raise ValueError("positive hopping coefficients are required")
    return intra, inter


def gamma_matrices():
    """Six Hermitian Clifford generators on eight retained components.

    These are a declared regulator representation, not a proof of physical
    spin, chirality, statistics or a minimal particle-species count.
    """
    identity = np.eye(2, dtype=complex)
    x = np.array(((0, 1), (1, 0)), dtype=complex)
    y = np.array(((0, -1j), (1j, 0)), dtype=complex)
    z = np.diag((1, -1)).astype(complex)
    return tuple(
        np.kron(np.kron(left, middle), right)
        for left, middle, right in (
            (x, identity, identity), (y, identity, identity),
            (z, x, identity), (z, y, identity),
            (z, z, x), (z, z, y),
        )
    )


def bloch_hamiltonian(momentum, *, intra=SOURCE_A, inter=SOURCE_B):
    """Free three-axis dimer operator; momentum is dimensionless Bloch k.

    The intercell convention is a minus sign: the minimum is at k=0 on
    every periodic torus, including odd side lengths.  All eight components
    remain present.
    """
    k = _momentum(momentum)
    intra, inter = _hoppings(intra, inter)
    gamma = gamma_matrices()
    result = np.zeros((8, 8), dtype=complex)
    for axis, value in enumerate(k):
        result += (intra - inter * math.cos(value)) * gamma[2 * axis]
        result += inter * math.sin(value) * gamma[2 * axis + 1]
    return result


def dispersion_squared(momentum, *, intra=SOURCE_A, inter=SOURCE_B):
    """Exact-form E(k)^2 = 3(a-b)^2 + 4ab sum sin(k_j/2)^2."""
    k = _momentum(momentum)
    intra, inter = _hoppings(intra, inter)
    return 3 * (intra - inter) ** 2 + 4 * intra * inter * math.fsum(
        math.sin(value / 2) ** 2 for value in k
    )


def finite_torus_gap(side, *, intra=SOURCE_A, inter=SOURCE_B):
    """Structural positive-band gap, exactly independent of torus volume.

    Every admitted periodic side contains k=0.  Positivity of all sine
    terms proves the minimum without enumerating side**3 momenta.
    """
    _level(side)
    intra, inter = _hoppings(intra, inter)
    return math.sqrt(3) * abs(intra - inter)


def refinement_coefficients(level, *, omega=SOURCE_OMEGA, imbalance=SOURCE_J):
    """a_n=(n Omega+J)/2, b_n=(n Omega-J)/2; spacing is 1/n.

    The growing hopping strengths/refinement prescription are construction
    choices.  This function supplies no physical drive or timing bound.
    """
    level = _level(level)
    omega, imbalance = _real(omega, "clock frequency"), _real(imbalance, "imbalance")
    if omega <= 0 or level * omega <= abs(imbalance):
        raise ValueError("positive omega and n*omega > abs(imbalance) are required")
    return ((level * omega + imbalance) / 2,
            (level * omega - imbalance) / 2)


def refined_dispersion_squared(momentum, level, *, omega=SOURCE_OMEGA,
                               imbalance=SOURCE_J):
    """Refined dispersion at a fixed physical-coordinate momentum p.

    'Physical-coordinate' distinguishes p from k=p/n; it does not certify
    that the chosen coordinate dictionary is physical spacetime.
    """
    p = _momentum(momentum)
    intra, inter = refinement_coefficients(level, omega=omega, imbalance=imbalance)
    # Use J directly to avoid subtracting two growing, nearly equal floats.
    return 3 * float(imbalance) ** 2 + 4 * intra * inter * math.fsum(
        math.sin(value / (2 * level)) ** 2 for value in p
    )


def continuum_dispersion_squared(momentum, *, omega=SOURCE_OMEGA,
                                 imbalance=SOURCE_J):
    p = _momentum(momentum)
    omega, imbalance = _real(omega, "clock frequency"), _real(imbalance, "imbalance")
    if omega <= 0:
        raise ValueError("a positive clock frequency is required")
    return 3 * imbalance ** 2 + omega ** 2 * math.fsum(x * x for x in p) / 4


def refinement_remainder_bound(momentum, level, *, omega=SOURCE_OMEGA,
                               imbalance=SOURCE_J):
    """Exact analytic upper formula for E_cont(p)^2-E_n(p)^2 >= 0.

    The mathematical bound is J²|p|²/(4n²)+Omega² sum(p_j^4)/(48n²).
    Its floating evaluation is not an outward numerical enclosure.
    """
    p = _momentum(momentum)
    refinement_coefficients(level, omega=omega, imbalance=imbalance)
    return (float(imbalance) ** 2 * math.fsum(x * x for x in p) / (4 * level ** 2)
            + float(omega) ** 2 * math.fsum(x ** 4 for x in p) / (48 * level ** 2))


def ring_sector_energies(*, intra=SOURCE_A, inter=SOURCE_B,
                         coupling=SOURCE_J, chemical_potential=0.0):
    """Independent exact-form spectra of the PRIMARY four-mode SSH ring.

    Hopping is -a on (0,1),(2,3), -b on (1,2),(3,0).  The density
    interaction is g sum_cycle (n_i-1/2)(n_j-1/2), with g>=0.
    A common chemical-potential term is +mu*N.  No CAR runtime or matrix
    diagonalization is consulted by these analytic formulas.
    """
    intra, inter = _hoppings(intra, inter)
    coupling = _real(coupling, "density interaction")
    chemical_potential = _real(chemical_potential, "chemical potential")
    if not intra > inter or coupling < 0:
        raise ValueError("the analytic ring branch requires a > b > 0 and g >= 0")
    omega, imbalance = intra + inter, intra - inter
    root_a = math.hypot(coupling, 4 * intra)
    root_b = math.hypot(coupling, 4 * inter)
    spectra = {
        0: (coupling,),
        1: (-omega, -imbalance, imbalance, omega),
        2: tuple(sorted((0., 0., (-coupling - root_a) / 2,
                         (-coupling + root_a) / 2, (-coupling - root_b) / 2,
                         (-coupling + root_b) / 2))),
        3: (-omega, -imbalance, imbalance, omega),
        4: (coupling,),
    }
    return {number: tuple(value + chemical_potential * number for value in values)
            for number, values in spectra.items()}


def ring_gap_predictions(*, intra=SOURCE_A, inter=SOURCE_B, coupling=SOURCE_J):
    """Finite charged gaps and a charge-conserving two-copy pair frequency.

    The neutral transfer compares N=(2,2) with N=(3,1) in two identical
    copies.  Its frequency is 2*Delta_4, distinct from the lowest neutral
    excitation within a single half-filled ring.  These are finite-cell
    energy differences, not all-volume interacting particle poles.
    """
    intra, inter = _hoppings(intra, inter)
    spectra = ring_sector_energies(intra=intra, inter=inter, coupling=coupling)
    ground = spectra[2][0]
    omega = intra + inter
    charged = (float(coupling) + math.hypot(float(coupling), 4 * intra)) / 2 - omega
    return {
        "half_filled_ground_energy": ground,
        "charged_addition_gap": charged,
        "charged_removal_gap": charged,
        "charge_transfer_pair_gap": 2 * charged,
        "single_ring_neutral_gap": spectra[2][1] - ground,
        "clock_frequency": omega,
        "charge_transfer_pair_gap_over_clock": 2 * charged / omega,
        "scope": "exact finite four-mode ring; no interacting infinite-volume pole claim",
    }
