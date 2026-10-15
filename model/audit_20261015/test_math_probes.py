"""Independent algebra probes for the 15 October audit; no packet is changed."""
from fractions import Fraction
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.anharmonic_transfer.environment import quadrature_bound
from model.observer.source_scalar import compiler
from model.scalar_observations.experimental_readout import canonical_source


def _number_phase_word(dimension, angle):
    bits = dimension.bit_length()-1
    return [compiler.pulse('I', -angle*(dimension-1)/2, [0])] + [
        compiler.pulse('Z', angle*2**bit/2, [bit]) for bit in range(bits)]


@pytest.mark.parametrize('dimension', [4, 8, 32, 64])
def test_number_drift_needs_only_one_scalar_and_one_phase_per_bit(dimension):
    """Check the full matrix, including global phase and high occupations."""
    angle = .123
    context = canonical_source._context()
    generic = canonical_source._diagonal_word(
        context, [-canonical_source._parameter(context, angle, 'angle')*n
                  for n in range(dimension)])
    short = _number_phase_word(dimension, angle)
    matrix = compiler.execute(short, dimension.bit_length()-1)
    target = np.diag(np.exp(-1j*angle*np.arange(dimension)))
    assert np.linalg.norm(matrix-target, 2) < 1e-13
    assert len(short) == dimension.bit_length()
    assert len(short) < len(generic)
    if dimension == 64:
        assert len(generic) == 838
        assert len(short) == 7


@pytest.mark.parametrize('dimension', [4, 8, 32, 64])
def test_controlled_occupation_parity_is_a_two_native_pulse_word(dimension):
    """The oscillator parity depends only on its least significant bit."""
    spin_bit = dimension.bit_length()-1
    phases = np.r_[np.zeros(dimension), math.pi*(np.arange(dimension) % 2)]
    generic = compiler.compile_diagonal(phases)
    short = [compiler.pulse('J', -math.pi, [spin_bit, 0]),
             compiler.pulse('I', math.pi, [0])]
    matrix = compiler.execute(short, spin_bit+1)
    assert np.linalg.norm(matrix-np.diag(np.exp(1j*phases)), 2) < 1e-13
    if dimension == 64:
        assert len(generic) == 2054
        assert len(short) == 2


@pytest.mark.parametrize('dimension,reference', [(3, 1), (4, 2), (8, 3)])
def test_bosonic_lindblad_truncation_residual_is_only_a_boundary_row(dimension, reference):
    """Embed the truncated state into exactly one extra oscillator level."""
    rng = np.random.default_rng(519+dimension+reference)
    spin_reference = 2*reference
    d = dimension
    raw = rng.normal(size=(spin_reference*d, 5))+1j*rng.normal(size=(spin_reference*d, 5))
    rho = raw@raw.conj().T
    rho /= np.trace(rho)
    embed = np.kron(np.eye(spin_reference), np.eye(d+1)[:, :d])
    embedded = embed@rho@embed.conj().T
    x = np.kron(np.eye(reference), [[0., 1.], [1., 0.]])
    delta, omega, gamma = 1.2, 2.7, .8

    def generator(dim, state):
        a = np.diag(np.sqrt(np.arange(1, dim)), 1)
        n = np.kron(np.eye(spin_reference), np.diag(np.arange(dim)))
        h = delta*n-omega/2*np.kron(x, a+a.T)
        commutator = n@state-state@n
        return -1j*(h@state-state@h)-gamma/2*(n@commutator-commutator@n)

    residual = generator(d+1, embedded)-embed@generator(d, rho)@embed.conj().T
    top_row = np.kron(np.eye(spin_reference), np.eye(d)[-1:])
    measured = np.linalg.norm(residual, ord='nuc')/2
    boundary = omega*math.sqrt(d)/2*np.linalg.norm(top_row@rho, ord='nuc')
    top_mass = np.trace(top_row@rho@top_row.conj().T).real
    upper = omega*math.sqrt(d*top_mass)/2
    assert measured == pytest.approx(boundary, rel=2e-13, abs=1e-13)
    assert boundary <= upper+1e-13


def test_infinite_oscillator_displacement_identity_fails_at_finite_cutoff():
    """Counterexample to an attractive but invalid finite-matrix replacement."""
    dimension, delta, force, duration = 8, 1.2, 1.1, .31
    a = np.diag(np.sqrt(np.arange(1, dimension)), 1)
    n = np.diag(np.arange(dimension))
    beta = force/(2*delta)
    direct = expm(-1j*duration*(delta*n-force/2*(a+a.T)))
    z = beta*(1-np.exp(-1j*delta*duration))
    candidate = (np.exp(1j*beta**2*(delta*duration-np.sin(delta*duration)))
                 * expm(z*a.T-z.conjugate()*a) @ expm(-1j*delta*duration*n))
    # The infinite-CCR identity is excellent on a low-energy state here but
    # cannot identify the two finite operators uniformly.
    assert np.linalg.norm((candidate-direct)[:, 0]) < 1e-8
    assert np.linalg.norm(candidate-direct, 2) > 1e-2


def test_smaller_gaussian_rules_have_an_a_priori_error_budget():
    """A numerical simplification, not a claim about changed data scores."""
    d, subdivisions, steps = 64, 256, 50
    half_duration = 50e-6/(2*subdivisions)
    maximum_gap = (math.pi/10)/(math.pi/35e-6)
    drive = Fraction(quadrature_bound(d, 18., half_duration, 4)
                     ['exact_half_diamond_error_upper'])
    gap = Fraction(quadrature_bound(d, 18., maximum_gap, 6)
                   ['exact_half_diamond_error_upper'])
    complete = 8*steps*subdivisions*drive+4*steps*gap
    assert 6.9e-6 < float(complete) < 7.1e-6
