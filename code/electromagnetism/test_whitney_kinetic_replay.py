"""Original-input controls for the independently assembled Whitney replay.

The verifier may not validate a rounded Schur subtraction by repeating it.
The controls below use the uniform-field scalar action and the analytic
icosahedral mode decomposition, not producer geometry or coefficient code.
"""
from copy import deepcopy
from decimal import Decimal
from fractions import Fraction
import hashlib
from pathlib import Path
import sys

import mpmath
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_whitney_quantum_state as state
import verify_whitney_quantum_packet as packet


def uniform_control(amplitude):
    """Original cone/Beta moments and the icosahedron adjacency spectrum.

    The absolute zero-field log determinant is the 80-digit evaluation of
      det(M) det(C.T M C) / det((M C).T M C) * det(2N)^2,
    where C is the vertex-grounded incidence matrix and M,N are the exact
    degree-two simplex masses on the golden cone. This independent constant
    is also derived from the original vertices in test_whitney_kinetic_precision.
    The field-dependent correction uses the radial mode and adjacency
    eigenvalues sqrt(5), -sqrt(5), -1 with multiplicities 3,3,5.
    """
    mp = mpmath.mp.clone()
    mp.dps = 80
    root = mp.sqrt(5)
    volume, kappa = 10+10*root/3, 30-10*root
    t = mp.mpf(amplitude)
    coupling = t/4
    a, b = 169*kappa, 114*volume/5
    remainder = a/(a+b*coupling**2)
    logdet = mp.mpf('-59.98181822720102208420414571293429952844')
    logdet -= mp.log1p(coupling**2*(19*volume/130)/(13*kappa/12))
    for adjacency, multiplicity in ((root, 3), (-root, 3), (-1, 5)):
        scalar_mass = volume*(5+adjacency)/100
        stiffness = mp.mpf(5)/3+(1-root)*adjacency/6
        logdet -= multiplicity*mp.log1p(coupling**2*scalar_mass/stiffness)
    return {
        'energy': float(b*remainder),
        'momentum': np.array([float(21*volume*remainder/10)]+[float(volume*remainder/5)]*12),
        'eta': float(-b*coupling/(a+b*coupling**2)),
        'potential': float(volume*(t*t/2+t**4/8)),
        'log_rho': float(logdet/2),
    }


@pytest.fixture(scope='module')
def geometry():
    _, vertices, edges, tets = state.exact_geometry()
    frame = np.asarray(packet.load()['orthonormal_coulomb_frame'])
    return vertices, edges, tets, frame


@pytest.mark.parametrize('amplitude', [0., 1., 1e6, 1e9])
@pytest.mark.parametrize('phase', [1, 1j])
def test_density_from_original_uniform_field(geometry, amplitude, phase):
    vertices, edges, tets, _ = geometry
    actual = state.independent_log_rho(np.full(13, phase*amplitude), vertices, edges, tets)
    assert actual == pytest.approx(uniform_control(amplitude)['log_rho'], rel=0, abs=2e-11)


@pytest.mark.parametrize('amplitude', [1., 1e6, 1e9])
@pytest.mark.parametrize('phase', [1, 1j])
def test_public_replay_retains_each_small_kinetic_component(geometry, amplitude, phase):
    vertices, edges, tets, frame = geometry
    matter = np.full(13, amplitude*phase)
    x = np.r_[12., -np.ones(12)]
    scalar_velocity = 1j*phase*x
    q = np.r_[np.zeros(42), matter.real, matter.imag]
    velocity = np.r_[np.zeros(42), scalar_velocity.real, scalar_velocity.imag]
    actual = packet.replay_phase_space(q, velocity, vertices, edges, tets, frame)
    expected = uniform_control(amplitude)
    energy = float(actual['momentum']@actual['velocity'])
    assert energy > 0
    assert energy == pytest.approx(expected['energy'], rel=2e-11, abs=0)
    block = actual['momentum'][43:] if phase == 1 else -actual['momentum'][30:43]
    np.testing.assert_allclose(block, expected['momentum'], rtol=2e-11, atol=0)
    np.testing.assert_allclose(actual['schur_scalar_potential'], expected['eta']*x, rtol=2e-11, atol=0)
    assert actual['potential'] == pytest.approx(expected['potential'], rel=2e-12, abs=0)
    assert actual['log_rho'] == pytest.approx(expected['log_rho'], rel=0, abs=2e-11)
    assert actual['gamma_min_eigenvalue'] > 0
    assert actual['cotangent_identity_defect'] < expected['energy']*1e-30


@pytest.mark.parametrize('phase', [1, 1j])
@pytest.mark.parametrize('delta', [0., 1e-100])
def test_rechart_retains_original_tiny_physical_velocity(geometry, phase, delta):
    vertices, edges, tets, frame = geometry
    incidence = np.zeros((42, 13))
    for row, (left, right) in enumerate(edges):
        incidence[row, left], incidence[row, right] = -1, 1
    xi = np.r_[0., 1., -1., np.zeros(10)]
    matter = np.full(13, phase)
    scalar_velocity = 1j*phase*xi/4
    scalar_velocity[0] = 1j*phase*delta
    q = np.r_[np.zeros(42), matter.real, matter.imag]
    velocity = np.r_[incidence@xi, scalar_velocity.real, scalar_velocity.imag]
    actual = packet.replay_phase_space(q, velocity, vertices, edges, tets, frame)
    # For every SPD edge mass, recharting the exact gradient D*xi returns
    # -xi. Thus its matter contribution cancels exactly, leaving the original
    # center perturbation. No computed mass/solve supplies this oracle.
    scalar = actual['full_velocity'][55:] if phase == 1 else -actual['full_velocity'][42:55]
    if delta:
        assert scalar[0] == pytest.approx(delta, rel=2e-12, abs=0)
        assert np.max(abs(scalar[1:])) < delta*1e-12
    else:
        # The exact-zero counterpart gets a finite absolute numerical bound;
        # no floating gauge solve is being certified as an exact zero.
        assert np.max(abs(scalar)) < 1e-70
    assert np.max(abs(actual['full_velocity'][:42])) < (delta*1e-12 if delta else 1e-70)


@pytest.mark.parametrize('phase', [1, 1j])
@pytest.mark.parametrize('representation', [Fraction, Decimal])
def test_rechart_preserves_tiny_exact_difference_near_one(geometry, phase, representation):
    vertices, edges, tets, frame = geometry
    incidence = np.zeros((42, 13))
    for row, (left, right) in enumerate(edges):
        incidence[row, left], incidence[row, right] = -1, 1
    xi = np.r_[0., 1., -1., np.zeros(10)]
    amplitude = (Fraction(10**100+1, 10**100) if representation is Fraction else
                 Decimal('1.'+'0'*99+'1'))
    q = [0]*42+([amplitude]*13+[0]*13 if phase == 1 else [0]*13+[amplitude]*13)
    velocity = np.r_[incidence@xi, np.zeros(13), xi/4] if phase == 1 else np.r_[incidence@xi, -xi/4, np.zeros(13)]
    actual = packet.replay_phase_space(q, velocity, vertices, edges, tets, frame)
    scalar = actual['full_velocity'][55:] if phase == 1 else -actual['full_velocity'][42:55]
    # i*xi/4 - i*(1+epsilon)*xi/4 = -i*epsilon*xi/4.
    # Magnitudes alone cannot reveal the 100 significant decimal places.
    assert scalar[1] == pytest.approx(-2.5e-101, rel=2e-12, abs=0)
    assert scalar[2] == pytest.approx(2.5e-101, rel=2e-12, abs=0)
    assert np.max(abs(scalar[[0, *range(3, 13)]])) < 1e-112


def test_moment_replay_accepts_exact_real_presentations(geometry):
    vertices, edges, tets, _ = geometry
    expected = uniform_control(1)['log_rho']
    for matter in ([1]*13, [Fraction(1)]*13, [Decimal('1')]*13, np.ones(13, np.float32)):
        assert state.independent_log_rho(matter, vertices, edges, tets) == pytest.approx(expected, abs=2e-11)


@pytest.mark.parametrize('phase', [1, 1j])
def test_nonuniform_complex_dressing_matches_individual_simplex_moments(geometry, phase):
    vertices, edges, tets, _ = geometry
    amplitude = 8.
    matter = np.zeros(13, complex)
    matter[0] = phase*amplitude
    metric, mass, _ = packet.metric_and_potential(matter, vertices, edges, tets)
    # Every center-to-boundary edge meets five tetrahedra. Its dressing
    # column is i*T/4 * lambda_0 lambda_i (real -T/4 for imaginary psi_0).
    # E[lambda_0^2 lambda_i]=1/60 and E[lambda_0^2 lambda_i^2]=1/210.
    element_volume = (3+np.sqrt(5))/6
    column = 55 if phase == 1 else 42
    sign = 1 if phase == 1 else -1
    for edge, (left, right) in enumerate(edges):
        if left == 0:
            assert metric[edge, column] == pytest.approx(sign*5*element_volume*amplitude/120,
                                                        rel=2e-13, abs=0)
            assert metric[edge, edge]-mass[edge, edge] == pytest.approx(
                5*element_volume*amplitude**2/1680, rel=2e-13, abs=0)


@pytest.mark.parametrize('bad', [True, '1', float('nan'), complex(1, float('inf')),
                                Decimal('NaN'), np.ma.masked, np.ma.array(1., mask=True)])
def test_state_replay_refuses_malformed_original_components(geometry, bad):
    vertices, edges, tets, _ = geometry
    matter = [1]*13
    matter[7] = bad
    with pytest.raises(ValueError, match='replay'):
        state.independent_log_rho(matter, vertices, edges, tets)


@pytest.mark.parametrize('bad', [True, '1', float('inf'), 1j, Decimal('Infinity'),
                                np.ma.masked, np.ma.array(1., mask=True)])
def test_packet_replay_refuses_malformed_original_components(geometry, bad):
    vertices, edges, tets, frame = geometry
    q = [0]*42+[1]*13+[0]*13
    q[43] = bad
    with pytest.raises(ValueError, match='replay'):
        packet.replay_phase_space(q, [0]*68, vertices, edges, tets, frame)


@pytest.mark.parametrize('diagonal', [(-1, -1), (0, 1)])
def test_constrained_replay_requires_spd_even_with_nonnegative_determinant(diagonal):
    mp = mpmath.mp.clone()
    mp.dps = 80
    # Two negative reduced directions have a positive determinant. The
    # unrelated gauge direction is positive and fully constrained.
    metric = mp.diag([*diagonal, 1])
    vertical = mp.matrix([[0], [0], [1]])
    section = mp.matrix([[1, 0], [0, 1], [0, 0]])
    with pytest.raises(ValueError, match='reduced kinetic metric.*positive'):
        state.reduced_moments(mp, metric, vertical, section)


def test_constrained_replay_accepts_positive_metric_and_rejects_asymmetry():
    mp = mpmath.mp.clone()
    mp.dps = 80
    vertical = mp.matrix([[0], [0], [1]])
    section = mp.matrix([[1, 0], [0, 1], [0, 0]])
    metric = mp.matrix([[2, 0, 1], [0, 3, 0], [1, 0, 4]])
    gamma, eta, log_rho = state.reduced_moments(mp, metric, vertical, section)
    assert gamma == mp.diag([mp.mpf(7)/4, 3])
    assert eta == mp.matrix([[-mp.mpf(1)/4, 0]])
    assert float(log_rho) == pytest.approx(float(mp.log(mp.mpf(21)/4)/2), rel=1e-15)
    metric[0, 1] = 1
    with pytest.raises(ValueError, match='symmetry'):
        state.reduced_moments(mp, metric, vertical, section)


def test_packet_replay_rejects_invalid_coulomb_frame(geometry):
    vertices, edges, tets, frame = geometry
    changed = frame.copy()
    changed[:, 0] *= 2
    with pytest.raises(ValueError, match='orthonormal Coulomb frame'):
        packet.replay_phase_space([0]*68, [0]*68, vertices, edges, tets, changed)


@pytest.mark.parametrize('verifier', [state, packet])
def test_existing_numerical_receipt_values_pass_fresh_independent_replay(verifier):
    # Keep the retained values unchanged. Only authenticate this working
    # source revision in the in-memory test copy; no artifact is regenerated.
    receipt = deepcopy(verifier.load())
    for name in receipt['source_pins']:
        receipt['source_pins'][name] = hashlib.sha256((verifier.ROOT/name).read_bytes()).hexdigest()
    assert verifier.verify(receipt)['accepted'] is True
