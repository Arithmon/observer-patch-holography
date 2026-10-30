"""Independent full matrices, translation, thermal and spectral controls."""
import json
import numpy as np
import pytest
from scipy.linalg import expm

from .dispersion import PeriodicSpectrum, ROOTS, build
from .test_response import jw


@pytest.fixture(scope='module')
def ring():
    return PeriodicSpectrum()


def independent_action():
    a = jw(8)
    n = [operator.conj().T@operator for operator in a]
    h = sum(n[2*x]@n[2*x+1] for x in range(4))
    for x in range(4):
        for spin, sign in enumerate((1, -1)):
            term = a[2*x+spin].conj().T@a[2*((x+1) % 4)+spin]
            h += 1j*sign*(term-term.conj().T)/16
    return h, h-sum(n)/2, a


def test_complete_periodic_action_and_translation_against_independent_jw(ring):
    h, k, a = independent_action()
    assert np.array_equal(h, ring.hamiltonian)
    assert np.array_equal(k, ring.generator)
    assert np.array_equal(np.linalg.matrix_power(ring.shift, 4), np.eye(256))
    for j in range(8):
        assert np.array_equal(ring.shift@a[j]@ring.shift.conj().T, a[(j+2) % 8])
    # Dropping occupation permutation signs would falsely implement bosonic translation.
    wrong = abs(ring.shift)
    assert np.linalg.norm(wrong@a[0]@wrong.T-a[2], 2) > 1
    eigenvalues = np.linalg.eigvalsh(k)
    assert np.max(abs(eigenvalues-ring.energies)) < 2e-14
    assert ring.ground_degeneracy == 1
    assert ring.labels[0] == (2, 2, 2)


def test_all_translation_blocks_and_spin_selection_keep_the_complete_field(ring):
    assert ring.vectors.shape == (256, 256)
    for j, root in enumerate(ROOTS):
        columns = ring.vectors[:, [label[2] == j for label in ring.labels]]
        assert np.linalg.norm(ring.shift@columns-root*columns, 2) < 1e-14
    sx, sy, sz = ring.spin
    assert np.array_equal(sx@sy-sy@sx, 1j*sz)
    assert np.array_equal(ring.shift@sx@ring.shift.conj().T, -sx)
    assert np.array_equal(ring.shift@sz@ring.shift.conj().T, sz)
    # The internal spin rotations are not a Lorentz-spin map commuting with translations.
    for number in (3, 5):
        indices = np.flatnonzero(ring.up+ring.down == number)
        p = ring.odd_spin_half[np.ix_(indices, indices)]
        assert np.linalg.norm(p@p-p, 2) < 1e-14


@pytest.mark.parametrize('beta', [.5, 1.])
def test_gibbs_control_uses_all_number_and_momentum_sectors(ring, beta):
    _, k, _ = independent_action()
    independent = expm(-beta*k); independent /= np.trace(independent)
    actual = ring.thermal_density(beta)
    assert np.linalg.norm(actual-independent, 2) < 1e-14
    assert np.linalg.eigvalsh(actual)[0] > 1e-4
    assert abs(np.trace(actual@ring.number)-4) < 1e-13
    assert np.linalg.norm(ring.shift@actual-actual@ring.shift, 2) < 1e-14


def test_momentum_probes_full_lehmann_weights_moments_and_time_domain(ring):
    _, k, a = independent_action()
    ground = ring.ground(); unitary = expm(-.31j*(k-ring.energies[0]*np.eye(256)))
    for momentum, root in enumerate(ROOTS):
        for spin in range(2):
            sums = []
            for charge in (1, -1):
                creation = sum(root**(-x)*a[2*x+spin].conj().T for x in range(4))/2
                op = creation if charge == 1 else creation.conj().T
                assert np.linalg.norm(ring.shift@op@ring.shift.conj().T-root**charge*op, 2) == 0
                response = ring.response(momentum, spin, charge)
                state = op@ground
                assert abs(np.vdot(state, unitary@state)
                           -response['full_weights']@np.exp(-.31j*response['gaps'])) < 2e-14
                assert abs(sum(row['weight'] for row in response['poles'])
                           -response['full_weights'].sum()) < 1e-14
                assert response['lowest_probe_residue'] > .32
                sums.append(response['zeroth_moment'])
            assert abs(sum(sums)-1) < 1e-14
    response = ring.response(1, 1, 1)
    # The lower quartet stays in the complete catalogue, with a proved spin-selection map.
    assert response['poles'][0]['minimum_gap'] < response['lowest_spin_half_gap']
    assert abs(response['poles'][0]['spin_half_multiplicity']) < 1e-14
    assert response['poles'][0]['weight'] < 1e-25


def test_fixed_velocity_failure_and_underconstrained_renormalized_fit_are_distinct():
    report = build(); json.dumps(report, allow_nan=False)
    b = report['benchmark']
    assert b['declared_dimer_anchor_and_speed_benchmark_satisfied'] is False
    assert b['same_volume_anchor_and_declared_speed_benchmark_satisfied'] is False
    assert b['lowest_envelope_same_volume_benchmark_satisfied'] is False
    assert b['lowest_envelope_same_volume_anchor_max_squared_residual'] > .03
    assert b['all_spin_resolved_best_fit_max_squared_residual'] > .06
    assert abs(b['envelope_fitted_slope_over_supplied_c_squared']-3.1217092863) < 1e-9
    assert b['distinct_sin_squared_coordinates'] == 2
    assert b['envelope_fit_is_underdetermined_as_a_dispersion_test'] is True
    assert b['emergent_Lorentz_with_renormalized_velocity_excluded'] is False
    assert report['whole_OPH_or_continuum_no_go'] is False


def test_free_bloch_control_retains_doublers_and_rejects_an_invented_unique_vacuum():
    free = PeriodicSpectrum(interaction=0.)
    assert free.ground_degeneracy == 16
    assert free.half_filled_ground_degeneracy == 6
    with pytest.raises(ValueError, match='unique'):
        free.ground()
    for momentum, root in enumerate(ROOTS):
        basis = np.vstack([root**(-x)*np.eye(2)/2 for x in range(4)])
        expected = .125*complex(root).imag*np.diag([1, -1])
        assert np.array_equal(basis.conj().T@free.single_h@basis, expected)


@pytest.mark.parametrize('args', [(True, 0, 1), (4, 0, 1), (0, -1, 1), (0, 0, True), (0, 0, 0)])
def test_malformed_momentum_probes_fail_closed(ring, args):
    with pytest.raises(ValueError):
        ring.probe(*args)
