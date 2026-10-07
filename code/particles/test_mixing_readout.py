"""Independent matrix controls for all existing CKM/PMNS readout entry points."""
import importlib
import math
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for directory in (ROOT, ROOT/'particles/flavor', ROOT/'particles/neutrino'):
    sys.path.insert(0, str(directory))

READERS = [
    ('neutrino.build_pmns_from_shared_flavor_basis', '_standard_pmns_parameters', 'signed'),
    ('neutrino.derive_neutrino_dimensionless_law_candidate_audit', '_pmns_parameters', 'degrees'),
    ('neutrino.derive_neutrino_two_parameter_exact_adapter', '_pmns_parameters', 'degrees'),
    ('neutrino.derive_neutrino_weighted_cycle_repair', '_pmns_parameters', 'canonical'),
    ('neutrino.derive_neutrino_physical_majorana_phase_theorem', '_standard_pmns_parameters', 'canonical'),
    ('neutrino.derive_neutrino_weighted_cycle_shared_basis_representation', '_standard_pmns_parameters', 'canonical'),
    ('flavor.derive_quark_d12_mass_branch_and_ckm_residual', '_standard_ckm_parameters', 'ckm'),
    ('flavor.sigma_ud_orbit_provider', '_standard_ckm_parameters', 'ckm'),
    ('flavor.enumerate_quark_local_basis_orbit_diagnostic', '_ckm_tuple', 'ckm'),
    ('flavor.derive_quark_transport_frame_diagnostic_orbit', '_ckm_tuple', 'ckm'),
]


@pytest.mark.parametrize('module,name', [(module, name) for module, name, _ in READERS])
@pytest.mark.parametrize('matrix', [np.full((3, 3), .2), np.full((3, 3), .2) + .5*np.eye(3)])
def test_every_reader_rejects_a_nonunitary_matrix(module, name, matrix):
    reader = getattr(importlib.import_module('particles.'+module), name)
    with pytest.raises(ValueError):
        reader(matrix)


from particles.mixing import mixing_parameters, pmns_signed


def rotations(a=.59, b=.81, c=.15, delta=1.2):
    """Independent forward PDG construction as three elementary rotations."""
    ca, cb, cc = map(math.cos, (a, b, c))
    sa, sb, sc = map(math.sin, (a, b, c))
    z = complex(math.cos(delta), math.sin(delta))
    r12 = np.array([[ca, sa, 0], [-sa, ca, 0], [0, 0, 1]])
    r13 = np.array([[cc, 0, sc*z.conjugate()], [0, 1, 0], [-sc*z, 0, cc]])
    r23 = np.array([[1, 0, 0], [0, cb, sb], [0, -sb, cb]])
    return r23 @ r13 @ r12, ca*sa*cb*sb*cc**2*sc*math.sin(delta)


@pytest.mark.parametrize('module,name,schema', READERS)
@pytest.mark.parametrize('delta', [0., .000001, 1.2, math.pi, 4.7, math.tau-.000001])
def test_legacy_schemas_against_independent_rotations(module, name, schema, delta):
    u, j = rotations(delta=delta)
    # Includes arbitrary row phases and Majorana-like column phases.
    u = np.exp(1j*np.array([.4, -1.7, 2.1]))[:, None]*u*np.exp(1j*np.array([-.8, 1.1, .2]))
    result = getattr(importlib.import_module('particles.'+module), name)(u)
    if schema in ('signed', 'ckm'):
        key = 'delta_pmns' if schema == 'signed' else 'delta_ckm'
        expected = dict(zip(('theta_12', 'theta_23', 'theta_13', key, 'jarlskog'), (.59, .81, .15, delta, j)))
        if key == 'delta_pmns' and delta > math.pi:
            expected[key] -= math.tau
    else:
        expected = {k+'_deg': math.degrees(v) for k, v in zip(('theta12', 'theta23', 'theta13', 'delta'), (.59, .81, .15, delta))}
        expected['J'] = j
        if schema == 'canonical':
            expected.update({k+'_rad': v for k, v in zip(('theta12', 'theta23', 'theta13', 'delta'), (.59, .81, .15, delta))})
    assert result.keys() == expected.keys()
    for key, value in expected.items():
        difference = result[key] - value
        if key.startswith('delta'):
            period = 360 if key.endswith('_deg') else math.tau
            assert (-math.pi <= result[key] <= math.pi) if schema == 'signed' else (0 <= result[key] < period)
            difference = math.remainder(difference, period)
        assert abs(difference) < 1e-11, (key, result[key], value)


@pytest.mark.parametrize('c,delta', [(1e-40, 1.2), (.0035, 1e-20), (1e-300, -1.2)])
def test_tiny_mixing_and_phase_are_not_floored(c, delta):
    u, j = rotations(c=c, delta=delta)
    result = pmns_signed(u)
    assert result['theta_13'] == pytest.approx(c, rel=2e-14, abs=0)
    assert result['delta_pmns'] == pytest.approx(delta, rel=2e-14, abs=0)
    assert result['jarlskog'] == pytest.approx(j, rel=2e-14, abs=0)


@pytest.mark.parametrize('bad', [np.eye(2), np.ones((3, 3)), np.eye(3)*1.01,
    rotations()[0]*1.01,
    np.full((3, 3), np.nan), np.full((3, 3), np.inf), np.full((3, 3), 1e308),
    np.eye(3, dtype=bool), [[True, 0, 0], [0, 1, 0], [0, 0, 1]],
    np.ma.array(np.eye(3), mask=np.eye(3)), np.eye(3).astype(str)])
def test_rejects_invalid_inputs_without_repair(bad):
    with pytest.raises(ValueError):
        mixing_parameters(bad)


@pytest.mark.parametrize('u', [np.eye(3), np.eye(3)[[1, 2, 0]], rotations(c=0)[0], rotations(a=0)[0], rotations(b=0)[0]])
def test_unidentifiable_phase_is_refused(u):
    with pytest.raises(ValueError, match='undefined'):
        mixing_parameters(u)


def test_cp_conjugation_and_random_roundtrips():
    rng = np.random.default_rng(1033)
    for _ in range(100):
        angles = rng.uniform(.01, 1.56, 3)
        delta = rng.uniform(.01, math.tau-.01)
        u, j = rotations(*angles, delta)
        result = mixing_parameters(u.conjugate())
        rebuilt, _ = rotations(*(result[k+'_rad'] for k in ('theta12', 'theta23', 'theta13', 'delta')))
        np.testing.assert_allclose(abs(rebuilt), abs(u), atol=2e-14, rtol=0)
        assert result['J'] == pytest.approx(-j, abs=2e-16)
        assert result['delta_rad'] == pytest.approx(math.tau-delta, abs=2e-14)


def test_nonzero_unrepresentable_jarlskog_is_not_reported_as_zero():
    u, _ = rotations(c=math.ulp(0.0), delta=1.2)
    with pytest.raises(ValueError, match='underflows'):
        mixing_parameters(u)


def test_large_theta13_does_not_erase_the_other_angles():
    u, _ = rotations(c=math.nextafter(math.pi/2, 0))
    result = mixing_parameters(u)
    assert result['theta12_rad'] == pytest.approx(.59, abs=1e-14)
    assert result['theta23_rad'] == pytest.approx(.81, abs=1e-14)
    assert result['delta_rad'] == pytest.approx(1.2, abs=1e-14)


def test_subnormal_anchor_with_zero_cp_phase_remains_readable():
    result = pmns_signed(rotations(c=math.ulp(0.), delta=0.)[0])
    assert result['theta_13'] == math.ulp(0.)
    assert result['delta_pmns'] == result['jarlskog'] == 0.


def test_jarlskog_is_correctly_rounded_on_supplied_entries():
    import mpmath as mp

    rng = np.random.default_rng(1056)
    # 8192 bits suffice for the exact four-factor product of binary64 entries.
    with mp.workprec(8192):
        for _ in range(30):
            u, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))
            v = [[mp.mpc(float(z.real), float(z.imag)) for z in row] for row in u]
            expected = mp.im(v[0][0]*v[1][1]*mp.conj(v[0][1])*mp.conj(v[1][0]))
            assert mixing_parameters(u)['J'] == float(expected)
