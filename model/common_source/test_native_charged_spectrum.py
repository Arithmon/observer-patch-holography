"""Independent CAR, Poisson, pole and finite-volume controls."""
from copy import deepcopy
from fractions import Fraction as F
import math

import numpy as np
import pytest
import sympy as s

from . import native_charged_spectrum as c


def _annihilator(site, cells):
    z, ladder = s.diag(1, -1), s.Matrix([[0, 1], [0, 0]])
    factors = [ladder if j == site else z if j < site else s.eye(2)
               for j in reversed(range(cells))]
    return s.kronecker_product(*factors)


def test_full_CAR_action_includes_density_interaction_in_every_charge_sector():
    cells = 4
    a = [_annihilator(x, cells) for x in range(cells)]
    n = [x.T*x for x in a]
    J, masses = c.native_parameters()
    m = masses[1]
    h = sum((s.Rational(m)*(-1)**x*n[x] for x in range(cells)), s.zeros(1 << cells))
    for x in range(cells-1):
        h += s.Rational(J)*(n[x]+n[x+1]-2*n[x]*n[x+1]-a[x].T*a[x+1]-a[x+1].T*a[x])
    for charge in range(cells+1):
        basis, actual = c.sector_matrix(cells, charge, m, J)
        assert h.extract(basis, basis) == actual


@pytest.mark.parametrize('sector', (0, 1, 2))
def test_exact_dimer_ground_charge_and_both_CAR_poles(sector):
    J, masses = c.native_parameters()
    m = masses[sector]
    cert = c.dimer_certificate(m, J)
    lo, hi = map(F, cert['charged_gap_interval'])
    assert lo > 0 and lo*(lo+2*J) <= m*m <= hi*(hi+2*J)
    assert hi-lo < F(1, 10**16)
    report, arrays = c.dimer_control(m, J)
    H, P = arrays['dimer_hamiltonian'], arrays['dimer_ground_projector']
    assert np.linalg.norm(P@P-P) < 1e-14
    assert np.trace(P).real == pytest.approx(1)
    assert np.linalg.norm(H@P+float((lo+hi)/2)*P) < 1e-14
    poles = report['local_CAR_poles']
    for site in (0, 1):
        assert sum(x['residue'] for x in poles if x['site'] == site) == pytest.approx(1)
    assert arrays['dimer_free_lift'][3, 3] == pytest.approx(2*float(J))
    assert H[3, 3] == 0


@pytest.mark.parametrize('cells', (2, 4, 6))
def test_source_exchange_is_exclusion_and_poisson_formula_is_exact(cells):
    lap = s.zeros(cells)
    for x in range(cells-1):
        lap[x, x] += 1; lap[x+1, x+1] += 1
        lap[x, x+1] -= 1; lap[x+1, x] -= 1
    w = s.Matrix([(-1)**x for x in range(cells)])
    f = (lap+s.ones(cells)/cells).inv()*w
    for charge in range(1, cells):
        basis, h = c.sector_matrix(cells, charge, F(1), F(1))
        diagonal = [sum(w[x]*((mask >> x)&1) for x in range(cells)) for mask in basis]
        H0 = h-s.diag(*diagonal)
        u = s.ones(len(basis), 1)
        v = s.Matrix(diagonal)
        z = s.Matrix([sum(f[x]*((mask >> x)&1) for x in range(cells)) for mask in basis])
        assert H0*u == s.zeros(len(basis), 1)
        assert H0*z == v
        assert (v.T*z)[0]/len(basis) == c.poisson_coefficient(cells, charge)
        assert c.poisson_coefficient(cells, charge) == s.Rational(charge*(cells-charge), 2*(cells-1))


def test_reflection_evenness_and_admitted_remainder():
    J, masses = c.native_parameters()
    m = masses[-1]
    basis, h = c.sector_matrix(4, 2, m, J)
    reflected = [basis.index(int(f'{mask:04b}'[::-1], 2)) for mask in basis]
    w = s.diag(*(sum((-1)**x*((mask >> x)&1) for x in range(4)) for mask in basis))
    assert h.extract(reflected, reflected) == h-2*s.Rational(m)*w
    bound = c.perturbative_bound(4, 2, m, J)
    assert bound['admitted']
    result, _ = c.finite_spectrum(4, m, J)
    actual = result['half_filled_ground_energy']
    assert abs(actual-float(F(bound['second_order_energy']))) <= float(F(bound['remainder_upper']))+c.EIGENPAIR_ADMISSION
    assert not c.perturbative_bound(10, 5, masses[0], J)['admitted']


def test_nonuniform_limit_is_explicit_in_the_gap_coefficient():
    J, masses = c.native_parameters()
    m = masses[-1]
    coefficients = []
    for L in c.SIZES:
        half = c.poisson_coefficient(L, L//2)
        adjacent = c.poisson_coefficient(L, L//2+1)
        coefficients.append(half-adjacent)
        assert half-adjacent == s.Rational(1, 2*(L-1))
    assert all(a > b for a, b in zip(coefficients, coefficients[1:]))
    assert c.perturbative_bound(10, 5, m, J)['isolation_ratio'] != '0'


def test_exact_dimer_certificate_preserves_tiny_positive_scales():
    scale = F(1, 10**100)
    ordinary = c.dimer_certificate(F(1), F(1))
    tiny = c.dimer_certificate(scale, scale)
    assert list(map(F, tiny['charged_gap_interval'])) == [scale*F(x) for x in ordinary['charged_gap_interval']]
    assert tiny['site_zero_removal_residue_interval'] == ordinary['site_zero_removal_residue_interval']
    assert tiny['site_zero_addition_residue_interval'] == ordinary['site_zero_addition_residue_interval']


@pytest.mark.parametrize('bad', (True, 0, -1, 0.1, float('nan'), '1/3'))
def test_original_exact_parameter_domain_refuses(bad):
    with pytest.raises(ValueError):
        c.dimer_certificate(bad, F(1))


@pytest.mark.parametrize('cells', (True, 1, 3, 12, 4.0))
def test_geometry_domain_refuses(cells):
    with pytest.raises(ValueError):
        c.sector_matrix(cells, 1, F(1), F(1))


def test_complete_prediction_rejects_mass_and_admission_forgery():
    report, _ = c.build()
    assert report['external_comparison']['ratio_overlaps'] == [False, False]
    assert report['empirical_inputs_to_prediction'] == []
    assert not report['scope']['observed_particle_spectrum_derived']
    for path in ('native_exchange_exact', 'numerical_eigenpair_admission'):
        forged = deepcopy(report)
        forged[path] = '1/2' if path.endswith('exact') else 2*c.EIGENPAIR_ADMISSION
        with pytest.raises(ValueError):
            c._compare(report, forged)
    forged = deepcopy(report)
    forged['finite_charge_spectra'][0]['addition_gap_diagnostic'] = float('nan')
    with pytest.raises(ValueError):
        c._compare(report, forged)
