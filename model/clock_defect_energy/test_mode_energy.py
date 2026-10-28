"""Independent native-kernel, source-word and dispersion checks."""
from fractions import Fraction as F

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

from . import mode_energy as m
from ..completion import observer


def test_complete_native_quantum_kernel_gives_the_current_form():
    value = m.current_forms()
    basis = m.native.states(binary=True)
    index = {state: i for i, state in enumerate(basis)}
    z = sp.Matrix([[sp.Rational(n)-sp.Rational(1, 6) for n in state] for state in basis])
    for law in ('C', 'E'):
        counts = sp.zeros(66)
        # Separate complete Kraus-map producer on H78, explicitly restricted.
        for seam in observer.SEAMS:
            for _, _, _, pairs in observer._raw_maps(law, seam):
                for i, j in pairs:
                    if observer.STATES[i] in index:
                        counts[index[observer.STATES[i]], index[observer.STATES[j]]] += 1
        assert counts*sp.ones(66, 1) == 60*sp.ones(66, 1)
        assert counts == counts.T
        conditional = (sp.eye(66)-counts/60)
        assert z.T*conditional*z/33 == value['increment_second_moment']
        assert z.T*z/66 == value['covariance']
    assert value['mode_variance'] == sp.Rational(20, 33)
    assert value['mode_increment_second_moment'] == sp.Rational(4, 33)


def test_form_is_positive_with_only_constant_kernel_and_rescaling_cancels():
    value = m.current_forms()
    form, cov = value['increment_second_moment']/2, value['covariance']
    assert form*sp.ones(12, 1) == sp.zeros(12, 1)
    assert form.rank() == cov.rank() == 11
    lap = sp.zeros(12)
    for a, b in observer.SEAMS:
        edge = sp.zeros(12, 1); edge[a], edge[b] = -1, 1
        lap += edge*edge.T
    assert form == lap/396  # Sum of squared real edge differences.
    v = sp.Matrix(value['mode_profile'])
    for factor in (sp.Rational(1, 7), -3, 11):
        w = factor*v
        assert (w.T*form*w)[0]/(w.T*cov*w)[0] == sp.Rational(1, 10)


def test_lean_attachment_transports_both_ports_and_signed_currents():
    value = m.registered_incidence_map()
    permutation = value['source_to_registered_lean_port']
    assert permutation != list(range(12))
    assert len(value['source_seam_maps']) == 30
    assert {row['registered_seam_index'] for row in value['source_seam_maps']} == set(range(30))
    assert {row['current_orientation'] for row in value['source_seam_maps']} == {-1, 1}
    for row in value['source_seam_maps']:
        a, b = row['source_seam']
        assert row['current_orientation']*(permutation[b]-permutation[a]) > 0


def test_existing_full_port_mode_and_product_word_are_reused():
    witness = m.native_mode()
    assert witness['ports'] == 240
    assert not witness['internal_fixture_survives_as_a_bounded_nonconstant_continuum_mode']
    # q=2,L=6,k=1 gives exactly radial frequency1/2.
    original = m.continuum.product_mode_source_control(
        radial_points=2, length=F(6), steps=1)
    assert m.source_word(F(1, 2), 1) == original['source_step_word']


@pytest.mark.parametrize('n', (1, 2, 4, 8))
def test_exact_ratio_enclosure_against_independent_high_precision_derivative(n):
    value = m.ratio_bounds(n)
    with mp.workdps(75):
        theta = mp.acos(mp.mpf(209)/210)
        q = mp.mpf(630)
        gap = lambda p: 2*n*mp.sqrt(q)*mp.acos(mp.cos(theta/(2*n))*mp.cos(p/(2*n*mp.sqrt(q))))
        rest = gap(0)
        inverse_mass = mp.diff(gap, 0, 2)
        ratio = rest*inverse_mass
        def exact(text):
            r = F(text)
            return mp.mpf(r.numerator)/r.denominator
        low, high = map(exact, value['rest_to_independent_inertia_ratio'])
        assert low <= ratio <= high < 1
        assert 1-ratio <= exact(value['mass_shell_defect_upper'])
        assert abs(rest-mp.sqrt(q)*theta) < mp.mpf('1e-65')


def test_subdivision_reduces_error_without_relabelling_carrier_refinement():
    bounds = [F(m.ratio_bounds(n)['mass_shell_defect_upper']) for n in (1, 2, 4, 8)]
    assert all(after < before/3 for before, after in zip(bounds, bounds[1:]))
    response = m.symbolic_response()
    assert not response['pulse_subdivision_is_carrier_refinement']
    assert not response['radial_symbol_is_native_physical_momentum']


def test_massless_symbol_measures_the_same_speed_independently():
    # Native constant angular mode has zero Z phase. Keep that zero pulse.
    for p in (F(1, 5), F(2, 5), F(4, 5)):
        word = [m.compiler.pulse('Z', 0., [0]),
                m.compiler.pulse('X', float(p)/(2*np.sqrt(m.Q)), [0])]
        unitary = m.compiler.execute(word, 1)
        gap = 2*np.sqrt(m.Q)*np.arccos(np.trace(unitary).real/2)
        assert abs(gap-float(p)) < 1e-10


def test_original_source_execution_keeps_all_reference_outputs():
    value = m.pulse_checks()
    assert len(value['rows']) == 12
    for row in value['rows']:
        assert row['source_pulses'] == 2*row['subdivisions']
        probabilities = row['complete_system_reference_probability_intervals']
        assert len(probabilities) == 4
        assert sum(F(p[0]) for p in probabilities) <= 1 <= sum(F(p[1]) for p in probabilities)
    assert not value['new_hardware_fidelity_certificate']


@pytest.mark.parametrize('bad', (True, 0, 3, F(2), 2.0))
def test_only_declared_subdivisions_are_admitted(bad):
    with pytest.raises(ValueError):
        m.ratio_bounds(bad)


def test_retained_receipt_replays():
    m.check_receipt(m.HERE/'mode_energy_results.json', m.report())
