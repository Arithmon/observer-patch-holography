"""Exact local-completion ambiguity and complete source-discriminator checks."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import math

import mpmath as mp
import numpy as np
import pytest
import sympy as s

from . import dispersion_identifiability as d, kinetic


def test_same_principal_speed_and_every_retained_Nambu_operator():
    certificate = d.exact_certificate()
    assert certificate['local_lattice_CAR_modes'] == 19440
    assert certificate['unretained_local_lattice_CAR_modes'] == 18810
    assert not certificate['whole_lattice_quantum_outputs_or_source_resources_equal']
    assert not certificate['retained_continuation_or_stop_uses_exterior_records_or_failures']
    assert sorted(row['kinetic_Weyl_chirality'] for row in certificate['kinetic_Weyl_nodes']) == [-1]*4+[1]*4
    for raw in certificate['retained_momenta']:
        momentum = tuple(s.sympify(x) for x in raw)
        effective = d.effective_momentum(momentum, 'two_hop_lattice')
        assert effective == momentum
        # The full spin/internal/Nambu comparison is independent of the
        # symbolic scalar identity, and retains all three families.
        actual = kinetic.bdg_symbol(tuple(map(float, effective)))
        expected = kinetic.bdg_symbol(tuple(map(float, momentum)))
        np.testing.assert_array_equal(actual, expected)
    lower, upper = map(F, certificate['node_sign_intervals_exact']['B_over_c'])
    assert lower <= upper < 0
    assert all(F(x) > 0 for x in certificate['node_sign_intervals_exact']['minimum_positive_sine_factor_over_c'])


def test_second_shell_is_not_another_first_shell_axis_or_a_speed_refit():
    c, first, second = d.completion_coefficients()
    p = s.Symbol('p', real=True)
    assert s.simplify(s.diff(d.lattice_component(p), p).subs(p, 0)-c) == 0
    assert s.simplify(d.lattice_component(d.SECOND_SHELL)/c-(1-s.sqrt(3)/s.pi)) == 0
    assert d.effective_momentum((0, 0, d.SECOND_SHELL), 'two_hop_lattice') != (0, 0, d.SECOND_SHELL)
    assert first != 0 and second != 0
    for bad in ((0., 0, 1), (False, 0, 1)):
        with pytest.raises(ValueError, match='exact real'):
            d.effective_momentum(bad, 'continuum')


@pytest.mark.parametrize('law', ['continuum', 'two_hop_lattice'])
def test_direct_high_precision_pair_evolution_lies_in_exact_joint_event_interval(law):
    ctx = mp.mp.clone()
    ctx.dps = 80
    number = lambda x: ctx.mpf(str(s.N(x, 85)))
    omega, mass = number(d.action.source_coefficients()['Omega']), number(d._selected_mass())
    effective = d.effective_momentum((0, 0, d.SECOND_SHELL), law)[2]
    frequency, time = number(d.completion_coefficients()[0]*effective), ctx.mpf(1)/4
    unitary = ctx.expm(-1j*time*ctx.matrix([[0, mass], [mass, -2*frequency]]))
    first = unitary[:, 0]
    probability = abs(first[1])**2
    pair_zero = ctx.mpf(0)
    for bit in (0, 1):
        changed = ctx.matrix(first)
        if bit:
            changed[1] *= 1j
        after = unitary*changed
        expectation = 2*ctx.re(ctx.conj(after[0])*after[1])
        pair_zero += (probability if bit else 1-probability)*(1+expectation)/2
    score = (ctx.mpf(15)/16)**2*(1+ctx.mpf(2)/3*ctx.cos(2*omega*time))/2*pair_zero
    lo, hi = d.ideal_event_interval(law)
    exact = lambda x: ctx.mpf(x.numerator)/x.denominator
    assert exact(lo) <= score <= exact(hi)
    assert hi-lo < F(1, 10**17)


@pytest.fixture(scope='module')
def executed():
    return d.execute(steps=16)


def test_two_complete_source_laws_use_one_clock_detector_and_all_Fock_factors(executed):
    report, arrays = executed
    assert len(arrays) == 52
    rows = report['second_shell_source_comparison']
    assert rows['continuum']['same_reference_clock_calibration'] == rows['two_hop_lattice']['same_reference_clock_calibration']
    for law, row in rows.items():
        assert row['CAR_modes_in_executed_momentum_pair'] == 180
        assert len(row['complete_law']['terminal']) == 80
        assert abs(row['complete_law']['total_probability']-1) < 2e-7
        assert row['source_joint_event_difference_diagnostic'] < 1e-5
        assert row['source_field_all_input_operator_error_diagnostic'] < .01
        assert arrays[law+'_probe_outputs'].shape == (80, 256)
        assert arrays[law+'_clock_reference_columns'].shape == (80, 4, 4)
        assert {r['flag'] for r in row['complete_law']['terminal']} == {'completed', 'first_selector_timeout', 'second_selector_timeout'}
        assert {r['pair_X'] for r in row['complete_law']['terminal']} == {0, 1}
    gap = (rows['continuum']['complete_law']['scored_joint_event_probability']
           -rows['two_hop_lattice']['complete_law']['scored_joint_event_probability'])
    assert gap > .004
    assert all(a.dtype.kind == 'c' and np.isfinite(a).all() for a in arrays.values())
    assert report['scope']['full810_mode_native_source_error_budget_reexecuted'] is False


def test_exact_all_attempt_sampling_design_does_not_claim_a_physical_calibration(executed):
    report, _ = executed
    design = report['fixed_decision_design']
    assert design['attempted_independent_trials'] == 1200000
    assert F(design['sampling_certificate']['familywise_error_upper']) < F(1, 100)
    assert F(design['sampling_certificate']['width_exact']) <= 2*(F(1, 500)-F(1, 2000))
    assert F(design['Bernoulli_probability_decision_margin_exact']) == F(3, 2000)
    assert 'Y=2*event-1' in design['sampling_bounded_observable']
    assert F(report['ideal_joint_event_separation_lower_exact']) > F(1, 250)
    assert design['source_numeric_or_physical_service_error_certified_below_assumed_error'] is False
    assert report['scope']['laboratory_observations_or_empirical_fit_used'] is False


def test_pair_detector_metadata_cannot_poison_cached_source():
    first = d.pair_detector()
    expected = dict(first[3])
    first[3]['native_pulses'] = -10
    first[2][0]['targets'][0] = -11
    second = d.pair_detector()
    assert second[3] == expected and all(i >= 0 for event in second[2] for i in event['targets'])
    with pytest.raises(ValueError):
        second[0].flags.writeable = True


def _write(path, report, arrays):
    document = d.canonical_report(report)
    document['source_bindings'] = d.source_bindings()
    artifact = path.with_suffix('.npz')
    np.savez_compressed(artifact, **arrays)
    document['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    path.write_text(json.dumps(document, allow_nan=False))
    return document


def test_owned_replay_checks_complex_phases_and_fixed_admissions(tmp_path, executed, monkeypatch):
    report, arrays = executed
    path = tmp_path/'comparison.json'
    document = _write(path, report, arrays)
    calls = []
    def fresh(*, steps):
        calls.append(steps)
        return executed
    monkeypatch.setattr(d, 'execute', fresh)
    d.check(path)
    assert calls == [16]
    changed = dict(arrays)
    changed['continuum_factor1'] = changed['continuum_factor1'].copy()
    changed['continuum_factor1'][2] *= 1j
    document = _write(path, report, changed)
    with pytest.raises(ValueError, match='complex'):
        d.check(path)
    document = _write(path, report, arrays)
    document['neutral_pair_detector']['source_CAR_operator_error_admission'] = 0.
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match='admission'):
        d.check(path)


def test_actual_retained_receipt_replays():
    d.check(d.ROOT/'model/common_source/chiral_candidate/dispersion_identifiability_results.json')
