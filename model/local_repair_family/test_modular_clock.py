from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import modular_clock as m


def _density():
    return np.diag(list(map(float, m.law.reference(m.parameters())))).astype(complex)


def _mean(outputs):
    return sum(key[2][3]*float(np.trace(value).real)
               for key, value in outputs.items() if key[0] == 'read')


def _distance(a, b):
    assert a.keys() == b.keys()
    return sum(np.linalg.norm(a[k]-b[k], 'nuc') for k in a)/2


def _independent_channel(density, parameters):
    """Endpoint oracle: no producer pair-transfer kernel or raw_maps."""
    states = m.law.native.states(binary=True); index = {x: i for i, x in enumerate(states)}
    result = float(1-parameters.repair_probability)*density
    for name, mixture in (('C', parameters.c_weight), ('E', 1-parameters.c_weight)):
        for edge, weight in zip(m.law.SEAMS, parameters.seam_weights, strict=True):
            a, b = edge
            for coin in (0, 1):
                groups = {}
                for i, state in enumerate(states):
                    target, rec = m.law.records.recorded_step(name, state, edge, coin)
                    groups.setdefault(rec, []).append((i, index[target]))
                probability = (parameters.repair_probability*mixture*F(weight, sum(parameters.seam_weights))
                               *F(parameters.port_weights[a if coin else b],
                                  parameters.port_weights[a]+parameters.port_weights[b]))
                for pairs in groups.values():
                    source, target = map(list, zip(*pairs))
                    result[np.ix_(target, target)] += float(probability)*density[np.ix_(source, source)]
    return result


def test_exact_regional_and_record_covariance_on_original_matrices():
    p = m.parameters(); sigma = m.law.reference(p)
    cert = m.exact_certificate(p)
    assert cert['complete_native_columns_checked'] == 7920
    assert cert['all_66_bare_two_port_algebras_invariant']
    # Exact mixed-tag histories, with two different seam orientations and
    # actual support changes. Ratios multiply without commuting K with sigma.
    tested = moving = 0
    for name1, name2 in product(('C', 'E'), repeat=2):
        for c1, l1, map1 in m.law.raw_maps(name1, (0, 5)):
            for c2, l2, map2 in m.law.raw_maps(name2, (1, 8)):
                next_map = dict(map2)
                composed = [(i, next_map[j]) for i, j in map1 if j in next_map]
                if not composed:
                    continue
                current1 = m.law.records.record_current(name1, (0, 5, c1, l1))
                current2 = m.law.records.record_current(name2, (1, 8, c2, l2))
                ratio = F(8, 25)**(current1+current2)
                assert all(sigma[j] == ratio*sigma[i] for i, j in composed)
                tested += 1; moving += any(i != j for i, j in composed)
    assert tested > 20 and moving > 10


@pytest.mark.parametrize('mixture', (F(), F(1, 2), F(1)))
def test_full_scheduler_response_against_independent_density_oracle(mixture):
    p = m.law.Parameters(m.parameters().port_weights, c_weight=mixture)
    oracle = m._ideal_arrays(); codes = np.asarray(m.common.full.support()[0])
    R = oracle['probe'][np.ix_(codes, codes)]
    rho = R@_density()@R.conj().T
    rho = _independent_channel(_independent_channel(rho, p), p)
    exact = m.scheduler_response(p)
    assert F(exact['cosine']) >= F(80801, 1501000)
    assert F(exact['constant']) == -F(69979084801, 647297244000)
    for s in (0., .5, 1.):
        U = np.diag([complex(float(x))**(-1j*s) for x in m.law.reference(p)])
        out = R.conj().T@U@rho@U.conj().T@R
        # Expected signed placement current on the actual biased seam.
        a, b = m.PROBE; expected = 0.
        for i, state in enumerate(m.law.basis()):
            if state[a]+state[b] == 1:
                score = F(p.port_weights[b], p.port_weights[a]+p.port_weights[b]) if state[a] else -F(p.port_weights[a], p.port_weights[a]+p.port_weights[b])
                expected += float(score)*out[i, i].real
        predicted = float(F(exact['constant']))+float(F(exact['cosine']))*math.cos(s*math.log(25/8))
        assert abs(expected-predicted) < 2e-14


def test_full_scheduler_holds_and_nonnative_weights_are_not_omitted():
    p = m.law.Parameters(m.parameters().port_weights, (1, 2)*15, F(3, 7), F(1, 16))
    response = m.scheduler_response(p)
    assert F(response['cosine']) > F(m.scheduler_response()['cosine'])
    assert m.scheduler_frame_error(p) < m.NATIVE_GATE
    assert m.scheduler_response(m.law.Parameters())['cosine'] == '0'
    assert m.scheduler_response(attempts=0)['cosine'] == '1649/15010'


@pytest.fixture(scope='module')
def source_report():
    return m.report()


def test_actual_source_report_has_both_complete_laws_and_full_scheduler(source_report):
    assert set(source_report['requested_complete_source_execution']) == {'C', 'E'}
    assert set(source_report['full_scheduler_response']) == {'C', 'E', 'mixed'}
    for value in source_report['full_scheduler_response'].values():
        assert F(value['mean_current_s1_minus_s0_interval'][1]) < -F(3, 100)
        assert value['inverse_mean_Lipschitz_upper'] == '32'
    assert F(source_report['source_budget']['complete_instrument_halfdiamond_upper']) < F(1, 10**7)
    assert source_report['source_budget']['hardware_and_physical_calibration_are_external_obligations']
    assert source_report['source_budget']['maximum_source_pulses'] == 24876
    assert source_report['source_sha256'] == m.source_sha256()


@pytest.mark.parametrize('name', ('C', 'E'))
def test_complete_source_covariance_with_complex_entangled_reference(source_report, name):
    i, j = m._pairs()[0]
    vector = np.zeros((66, 2), complex)
    vector[i, 0] = 1/2; vector[j, 1] = 1j/2
    vector[20, 0] += .5; vector[20, 1] += .5j
    flat = vector.reshape(-1); rho = np.outer(flat, flat.conj())
    before = m.execute(rho, name, order='endpoint')
    between = m.execute(rho, name)
    ideal = m.execute(rho, name, raw=False)
    budget = float(F(source_report['source_budget']['complete_instrument_halfdiamond_upper']))
    assert _distance(between, before) <= 2*budget
    assert _distance(between, ideal) <= budget
    marginal = np.zeros((2, 2), complex)
    for key, value in between.items():
        marginal += value if key[0] == 'abort' else np.einsum('aras->rs', value.reshape(66, 2, 66, 2))
    assert np.linalg.norm(marginal-vector.T@vector.conj()) < 2*budget
    raw_zero = ((0, 3), (0, 3), (0, 3))
    assert not np.any(m.raw_output(between, name, raw_zero))
    read_key = next(k for k in between if k[0] == 'read')
    raw_live = tuple(x[:2] for x in (*read_key[1], read_key[2]))
    assert np.array_equal(m.raw_output(between, name, raw_live), between[read_key])
    del between[read_key]
    with pytest.raises(KeyError):
        m.raw_output(between, name, raw_live)


def test_overlapping_E_repair_kills_visibility_without_breaking_covariance():
    sigma = _density(); seams = ((1, 8), (0, 5))
    phased = m.execute(sigma, 'E', seams=seams, raw=False)
    unphased = m.execute(sigma, 'E', seams=seams, raw=False, order='zero')
    assert _distance(phased, unphased) < 1e-14
    ordinary = m.execute(sigma, 'E', raw=False)
    zero = m.execute(sigma, 'E', raw=False, order='zero')
    assert abs(_mean(ordinary)-_mean(zero)) > .06
    assert any(k[0] == 'read' and any(r[3] for r in k[1]) and np.trace(v).real > 0
               for k, v in ordinary.items())


def test_padding_export_retains_complex_reference_against_explicit_contraction():
    from ..native_modular_reference import history_clock
    codes = m.common.full.support()[0]; outside = next(i for i in range(512) if i not in codes)
    U = np.eye(512, dtype=complex); a = codes[0]
    U[np.ix_((a, outside), (a, outside))] = [[.6, .8j], [.8j, .6]]
    vector = np.zeros((66, 2), complex); vector[0] = [1/math.sqrt(2), 1j/math.sqrt(2)]
    rho = np.outer(vector.reshape(-1), vector.reshape(-1).conj())
    fast = m._export(U, rho); explicit = history_clock._export(U, rho)
    assert all(np.allclose(x, y, atol=1e-15, rtol=0) for x, y in zip(fast, explicit))
    assert abs(fast[1][0, 1].imag) > .3
    assert abs(np.trace(fast[1]).real-.64) < 1e-15


def test_detached_outputs_and_failed_original_inputs(source_report):
    arrays = m.source_arrays(); arrays['clock'][:] = 0
    assert np.linalg.norm(m.source_arrays()['clock']) > 20
    branches = m.native_branches('C', (0, 5)); p = branches[0][3]
    object.__setattr__(p, '_numerator', 999)
    assert m.native_branches('C', (0, 5))[0][3] == F(8, 33)
    prepared = m.parameters()
    object.__setattr__(prepared.c_weight, '_denominator', 3)
    object.__setattr__(prepared.repair_probability, '_denominator', 2)
    assert m.parameters().c_weight == F(1, 2)
    assert m.parameters().repair_probability == F(1)
    for bad in (1, 2, -1):
        with pytest.raises(ValueError):
            m.execute(_density(), 'C', raw=bad)
    rho = _density(); rho[0, 0] = -1e-100
    with pytest.raises(ValueError):
        m.execute(rho, 'C', raw=False)
    tiny = np.zeros((66, 66), complex); tiny[0, 0] = 1e-300
    assert m._validate_input(tiny)[0, 0] == 1e-300
    with pytest.raises(ValueError):
        m.raw_output({}, 'C', ((False, 0), (0, 0), (0, 0)))


def test_report_refuses_nonfinite_or_disconnected_clock(source_report, monkeypatch):
    original = m.scheduler_execute
    monkeypatch.setattr(m, 'scheduler_execute', lambda *a, **k: {'mean_current': float('nan'), 'complete_mass': 1., 'export_abort_mass': 0.})
    with pytest.raises(ArithmeticError, match='readout'):
        m.report()
    monkeypatch.setattr(m, 'scheduler_execute', original)
    source = m.source_arrays()
    source['clock'] = source['half_clock'] = np.eye(512, dtype=complex)
    monkeypatch.setattr(m, 'source_arrays', lambda: source)
    with pytest.raises(ArithmeticError, match='readout'):
        m.report()


def test_full_source_certificate_checks_unoccupied_padding(source_report):
    array = m.source_arrays()['clock']; codes = m.common.full.support()[0]
    outside = next(i for i in range(512) if i not in codes)
    array[outside, outside] *= -1
    assert m.common.certify_diagonal(F(1), m._ratios(), array) > m.WORD_GATE


def test_physical_allowance_is_conditional_exact_operator_composition():
    assert m.physical_allowance()['halfdiamond_upper'] == '0'
    result = m.physical_allowance(F(1, 10**10), F(1, 10**11), F(1, 10**9))
    total = 24876*F(1, 10**10)+49752*F(1, 10**11)+3*F(1, 10**9)
    delta = total/(1-total)
    assert F(result['operator_upper']) == delta
    assert F(result['halfdiamond_upper']) == delta*(2+delta)/2
    with pytest.raises(ValueError):
        m.physical_allowance(F(1, 100))
    damaged = F(1, 2); object.__setattr__(damaged, '_denominator', -2)
    with pytest.raises(ValueError):
        m.physical_allowance(damaged)


def test_first_failure_gain_bound_requires_its_denominator_domain():
    with pytest.raises(ArithmeticError, match='first-failure gain domain'):
        m._processor_failure_bound(10**11)


def test_exported_literal_program_controls_cached_execution_and_budget(source_report, monkeypatch):
    original = m.source_program
    before = original()
    expected_menu = {(p['kind'], p['angle'], len(p['targets']))
                     for word in before.values() for p in word}
    assert source_report['source_budget']['actual_unique_processor_pulses_checked'] == len(expected_menu)
    # Mutating a returned object cannot poison any future program or array.
    damaged = original()
    damaged['clock'][0]['targets'][0] = (damaged['clock'][0]['targets'][0]+1) % 9
    assert original() == before

    def wrong_program():
        result = original()
        for pulse in result['clock']:
            pulse['targets'] = [(bit+1) % 9 for bit in pulse['targets']]
        return result

    wrong = wrong_program()
    assert {name: len(word) for name, word in wrong.items()} == {name: len(word) for name, word in before.items()}
    assert {(p['kind'], p['angle'], len(p['targets']))
            for word in wrong.values() for p in word} == expected_menu
    monkeypatch.setattr(m, 'source_program', wrong_program)
    # Old arrays are already cached. Same counts and pulse menu must not reuse
    # them: the actual wrong word has to execute and fail its target check.
    with pytest.raises(ArithmeticError, match='literal full source phase'):
        m.source_budget()
    monkeypatch.setattr(m, 'source_program', original)
    assert m.source_budget() == source_report['source_budget']


@pytest.mark.parametrize('damage', ('boolean_target', 'nonfinite_angle', 'integer_angle', 'pulse_subclass', 'word_subclass'))
def test_literal_program_rejects_noncanonical_objects_before_cache_lookup(monkeypatch, damage):
    program = m.source_program()
    class Pulse(dict):
        pass
    class Word(list):
        pass
    if damage == 'boolean_target': program['clock'][0]['targets'][0] = False
    elif damage == 'nonfinite_angle': program['clock'][0]['angle'] = float('nan')
    elif damage == 'integer_angle': program['clock'][0]['angle'] = 0
    elif damage == 'pulse_subclass': program['clock'][0] = Pulse(program['clock'][0])
    elif damage == 'word_subclass': program['clock'] = Word(program['clock'])
    monkeypatch.setattr(m, 'source_program', lambda: program)
    with pytest.raises(ValueError):
        m.source_arrays()


@pytest.mark.parametrize('name', ('C', 'E'))
@pytest.mark.parametrize('slots', (((0, 3), ('bad',), (False, 0)), ((0, 0), (0, 3), ('bad',))))
def test_zero_raw_prefix_never_hides_malformed_later_slots(name, slots):
    with pytest.raises(ValueError, match='integer coin/raw'):
        m.raw_output({('abort', 'probe'): np.zeros((1, 1))}, name, slots)


def test_clock_order_requires_original_string_before_source_execution():
    class Order:
        def __eq__(self, other):
            return other == 'zero'
    for order in (Order(), np.str_('zero')):
        with pytest.raises(ValueError, match='clock order'):
            m.execute(_density(), 'C', raw=False, order=order)
        with pytest.raises(ValueError, match='clock order'):
            m.scheduler_execute(raw=False, order=order)


def test_processor_bounds_resolve_tiny_original_components_before_norms():
    array = np.eye(6, dtype=complex)
    array[2, 0] = 1j*2.**-540
    array[0, 1] = 2.**-540
    leakage, gain = m._processor_bounds(array, 2)
    # Squaring either original component in binary64 erases it. Its exact
    # dyadic square, and the off-diagonal Gram defect, remain nonzero here.
    assert float(abs(array[2, 0]))**2 == 0
    assert leakage == F(1, 2**540)
    assert gain > F(1, 2**540)
    assert m._processor_bounds(np.eye(6, dtype=complex), 4) == (F(), F())


@pytest.mark.parametrize('damage', ('complement', 'gain'))
def test_complete_budget_rejects_small_processor_failures_or_gain(source_report, monkeypatch, damage):
    original = m.common.full.compiler.native_pulse
    def wrong_pulse(kind, angle):
        array = original(kind, angle)
        if kind == 'I':
            if damage == 'complement': array[5, 0] += 2.**-34
            else: array[:2, :2] *= 1+2.**-34
        return array
    monkeypatch.setattr(m.common.full.compiler, 'native_pulse', wrong_pulse)
    # Warm word caches must not bypass the fresh literal pulse-menu check.
    with pytest.raises(ArithmeticError, match='processor complement or gain'):
        m.source_budget()


def test_processor_certificate_rejects_array_subclass_before_overridden_copy():
    class Hidden(np.ndarray):
        def copy(self):
            return np.eye(6, dtype=complex)
    array = np.eye(6, dtype=complex).view(Hidden)
    array[5, 0] = 1
    with pytest.raises(ValueError, match='original six-level pulse'):
        m._processor_bounds(array, 2)


@pytest.mark.parametrize('method', (m.scheduler_frame_error, m.detector_frame_error))
def test_public_frame_error_cannot_poison_cached_outward_root(method):
    result = method(m.parameters())
    expected = F(result.numerator, result.denominator)
    assert expected > 0
    object.__setattr__(result, '_numerator', -1)
    repeated = method(m.parameters())
    assert repeated == expected and repeated is not result


def test_processor_bounds_detach_even_exact_zero_root_cache_entries():
    array = np.eye(6, dtype=complex)
    result = m._processor_bounds(array, 2)
    assert result == (F(), F())
    object.__setattr__(result[0], '_numerator', -1)
    object.__setattr__(result[1], '_denominator', -1)
    assert m._processor_bounds(array, 2) == (F(), F())
