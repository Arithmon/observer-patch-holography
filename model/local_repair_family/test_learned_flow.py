from fractions import Fraction as F
from itertools import product
import math

import numpy as np
import pytest

from . import learned_flow as m


def _weights():
    # Twelve distinct, noncanonical values; no true-weight table is the input.
    return (F(1), F(999, 1000))+tuple(F(8, 25)*F(1000+i, 1000) for i in range(2, 12))


def _density():
    # Complex work/reference coherence spanning several spectator sectors.
    i, j = m.clock._pairs()[0]
    vector = np.zeros((66, 2), complex)
    vector[i, 0] = .5; vector[j, 1] = .5j
    vector[20, 0] += .5; vector[20, 1] += .5j
    flat = vector.reshape(-1)
    return np.outer(flat, flat.conj()), vector.T@vector.conj()


def _distance(a, b):
    assert a.keys() == b.keys()
    return sum(np.linalg.norm(a[k]-b[k], 'nuc') for k in a)/2


def test_exact_inferred_covariance_includes_mixed_tags_and_changed_supports():
    weights = _weights(); states = m.law.basis()
    masses = [math.prod(w**n for w, n in zip(weights, x, strict=True)) for x in states]
    cert = m.exact_certificate(weights)
    assert cert['complete_native_columns_checked'] == 7920
    assert not cert['stationarity_of_misspecified_reference_claimed']
    changed = 0
    for name1, name2 in product(('C', 'E'), repeat=2):
        for coin1, label1, pairs1 in m.law.raw_maps(name1, (0, 5)):
            for coin2, label2, pairs2 in m.law.raw_maps(name2, (5, 9)):
                mapping = dict(pairs2)
                current1 = m.law.records.record_current(name1, (0, 5, coin1, label1))
                current2 = m.law.records.record_current(name2, (5, 9, coin2, label2))
                ratio = (weights[5]/weights[0])**current1*(weights[9]/weights[5])**current2
                for i, j in pairs1:
                    if j in mapping:
                        assert masses[mapping[j]] == ratio*masses[i]
                        changed += int(i != mapping[j])
    assert changed > 100


def test_only_supplied_weights_determine_every_source_and_padding_phase(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('hidden true weights were read')
    monkeypatch.setattr(m.clock, 'parameters', forbidden)
    monkeypatch.setattr(m.law, 'reference', forbidden)
    weights = _weights(); program = m.source_program(weights)
    assert program == m.source_program(tuple(13*x for x in weights))
    diagonal = m.ratios(weights); codes = m.common.full.support()[0]
    for address, state in zip(codes, m.law.basis(), strict=True):
        assert diagonal[address] == math.prod(w**n for w, n in zip(weights, state, strict=True))
    assert all(diagonal[i] == 1 for i in set(range(512))-set(codes))
    changed = list(weights); changed[10] *= F(1001, 1000)
    assert m.source_program(tuple(changed))['clock'] != program['clock']
    assert m.source_program(tuple(changed))['probe'] == program['probe']


@pytest.mark.parametrize('bad', [None, [1]*12, (1,)*11, (True,)+(1,)*11,
                                (0,)+(1,)*11, (-1,)+(1,)*11,
                                (1.,)+(1,)*11, (complex(1),)+(1,)*11])
def test_original_positive_rational_domain(bad):
    with pytest.raises(ValueError):
        m.source_program(bad)


@pytest.mark.parametrize('slots', [(True, 1), (1, True), (2, 2), (-1, -1), (1, 0)])
def test_malformed_fraction_slots_rejected_before_cache_alias(slots):
    bad = F(1)
    bad._numerator, bad._denominator = slots
    with pytest.raises(ValueError):
        m.source_arrays((bad,)+(F(1),)*11)


def test_returned_ratios_and_programs_do_not_own_input_or_cache_objects():
    weights = _weights(); original = m.source_program(weights)
    exposed = m.ratios(weights)
    address = m.common.full.support()[0][0]
    exposed[address]._numerator = 7
    original['clock'][0]['angle'] += .7
    assert m.source_program(weights)['clock'][0]['angle'] != original['clock'][0]['angle']
    assert m.ratios(weights)[address] != exposed[address]
    before = m.ratios(weights)
    weights[1]._numerator = 997
    assert m.ratios(_weights()) == before
    assert m.ratios(weights) != before


def test_dyadic_phase_enclosure_against_independent_high_precision():
    import mpmath as mp
    ratios = (F(1), F(25, 8), F(8, 25),
              F(1048573, 262147)**6, F(262147, 1048573)**6,
              F(2**4096+1, 2**4096-1), F(1, 2**1024), F(2**1024))
    with mp.workdps(160):
        for ratio, parameter in product(ratios, (F(1, 2), F(1))):
            re, im, radius = m._phase_enclosure(parameter, ratio)
            def exact(x):
                return mp.mpf(x.numerator)/x.denominator
            target = mp.exp(-1j*exact(parameter)*mp.log(exact(ratio)))
            error = abs(target-(exact(re)+1j*exact(im)))
            assert error <= exact(radius)+mp.mpf('1e-150')
            assert 0 <= radius < F(1, 10**25)
            assert re.denominator <= 2**128 and im.denominator <= 2**128
    # Publicly returned numerical enclosures cannot modify the cache payload.
    re._numerator = 123
    assert m._phase_enclosure(F(1), F(2**1024))[0] != re


def test_mutable_inherited_log_cache_cannot_define_the_proof_target():
    ratio = F(19, 7); weights = (F(1), ratio)+(F(1),)*10
    endpoints = m.common.logarithm_interval(ratio)
    original = [(x, x.numerator, x.denominator) for x in endpoints]
    try:
        for x in endpoints:
            x._numerator, x._denominator = 0, 1
        # The inherited compiler would produce an identity phase. The local
        # proof authority must independently reject it, including cold caches.
        assert m.certify_phase(F(1), weights, np.eye(512, dtype=complex)) > 1
    finally:
        for x, n, d in original:
            x._numerator, x._denominator = n, d


def test_mutable_inherited_sqrt_cache_cannot_erase_a_wrong_unitary():
    endpoints = m.common.full.exact._sqrt_interval(F(4))
    original = [(x, x.numerator, x.denominator) for x in endpoints]
    wrong = np.eye(512, dtype=complex)
    wrong[m.common.full.support()[0][0], m.common.full.support()[0][0]] = -1
    try:
        for x in endpoints:
            x._numerator, x._denominator = 0, 1
        assert m.certify_phase(F(1), (F(1),)*12, wrong) >= 2
        pulse = np.eye(6, dtype=complex); pulse[5, 0] = 2
        assert m._processor_bounds(pulse, 2)[0] == 2
    finally:
        for x, n, d in original:
            x._numerator, x._denominator = n, d


@pytest.fixture(scope='module')
def source_report():
    return m.report(_weights())


def test_full512_executed_words_against_independent_diagonal(source_report, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('compiler consulted the apparatus reference')
    monkeypatch.setattr(m.clock, 'parameters', forbidden)
    monkeypatch.setattr(m.law, 'reference', forbidden)
    arrays = m.source_arrays(_weights())
    assert m.source_budget(_weights()) == source_report['source_budget']
    diagonal = np.ones(512, complex)
    weights = _weights()
    for address, state in zip(m.common.full.support()[0], m.law.basis(), strict=True):
        energy = sum(n*(math.log(w.numerator)-math.log(w.denominator))
                     for w, n in zip(weights, state, strict=True))
        diagonal[address] = np.exp(-1j*energy)
    assert np.linalg.norm(arrays['clock']-np.diag(diagonal), 'fro') < 1e-9
    assert np.linalg.norm(arrays['half_clock']@arrays['half_clock']-arrays['clock'], 'fro') < 2e-9
    arrays['clock'][:] = 0
    assert np.linalg.norm(m.source_arrays(_weights())['clock']) > 20
    budget = source_report['source_budget']
    assert budget['maximum_source_pulses'] == 24876
    assert budget['coherent_load_export_calls'] == 49752
    assert F(budget['complete_instrument_halfdiamond_upper']) < F(1, 10**7)
    assert source_report['source_sha256'] == m.source_sha256()


@pytest.mark.parametrize('name', ('C', 'E'))
def test_complete_source_history_with_untouched_complex_reference(source_report, name):
    rho, reference = _density(); weights = _weights()
    actual = m.execute(rho, weights, name)
    endpoint = m.execute(rho, weights, name, order='endpoint')
    ideal = m.execute(rho, weights, name, raw=False)
    allowance = float(F(source_report['source_budget']['complete_instrument_halfdiamond_upper']))
    assert _distance(actual, ideal) < allowance
    assert _distance(actual, endpoint) < 2*allowance
    marginal = np.zeros((2, 2), complex)
    for key, value in actual.items():
        marginal += value if key[0] == 'abort' else np.einsum('aras->rs', value.reshape(66, 2, 66, 2))
    assert np.linalg.norm(marginal-reference, 'nuc') < 2*allowance
    assert not np.any(m.raw_output(actual, name, ((0, 3),)*3))
    live = next(key for key in actual if key[0] == 'read')
    slots = tuple(row[:2] for row in (*live[1], live[2]))
    assert np.array_equal(m.raw_output(actual, name, slots), actual[live])
    del actual[live]
    with pytest.raises(KeyError):
        m.raw_output(actual, name, slots)


def test_inferred_phase_has_a_nonzero_complete_readout_without_branch_selection():
    rho = np.diag(list(map(float, m.law.reference(m.clock.parameters())))).astype(complex)
    for name in ('C', 'E'):
        phased = m.execute(rho, _weights(), name, raw=False)
        zero = m.execute(rho, _weights(), name, raw=False, order='zero')
        def mean(outputs):
            return sum(key[2][3]*float(np.trace(value).real)
                       for key, value in outputs.items() if key[0] == 'read')
        assert abs(mean(phased)-mean(zero)) > .03
        assert abs(sum(np.trace(x).real for x in phased.values())-1) < 1e-13


def test_complete_report_rejects_a_disconnected_or_nonfinite_clock(source_report, monkeypatch):
    arrays = m.source_arrays(_weights())
    arrays['clock'] = arrays['half_clock'] = np.eye(512, dtype=complex)
    monkeypatch.setattr(m, 'source_arrays', lambda weights: arrays)
    with pytest.raises(ArithmeticError, match='complete learned source'):
        m.report(_weights())


def test_report_refuses_nan_before_comparing_residuals(source_report, monkeypatch):
    original = m.execute
    def changed(*args, **kwargs):
        outputs = original(*args, **kwargs)
        outputs['abort', 'probe'][0, 0] = float('nan')
        return outputs
    monkeypatch.setattr(m, 'execute', changed)
    with pytest.raises(ArithmeticError, match='nonfinite'):
        m.report(_weights())


def test_literal_program_mutation_cannot_reuse_warm_target_cache(source_report, monkeypatch):
    original = m.source_program
    def changed(weights):
        program = original(weights)
        # A literal I-word with a wrong global phase is unitary on every column
        # but disagrees with the inferred source and declared padding target.
        program['clock'] = [m.common.full.compiler.pulse('I', .1, [0])]
        return program
    monkeypatch.setattr(m, 'source_program', changed)
    with pytest.raises(ArithmeticError, match='literal learned full source phase'):
        m.source_arrays(_weights())


@pytest.mark.parametrize('damage', ('boolean_target', 'nan', 'integer_angle', 'word_subclass'))
def test_literal_validation_precedes_cache_lookup(monkeypatch, damage):
    program = m.source_program(_weights())
    if damage == 'boolean_target': program['clock'][0]['targets'][0] = False
    elif damage == 'nan': program['clock'][0]['angle'] = float('nan')
    elif damage == 'integer_angle': program['clock'][0]['angle'] = 0
    else:
        class Word(list):
            pass
        program['clock'] = Word(program['clock'])
    monkeypatch.setattr(m, 'source_program', lambda weights: program)
    with pytest.raises(ValueError):
        m.source_arrays(_weights())


def test_padding_escape_preserves_reference_and_invalid_raw_tails_refuse():
    rho, reference = _density()
    identity = np.eye(512, dtype=complex); codes = m.common.full.support()[0]
    outside = next(i for i in range(512) if i not in codes)
    occupied = int(np.flatnonzero(np.diag(rho).real)[0])//2
    identity[[codes[occupied], outside]] = identity[[outside, codes[occupied]]]
    retained, abort = m.clock._export(identity, rho)
    marginal = np.einsum('aras->rs', retained.reshape(66, 2, 66, 2))+abort
    assert np.allclose(marginal, reference, atol=1e-15)
    assert np.trace(abort).real > 0
    with pytest.raises(ValueError):
        m.raw_output({('abort', 'probe'): np.zeros((2, 2))}, 'C', ((0, 3), (0, 0), (False, 0)))


def test_physical_errors_are_conditional_and_denominator_guarded(source_report):
    zero = m.physical_allowance(_weights())
    assert zero['halfdiamond_upper'] == '0'
    small = m.physical_allowance(_weights(), F(1, 10**12), F(1, 10**13), F(1, 10**10))
    total = F(24876, 10**12)+F(49752, 10**13)+F(3, 10**10)
    assert F(small['operator_upper']) == total/(1-total)
    with pytest.raises(ValueError, match='below one'):
        m.physical_allowance(_weights(), F(1, 24876))
