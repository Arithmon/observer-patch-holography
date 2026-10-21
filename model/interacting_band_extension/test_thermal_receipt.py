"""Receipt checking never trusts an uploaded bound or drops unexpected fields."""
from collections import Counter
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from . import thermal_certificate as certificate


@pytest.fixture
def reference(monkeypatch):
    # This tests the CLI acceptance boundary. Mathematical tests independently
    # exercise the producer; model.verify replays the full expensive receipt.
    result=json.loads(Path(__file__).with_name('thermal_results.json').read_text())
    monkeypatch.setattr(certificate,'run',lambda: copy.deepcopy(result))
    return result


@pytest.mark.parametrize('mutation',[
    lambda r:r['cold_protocol']['histories'].pop(),
    lambda r:r['cold_protocol']['histories'][0]['record_alphabet'].pop(),
    lambda r:r['cold_protocol']['histories'][0].update(retained_reference_dimension=True),
    lambda r:r.update(complete_source_to_infinite_output_error_upper_under_numerical_guards='0'),
    lambda r:r['boundary'].update(physical_seconds_or_universal_empirical_agreement_established=True),
    lambda r:r['source_sha256'].clear(),
    lambda r:r.update(extra_claim='all premises proved'),
    lambda r:r['source_execution']['histories'][0]['probabilities'].__setitem__(0,float('nan')),
])
def test_semantic_mutations_fail_closed(reference,mutation,tmp_path,monkeypatch):
    bad=copy.deepcopy(reference);mutation(bad)
    path=tmp_path/'bad.json';path.write_text(json.dumps(bad))
    monkeypatch.setattr(sys,'argv',['thermal','--check',str(path)])
    with pytest.raises((ValueError,ArithmeticError)):
        certificate.main()


def test_duplicate_fields_and_overwrite_fail(reference,tmp_path,monkeypatch):
    path=tmp_path/'bad.json';path.write_text('{"schema":"one","schema":"two"}')
    monkeypatch.setattr(sys,'argv',['thermal','--check',str(path)])
    with pytest.raises(ValueError):certificate.main()
    monkeypatch.setattr(sys,'argv',['thermal','--write',str(path)])
    with pytest.raises(FileExistsError):certificate.main()


def test_outward_display_does_not_round_down():
    from fractions import Fraction as F
    for value in (F(1,3),F(-1,3),F(1,10**30),F(17,8)):
        assert F(certificate.upper(value))>=value
    with pytest.raises(ValueError):certificate.upper(True)


def test_receipt_binds_inherited_controller_clock_and_source_inputs():
    from .construction import SOURCE_PATHS
    assert set(SOURCE_PATHS) <= set(certificate.SOURCE_PATHS)
    assert 'model/interacting_band_extension/thermal_controller.py' in certificate.SOURCE_PATHS
    assert len(certificate.SOURCE_PATHS) == len(set(certificate.SOURCE_PATHS))


@pytest.mark.parametrize('value', [0, 1, F(1, 3), F(1)])
def test_certified_upper_accepts_only_values_inside_exact_ceiling(value):
    assert certificate.certified_upper(value, F(1), 'test bound') == '1'


@pytest.mark.parametrize('helper', ['certified_upper', 'certified_lower'])
@pytest.mark.parametrize('bad', [
    None, True, False, 0.0, float('nan'), float('inf'), float('-inf'),
    '1/2', complex(1), np.float64(1), np.int64(1), np.bool_(True), F(-1),
])
def test_certified_bounds_reject_inexact_or_invalid_values(helper, bad):
    with pytest.raises(ValueError):
        getattr(certificate, helper)(bad, F(1), 'test bound')


@pytest.mark.parametrize('helper', ['certified_upper', 'certified_lower'])
@pytest.mark.parametrize('bad', [
    None, True, False, 1.0, float('nan'), float('inf'), float('-inf'),
    '1', np.float64(1), np.int64(1), np.bool_(True), 0, F(0), F(-1),
])
def test_certified_bounds_require_exact_positive_declared_limits(helper, bad):
    with pytest.raises(ValueError):
        getattr(certificate, helper)(F(1), bad, 'test bound')


def test_certified_bounds_preserve_exact_boundaries_without_float_tolerance():
    boundary = F(1, 10**9)
    tiny = F(1, 10**100)
    assert certificate.certified_upper(boundary, boundary, 'upper') == str(boundary)
    assert certificate.certified_lower(boundary, boundary, 'lower') == str(boundary)
    assert certificate.certified_upper(boundary-tiny, boundary, 'upper') == str(boundary)
    assert certificate.certified_lower(boundary+tiny, boundary, 'lower') == str(boundary)
    with pytest.raises(ArithmeticError):
        certificate.certified_upper(boundary+tiny, boundary, 'upper')
    with pytest.raises(ArithmeticError):
        certificate.certified_lower(boundary-tiny, boundary, 'lower')
    with pytest.raises(ArithmeticError):
        certificate.certified_lower(0, boundary, 'lower')


def test_seed_certificate_tolerates_eigensolver_rounding_without_changing_declared_bound(monkeypatch):
    # Exact dyadic residuals certify both candidates; an eigensolver's last-bit
    # choice is not a premise of the mathematical certificate or receipt.
    baseline, baseline_budget = certificate.finite_gibbs_seed()
    eigh = np.linalg.eigh

    def one_ulp_candidate(matrix):
        values, vectors = eigh(matrix)
        return np.nextafter(values, np.inf), np.nextafter(vectors, np.inf)

    monkeypatch.setattr(np.linalg, 'eigh', one_ulp_candidate)
    shifted, shifted_budget = certificate.finite_gibbs_seed()
    ceiling = certificate.CANONICAL_SEED_ERROR_CEILING
    for budget in (baseline_budget, shifted_budget):
        assert certificate.certified_upper(
            budget['canonical_seed_vector_error_upper'], ceiling, 'seed') == str(ceiling)
    # The independently certified errors bound their common exact seed.
    assert np.linalg.norm(baseline-shifted, 'fro') <= float(
        baseline_budget['canonical_seed_vector_error_upper']
        + shifted_budget['canonical_seed_vector_error_upper'])


def test_seed_certificate_rejects_candidate_whose_certified_error_exceeds_ceiling(monkeypatch):
    eigh = np.linalg.eigh

    def degraded_candidate(matrix):
        values, vectors = eigh(matrix)
        vectors[0, 0] += 1e-8
        return values, vectors

    monkeypatch.setattr(np.linalg, 'eigh', degraded_candidate)
    _, budget = certificate.finite_gibbs_seed()
    with pytest.raises(ArithmeticError):
        certificate.certified_upper(budget['canonical_seed_vector_error_upper'],
                                    certificate.CANONICAL_SEED_ERROR_CEILING, 'seed')


def test_run_rejects_excess_seed_allowance_before_composing_history(monkeypatch):
    # Bypass costly, unrelated preparation inputs, but exercise run's actual
    # seed check. A threshold comparison after float conversion would miss this.
    def preparation(occupation, shell, beta=F(1, 4)):
        lower, upper, error = certificate.BRACKET_ENVELOPES[
            'hot' if occupation == 24 else 'cold']
        return {'occupation': occupation, 'lower_form_shell': shell,
                'inverse_temperature': beta, 'coupling': F(1),
                'infinite_partition_lower': lower, 'infinite_partition_upper': upper,
                'canonical_preparation_trace_distance_upper': error,
                'harmonic_complement_upper': F(0)}

    monkeypatch.setattr(certificate, 'preparation_bracket', preparation)
    monkeypatch.setattr(certificate, 'certify_spectrum',
                        lambda *args, **kwargs: {'partition_function_interval': (F(1), F(1))})
    monkeypatch.setattr(certificate, 'finite_gibbs_seed', lambda: (
        np.zeros((64, 64), complex), {'canonical_seed_vector_error_upper':
            certificate.CANONICAL_SEED_ERROR_CEILING+F(1, 10**100)}))

    def unexpected_history(*args, **kwargs):
        pytest.fail('an uncertified seed reached history composition')

    monkeypatch.setattr(certificate, 'history', unexpected_history)
    with pytest.raises(ArithmeticError):
        certificate.run()


@pytest.mark.parametrize('direction', [-1, 1])
def test_receipt_rejects_lowered_or_raised_declared_seed_ceiling(
        reference, direction, tmp_path, monkeypatch):
    bad = copy.deepcopy(reference)
    saved = F(bad['cold_protocol']['canonical_seed_vector_error_upper'])
    bad['cold_protocol']['canonical_seed_vector_error_upper'] = str(
        saved+direction*F(1, 10**100))
    path = tmp_path/'changed-seed-ceiling.json'
    path.write_text(json.dumps(bad))
    monkeypatch.setattr(sys, 'argv', ['thermal', '--check', str(path)])
    with pytest.raises(ValueError):
        certificate.main()


@pytest.fixture
def small_preparation():
    # An asymmetric normalized target makes the work/reference flattening
    # convention observable without compiling the full twelve-qubit word.
    target = np.array([[.6, .8], [0., 0.]])
    word = certificate.compiler.compile_real_state(target.T.reshape(-1))
    executed = {'preparation': {'source_resources': certificate.compiler.resources(word)}}
    return executed, target


def test_preparation_resource_certificate_recompiles_actual_reference_order(
        small_preparation, monkeypatch):
    executed, target = small_preparation
    original = copy.deepcopy(executed)
    compile_state = certificate.compiler.compile_real_state
    arguments = []

    def capture(vector):
        arguments.append(vector.copy())
        return compile_state(vector)

    monkeypatch.setattr(certificate.compiler, 'compile_real_state', capture)
    certificate._certify_preparation_resources(executed, target)
    assert len(arguments) == 1
    np.testing.assert_array_equal(arguments[0], target.T.reshape(-1))
    original['preparation']['source_resources']['integrated_processor_norm_upper'] = str(
        certificate.PREPARATION_PROCESSOR_NORM_CEILING)
    assert executed == original


@pytest.mark.parametrize('field, replacement', [
    ('native_pulses', -1),
    ('proper_code_service_events', True),
    ('zero_angle_requests', 10**9),
    ('integrated_processor_norm_upper', 0.),
    ('integrated_processor_norm_upper', 277000.),
    ('integrated_processor_norm_upper', float('nan')),
    ('maximum_numerical_valid_code_leakage', 1.),
    ('service_failure_outcomes_retained', False),
])
def test_preparation_resource_certificate_rejects_unverified_resource_claims(
        small_preparation, field, replacement):
    executed, target = small_preparation
    executed['preparation']['source_resources'][field] = replacement
    with pytest.raises(ValueError):
        certificate._certify_preparation_resources(executed, target)


@pytest.mark.parametrize('angles', [
    [-277000.], [276999., 1.], [-276999., -1.],
    [263867.88048755], [float(np.nextafter(263867.88048755, np.inf))],
])
def test_preparation_resource_ceiling_accepts_exact_boundary_and_numerical_variation(
        angles, monkeypatch):
    word = [{'kind': 'Y', 'targets': [0], 'angle': angle} for angle in angles]
    executed = {'preparation': {'source_resources': certificate.compiler.resources(word)}}
    monkeypatch.setattr(certificate.compiler, 'compile_real_state', lambda vector: word)
    certificate._certify_preparation_resources(executed, np.eye(2)/np.sqrt(2))
    assert executed['preparation']['source_resources']['integrated_processor_norm_upper'] == '277000'


def test_preparation_resource_ceiling_uses_exact_absolute_angle_sum(monkeypatch):
    # The fresh floating resource sum equals the ceiling, but the exact sum
    # exceeds it. Summing signed angles or rounding first would accept this.
    word = [{'kind': 'Y', 'targets': [0], 'angle': -277000.},
            {'kind': 'Y', 'targets': [0], 'angle': 2.**-40}]
    resources = certificate.compiler.resources(word)
    assert resources['integrated_processor_norm_upper'] == 277000.
    executed = {'preparation': {'source_resources': resources}}
    monkeypatch.setattr(certificate.compiler, 'compile_real_state', lambda vector: word)
    with pytest.raises(ArithmeticError):
        certificate._certify_preparation_resources(executed, np.eye(2)/np.sqrt(2))
    assert executed['preparation']['source_resources']['integrated_processor_norm_upper'] == 277000.


@pytest.mark.parametrize('angle', [True, None, '1', complex(1),
                                  float('nan'), float('inf'), float('-inf')])
def test_preparation_resource_certificate_preserves_compiler_angle_validation(angle, monkeypatch):
    word = [{'kind': 'Y', 'targets': [0], 'angle': angle}]
    monkeypatch.setattr(certificate.compiler, 'compile_real_state', lambda vector: word)
    executed = {'preparation': {'source_resources': {}}}
    with pytest.raises(ValueError):
        certificate._certify_preparation_resources(executed, np.eye(2)/np.sqrt(2))


def test_twelve_qubit_preparation_template_fits_structural_resource_ceiling():
    # All 4096 amplitudes are exactly representable and normalized. Compile
    # every scheduled pulse, including zero-angle requests; do not execute.
    word = certificate.compiler.compile_real_state(np.full(4096, 1/64))
    counts = Counter(event['kind'] for event in word)
    assert counts == {'J': 40962, 'H': 81924, 'X': 8190, 'Z': 4095}
    # The template has angle bounds pi, pi/2, pi/4, pi, respectively.
    # Replacing pi by the strict rational upper bound 22/7 remains outward.
    caps = {'J': F(22, 7), 'H': F(11, 7), 'X': F(11, 14), 'Z': F(22, 7)}
    assert all(abs(F(event['angle'])) < caps[event['kind']] for event in word)
    structural_ceiling = sum((counts[kind]*cap for kind, cap in caps.items()), F())
    assert structural_ceiling == F(1937463, 7)
    assert structural_ceiling < certificate.PREPARATION_PROCESSOR_NORM_CEILING == F(277000)
