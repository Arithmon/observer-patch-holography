"""Independent operator and separating-cone tests for the mapping obstruction."""
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math
import sys

import numpy as np
import pytest

from model.identifiability_checks import dephasing_observability as d


@pytest.mark.parametrize('bits', list(product((0, 1), repeat=3)))
def test_cut_formula_matches_projective_recording_on_every_matrix_unit(bits):
    projectors = [np.diag([int(bit == outcome) for bit in bits]) for outcome in (0, 1)]
    rates = d.cut_rates(bits, F(3, 7))
    for a in range(3):
        for b in range(3):
            unit = np.zeros((3, 3))
            unit[a, b] = 1
            generator = (3/7)*(sum(p@unit@p for p in projectors)-unit)
            expected = -float(F(3, 7)*int(bits[a] != bits[b]))*unit
            np.testing.assert_array_equal(generator, expected)
    assert rates[0]+rates[1] >= rates[2]


@pytest.mark.parametrize('gamma', [F(0), F(1, 7), F(18), F(1000, 3)])
def test_separation_and_sharp_minimax_error_are_exact(gamma):
    result = d.generator_separation(gamma)
    assert F(result['target_functional_exact']) == -gamma
    assert F(result['sharp_maximum_entry_error_per_apparatus_second']) == gamma/3
    assert result['positive_recording_generators_can_equal_target'] == (gamma == 0)
    errors = tuple(map(F, result['attaining_errors']))
    assert errors == (gamma/3, gamma/3, -gamma/3)


def test_number_rates_are_the_independent_lindblad_matrix_unit_eigenvalues():
    number = np.diag([0., 1., 2.])
    for (a, b), rate in zip(d.PAIRS, d.number_rates(), strict=True):
        unit = np.zeros((3, 3))
        unit[a, b] = 1
        generator = 18*(number@unit@number-(number@number@unit+unit@number@number)/2)
        np.testing.assert_array_equal(generator, -float(rate)*unit)
    assert d.number_rates() == (F(9), F(9), F(36))


def test_positive_mixtures_obey_all_triangle_inequalities():
    # Several cuts and unequal rational rates test conic closure; the proof
    # itself is the nonnegative weighted sum of the eight exact inequalities.
    cuts = list(product((0, 1), repeat=3))
    total = tuple(sum((d.cut_rates(bits, F(k+1, k+2))[axis]
                       for k, bits in enumerate(cuts)), F(0)) for axis in range(3))
    assert all(2*value <= sum(total) for value in total)


@pytest.mark.parametrize('gamma,duration', [(F(18), F(1, 1000)), (F(3, 2), F(7, 3)),
                                           (F(0), F(1)), (F(18), F(0))])
def test_finite_time_gap_bound_and_identity_boundaries(gamma, duration):
    result = d.finite_time_witness(gamma, duration)
    exponent = 2*gamma*duration
    gap = -math.expm1(-float(exponent))/2
    assert gap+1e-15 >= float(F(result['probability_gap_lower_bound_exact']))
    assert result['strict_separation'] == bool(exponent)
    assert result['partial_parity_plus_probability_exact'] == '1'


@pytest.mark.parametrize('strength', [0., .125, .5, 1.])
def test_every_partial_parity_channel_preserves_same_parity_coherence(strength):
    rho = np.zeros((3, 3))
    rho[np.ix_((0, 2), (0, 2))] = .5
    projectors = (np.diag([1., 0., 1.]), np.diag([0., 1., 0.]))
    recorded = sum(p@rho@p for p in projectors)
    np.testing.assert_array_equal((1-strength)*rho+strength*recorded, rho)


def test_two_levels_do_not_obey_the_three_level_obstruction():
    # One bit at rate gamma/2 exactly matches the sole 0/1 coherence rate.
    target = d.number_rates()[0]
    assert d.cut_rates((0, 1, 1), target)[0] == target
    assert d.cut_rates((0, 1, 1), target)[2] != d.number_rates()[2]


@pytest.mark.parametrize('bad', [-1, F(-1, 2), True, 18., float('nan'), '18'])
def test_invalid_exact_rates_are_rejected(bad):
    with pytest.raises(ValueError):
        d.generator_separation(bad)
    with pytest.raises(ValueError):
        d.finite_time_witness(duration=bad)


@pytest.mark.parametrize('bits', [(0, 1), (0, 1, 2), (False, 1, 0), (0., 1, 0)])
def test_binary_partition_domain_is_explicit(bits):
    with pytest.raises(ValueError):
        d.cut_rates(bits)


def test_actual_C_E_control_preserves_even_coherence_and_distinguishes_odd():
    result = d.native_parity_control()
    proof = result['full_encoded_parity_refinement']
    assert proof['matrix_units_checked'] == 78**2
    assert proof['E_is_recorded_parity_after_C']
    even, odd = result['actual_native_instrument_matrix_unit_controls']
    assert (even['C_output_coefficient_sum_exact'], even['E_output_coefficient_sum_exact']) == ('1', '1')
    assert (odd['C_output_coefficient_sum_exact'], odd['E_output_coefficient_sum_exact']) == ('1/2', '0')
    assert not result['seam_loads_identified_with_laboratory_phonon_numbers']


def test_receipt_is_serializable_and_separates_math_from_empirical_verdict():
    result = d.run()
    assert json.loads(json.dumps(result)) == result
    assert not result['physical_mapping_established']
    assert not result['laboratory_candidate_likelihoods_computed']
    assert not result['C_or_E_rejected_by_data']
    assert result['primary_source']['paper'] == d.PAPER_URL
    assert result['primary_source']['data'] == d.DATA_URL
    assert result['source_bindings'] == {
        path: hashlib.sha256((d.ROOT/path).read_bytes()).hexdigest() for path in d.SOURCE_PATHS}


def test_cli_writes_once_replays_and_rejects_exact_receipt_tampering(tmp_path, monkeypatch, capsys):
    receipt = tmp_path/'receipt.json'
    monkeypatch.setattr(sys, 'argv', ['dephasing_observability', '--write', str(receipt)])
    d.main()
    with pytest.raises(FileExistsError):
        d.main()
    monkeypatch.setattr(sys, 'argv', ['dephasing_observability', '--check', str(receipt)])
    d.main()
    assert 'obstruction replayed' in capsys.readouterr().out
    value = json.loads(receipt.read_text(encoding='utf-8'))
    value['generator_obstruction']['target_functional_exact'] = '0'
    receipt.write_text(json.dumps(value), encoding='utf-8')
    with pytest.raises(ValueError, match='Receipt mismatch'):
        d.main()


def test_changed_laboratory_premise_is_rejected(monkeypatch):
    protocol = d.strict_json.load(d.ROOT/'model/anharmonic_observations/protocol.json')
    protocol['apparatus']['dephasing_per_second'] = 19
    monkeypatch.setattr(d.strict_json, 'load', lambda path: protocol)
    with pytest.raises(AssertionError, match='retained laboratory premise'):
        d.run()
