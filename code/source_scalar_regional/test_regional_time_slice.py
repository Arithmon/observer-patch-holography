"""False-green controls for exact minimal collar and regional algebra scope."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

import regional_time_slice as producer
import verify_regional_time_slice as verifier


@pytest.fixture(scope='module')
def receipt():
    return verifier.load(verifier.OUTPUT)


def test_frozen_receipt_has_exact_reproduction(receipt):
    assert producer.canonical(receipt) == producer.canonical(producer.build())
    result = verifier.verify(receipt)
    assert result['region_count'] == 8
    assert result['missing_row_counterexamples'] == 62


def test_one_force_combination_is_not_one_coordinate_read(receipt):
    singleton = receipt['regions'][2]
    assert singleton['minimal_linear_collar_dimension'] == 1
    assert len(singleton['coordinate_collar_sites']) == 6
    row = singleton['collar_rows_Qphi'][0]
    assert sum(v != ['0', '0'] for v in row) == 6


def test_missing_collar_counterexamples_are_actual_equal_data(receipt):
    # Independently evaluate the linear data map for all 62 witnesses.
    masses, rows, _, _, _ = verifier.model()
    for region in receipt['regions']:
        for witness in region['missing_row_controls']:
            site = witness['exterior_unit_site']
            velocity = [verifier.parse(x) for x in witness['local_velocity_shift_Qphi']]
            assert any(v != verifier.R() for v in velocity)
            for i, p in zip(region['region_sites'], velocity):
                next_field = verifier.TAU*p - verifier.TAU**2/2*dict(rows[i]).get(site, verifier.R())
                assert next_field == verifier.R()


def mutate(receipt, name):
    r = deepcopy(receipt)
    cube = r['regions'][3]
    if name == 'attested_outcomes':
        r['scope']['quantum_outcomes_attested'] = True
    elif name == 'physical_timeslice':
        r['scope']['physical_time_slice_or_clock'] = True
    elif name == 'classical_records_as_weyl':
        r['scope']['classical_records_are_not_weyl_access'] = False
    elif name == 'missing_region':
        r['regions'].pop(2)
    elif name == 'zero_collar_claim':
        cube['no_collar_recovers_original_regional_algebra'] = True
    elif name == 'lower_rank':
        cube['minimal_linear_collar_dimension'] -= 1
    elif name == 'boolean_rank':
        r['regions'][1]['minimal_linear_collar_dimension'] = True
    elif name == 'force_factor':
        cube['force_from_collar_rows_Qphi'][0][0] = ['0', '0']
    elif name == 'span_factor':
        cube['collar_from_force_rows_Qphi'][0][0] = ['0', '0']
    elif name == 'identity_minor':
        cube['collar_rows_Qphi'][0][cube['pivot_columns'][0]] = ['2', '0']
    elif name == 'dropped_control':
        cube['missing_row_controls'].pop()
    elif name == 'zero_perturbation':
        cube['missing_row_controls'][0]['local_velocity_shift_Qphi'] = [['0', '0']]*8
    elif name == 'lost_external_site':
        cube['exterior_sites'].pop()
    elif name == 'coordinate_access_cost':
        r['regions'][2]['coordinate_collar_sites'].pop()
    elif name == 'commuting_cross_time':
        cube['q_next_cross_bracket_over_i_hbar_Qphi'][0] = ['0', '0']
    elif name == 'changed_mass':
        r['mass_Qphi'][0] = ['1', '0']
    elif name == 'changed_clock':
        r['tau_Qphi'] = ['1', '0']
    elif name == 'extra_scope':
        r['scope']['native_access_derived'] = True
    elif name == 'rank_inflation':
        cube['two_layer_plus_collar_real_dimension'] += 1
    else:
        raise AssertionError(name)
    return r


@pytest.mark.parametrize('name', [
    'attested_outcomes', 'physical_timeslice', 'classical_records_as_weyl',
    'missing_region', 'zero_collar_claim', 'lower_rank', 'boolean_rank',
    'force_factor', 'span_factor', 'identity_minor', 'dropped_control',
    'zero_perturbation', 'lost_external_site', 'coordinate_access_cost',
    'commuting_cross_time', 'changed_mass', 'changed_clock', 'extra_scope', 'rank_inflation',
])
def test_verifier_rejects_false_green(receipt, name):
    with pytest.raises(ValueError):
        verifier.verify(mutate(receipt, name))


def test_coherently_erased_coupling_still_fails(receipt):
    r = deepcopy(receipt)
    singleton = r['regions'][2]
    singleton.update({
        'minimal_linear_collar_dimension': 0,
        'two_layer_plus_collar_real_dimension': 2,
        'no_collar_recovers_original_regional_algebra': True,
        'collar_rows_Qphi': [], 'collar_from_force_rows_Qphi': [],
        'force_from_collar_rows_Qphi': [[]], 'pivot_columns': [],
        'missing_row_controls': [], 'coordinate_collar_sites': [],
    })
    with pytest.raises(ValueError, match='full exterior force not reconstructed'):
        verifier.verify(r)


def test_noncanonical_and_duplicate_json_rejected(tmp_path):
    p = tmp_path/'bad.json'
    for content in ('{"a":1,"a":2}', '{"a":0.5}', '{"a":NaN}'):
        p.write_text(content)
        with pytest.raises(ValueError):
            verifier.load(p)
    with pytest.raises(ValueError):
        verifier.parse(['2/2', '0'])
