"""Gate the complete source-word relational charge protocol."""
from fractions import Fraction as F

import numpy as np
import pytest

from . import charged_reference_source as protocol
from .construction import source_kernel


def test_compiled_preparation_and_full_analyzer_effect():
    prepared = protocol._preparation()
    target = np.zeros(16, complex)
    target[2] = target[8] = 1j/np.sqrt(2)
    assert np.linalg.norm(prepared-target) < protocol.TOLERANCE

    columns = protocol._analyzer_columns()
    plus, minus = columns[:16], columns[16:]
    expected = np.eye(16, dtype=complex)/2
    expected[2, 8] = expected[8, 2] = 1/2
    assert np.linalg.norm(plus.conj().T@plus-expected, 2) < protocol.TOLERANCE
    assert np.linalg.norm(minus.conj().T@minus-(np.eye(16)-expected), 2) < protocol.TOLERANCE
    assert abs(np.linalg.norm(plus[:, 2])**2-F(1, 2)) < protocol.TOLERANCE
    assert abs(np.linalg.norm(plus[:, 1])**2-F(1, 2)) < protocol.TOLERANCE


@pytest.mark.parametrize('law', ['C', 'E'])
@pytest.mark.parametrize('coin', [0, 1])
def test_compiled_local_repair_all_binary_columns_and_raw_codes(law, coin):
    columns = protocol._repair_columns(law, coin)
    assert columns.shape == (1 << (3+(2 if law == 'C' else 3)), 4)
    assert np.linalg.norm(columns.conj().T@columns-np.eye(4), 2) < protocol.TOLERANCE
    for mask in range(4):
        target, current, imbalance = protocol._local_image(mask, coin)
        address = protocol.LOCAL_CODES[target]+(
            protocol._record_code(law, current, imbalance) << 3)
        assert protocol._permutation_image(law, coin, protocol.LOCAL_CODES[mask]) == address
        expected = np.zeros(len(columns), complex)
        expected[address] = 1
        assert np.linalg.norm(columns[:, mask]-expected) < protocol.TOLERANCE


def test_every_coin_record_pointer_branch_and_exact_uniform_aggregate():
    result = protocol.source_word_readout(source_kernel().seams)
    assert result['source_preparation_pulses'] == 5
    assert result['source_analyzer_pulses'] == 598
    assert result['source_repair_pulses'] == {'C': 5370, 'E': 5373}
    assert result['exact_repair_column_maps_checked'] == 16
    assert result['rows']['C']['all_coin_raw_record_pointer_branches'] == 16
    assert result['rows']['E']['all_coin_raw_record_pointer_branches'] == 32
    assert result['rows']['C']['zero_probability_abort_branches'] == 4
    assert result['rows']['E']['zero_probability_abort_branches'] == 12
    assert result['rows']['C']['uniform_thirty_seam_plus_probability_exact'] == '23/24'
    assert result['rows']['E']['uniform_thirty_seam_plus_probability_exact'] == '11/12'
    assert result['remote_seam_source_word_vacuum_columns_checked'] == 4
    assert result['remote_seam_exact_coherence_units_checked'] == 50
    assert result['maximum_floating_branch_or_probability_error'] < protocol.TOLERANCE
    assert not result['proper_code_service_failures_physically_bounded']
    assert not result['preparation_analyzer_or_C_E_law_selected_by_native_records']
    assert not result['particle_pole_or_rest_mass_identified']


def test_source_protocol_rejects_different_carrier():
    seams = source_kernel().seams
    with pytest.raises(ValueError):
        protocol.source_word_readout(seams[:-1])
    with pytest.raises(ValueError):
        protocol.source_word_readout(tuple(seams[:-1])+(seams[0],))
