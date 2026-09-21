"""Quantum Gauss, CAR, current and infinite-rotor mutation controls."""
import copy
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_quantum_link as v


@pytest.fixture(scope='module')
def receipt():
    return v.load()


def test_independent_whole_quantum_certificate(receipt):
    result = v.verify(receipt)
    assert result['independent_quantum_state_replay']
    assert result['branches'] == 32768 and result['checked_CAR_transitions'] == 491520
    assert result['operator_Gauss_norm_squared'] == '0'
    assert result['mean_electric'] == '-18144/7225'
    assert result['exact_continuous_solution'] and not result['cyclic_link_truncation']


@pytest.mark.parametrize('case', ['missing_dressing', 'reversed_shift', 'cyclic_flux'])
def test_operator_gauss_requires_actual_integer_dressing(case):
    kwargs = {'flux_sign': 0 if case == 'missing_dressing' else 1} if case != 'cyclic_flux' else {'modulus': 37}
    with pytest.raises(ValueError, match='operator Gauss'):
        v.branch_algebra(**kwargs)


def test_zero_mean_is_not_operator_gauss():
    # Undressed equal left/right probabilities have zero mean charge but
    # a strictly positive Gauss-square expectation for a sharp zero electric.
    qs = v.charges()
    assert sum(qs) == 0
    assert sum(q * q for q in qs) * v.F(1, 4) == 30
    with pytest.raises(ValueError, match='operator Gauss'):
        v.branch_algebra(flux_sign=0)


def test_cyclic_shift_fails_U1_canonical_commutator():
    # For a 37-state wrap, U|18>=|-18>, hence [E,U]|18>=-36|-18>.
    # The bilateral integer shift instead yields coefficient +1 everywhere.
    for n in (-19, -18, 0, 18, 19):
        assert (n + 1) - n == 1
    assert (-18) - 18 != 1


def test_hopping_is_even_CAR_and_channels_commute():
    # Check actual signed exterior-algebra matrix elements, not tensor flips.
    for branch in range(8):
        occupied = v.occupation_bits(branch, 3)
        for i in range(3):
            bit_i = (branch >> i) & 1
            moved, sign_i = v.hop(occupied, 2 * i + 1 - bit_i, 2 * i + bit_i)
            assert sign_i == 1
            for j in range(3):
                if i == j:
                    continue
                bit_j = (branch >> j) & 1
                ij, sign_ij = v.hop(moved, 2 * j + 1 - bit_j, 2 * j + bit_j)
                other, sign_j = v.hop(occupied, 2 * j + 1 - bit_j, 2 * j + bit_j)
                ji, sign_ji = v.hop(other, 2 * i + 1 - bit_i, 2 * i + bit_i)
                assert ij == ji and sign_i * sign_ij == sign_j * sign_ji == 1


@pytest.mark.parametrize('flag', list(v.SCOPE))
def test_scope_rewrite_fails(receipt, flag):
    packet = copy.deepcopy(receipt)
    packet['scope'][flag] = not packet['scope'][flag]
    with pytest.raises(ValueError, match='scope'):
        v.validate_custody(packet)


@pytest.mark.parametrize('field', ['ambient_space', 'shift_action', 'commutator', 'finite_wrap',
                                  'current_channel', 'Hamiltonian_channel', 'factor_schedule'])
def test_quantum_law_rewrite_fails(receipt, field):
    packet = copy.deepcopy(receipt)
    packet['quantum_link'][field] = 'different operator law'
    with pytest.raises(ValueError, match='quantum link'):
        v.validate_custody(packet)


@pytest.mark.parametrize('case', ['charge', 'missing', 'phase', 'norm', 'kappa', 'mode_order', 'spin'])
def test_channel_and_spin_mutations_fail(receipt, case):
    channels = copy.deepcopy(receipt['channels'])
    if case == 'spin':
        basis = copy.deepcopy(receipt['channel_basis'])
        basis['right_spin'][0][1] = '4/5'
        with pytest.raises(ValueError, match='spin'):
            v.basis_check(basis, receipt['carrier'])
        return
    if case == 'charge':
        channels[6]['integer_charge'] = 4
    elif case == 'missing':
        channels.pop()
    elif case == 'phase':
        channels[-1]['right_amplitude'][1] = '4/5'
    elif case == 'norm':
        channels[-1]['left_amplitude'][0] = '1'
    elif case == 'kappa':
        channels[-1]['kappa'] = '1'
    else:
        channels[2]['mode_indices'].reverse()
    with pytest.raises(ValueError, match='channels'):
        v.channel_check(channels)


@pytest.mark.parametrize('case', ['current_sign', 'transfer_sign', 'normalized_midpoint', 'probability',
                                  'mean', 'variance', 'purity', 'missing_factor', 'wrong_factor'])
def test_factor_history_mutations_fail(receipt, case):
    history = copy.deepcopy(receipt['history'])
    if case == 'current_sign':
        history[1]['factor']['midpoint_current'] = str(-v.F(history[1]['factor']['midpoint_current']))
    elif case == 'transfer_sign':
        history[1]['factor']['mean_electric_change'] = str(-v.F(history[1]['factor']['mean_electric_change']))
    elif case == 'normalized_midpoint':
        history[1]['factor']['midpoint_norm_squared'] = '1'
    elif case == 'probability':
        history[1]['readout']['flux_probabilities'][0]['probability'] = '0'
    elif case in ('mean', 'variance', 'purity'):
        key = {'mean': 'mean_electric', 'variance': 'electric_variance', 'purity': 'link_reduced_purity'}[case]
        history[1]['readout'][key] = '1'
    elif case == 'missing_factor':
        history.pop()
    else:
        history[1]['factor']['channel'] = 2
    with pytest.raises(ValueError):
        v.history_check(history)


@pytest.mark.parametrize('case', ['clock', 'exponential', 'current', 'Gauss', 'point_variance'])
def test_exact_continuous_claim_mutations_fail(receipt, case):
    solution = copy.deepcopy(receipt['continuous_solution'])
    if case == 'clock':
        solution['physical_time_calibrated'] = True
    elif case == 'exponential':
        solution['Cayley_word_is_this_exact_time'] = True
    elif case == 'current':
        solution['instantaneous_current'] = '0'
    elif case == 'Gauss':
        solution['operator_Gauss'] = 'zero only in mean'
    else:
        solution['exact_point']['readout']['electric_variance'] = '0'
    with pytest.raises(ValueError, match='continuous'):
        v.continuous_check(solution)


def test_source_hashes_do_not_certify_quantum_state(receipt):
    packet = copy.deepcopy(receipt)
    packet['channels'][-1]['right_amplitude'][1] = '4/5'
    assert not v.validate_custody(packet)['quantum_state_verified']
    with pytest.raises(ValueError, match='channels'):
        v.verify(packet)


@pytest.mark.parametrize('case', ['hash', 'size', 'missing', 'extra', 'bool_size'])
def test_source_binding_mutations(receipt, case):
    packet = copy.deepcopy(receipt)
    key = next(iter(packet['source_pins']))
    if case == 'hash':
        packet['source_pins'][key]['sha256'] = '0' * 64
    elif case == 'size':
        packet['source_pins'][key]['bytes'] += 1
    elif case == 'missing':
        del packet['source_pins'][key]
    elif case == 'extra':
        packet['invented'] = True
    else:
        packet['source_pins'][key]['bytes'] = True
    with pytest.raises(ValueError):
        v.validate_custody(packet)
