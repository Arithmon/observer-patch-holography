"""Complete source support, memory dynamics and separated hierarchy comparison."""
import numpy as np
import pytest
import sympy as s

from . import native_mass_record_cycles as m


@pytest.fixture(scope='module')
def data(): return m.architecture()


@pytest.fixture(scope='module')
def execution(): return m.execute()


def test_event_is_selected_from_source_support_and_moves_work(data):
    supported, moving = m.moving_E_labels()
    assert data['E_event'] == (0, (moving or supported)[0])
    assert data['q'] > 0 and data['coincidence'] == data['p']*data['q']
    assert data['zero_E_control_weight'] == 0
    assert 1 not in supported


def test_actual_edge_rotations_and_complete_stationary_memory(data):
    C, E = map(m.objects.permutation_matrix, (data['gC'], data['gE']))
    assert C*C == E*E == s.eye(12) and C*E != E*C
    assert data['T']*data['stationary'] == data['stationary']
    assert sum(data['stationary']) == 1 and all(x >= 0 for x in data['stationary'])
    assert len(data['labels']) == 1536
    assert {row[4] for row in data['labels']} == set(range(12))
    assert len({row[:4] for row in data['labels']}) == 128


def test_zero_sum_contraction_including_transient_ports(data):
    T = data['T']
    for i in range(12):
        for j in range(12):
            assert sum(abs(x) for x in T[:, i]-T[:, j]) <= 2*data['p']


def test_positive_kernel_and_exact_root_certificate(data):
    Y = data['Y']; certificate = m.spectral_certificate()
    assert Y == Y.T and s.trace(Y) == 1
    assert certificate['rank'] == Y.rank()
    if Y.rank() == 3:
        bounds = [[s.Rational(x) for x in pair] for pair in certificate['root_intervals_exact']]
        assert all(0 < lo < hi for lo, hi in bounds)
        assert sum(lo for lo, _ in bounds) < 1 < sum(hi for _, hi in bounds)
        assert bounds[0][1] < bounds[1][0] and bounds[1][1] < bounds[2][0]
        assert certificate['three_positive_couplings']


def test_complete_finite_channel_is_cp_tp_and_has_declared_scope(execution):
    report, arrays = execution
    channel = arrays['cycle_finite_channel']
    choi = channel.reshape(6, 6, 6, 6).transpose(0, 2, 1, 3).reshape(36, 36)
    assert np.linalg.norm(choi-choi.conj().T) < 1e-12
    assert np.linalg.eigvalsh(choi).min() > -1e-12
    assert np.linalg.norm(np.eye(6).reshape(36)@channel-np.eye(6).reshape(36)) < 1e-12
    assert not report['scope']['source_native_pulse_history_executed']
    assert not report['scope']['physical_fermion_mass_map_derived']


def test_all_raw_branches_preserve_every_old_memory_input(execution):
    _, arrays = execution
    grams = np.zeros((12, 6, 6), complex)
    for index, row in enumerate(arrays['cycle_labels']):
        U = arrays['cycle_pulses'][row[5]]
        grams[row[4]] += arrays['cycle_raw_weights'][index//12]*(U.conj().T@U)
    assert np.max(abs(grams-np.eye(6))) < 1e-12


@pytest.mark.parametrize('key', ['cycle_labels', 'cycle_raw_weights', 'cycle_finite_channel'])
def test_memory_record_weight_and_output_tampering_rejected(execution, key):
    _, arrays = execution; changed = dict(arrays); changed[key] = arrays[key].copy()
    changed[key].flat[0] += 1 if changed[key].dtype.kind in 'iu' else 1e-5
    with pytest.raises(ValueError): m.validate_arrays(changed)


def test_source_and_root_producers_do_not_read_empirical_tables(monkeypatch):
    def forbidden(*args, **kwargs): raise AssertionError('empirical packet entered source producer')
    monkeypatch.setattr(m.strict_json, 'load', forbidden)
    m.architecture.cache_clear()
    assert m.architecture()['p'] == s.Rational(1, 12)
    assert m.spectral_certificate()['rank'] >= 1
