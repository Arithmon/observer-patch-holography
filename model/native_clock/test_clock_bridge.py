"""Independent controls for count conversion, waiting tails and host replay."""
from fractions import Fraction as F
import copy

import pytest
import sympy as sp

from . import clock_bridge as b


@pytest.fixture(scope='module')
def source():
    kernel = b.transfer.native_counts()
    marks = b.marked_source(kernel)
    return kernel, marks, b.conversion_certificate(kernel, marks)


def test_complete_incident_resolvent_has_certified_unretained_waiting_tail(source):
    kernel, marks, _ = source
    # Independent finite native execution: stop at the first incident attempt.
    # The unresolved branch is retained, never conditionally normalized away.
    state_index = {state: i for i, state in enumerate(kernel.states)}
    waiting = {0: F(1)}
    stopped = [[F(0)]*66 for _ in range(2)]
    for _ in range(12):
        future = {}
        for i, mass in waiting.items():
            for seam in kernel.seams:
                for coin in (0, 1):
                    output, (_, current) = b.native.native(kernel.states[i], seam, coin)
                    j = state_index[output]
                    if b.PORT in seam:
                        stopped[abs(current)][j] += mass/60
                    else:
                        future[j] = future.get(j, F(0))+mass/60
        waiting = future
    tail = sum(waiting.values())
    assert tail == F(5, 6)**12
    assert sum(map(sum, stopped))+tail == 1
    for tick, Q in enumerate((marks['Q0'], marks['Q1'])):
        for j in range(66):
            assert stopped[tick][j] <= Q[0, j] <= stopped[tick][j]+tail
    assert sum(F(marks['Q1'][0, j]) for j in range(66)) > 0


def test_steady_rate_does_not_replace_joint_clock_law(source):
    _, marks, _ = source
    stationary, one = sp.ones(1, 66)/66, sp.ones(66, 1)
    p1 = F((stationary*marks['Q1']*one)[0])
    p11 = F((stationary*marks['Q1']*marks['Q1']*one)[0])
    assert p1 == F(5, 33)
    assert p11 == F(325, 7326)
    assert p11 != p1*p1
    assert p11/p1 == F(65, 222)


def test_joint_movement_generating_function_retains_geometric_incident_holds(source):
    kernel, marks, cert = source
    M = b.movement_generating_kernel(kernel, marks, F(1), F(1))
    for i, row in enumerate(cert['jump']):
        assert list(M.row(i)) == [dict(row).get(j, F(0)) for j in range(66)]
    weighted = b.movement_generating_kernel(kernel, marks, F(1, 2), F(2, 3))
    # Direct partial sum over the number g of incident holds before the move.
    # Remote zero holds have already been summed; retain the exact g>8 tail.
    for i, row in enumerate(cert['jump']):
        p = F(sum(kernel.marked_counts[1][i]), 60)
        q = marks['nonmoving_incident'][i]
        for j, jump_p in row:
            tick = abs(kernel.states[j][b.PORT]-kernel.states[i][b.PORT])
            partial = sum(jump_p*p/(p+q)*(q/(p+q))**g
                          *F(1, 2)**(g+tick)*F(2, 3)**tick for g in range(9))
            omitted_mass = jump_p*(q/(p+q))**9
            assert partial <= weighted[i, j] <= partial+omitted_mass
    assert all(sum(weighted.row(i)) < 1 for i in range(66))


@pytest.mark.parametrize('initial', [0, 10, 65, 'stationary'])
@pytest.mark.parametrize('feedback', [False, True])
def test_joint_field_clock_feedback_law_keeps_zero_hits(source, initial, feedback):
    kernel, marks, _ = source
    law = b.count_field_law(kernel, marks, 2, initial, feedback=feedback)
    assert sum(law.values()) == 1
    assert any(count == 0 for count, _, _, _ in law)
    for (count, first, final, incident), mass in law.items():
        assert mass > 0
        assert incident == 2+int(feedback and first)
        assert 0 <= count <= incident
        if initial != 'stationary':
            assert final == (kernel.states[initial][b.PORT]+count) % 2
    if feedback:
        assert any(incident == 3 for _, _, _, incident in law)


def _direct_moments(kernel, marks, prior, n):
    """Propagate uncorrected clock moments using direct native branch counts."""
    mass, first, second = list(prior), [F(0)]*66, [F(0)]*66
    for _ in range(n):
        future = [[F(0)]*66 for _ in range(3)]
        for i, state in enumerate(kernel.states):
            escape = F(sum(kernel.marked_counts[1][i]), 60)
            q = marks['nonmoving_incident'][i]
            mean_g = q/escape
            var_g = q*(escape+q)/escape**2
            for j, choices in enumerate(kernel.marked_counts[1][i]):
                if not choices:
                    continue
                probability = F(choices, 60)/escape
                tick = abs(kernel.states[j][b.PORT]-state[b.PORT])
                mean = mean_g-(b.RATIO-1)*tick
                moment2 = var_g+mean**2
                future[0][j] += probability*mass[i]
                future[1][j] += probability*(first[i]+mean*mass[i])
                future[2][j] += probability*(second[i]+2*mean*first[i]+moment2*mass[i])
        mass, first, second = future
    return sum(first), sum(second)


@pytest.mark.parametrize('initial', [0, 10, 65, 'palm'])
def test_nonstationary_conversion_bound_against_direct_moment_propagation(source, initial):
    kernel, marks, cert = source
    prior = cert['palm'] if initial == 'palm' else tuple(F(i == initial) for i in range(66))
    mean, second = _direct_moments(kernel, marks, prior, 8)
    assert abs(mean) <= cert['range']
    assert second <= 2*8*cert['variance_ceiling']+2*cert['range']**2
    if initial == 'palm':
        assert mean == 0
    else:
        assert mean != 0
    assert cert['range'] == F(1254, 125)
    assert cert['variance_ceiling'] == F(29473097353, 1176125000)
    assert cert['asymptotic_variance_per_global_movement'] == F(43736681, 3637500)


def test_error_radius_is_conservative_and_strictly_unit_independent(source):
    _, _, cert = source
    assert b.error_radius(cert, 10000, F(1, 60)) == 3889
    for n in (1, 10000, 10**12):
        radius = b.error_radius(cert, n, F(1, 60))
        assert radius > cert['range']
        assert n*cert['variance_ceiling']/(radius-cert['range'])**2 <= F(1, 60)


def test_retained_host_comparison_preserves_calibration_and_rejection(source):
    kernel, _, cert = source
    actual = b.retained_comparison(kernel, cert)
    assert not actual['new_timestamps_or_fitted_parameters']
    assert actual['calibration']['clock_port'] == 2
    assert actual['comparison_port'] == 3
    rows = actual['held_out_comparisons']
    assert [row['incident_hits'] for row in rows[:3]] == [11168, 11060, 11160]
    assert [row['zero_current_incident_hits'] for row in rows[:3]] == [9564, 9366, 9568]
    assert all(row['within_native_count_bound'] for row in rows)
    assert [row['host_interval_verdict'] for row in rows] == ['within_tolerance']*3+['reject']*3
    for a, c in zip(rows[:3], rows[3:]):
        assert a['incident_hits'] == c['incident_hits']
        assert a['local_movement_ticks'] == c['local_movement_ticks']
        assert a['host_duration_ns'] != c['host_duration_ns']


def test_corrupted_waiting_marks_or_raw_currents_cannot_pass(source, tmp_path):
    kernel, _, cert = source
    raw = b.strict_json.loads(b.host_clock.RAW_PATH.read_text())
    bad = copy.deepcopy(raw)
    records = bytearray.fromhex(bad['runs'][0]['records_hex'])
    records[0] += 1  # do not silently accept a different current in a held record
    bad['runs'][0]['records_hex'] = records.hex()
    path = tmp_path/'bad.json'
    path.write_text(b.canonical(bad))
    with pytest.raises(ValueError, match='native current mismatch'):
        b.retained_comparison(kernel, cert, path)


@pytest.mark.parametrize('n,p', [(True, F(1, 60)), (0, F(1, 60)), (1, 0.05), (1, F(0)), (1, F(1))])
def test_conversion_bound_rejects_undeclared_domains(source, n, p):
    with pytest.raises(ValueError):
        b.error_radius(source[2], n, p)
