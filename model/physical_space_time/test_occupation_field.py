"""Independent native-history, random-clock and spatial-limit controls."""
from fractions import Fraction as F
from itertools import combinations, product

import numpy as np
import pytest
import sympy as sp

from . import occupation_field as f


def test_original_native_branches_and_canonical_covariance():
    states, a, b, z, checked = f.native_counts()
    assert checked == 3960 and len(states) == 66
    assert z*sp.ones(66, 1) == sp.zeros(12, 1)
    assert z*z.T/66 == f.CHI*f.P0
    assert a == a.T and b == b.T
    assert min((a+b)[i, i] for i in range(66)) > 0
    # Every occupation state is reached in the original exchange graph.
    reached, pending = {0}, [0]
    while pending:
        old = pending.pop()
        for new in range(66):
            if (a+b)[new, old] and new not in reached:
                reached.add(new); pending.append(new)
    assert len(reached) == 66


@pytest.mark.parametrize('clock', ['attempt', 'incident'])
def test_complete_configuration_law_closes_the_same_field(clock):
    _, _, _, z, _ = f.native_counts()
    full = f.configuration_kernel(clock)
    reduced = f.field(clock)['kernel']
    assert z*full == reduced*z
    assert full*sp.ones(66, 1) == sp.ones(66, 1)
    assert sp.ones(1, 66)*full == sp.ones(1, 66)
    assert min(full[i, i] for i in range(66)) > 0
    # Poisson equation and exact martingale covariance on the native state space.
    inverse = (sp.eye(12)-reduced+sp.ones(12)/12).inv()-sp.ones(12)/12
    potential = inverse*z
    assert potential-potential*full == z
    noise = f.CHI*(f.P0-reduced*reduced.T)
    assert inverse*noise*inverse.T == f.field(clock)['covariance']


def test_stopping_retains_remote_words_and_zero_current_incidents():
    _, a, b, _, _ = f.native_counts()
    complete = f.configuration_kernel('incident')
    partial, power = sp.zeros(66), sp.eye(66)
    for _ in range(3):
        partial += b*power/60
        power = a*power/60
    remainder = complete-partial
    assert all(x >= 0 for x in remainder)
    assert sp.ones(1, 66)*remainder == sp.Rational(5, 6)**3*sp.ones(1, 66)
    assert complete != b/10  # Replacing remote histories by identity fails.


@pytest.mark.parametrize('clock', ['attempt', 'incident'])
def test_precision_is_the_actual_Gaussian_action(clock):
    data = f.field(clock)
    assert data['covariance']*data['precision'] == f.P0
    assert data['covariance']*sp.ones(12, 1) == sp.zeros(12, 1)
    assert all(p > 0 for p in data['positive_LDL_pivots'])
    assert data['precision'] != data['laplacian']/(4*data['m']*f.CHI)


def test_mean_clock_conversion_misses_anisotropic_sampling_correction():
    attempt, incident = f.field('attempt'), f.field('incident')
    assert incident['kernel'] != incident['kernel'].T
    assert incident['covariance'] != attempt['covariance']/6
    delta = incident['covariance']/5-attempt['covariance']/30
    assert delta[0, 0] != delta[1, 1]
    # Neither a scalar amplitude nor a scalar clock can remove the defect.
    assert incident['covariance'][0, 0]*attempt['covariance'][1, 1] != (
        incident['covariance'][1, 1]*attempt['covariance'][0, 0])


def test_actual_rank_three_field_uses_the_same_gram_and_history():
    band = f.carrier_band()
    p, radial, tangent = band['projector'], band['radial'], band['tangent']
    assert f._expanded_product(p, p) == p
    assert p.trace() == 3 and set(p.diagonal()) == {sp.Rational(1, 4)}
    assert f._expanded_product(radial, tangent) == sp.zeros(12)
    assert radial.trace() == 1 and tangent.trace() == 2
    assert f._expanded_product(f.field()['laplacian'], p) == ((5-sp.sqrt(5))*p).applyfunc(sp.expand)
    assert band['attempt']['radial_variance'] == band['attempt']['tangent_variance']
    assert sp.simplify(band['incident']['radial_variance']-(5+sp.sqrt(5))/6) == 0
    assert sp.simplify(band['incident']['tangent_variance']-(115+19*sp.sqrt(5))/132) == 0
    assert band['incident']['radial_variance'] != band['incident']['tangent_variance']


@pytest.mark.parametrize('clock', ['attempt', 'incident'])
def test_projected_lags_against_full_configuration_process(clock):
    p = f.carrier_band()['projector']
    _, _, _, z, _ = f.native_counts()
    kernel = f.configuration_kernel(clock)
    for lag in (1, 2):
        full = f._expanded_product(p, z*(kernel**lag)*z.T/66, p)
        predicted = f._expanded_product(p, f.field(clock)['kernel']**lag, p)*f.CHI
        assert full == predicted.applyfunc(sp.expand)


def test_incident_projection_needs_eliminated_modes_and_marginal_precision():
    data = f.carrier_band()
    p, k = data['projector'], f.field('incident')['kernel']
    actual = f._expanded_product(p, k, k, p)
    reduced_twice = f._expanded_product(p, k, p, k, p)
    assert (actual-reduced_twice).applyfunc(sp.expand) == data['incident']['two_tick_closure_defect']
    assert data['incident']['two_tick_closure_defect'][0, 0] != 0
    assert data['attempt']['two_tick_closure_defect'] == sp.zeros(12)
    conditioned = f._expanded_product(p, f.field('incident')['precision'], p)
    assert conditioned != data['incident']['precision']


def test_local_clock_tape_does_not_disclose_the_remote_field():
    control = f.incomplete_tape_control()
    assert control['identical_incident_record'] == [0, 0]
    assert control['terminal_occupied_ports'] == [[2, 5], [1, 2]]
    assert not control['port_only_tape_determines_the_full_field']


@pytest.mark.parametrize('clock', ['attempt', 'incident'])
def test_finite_history_moments_against_independent_two_step_enumeration(clock):
    states, _, _, _, _ = f.native_counts()
    kernel = f.configuration_kernel(clock)
    values = [sp.Rational(state[0])-sp.Rational(1, 6) for state in states]
    moments = [sum(kernel[j, i]*(values[i]+values[j])**r for i, j in product(range(66), repeat=2))/66
               for r in range(5)]
    result = f.scalar_moments(clock, 2)
    assert result['sum_moments'] == [str(x) for x in moments]
    assert F(result['normalized_fourth_cumulant']) != 0
    assert f.covariance_at_count(clock, 1) == f.CHI*f.P0


@pytest.mark.parametrize('clock', ['attempt', 'incident'])
def test_finite_count_covariance_approaches_limit_with_certified_bound(clock):
    limit = np.asarray(f.field(clock)['covariance'], dtype=float)
    errors = []
    for n in (4, 16, 64):
        actual = np.asarray(f.covariance_at_count(clock, n), dtype=float)
        error = float(np.linalg.norm(actual-limit, 2))
        assert error <= float(f.covariance_error_bound(clock, n))
        errors.append(error)
    assert errors[-1] < errors[0]
    assert f.covariance_error_bound(clock, 8192) == f.covariance_error_bound(clock, 4096)/2


def test_conserved_mode_must_be_removed_before_time_average_limit():
    # A grand-canonical preparation has a random conserved mass. Its normalized
    # temporal sum's constant-mode variance grows with count instead of settling.
    rho, ports = F(1, 3), 12
    variance_at_one = rho*(1-rho)
    assert 64*variance_at_one > variance_at_one
    canonical = f.field('attempt')['covariance']
    assert (sp.ones(1, ports)*canonical*sp.ones(ports, 1))[0] == 0


def test_torus_is_a_declared_architecture_using_the_same_native_law():
    control = f.torus_native_control(4)
    assert control['native_branch_checks'] == 1024
    assert not control['this_is_the_existing_spherical_source_screen']
    edges = f.torus_edges(8)
    assert len(edges) == 128
    assert all(sum(vertex in edge for edge in edges) == 4 for vertex in range(64))


def test_mean_preserving_binary_leakage_is_rejected(monkeypatch):
    original = f.native.native

    def mutation(loads, seam, coin):
        if loads == (1, 1):
            return ((2, 0) if coin == 0 else (0, 2)), (coin, 2*coin-1)
        return original(loads, seam, coin)

    assert tuple(sum(mutation((1, 1), (0, 1), coin)[0][i] for coin in (0, 1))
                 for i in (0, 1)) == (2, 2)
    monkeypatch.setattr(f.native, 'native', mutation)
    with pytest.raises(ArithmeticError, match='binary exchange'):
        f.torus_native_control(4)


def test_fixed_band_Green_and_actual_local_clock_converge_together():
    rows = [f.torus_diagnostic(n) for n in (8, 16, 32, 64)]
    for key in ('green_operator_error_bound', 'attempt_total_covariance_error_bound',
                'incident_total_covariance_error_bound', 'clock_to_clock_covariance_error_bound'):
        assert all(F(b[key]) < F(a[key])/3 for a, b in zip(rows, rows[1:]))
    assert rows[0]['incident_cross_mode_diagnostic'] != 0
    assert F(rows[-1]['incident_total_covariance_error_bound']) < F(1, 500)


def test_native_time_response_stays_diffusive_in_the_torus_limit():
    # The same mode's density generator tends to -|k|^2/2 in per-edge attempt
    # time times a, whereas the spatial fluctuation action tends to |k|^2.
    # An inverse-Laplacian static covariance is not a propagating wave proof.
    for n in (16, 32, 64):
        energy = [(n/np.pi)**2*np.sin(np.pi*k/n)**2 for k in (1, 2)]
        assert energy[1]/energy[0] > 3.8
        assert abs(energy[1]/energy[0]-4) < 0.16


@pytest.mark.parametrize('ports,edges,total,clock_edges', [
    (3, ((0, 1), (1, 2)), 1, ((0, 1),)),
    (3, ((0, 1), (0, 2), (1, 2)), 2, ((0, 1), (0, 2))),
    (4, ((0, 1), (0, 3), (1, 2), (2, 3)), 2, ((0, 1),)),
])
def test_general_graph_formula_against_full_configuration_Poisson_solve(ports, edges, total, clock_edges):
    states = tuple(tuple(int(i in occupied) for i in range(ports))
                   for occupied in combinations(range(ports), total))
    index, dimension = {x: i for i, x in enumerate(states)}, len(states)
    a, b = sp.zeros(dimension), sp.zeros(dimension)
    for column, state in enumerate(states):
        for edge, coin in product(edges, (0, 1)):
            after, _ = f.native.native(state, edge, coin)
            (b if edge in clock_edges else a)[index[after], column] += sp.Rational(1, 2*len(edges))
    k = b*(sp.eye(dimension)-a).inv()
    state_projection = sp.ones(dimension)/dimension
    h = (sp.eye(dimension)-k+state_projection).inv()-state_projection
    z = sp.Matrix(ports, dimension, lambda i, j: sp.Rational(states[j][i])-sp.Rational(total, ports))
    actual = z*(h+h.T-sp.eye(dimension))*z.T/dimension
    l = f._laplacian(ports, edges)
    lr = f._laplacian(ports, tuple(e for e in edges if e not in clock_edges))
    projection = sp.ones(ports)/ports
    green = (l+projection).inv()-projection
    chi = sp.Rational(total*(ports-total), ports*(ports-1))
    expected = chi*(4*len(clock_edges)*green+lr*green+green*lr-(sp.eye(ports)-projection))
    assert actual == expected


def test_torus_incident_formula_against_full_density_resolvent():
    size = 8
    edges = f.torus_edges(size)
    local = tuple(e for e in edges if 0 in e)
    laplacian = np.asarray(f._laplacian(size**2, edges), dtype=float)
    li = np.asarray(f._laplacian(size**2, local), dtype=float)
    d = (8*np.eye(size**2)+laplacian-li)/(2*len(edges))
    b = (8*np.eye(size**2)-li)/(2*len(edges))
    k = np.linalg.solve(d.T, b.T).T
    pi = np.ones((size**2, size**2))/size**2
    h = np.linalg.inv(np.eye(size**2)-k+pi)-pi
    covariance = h+h.T-(np.eye(size**2)-pi)
    points = np.asarray(tuple(product(range(size), repeat=2)))
    columns, energies = [], []
    for wave in f.MODES:
        phase = 2*np.pi*(points@np.asarray(wave))/size
        columns.extend((np.sqrt(2)/size*np.cos(phase), np.sqrt(2)/size*np.sin(phase)))
        energy = (size/np.pi)**2*sum(np.sin(np.pi*x/size)**2 for x in wave)
        energies.extend((energy, energy))
    w = np.asarray(columns).T
    a_scale = size**2/(4*np.pi**2)
    actual = w.T@covariance@w/(16*a_scale)
    green = np.diag(1/np.asarray(energies))
    projected = w.T@li@w
    expected = green+np.eye(8)/(16*a_scale)-(projected@green+green@projected)/16
    assert np.linalg.norm(actual-expected, 2) < 1e-12


@pytest.mark.parametrize('value', [True, 0, -1, 1.5, '8'])
def test_invalid_count_rejected(value):
    with pytest.raises(ValueError):
        f.covariance_at_count('attempt', value)


def test_invalid_clock_or_aliased_torus_rejected():
    with pytest.raises(ValueError):
        f.field('movement')
    with pytest.raises(ValueError):
        f.torus_bounds(4)


def test_cached_theorem_data_cannot_be_replaced_by_a_caller():
    data = f.field('incident')
    with pytest.raises(TypeError):
        data['covariance'] = sp.eye(12)
    with pytest.raises(TypeError):
        data['kernel'][0, 0] = 1
