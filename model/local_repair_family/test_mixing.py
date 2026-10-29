"""Independent endpoint, matrix-power and adversarial custody controls."""
from dataclasses import replace
from fractions import Fraction as F

import numpy as np
import pytest

from . import local_bias as local, mixing as m
from ..observer import native


def independent_branches(parameters):
    states = native.states(binary=True)
    index = {x: i for i, x in enumerate(states)}
    branches = []
    for edge, seam_weight in zip(m.records.seams(), parameters.seam_weights, strict=True):
        a, b = edge
        for law, law_weight in (('C', parameters.c_weight), ('E', 1-parameters.c_weight)):
            for coin in (0, 1):
                groups = {}
                for i, state in enumerate(states):
                    total = state[a]+state[b]
                    output = list(state)
                    output[a], output[b] = (total+coin)//2, (total+1-coin)//2
                    current = state[a]-output[a]
                    label = current if law == 'C' else state[a]-state[b]
                    groups.setdefault(label, []).append((i, index[tuple(output)]))
                weight = (parameters.repair_probability*F(seam_weight, sum(parameters.seam_weights))
                          * law_weight*F(parameters.port_weights[a if coin else b],
                                        parameters.port_weights[a]+parameters.port_weights[b]))
                branches.extend((weight, tuple(pairs)) for pairs in groups.values())
    branches.append((1-parameters.repair_probability, tuple((i, i) for i in range(66))))
    return branches


@pytest.mark.parametrize('c,repair', [(F(1), F(1)), (F(0), F(1)), (F(1, 3), F(3, 4))])
def test_kernel_and_common_record_rate_from_independent_endpoint_arithmetic(c, repair):
    parameters = replace(m.canonical_parameters(), c_weight=c, repair_probability=repair,
                         seam_weights=tuple(1+i % 3 for i in range(30)))
    kernel = [[F() for _ in range(66)] for _ in range(66)]
    agreement = [[F() for _ in range(66)] for _ in range(66)]
    for weight, pairs in independent_branches(parameters):
        for i, j in pairs:
            kernel[i][j] += weight
            for other, _ in pairs:
                agreement[i][other] += weight
    assert m.probability_kernel(parameters) == tuple(map(tuple, kernel))
    assert m.coherence_rate(parameters) == max(agreement[i][j] for i in range(66) for j in range(66) if i != j)
    pi = local.reference(parameters)
    assert tuple(sum(pi[i]*kernel[i][j] for i in range(66)) for j in range(66)) == pi


def test_twenty_step_minorization_reconstructed_by_sequential_integer_powers():
    parameters = m.canonical_parameters()
    kernel = m.probability_kernel(parameters)
    base = np.array([[int(x*1980) for x in row] for row in kernel], dtype=object)
    actual = np.eye(66, dtype=object)
    for _ in range(20):
        actual = actual@base
    rows, scale = m.population_power(parameters, 20)
    assert scale == 1980**20
    assert rows == tuple(map(tuple, actual))
    minor = m.population_minorization(parameters)
    minima = tuple(min(row[j] for row in rows) for j in range(66))
    assert minor['column_minima'] == minima
    assert minor['beta'] == F(sum(minima), scale) > F(33, 200)
    # The residual is an actual stochastic matrix, not just a scalar claim.
    residual_mass = scale-sum(minima)
    for row in rows:
        assert all(row[j] >= minima[j] for j in range(66))
        assert sum(row[j]-minima[j] for j in range(66)) == residual_mass


def test_general_geometry_and_canonical_exact_finite_threshold():
    general = m.general_certificate()
    assert general['configuration_graph_diameter'] == 4
    assert F(general['population_minorization_lower']) == F(66, 38400**4)
    assert F(general['off_diagonal_rate_upper']) == F(3839, 3840)
    assert F(general['E_only_off_diagonal_rate_upper']) == F(6397, 6400)
    canonical = m.canonical_certificate()
    assert canonical['C_off_diagonal_rate'] == '893/990'
    assert canonical['E_off_diagonal_rate'] == '7/10'
    loose = F(65, 2)*F(893, 990)**1600+F(167, 200)**80
    assert loose <= F(canonical['half_diamond_upper']) < F(1, 10**6)
    assert m.mixing_bound(m.canonical_parameters(), 1600) <= loose
    assert m.mixing_bound(m.canonical_parameters(), 0) == 1


def test_canonical_fixed_tagged_mixtures_share_kernel_and_obey_convex_coherence_bound():
    base = m.canonical_parameters()
    kernel = m.probability_kernel(base)
    for weight in (F(), F(1, 256), F(1, 3), F(1, 2), F(255, 256), F(1)):
        mixture = replace(base, c_weight=weight)
        assert m.probability_kernel(mixture) == kernel
        assert m.coherence_rate(mixture) <= weight*F(893, 990)+(1-weight)*F(7, 10)


def test_canonical_certificate_refuses_a_different_but_stochastic_E_kernel(monkeypatch):
    original = m.probability_kernel

    def wrong_endpoint(parameters):
        rows = original(parameters)
        if parameters.c_weight == 0:
            changed = tuple((row[1], row[0], *row[2:]) for row in rows)
            assert changed != rows and all(sum(row) == 1 for row in changed)
            return changed
        return rows

    monkeypatch.setattr(m, 'probability_kernel', wrong_endpoint)
    with pytest.raises(ArithmeticError, match='endpoint population kernels disagree'):
        m.canonical_certificate()


def test_extremal_family_member_contains_the_complete_fair_uniform_instrument():
    parameters = local.Parameters(port_weights=(4, 4)+(1,)*10,
        seam_weights=(1,)+(16,)*29, c_weight=F(1, 256), repair_probability=F(1, 16))
    comparison = local.Parameters(c_weight=F(1, 256), repair_probability=F(1))
    alpha = F(1, 640)
    remainder_effect = [F()]*66
    for (weight, pairs), (fair_weight, fair_pairs) in zip(
            independent_branches(parameters), independent_branches(comparison), strict=True):
        assert pairs == fair_pairs
        remainder = weight-alpha*fair_weight
        assert remainder >= 0
        for i, _ in pairs:
            remainder_effect[i] += remainder
    assert remainder_effect == [1-alpha]*66
    kernel = m.probability_kernel(parameters)
    assert min(kernel[i][i] for i in range(66)) >= F(1, 5)
    assert min(kernel[i][j] for i in range(66) for j in range(66)
               if i != j and kernel[i][j]) >= F(1, 38400)
    assert m.coherence_rate(parameters) <= F(3839, 3840)


def test_wrong_uniform_reference_is_not_substituted_for_the_biased_target():
    kernel = m.probability_kernel(m.canonical_parameters())
    after = tuple(sum(kernel[i][j] for i in range(66))/66 for j in range(66))
    assert after != (F(1, 66),)*66
    assert local.reference(m.canonical_parameters()) != (F(1, 66),)*66


def test_finite_quantum_coherence_survives_and_reference_marginal_is_retained():
    parameters = m.canonical_parameters()
    states = native.states(binary=True)
    i = states.index((0, 0, 1, 1)+(0,)*8)
    j = states.index((0, 0, 1, 0, 1)+(0,)*7)
    # |i,0>+i|j,1>, normalized; apply every complete branch independently.
    density = np.zeros((66, 2, 66, 2), complex)
    density[i, 0, i, 0] = density[j, 1, j, 1] = .5
    density[i, 0, j, 1], density[j, 1, i, 0] = -.5j, .5j
    output = np.zeros_like(density)
    for weight, pairs in independent_branches(parameters):
        mapping = dict(pairs)
        for a, b in ((i, i), (i, j), (j, i), (j, j)):
            if a in mapping and b in mapping:
                output[mapping[a], :, mapping[b], :] += float(weight)*density[a, :, b, :]
    assert np.linalg.norm(output[:, 0, :, 1]) > 0
    assert np.allclose(np.einsum('iaib->ab', output), np.eye(2)/2, atol=1e-14, rtol=0)
    replacement = np.einsum('ij,ab->iajb', np.diag(list(map(float, local.reference(parameters)))), np.eye(2)/2)
    distance = np.linalg.svd((output-replacement).reshape(132, 132), compute_uv=False).sum()/2
    assert distance > .1
    assert distance <= float(m.mixing_bound(parameters, 1))


@pytest.mark.parametrize('field,bad', [
    ('port_weights', (True,)+(1,)*11), ('port_weights', (2,)*12),
    ('port_weights', (1,)*11+(5,)), ('port_weights', [1]*12),
    ('seam_weights', (1,)*29+(17,)), ('seam_weights', (1,)*29+(1.0,)),
    ('c_weight', .5), ('c_weight', True), ('repair_probability', F(0)),
    ('c_weight', F(1, 257)),
])
def test_mutated_original_parameters_are_revalidated_before_cached_results(field, bad):
    parameters = m.canonical_parameters()
    m.probability_kernel(parameters)
    object.__setattr__(parameters, field, bad)
    with pytest.raises(ValueError):
        m.mixing_bound(parameters, 1600)


@pytest.mark.parametrize('slots', [(True, 1), (1, True), (2, 2), (0, 0), (-1, -2), (1, -1)])
def test_noncanonical_fraction_internal_slots_are_rejected(slots):
    parameters = m.canonical_parameters()
    value = F(1)
    value._numerator, value._denominator = slots
    object.__setattr__(parameters, 'c_weight', value)
    with pytest.raises(ValueError):
        m.coherence_rate(parameters)


def test_public_fraction_and_parameter_aliases_cannot_poison_owned_cache():
    parameters = m.canonical_parameters()
    first = m.probability_kernel(parameters)
    expected = m.probability_kernel(parameters)
    first[0][0]._numerator = 999
    rate = m.coherence_rate(parameters)
    rate._denominator = 3
    minor = m.population_minorization(parameters)
    minor['beta']._numerator = 999
    minor['column_minima'] = ()
    parameters.c_weight._numerator = 0
    fresh = m.canonical_parameters()
    assert fresh.c_weight == 1
    assert m.probability_kernel(fresh) == expected
    assert m.coherence_rate(fresh) == F(893, 990)
    assert m.population_minorization(fresh)['beta'] > F(33, 200)


@pytest.mark.parametrize('bad', [True, 1.0, F(1), -1, m.MAX_POWER+1])
def test_original_matrix_horizons_rejected(bad):
    with pytest.raises(ValueError):
        m.population_power(m.canonical_parameters(), bad)


@pytest.mark.parametrize('bad', [True, 1.0, F(1), -1, m.MAX_ATTEMPTS+1])
def test_original_mixing_horizons_rejected(bad):
    with pytest.raises(ValueError):
        m.mixing_bound(m.canonical_parameters(), bad)


def test_zero_minorization_block_and_foreign_parameter_types_rejected():
    with pytest.raises(ValueError):
        m.population_minorization(m.canonical_parameters(), 0)
    with pytest.raises(ValueError):
        m.probability_kernel({})


def test_changed_native_map_does_not_reuse_warm_numeric_cache(monkeypatch):
    parameters = m.canonical_parameters()
    m.probability_kernel(parameters)
    original = native.native

    def wrong(state, edge, coin):
        output, record = original(state, edge, coin)
        return output, (record[0], record[1]+1)

    monkeypatch.setattr(native, 'native', wrong)
    with pytest.raises(ArithmeticError, match='endpoint output or emitted record'):
        m.probability_kernel(parameters)


def test_report_is_exact_compact_and_preserves_archive_boundary():
    value = m.report()
    assert value['general_family']['unique_full_quantum_source_stationary_state']
    assert not value['scope']['new_preparation_archive_in_replacement_target']
    assert not value['scope']['adaptive_feedback_stationarity_claimed']
    assert not value['scope']['native_source_pulses_executed']
    assert len(value['canonical']['population_column_minima']) == 66
