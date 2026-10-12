"""Quantum region restrictions, actual interaction support and cutoff controls."""
import numpy as np
import pytest
from scipy.linalg import expm

from . import scalar_coupling as scalar
from . import spatial as s


def hermitian(size, seed=716):
    rng = np.random.default_rng(seed)
    value = rng.normal(size=(size, size)) + 1j*rng.normal(size=(size, size))
    return (value + value.conj().T)/2


@pytest.mark.parametrize("d", (2, 4))
def test_scalar_regions_commute_on_the_actual_shared_register(d):
    left, right = hermitian(d), hermitian(d, 911)
    a, b = s.scalar_operator(d, 0, left), s.scalar_operator(d, 1, right)
    assert np.allclose(a, np.kron(np.kron(np.eye(d), left), np.eye(4)))
    assert np.allclose(b, np.kron(np.kron(right, np.eye(d)), np.eye(4)))
    assert np.linalg.norm(s.commutator(a, b), 2) < 1e-13
    # Arbitrary scalar words are neutral and leave the entire physical gauge
    # sector alone; neither occupation bit is relabeled as a local scalar site.
    gauge = s.embed_operator(d, "gauge_matter_link", hermitian(4))
    assert np.linalg.norm(s.commutator(a, gauge), 2) < 1e-13


@pytest.mark.parametrize("region", [(), ("scalar_site0",),
    ("scalar_site0", "gauge_matter_link"), s.REGIONS])
def test_cp_restriction_kraus_and_reference_positivity(region):
    d, size = 2, 16
    maps = s.restriction_kraus(d, region)
    out_size = maps[0].shape[0]
    assert np.allclose(sum(k @ k.conj().T for k in maps), np.eye(out_size))
    # Positive entangled input on a spectator qubit; retain every Kraus arm.
    rng = np.random.default_rng(267)
    v = rng.normal(size=2*size) + 1j*rng.normal(size=2*size)
    v /= np.linalg.norm(v)
    rho = np.outer(v, v.conj())
    restricted = sum(np.kron(np.eye(2), k) @ rho @ np.kron(np.eye(2), k.conj().T) for k in maps)
    assert np.linalg.eigvalsh(restricted)[0] >= -2e-15
    observable = hermitian(size)
    by_kraus = sum(k @ observable @ k.conj().T for k in maps)
    assert np.allclose(s.embed_operator(d, region, by_kraus), s.conditional_expectation(d, observable, region))


@pytest.mark.parametrize("d", (2, 4))
def test_conditional_expectation_is_bimodule_projection_and_not_homomorphism(d):
    x = hermitian(4*d*d)
    a, b = s.scalar_operator(d, 0, hermitian(d)), s.scalar_operator(d, 0, hermitian(d, 911))
    e = lambda value: s.conditional_expectation(d, value, "scalar_site0")
    assert np.allclose(e(np.eye(4*d*d)), np.eye(4*d*d))
    assert np.allclose(e(e(x)), e(x))
    assert np.allclose(e(a), a)
    assert np.allclose(e(a @ x @ b), a @ e(x) @ b)
    assert np.trace(e(x)) == pytest.approx(np.trace(x))
    remote = np.zeros((d, d)); remote[0, 1] = remote[1, 0] = 1
    y = s.scalar_operator(d, 1, remote)
    assert np.allclose(e(y), 0)
    assert np.allclose(e(y @ y), (2/d)*np.eye(4*d*d))
    assert np.linalg.norm(e(y @ y)-e(y) @ e(y), 2) == pytest.approx(2/d)
    # Intersecting region restrictions compose without duplicating a factor.
    first = s.conditional_expectation(d, x, ("scalar_site0", "gauge_matter_link"))
    assert np.allclose(s.conditional_expectation(d, first, ("scalar_site0", "scalar_site1")), e(x))


@pytest.mark.parametrize("d", (2, 4))
def test_regional_states_are_actual_source_marginals_and_keep_branch_weight(d):
    state, _ = scalar.source_evolve(d, 16)
    for site in (0, 1):
        a = hermitian(d, 718+site)
        density = s.regional_density(d, state, f"scalar_site{site}")
        assert np.trace(density) == pytest.approx(1, abs=2e-12)
        assert np.linalg.eigvalsh(density)[0] >= -2e-15
        assert np.trace(density @ a) == pytest.approx(np.vdot(state, s.scalar_operator(d, site, a) @ state))
        branch = s.regional_density(d, np.sqrt(.17)*state, f"scalar_site{site}")
        assert np.allclose(branch, .17*density)
    joint = s.regional_density(d, state, s.REGIONS)
    assert np.allclose(joint, np.outer(state, state.conj()))
    gauge = s.regional_density(d, state, "gauge_matter_link")
    assert np.trace(gauge @ gauge).real < .96


@pytest.mark.parametrize("d", (2, 4))
def test_graph_recovers_only_actual_cross_region_terms_and_deletion_controls(d):
    expected = [["scalar_site0", "scalar_site1"], ["scalar_site0", "gauge_matter_link"]]
    graph = s.interaction_graph(d)
    assert graph["edges"] == expected
    assert graph["factorization_and_localization_supplied"]
    assert not graph["charged_endpoint_tensor_factorization_claimed"]
    assert not graph["metric_or_spatial_dimension_recovered"]
    witnesses = {tuple(row["regions"]): row["norm"] for row in graph["double_commutator_witnesses"]}
    assert witnesses[tuple(expected[0])] > .14
    assert witnesses[tuple(expected[1])] > .49
    assert witnesses[("scalar_site1", "gauge_matter_link")] < 1e-13
    terms = scalar.hamiltonian_terms(d)
    h = sum(terms.values())
    assert s.interaction_graph(d, coupling=0)["edges"] == expected[:1]
    assert s.interaction_graph(d, hamiltonian=h-terms["bond"])["edges"] == expected[1:]
    assert s.interaction_graph(d, hamiltonian=h-terms["bond"]-terms["interaction"])["edges"] == []
    pieces = s.interaction_components(d, h)
    assert np.allclose(sum(pieces.values()), h)
    assert np.allclose(pieces[("scalar_site0", "scalar_site1")], terms["bond"])
    assert np.allclose(pieces[("scalar_site0", "gauge_matter_link")], terms["interaction"])
    for group, component in pieces.items():
        for deleted in group:
            retained = tuple(name for name in s.REGIONS if name != deleted)
            assert np.linalg.norm(s.conditional_expectation(d, component, retained), 2) < 2e-14


def test_graph_support_is_invariant_under_local_coordinate_changes():
    h = scalar.ideal_hamiltonian(2)
    unitary = expm(1j*s.scalar_operator(2, 0, hermitian(2)))
    original = s.interaction_graph(2)
    changed = s.interaction_graph(2, hamiltonian=unitary @ h @ unitary.conj().T)
    assert changed["edges"] == original["edges"]
    for a, b in zip(original["centered_component_norms"], changed["centered_component_norms"], strict=True):
        assert a["region"] == b["region"]
        assert a["norm"] == pytest.approx(b["norm"], abs=1e-13)


def test_cutoff_compression_agrees_but_generator_and_unital_embeddings_fail():
    v = s.occupation_embedding()
    h2, h4 = scalar.ideal_hamiltonian(2), scalar.ideal_hamiltonian(4)
    assert np.array_equal(v.T @ v, np.eye(16))
    assert np.allclose(v.T @ h4 @ v, h2)
    assert np.linalg.norm(v @ v.T-np.eye(64), 2) == 1
    residual = h4 @ v-v @ h2
    blank = residual[:, 0]
    expected = np.zeros(64)
    expected[[8, 32]] = 9*np.sqrt(2)/128
    assert np.allclose(blank, expected)
    assert np.linalg.norm(blank) == pytest.approx(9/64)
    assert np.linalg.norm(residual, 2) > .6


def test_receipt_keeps_source_counts_and_physical_scope():
    result = s.spatial_receipt(2)
    assert result["source_execution"]["evolution"]["resources"]["native_pulses"] > 0
    assert len(result["regional_states"]) == 3
    assert not result["physical_clock_or_spatial_continuum_claimed"]


@pytest.mark.parametrize("region", [("bad",), ("scalar_site0", "scalar_site0"), {"scalar_site0"}, 0])
def test_rejects_ambiguous_regions(region):
    with pytest.raises(ValueError):
        s.conditional_expectation(2, np.eye(16), region)


def test_rejects_invalid_states_dimensions_and_operators():
    for density in (-np.eye(16)/16, np.eye(16), np.full((16, 16), np.nan)):
        with pytest.raises(ValueError):
            s.regional_density(2, density, "scalar_site0")
    for bad in (True, 1, 8):
        with pytest.raises(ValueError):
            s.interaction_graph(bad)
    with pytest.raises(ValueError):
        s.scalar_operator(2, True, np.eye(2))
    with pytest.raises(ValueError):
        s.interaction_components(2, np.triu(np.ones((16, 16))))
    with pytest.raises(ValueError):
        s.occupation_embedding(4, 2)
