from fractions import Fraction as F

import numpy as np
import pytest
from scipy.linalg import expm

from . import dispersion, model


def test_parameters_are_the_actual_native_reference_and_connected_odds():
    p = model.parameters()
    assert p.weights == (F(2, 225), F(1, 36), F(2, 45))
    assert abs(p.intra-p.inter-p.interaction) < 3e-16
    assert abs(p.intra+p.inter-p.clock) < 3e-16
    assert p.record()["empirical_mass_inputs"] == []


def test_full_fock_car_and_conserved_charge():
    identity = np.eye(16)
    for i, a in enumerate(model.annihilators()):
        for j, b in enumerate(model.annihilators()):
            assert np.array_equal(a@b+b@a, np.zeros((16, 16)))
            assert np.array_equal(a@b.T+b.T@a, identity if i == j else np.zeros((16, 16)))
    h = model.ring_hamiltonian()
    charge = np.diag(model.number_diagonal())
    assert np.max(abs(h@charge-charge@h)) == 0


def test_density_term_is_genuinely_quartic_not_a_one_body_offset():
    diagonal = np.diag(model.ring_terms()["interaction"])
    for i, j in model.RING_EDGES:
        difference = diagonal[(1<<i)+(1<<j)]-diagonal[1<<i]-diagonal[1<<j]+diagonal[0]
        assert abs(difference-model.parameters().interaction) < 1e-15


@pytest.mark.parametrize("interaction", [False, True])
@pytest.mark.parametrize("mu", [0., .37, -.81])
def test_sector_spectra_and_neutral_gap_match_independent_exact_formulas(interaction, mu):
    p = model.parameters()
    h = model.ring_hamiltonian(interaction=interaction, chemical_potential=mu)
    exact = dispersion.ring_sector_energies(coupling=p.interaction if interaction else 0., chemical_potential=mu)
    for number in range(5):
        indices = [i for i in range(16) if i.bit_count() == number]
        actual = np.linalg.eigvalsh(h[np.ix_(indices, indices)])
        np.testing.assert_allclose(actual, exact[number], atol=2e-14)


def test_charge_transfer_uses_a_real_reservoir_with_fixed_total_charge():
    vector, spec = model.charge_transfer_preparation()
    charges = np.array([i.bit_count() for i in range(256)])
    assert np.linalg.norm(vector[charges != 4]) == 0
    assert abs(np.vdot(vector, vector)-1) < 1e-14
    assert spec["state_preparation_is_an_architecture_service"]
    exact = dispersion.ring_gap_predictions()
    assert abs(spec["pair_gap"]-exact["charge_transfer_pair_gap"]) < 1e-14


def test_native_source_unitaries_converge_on_all_fock_columns():
    target = expm(-.25j*model.ring_hamiltonian())
    errors = []
    for steps in (2, 4, 8):
        actual, cost = model.source_cell_unitary(.25, steps)
        errors.append(np.linalg.norm(actual-target, 2))
        assert errors[-1] <= float(model.split_operator_bound(F(1, 4), steps))/2
        np.testing.assert_allclose(actual.conj().T@actual, np.eye(16), atol=3e-12)
        assert cost["native_pulses"] > 0
    assert errors[1] < errors[0]/3.9 and errors[2] < errors[1]/3.9


def test_source_path_does_not_consult_target_hamiltonian(monkeypatch):
    def reject(*args, **kwargs):
        raise AssertionError("target matrix used by source executor")
    monkeypatch.setattr(model, "ring_hamiltonian", reject)
    matrix, _ = model.source_cell_unitary(.1234567, 3)
    assert matrix.shape == (16, 16)


def test_cached_source_matrix_and_costs_cannot_be_corrupted():
    matrix, cost = model.source_cell_unitary(.1, 2)
    with pytest.raises(ValueError):
        matrix.setflags(write=True)
    cost["native_pulses"] = -1
    assert model.source_cell_unitary(.1, 2)[1]["native_pulses"] > 0
    with pytest.raises(ValueError):
        model.source_cell_unitary(.1, True)


@pytest.mark.parametrize("targets", [(0, 1), (1, 3), (3, 0), (2, 0)])
def test_tensor_contraction_preserves_order_and_entangled_reference(targets):
    rng = np.random.default_rng(241)
    state = rng.normal(size=(16, 3)) + 1j*rng.normal(size=(16, 3))
    gate = rng.normal(size=(4, 4)) + 1j*rng.normal(size=(4, 4))
    expected = np.zeros_like(state)
    for incoming in range(16):
        local = sum(((incoming>>bit)&1)<<j for j, bit in enumerate(targets))
        for outgoing in range(4):
            mask = incoming
            for j, bit in enumerate(targets):
                mask = (mask & ~(1<<bit)) | (((outgoing>>j)&1)<<bit)
            expected[mask] += gate[outgoing, local]*state[incoming]
    np.testing.assert_allclose(model.apply_subsystem(state, gate, targets, 4), expected, atol=2e-15)
