"""Independent CAR/flux reconstruction and complete native-channel controls."""
from fractions import Fraction as Q
import math

import numpy as np
import pytest
from scipy.linalg import expm

from model.observer.source_scalar import compiler
from model.common_source import gauge_matter as gm


def independent_annihilator(mode):
    # Jordan-Wigner tensor representation; no producer ladder helper.
    z = np.diag([1, -1])
    lower = np.array([[0, 1], [0, 0]])
    result = np.ones((1, 1), dtype=int)
    for j in reversed(range(4)):
        result = np.kron(result, lower if j == mode else z if j < mode else np.eye(2, dtype=int))
    return result


def independent_isometry_and_ambient():
    # Finite array for the three populated integer flux values only. This
    # isometry is NOT a cyclic rotor, and no full-link shift is inferred.
    fluxes = (-1, 0, 1)
    code = np.zeros((48, 4), dtype=int)
    for b1 in range(2):
        for b0 in range(2):
            occupation = (1 << b0) | (1 << (2+b1))
            code[16*fluxes.index(b1-b0)+occupation, b0+2*b1] = 1
    return code, fluxes


def test_car_signs_on_full_fock_space_and_nonadjacent_countercontrol():
    for mode in range(4):
        destroy = independent_annihilator(mode)
        for creation, operator in ((False, destroy), (True, destroy.T)):
            for mask in range(16):
                column = operator[:, mask]
                expected = [(int(column[k]), k) for k in np.flatnonzero(column)]
                actual = gm.car_ladder(mask, mode, creation)
                assert expected == ([] if actual is None else [actual])
    annihilators = [independent_annihilator(j) for j in range(4)]
    for i, a in enumerate(annihilators):
        for j, b in enumerate(annihilators):
            assert np.array_equal(a@b+b@a, np.zeros((16, 16), int))
            assert np.array_equal(a@b.T+b.T@a, np.eye(16, dtype=int)*(i == j))
    # A commuting-bit implementation misses this negative fermionic sign.
    assert (annihilators[2].T@annihilators[0])[6, 3] == -1


def test_embedding_both_operator_gauss_constraints_and_dressed_hopping():
    code, fluxes = independent_isometry_and_ambient()
    assert np.array_equal(code.T@code, np.eye(4, dtype=int))
    a = [independent_annihilator(j) for j in range(4)]
    number = [item.T@item for item in a]
    flux = np.kron(np.diag(fluxes), np.eye(16, dtype=int))
    charge_left = np.kron(np.eye(3, dtype=int), number[0]-number[2])
    charge_right = np.kron(np.eye(3, dtype=int), number[1]-number[3])
    assert np.array_equal((flux-charge_left)@code, np.zeros((48, 4), int))
    assert np.array_equal((-flux-charge_right)@code, np.zeros((48, 4), int))
    ambient_hopping = np.zeros((48, 48), dtype=int)
    for channel, charge in enumerate((1, -1)):
        # A finite compression of U^q is used ONLY to compute the invariant
        # code action. A wrong cyclic wrap is never permitted.
        up = np.array([[int(target == source+charge) for source in fluxes]
                       for target in fluxes])
        right_to_left = a[2*channel].T@a[2*channel+1]
        ambient_hopping += np.kron(up, right_to_left)+np.kron(up.T, right_to_left.T)
    expected = code.T@ambient_hopping@code
    assert np.array_equal(ambient_hopping@code, code@expected)
    assert np.array_equal(gm.numerical_matrix(gm.hopping_hamiltonian()), expected)
    assert np.array_equal(gm.numerical_matrix(gm.electric_operator()), code.T@flux@code)
    assert np.array_equal(gm.numerical_matrix(gm.electric_energy()), code.T@flux@flux@code)
    for row in gm.physical_basis():
        index = row["index"]
        nonzero = int(np.flatnonzero(code[:, index])[0])
        assert row["fermion_occupation_mask"] == nonzero % 16
        assert row["electric_flux"] == fluxes[nonzero//16]


def test_integer_link_no_wrap_in_far_away_ambient_inputs():
    for flux in (-1000000, 1000000):
        for mask in range(16):
            for channel in (0, 1):
                for _, output, moved in gm.dressed_hopping_images(mask, flux, channel):
                    assert abs(moved-flux) == 1
                    assert gm.gauss_values(output, moved) == gm.gauss_values(mask, flux)


def test_exact_echo_independent_closed_form_all_columns():
    # Closed echo polynomial in the independently represented ZZ observable.
    u = np.array([[3, -4j], [-4j, 3]], complex)/5
    hopping = np.kron(u, u)
    electric = np.diag([1, (3-4j)/5, (3-4j)/5, 1])
    echo = hopping.conj().T@electric@hopping
    assert np.linalg.norm(gm.numerical_matrix(gm.ideal_echo())-echo, 2) < 2e-15
    assert gm.matmul(gm.dagger(gm.ideal_echo()), gm.ideal_echo()) == gm.identity(4)
    assert gm.ideal_echo()[0][0] == (Q(2549, 3125), -Q(1152, 3125))
    assert gm.ideal_echo()[1][0] == (Q(336, 3125), -Q(168, 3125))
    assert gm.ideal_echo()[2][0] == (Q(336, 3125), -Q(168, 3125))
    assert gm.ideal_echo()[3][0] == (-Q(576, 3125), -Q(1152, 3125))


@pytest.mark.parametrize("bit", (0, 1))
@pytest.mark.parametrize("initial", range(4))
def test_complete_binary_instrument_on_every_input(bit, initial):
    probabilities = [gm.branch_probability(gm.ideal_echo(), initial, bit, outcome)
                     for outcome in (0, 1)]
    assert sum(probabilities, Q(0)) == 1
    assert probabilities[1-((initial >> bit) & 1)] == Q(576, 3125)
    for outcome in (0, 1):
        p = gm.occupation_projector(bit, outcome)
        assert gm.matmul(p, p) == p
        assert gm.matmul(p, gm.electric_operator()) == gm.matmul(gm.electric_operator(), p)


def test_interacting_electric_term_and_counterfactual_differ():
    h = gm.numerical_matrix(gm.hopping_hamiltonian())
    e = gm.numerical_matrix(gm.electric_operator())
    v = e@e
    assert np.linalg.norm(h@v-v@h, 2) == pytest.approx(2.)
    current = gm.numerical_matrix(gm.electric_current_operator())
    assert np.array_equal(1j*((h+v)@e-e@(h+v)), -current)
    no_electric = gm.matmul(gm.dagger(gm.ideal_hopping()), gm.ideal_hopping())
    assert no_electric == gm.identity(4)
    assert gm.branch_probability(no_electric, 0, 0, 1) == 0
    assert gm.branch_probability(gm.ideal_echo(), 0, 0, 1) == Q(576, 3125)
    # Electric operation cannot be factored as phases on independent channels.
    diagonal = np.diag(gm.numerical_matrix(gm.ideal_electric()))
    assert abs(diagonal[0]*diagonal[3]-diagonal[1]*diagonal[2]) > 1


@pytest.mark.parametrize("ideal,word", ((gm.ideal_hopping, gm.source_hopping_word),
                                      (gm.ideal_electric, gm.source_electric_word),
                                      (gm.ideal_echo, gm.source_echo_word)))
def test_actual_m6_word_matches_all_columns_and_channel_matrix_units(ideal, word):
    actual = compiler.execute(word(), 2)
    target = gm.numerical_matrix(ideal())
    assert np.linalg.norm(actual-target, 2) < 1e-13
    for row in range(4):
        for col in range(4):
            assert np.linalg.norm(np.outer(actual[:, row], actual[:, col].conj())
                                  - np.outer(target[:, row], target[:, col].conj()), 2) < 1e-13
    # An untouched quantum reference tests coherent full-code preservation.
    entangled = np.eye(4, dtype=complex)/2
    assert np.linalg.norm(compiler.execute(word(), 2, entangled)-target@entangled, 2) < 1e-13


def test_remapped_word_preserves_scalar_slots_and_bit_order():
    actual = compiler.execute(compiler.remap(gm.source_echo_word(), [0, 2]), 4)
    expected = np.zeros((16, 16), complex)
    u = gm.numerical_matrix(gm.ideal_echo())
    for source in range(16):
        logical = (source & 1)+2*((source >> 2) & 1)
        spectator = source & 0b1010
        for destination in range(4):
            output = spectator+(destination & 1)+4*((destination >> 1) & 1)
            expected[output, source] = u[destination, logical]
    assert np.linalg.norm(actual-expected, 2) < 1e-13


@pytest.mark.parametrize("duration", (0., .1, -.2))
def test_general_symmetric_step_uses_same_declared_hamiltonian(duration):
    h = gm.numerical_matrix(gm.hopping_hamiltonian())
    v = gm.numerical_matrix(gm.electric_energy())
    actual = compiler.execute(gm.source_step_word(duration), 2)
    half = expm(-.5j*duration*h)
    expected = half@expm(-1j*duration*v)@half
    assert np.linalg.norm(actual-expected, 2) < 1e-13
    assert np.array_equal(gm.numerical_matrix(gm.hamiltonian()), h+v)


def test_audit_full_channel_resources_and_source_pins():
    result = gm.audit()
    assert result["exact_checks"]["CAR_cases"] == 512
    assert result["exact_checks"]["physical_code_columns"] == 4
    assert result["native_all_16_matrix_units_channel_residual"] < 1e-13
    assert result["resources"]["native_pulses"] == 22
    assert result["resources"]["proper_code_service_events"] == 44
    assert result["resources"]["service_failure_outcomes_retained"]
    assert len(result["source_sha256"]) == 3
    assert not result["continuum_limit_or_physical_identification_supplied"]


@pytest.mark.parametrize("call", (lambda: gm.car_ladder(True, 0, False),
                                 lambda: gm.car_ladder(0, 4, False),
                                 lambda: gm.car_ladder(0, 0, 1),
                                 lambda: gm.occupation_projector(True, 0),
                                 lambda: gm.occupation_projector(0, 2),
                                 lambda: gm.gauss_values(0, 1.5),
                                 lambda: gm.dressed_hopping_images(0, 0, -1),
                                 lambda: gm.branch_probability(gm.ideal_echo(), 4, 0, 0),
                                 lambda: gm.source_step_word(True),
                                 lambda: gm.source_step_word(float("nan"))))
def test_invalid_physical_indices_rejected(call):
    with pytest.raises(ValueError):
        call()
