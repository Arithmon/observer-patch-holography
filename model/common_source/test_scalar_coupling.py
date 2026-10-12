"""Independent carrier, interaction, native execution and refinement checks."""
from fractions import Fraction as F
from decimal import Decimal, localcontext

import numpy as np
import pytest
from scipy.linalg import expm

from . import scalar_coupling as s
from . import gauge_matter as g
from ..observer.source_scalar import compiler as c


def independent_field(d):
    """Infinite-oscillator matrix elements projected before multiplying."""
    q = np.zeros((d, d))
    q2, p2 = np.diag(np.arange(d) + .5), np.diag(np.arange(d) + .5)
    for n in range(d - 1):
        q[n, n + 1] = q[n + 1, n] = np.sqrt((n + 1) / 2)
    for n in range(d - 2):
        q2[n, n + 2] = q2[n + 2, n] = np.sqrt((n + 1) * (n + 2)) / 2
        p2[n, n + 2] = p2[n + 2, n] = -q2[n, n + 2]
    return q, q2, p2


def independent_terms(d, coupling=.25):
    q, q2, p2 = independent_field(d)
    x, z, i = np.array([[0., 1.], [1., 0.]]), np.diag([1., -1.]), np.eye(2)
    hop = np.kron(i, x) + np.kron(x, i)
    electric = np.diag([0., 1., 1., 0.])
    onsite = (p2 + (1 + 9 / 32) * q2) / 2
    return {
        "hopping": np.kron(np.eye(d * d), hop),
        "electric": np.kron(np.eye(d * d), electric),
        "onsite": np.kron(np.kron(np.eye(d), onsite) + np.kron(onsite, np.eye(d)), np.eye(4)),
        "bond": -9 / 32 * np.kron(np.kron(q, q), np.eye(4)),
        "interaction": coupling * np.sqrt(2) * np.kron(np.kron(np.eye(d), q), np.kron(i, z)),
    }


@pytest.mark.parametrize("d", (2, 4))
def test_common_target_uses_projected_physical_squares_and_correct_charge_sign(d):
    independent = independent_terms(d)
    for name, actual in s.hamiltonian_terms(d).items():
        assert np.linalg.norm(actual - independent[name], 2) < 3e-15
    q, q2, _ = independent_field(d)
    assert np.linalg.norm(q2 - q @ q, 2) == pytest.approx(d / 2)
    # A scalar coordinate transition with gauge bit zero has +g coefficient;
    # gauge bit one has the opposite sign. This catches a left/right reversal.
    interaction = s.hamiltonian_terms(d)["interaction"]
    assert interaction[4, 0] == pytest.approx(.25)
    assert interaction[5, 1] == pytest.approx(-.25)


@pytest.mark.parametrize("d", (2, 4))
def test_scalar_charge_word_executes_complete_noncommuting_interaction(d):
    dt = .31
    word = s.scalar_subwords(d, dt)["interaction"]
    actual = c.execute(word, s.qubits(d))
    term = independent_terms(d)["interaction"]
    assert np.linalg.norm(actual - expm(-1j * dt * term), 2) < 5e-13
    hopping = independent_terms(d)["hopping"]
    assert np.linalg.norm(term @ hopping - hopping @ term, 2) > .4
    assert all(set(event) == {"kind", "targets", "angle"} for event in word)


@pytest.mark.parametrize("d", (2, 4))
def test_same_native_register_matches_independent_split_on_all_columns(d):
    actual, metadata = s.source_unitary(d, 16)
    expected = s.independent_split_unitary(d, 16)
    assert np.linalg.norm(actual - expected, 2) < 4e-12
    assert np.linalg.norm(actual.conj().T @ actual - np.eye(len(actual)), 2) < 3e-12
    assert metadata["resources"]["native_pulses"] == 16 * len(s.step_word(d, 1 / 16))
    assert metadata["resources"]["proper_code_service_events"] == 2 * metadata["resources"]["native_pulses"]
    # Arbitrary entangled reference columns use this same actual unitary.
    rng = np.random.default_rng(173 + d)
    columns = rng.normal(size=(4 * d * d, 3)) + 1j * rng.normal(size=(4 * d * d, 3))
    columns /= np.linalg.norm(columns)
    assert np.linalg.norm(actual @ columns - expected @ columns) < 4e-12


@pytest.mark.parametrize("d", (2, 4))
def test_source_preparation_and_coupling_make_joint_state_not_tensor_probabilities(d):
    initial, metadata = s.preparation(d)
    assert np.linalg.norm(initial) == pytest.approx(1, abs=2e-12)
    assert np.linalg.norm(initial.reshape(d * d, 4)[:, 1:]) < 2e-13
    # Independently project the omega=5/4 difference-mode squeezed vacuum
    # through total occupation two; d=2 retains only its blank component.
    target = np.zeros(d * d)
    target[0] = 1
    if d == 4:
        target[2] = target[2 * d] = -1 / (18 * np.sqrt(2))
        target[1 + d] = 1 / 18
    target /= np.linalg.norm(target)
    assert np.linalg.norm(initial - np.kron(target, [1, 0, 0, 0])) < 3e-12
    assert metadata["resources"]["native_pulses"] == len(s.initialpreparation_word(d))
    purities = []
    for coupling in (0, .25):
        evolution, _ = s.source_unitary(d, 16, coupling=coupling)
        state = (evolution @ initial).reshape(d * d, 4)
        reduced = state.conj().T @ state
        purities.append(float(np.trace(reduced @ reduced).real))
    assert purities[0] == pytest.approx(1, abs=3e-12)
    assert purities[1] < .96


def _ambient_gauge():
    """Independent Jordan-Wigner/CAR construction with five integer fluxes."""
    ladders = []
    for mode in range(4):
        a = np.zeros((16, 16))
        for incoming in range(16):
            if incoming >> mode & 1:
                a[incoming ^ (1 << mode), incoming] = (-1) ** ((incoming & ((1 << mode) - 1)).bit_count())
        ladders.append(a)
    up = np.diag(np.ones(4), -1)
    down = up.T
    hopping = np.zeros((80, 80))
    for channel, charge in enumerate((1, -1)):
        left, right = ladders[2 * channel:2 * channel + 2]
        shift = up if charge == 1 else down
        hopping += np.kron(shift, left.T @ right) + np.kron(shift.T, right.T @ left)
    n = [a.T @ a for a in ladders]
    e = np.kron(np.diag(np.arange(-2, 3)), np.eye(16))
    gl = e - np.kron(np.eye(5), n[0] - n[2])
    gr = -e - np.kron(np.eye(5), n[1] - n[3])
    charge_difference = np.kron(np.eye(5), n[0] - n[1])
    v = np.zeros((80, 4))
    for index in range(4):
        b0, b1 = index & 1, index >> 1
        occupation, flux = (1 << b0) + (1 << (2 + b1)), -b0 + b1
        v[(flux + 2) * 16 + occupation, index] = 1
    return v, hopping + e @ e, gl, gr, charge_difference


def test_common_interaction_preserves_operator_gauss_and_car_parity():
    v, ambient_h, gl, gr, difference = _ambient_gauge()
    assert np.array_equal(gl @ v, np.zeros_like(v))
    assert np.array_equal(gr @ v, np.zeros_like(v))
    assert np.array_equal(ambient_h @ v, v @ g.numerical_matrix(g.hamiltonian()))
    z0 = np.diag([1., -1., 1., -1.])
    assert np.array_equal(difference @ v, v @ z0)
    assert np.array_equal(gl @ difference - difference @ gl, np.zeros_like(gl))
    assert np.array_equal(gr @ difference - difference @ gr, np.zeros_like(gr))
    parity = np.kron(np.eye(5), np.diag([(-1) ** n.bit_count() for n in range(16)]))
    assert np.array_equal(parity @ v, v)
    state, _ = s.source_evolve(2, 16)
    # No branch in the fully evolved entangled state leaves the exact Gauss
    # sector or changes fermion parity. Flux values are not independently reset.
    ambient_columns = v @ state.reshape(4, 4).T
    assert np.linalg.norm(gl @ ambient_columns) == 0
    assert np.linalg.norm(gr @ ambient_columns) == 0
    assert np.linalg.norm(parity @ ambient_columns - ambient_columns) == 0


@pytest.mark.parametrize("d", (2, 4))
def test_time_refinement_compares_same_continuous_common_hamiltonian(d):
    target = expm(-1j * sum(independent_terms(d).values()))
    errors = [np.linalg.norm(s.independent_split_unitary(d, n) - target, 2) for n in (8, 16, 32)]
    assert 3.9 < errors[0] / errors[1] < 4.1
    assert 3.9 < errors[1] / errors[2] < 4.1
    for n, error in zip((8, 16, 32), errors, strict=True):
        assert error < float(F(s.time_error_bound(d, n)["ideal_split_operator_error_upper"]))
    bound = [F(s.time_error_bound(d, n)["unclipped_operator_error_upper"]) for n in (128, 256, 512)]
    assert bound[0] == 4 * bound[1] == 16 * bound[2]


@pytest.mark.parametrize("d", (2, 4))
def test_exact_norm_budget_and_rational_root_bound(d):
    certificate = s.time_error_bound(d, 512)
    bound = F(certificate["sum_generator_norms_upper"])
    assert sum(np.linalg.norm(h, 2) for h in independent_terms(d).values()) < float(bound)
    with localcontext() as ctx:
        ctx.prec = 80
        exact_s = Decimal(3) + Decimal(25) / 16 * (2 * d - 1) + Decimal(4 * d - 2).sqrt() / 4
        rational_s = Decimal(bound.numerator) / Decimal(bound.denominator)
        assert rational_s >= exact_s
        assert rational_s - exact_s < Decimal("1e-15")
    assert not certificate["spatial_or_occupation_error_included"]
    assert not certificate["floating_execution_or_physical_error_included"]


def test_high_precision_d2_palindrome_against_exact_rational_hamiltonian():
    import mpmath as mp
    with mp.workdps(60):
        terms = {key: mp.zeros(16) for key in ("hopping", "electric", "onsite", "bond", "interaction")}
        for k in range(16):
            b0, b1, n0, n1 = [(k >> bit) & 1 for bit in range(4)]
            terms["hopping"][k ^ 1, k] = terms["hopping"][k ^ 2, k] = 1
            terms["electric"][k, k] = (b0 - b1)**2
            terms["onsite"][k, k] = mp.mpf(73) * (n0 + n1 + 1) / 64
            terms["bond"][k ^ 12, k] = -mp.mpf(9) / 64
            terms["interaction"][k ^ 4, k] = mp.mpf((-1)**b0) / 4
        dt = mp.mpf(1) / 32
        h = sum(terms.values(), mp.zeros(16))
        quarter = mp.expm(-mp.j * dt * terms["hopping"] / 4)
        gauge = quarter * mp.expm(-mp.j * dt * terms["electric"] / 2) * quarter
        onsite = mp.expm(-mp.j * dt * terms["onsite"] / 2)
        bond = mp.expm(-mp.j * dt * terms["bond"] / 2)
        interaction = mp.expm(-mp.j * dt * terms["interaction"])
        evolution = (gauge * onsite * bond * interaction * bond * onsite * gauge)**32
        difference = evolution - mp.expm(-mp.j * h)
        frobenius_error = mp.sqrt(sum(abs(x)**2 for x in difference))
        bound = F(s.time_error_bound(2, 32)["ideal_split_operator_error_upper"])
        assert frobenius_error < mp.mpf(bound.numerator) / bound.denominator
        numeric = np.array([[complex(evolution[i, j]) for j in range(16)] for i in range(16)])
        assert np.linalg.norm(numeric - s.independent_split_unitary(2, 32), 2) < 2e-14


def test_decoupled_word_keeps_same_service_schedule_and_global_phases():
    for d in (2, 4):
        actual = s.step_word(d, .1, .25)
        zero = s.step_word(d, .1, 0)
        assert [(e["kind"], e["targets"]) for e in actual] == [(e["kind"], e["targets"]) for e in zero]
        assert any(e["angle"] == 0 for e in zero)
    actual, _ = s.source_unitary(2, 2, time=0)
    assert np.linalg.norm(actual - np.eye(16), 2) < 1e-12


@pytest.mark.parametrize("arguments", ((True, 2), (3, 2), (2, True), (2, 0), (2, 1.5)))
def test_invalid_common_evolution_dimensions_and_steps_rejected(arguments):
    with pytest.raises(ValueError):
        s.source_unitary(*arguments)


@pytest.mark.parametrize("field,value", (("time", -.1), ("time", float("nan")), ("coupling", True), ("coupling", float("inf"))))
def test_invalid_model_parameters_rejected(field, value):
    with pytest.raises(ValueError):
        s.source_unitary(2, 2, **{field: value})


@pytest.mark.parametrize("field,value", (("time", 1.), ("time", True), ("coupling", .25), ("coupling", float("nan"))))
def test_bound_requires_exact_rational_inputs(field, value):
    with pytest.raises(ValueError):
        s.time_error_bound(2, 2, **{field: value})
