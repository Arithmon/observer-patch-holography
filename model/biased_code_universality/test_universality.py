"""Independent endpoint effects, arbitrary complex codes and clock joins."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations

import numpy as np
import pytest

from . import universality as u


PAIRS = tuple(combinations(range(12), 2))
INDEX = {pair: i for i, pair in enumerate(PAIRS)}


def _generic_code(core, rng):
    # Balanced pairs of unrelated complex frames give unequal exterior
    # weights while satisfying sum w U* Z U=0 exactly in the construction.
    safe = [j for j in range(12) if j not in core and all(
        tuple(sorted((j, a))) not in u.transport.SEAMS for a in core)]
    weights = [F(1, 2)] * 2 if len(safe) == 2 else [F(1, 6)] * 2 + [F(1, 3)] * 2
    answer = np.zeros((66, 2), complex)
    swap = np.array([[0., 1.], [1., 0.]])
    for n in range(0, len(safe), 2):
        matrix = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        frame, _ = np.linalg.qr(matrix)
        for j, weight, current in zip(safe[n:n+2], weights[n:n+2],
                                      (frame, swap @ frame), strict=True):
            for a, row in zip(core, current, strict=True):
                answer[INDEX[tuple(sorted((a, j)))]] = np.sqrt(float(weight)) * row
    return answer


def _raw_domains(law, seam, coin):
    # Endpoint arithmetic only: no producer branch, code or probability map.
    result = {raw: [] for raw in range(4 if law == "C" else 8)}
    a, b = seam
    for i, pair in enumerate(PAIRS):
        left, right = int(a in pair), int(b in pair)
        parity = (left + right) % 2
        after_left = (left + right) // 2 + coin * parity
        current = left - after_left
        raw = {0: 0, 1: 1, -1: 2}[current]
        if law == "E":
            raw += 4 * parity
        result[raw].append(i)
    return result


@pytest.mark.parametrize("parameters", [
    u.transport.parameters(),
    u.transport.local.Parameters(),
    u.transport.local.Parameters((2, 3, 4) * 4, tuple(1 + i % 4 for i in range(30)),
                                 F(3, 7), F(5, 8)),
])
def test_universal_formula_matches_both_complete_native_laws(parameters):
    derived = u.kernel(parameters)
    for law in ("C", "E"):
        assert derived == u.transport.kernel(parameters, law)
    assert all(sum(row) == 1 and min(row) >= 0 for row in derived)


def test_generic_complex_codes_have_scalar_complete_effects_and_the_same_biased_walk():
    rng = np.random.default_rng(20261029)
    p = u.transport.parameters()
    expected = u.kernel(p)
    complete_slots = 0
    domains = {(law, seam, coin): _raw_domains(law, seam, coin)
               for law in ("C", "E") for seam in u.transport.SEAMS for coin in (0, 1)}
    for i, core in enumerate(u.transport.LIVE):
        code = _generic_code(core, rng)
        assert np.linalg.norm(code.conj().T @ code - np.eye(2)) < 2e-14
        for a in core:
            rows = code[[j for j, pair in enumerate(PAIRS) if a in pair]]
            assert np.linalg.norm(rows.conj().T @ rows - np.eye(2) / 2) < 2e-14
        for law in ("C", "E"):
            moved = Counter()
            for seam in u.transport.SEAMS:
                for coin in (0, 1):
                    for raw, indices in domains[law, seam, coin].items():
                        complete_slots += 1
                        rows = code[indices]
                        gram = rows.conj().T @ rows
                        probability = float(gram.trace().real / 2)
                        assert np.linalg.norm(gram - probability * np.eye(2)) < 2e-14
                        current = {0: 0, 1: 1, 2: -1}.get(raw % 4, 0)
                        a, b = seam
                        target = tuple(sorted(b if k == a else a if k == b else k
                                              for k in core)) if current else core
                        coin_weight = F(p.port_weights[a if coin else b],
                                        p.port_weights[a] + p.port_weights[b])
                        moved[target] += probability * float(coin_weight) / 30
            assert max(abs(moved[c] - float(value)) for c, value in
                       zip(u.transport.CORES, expected[i], strict=True)) < 3e-15
    assert complete_slots == 43200


def test_all66_cores_and_the_arbitrary_code_reset_bound():
    assert len(u.transport.LIVE) == 60 and len(u.transport.STOP) == 6
    for core in u.transport.STOP:
        assert not u.transport.screening.safe_exterior(core)
    certificate = u.reset_scope_certificate()
    assert certificate["maximizing_reset_cores"] == [[0, 1]]
    assert certificate["all_complete_C_E_core_kernel_entries_checked"] == 7920
    assert certificate["adaptive_arbitrary_exact_source_qubit_codes_allowed"]
    assert not certificate["all_finite_horizon_biased_optimality_proved"]
    assert not certificate["physical_clock_energy_mass_or_particle_identification"]
    _, _, _, _, h, maximum, _ = u.reset.exact_chain()
    assert all(after - before <= maximum for before in (*h, F()) for after in h)


def test_sharp_charge_without_exact_recovery_does_not_get_the_universal_kernel():
    # A sharp orthogonal code with unequal endpoint weights is inadmissible.
    rng = np.random.default_rng(1153)
    core = (0, 1)
    code = _generic_code(core, rng)
    code = code * np.sqrt([1.2 if core[0] in pair else .8 for pair in PAIRS])[:, None]
    # Normalizing the logical Gram still leaves at least one nonscalar flag.
    values, vectors = np.linalg.eigh(code.conj().T @ code)
    code = code @ vectors @ np.diag(values ** -.5) @ vectors.conj().T
    residual = 0.
    for seam in u.transport.SEAMS:
        for indices in _raw_domains("C", seam, 0).values():
            rows = code[indices]
            gram = rows.conj().T @ rows
            residual = max(residual, np.linalg.norm(gram - np.trace(gram) * np.eye(2) / 2))
    assert residual > .01


def test_complete_source_memory_branches_compose_with_a_separate_clock_factor():
    # A supplied nontrivial clock instrument, not a native-clock preparation
    # or a physical calibration. The input entangles logical, clock and R.
    rng = np.random.default_rng(786)
    logical = rng.normal(size=(2, 2, 3)) + 1j * rng.normal(size=(2, 2, 3))
    logical /= np.linalg.norm(logical)
    phase = np.diag(np.exp(1j * .73 * np.log([25 / 33, 8 / 33])))
    clocks = [np.diag(np.sqrt([.3, .7])) @ phase,
              np.diag(np.sqrt([.7, .3])) @ phase]
    assert np.linalg.norm(sum(op.conj().T @ op for op in clocks) - np.eye(2)) < 1e-14
    complete_joint_outputs = resets = 0
    for core, law, seam, _, _ in u.transport.witness_path():
        source = u.transport.encode(core, logical.reshape(2, 6))
        branches = u.reset.source_branches(core, law, seam, source)
        total_mass = 0.
        for branch in branches:
            resets += branch["reset_performed"]
            for op in clocks:
                complete_joint_outputs += 1
                actual = np.einsum("ij,sjr->sir", op, branch["amplitude"].reshape(512, 2, 3))
                target_logical = np.einsum("ij,ljr->lir", op, logical)
                if branch["probability"]:
                    expected = (np.sqrt(float(branch["probability"])) *
                                (u.transport.ideal_code(branch["core_after"]) @
                                 target_logical.reshape(2, 6))).reshape(512, 2, 3)
                else:
                    expected = np.zeros_like(actual)
                assert np.linalg.norm(actual - expected) < 1e-8
                total_mass += np.linalg.norm(actual) ** 2
        assert abs(total_mass - 1) < 1e-8
    assert complete_joint_outputs == 48
    assert resets > 0
