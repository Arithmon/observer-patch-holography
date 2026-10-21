"""Larger source words preserve signed amplitudes, phase and fixed schedules."""
import math

import numpy as np
import pytest

from . import compiler as c, field


@pytest.mark.parametrize("bits", (7, 8))
@pytest.mark.parametrize("kind", ("dense", "sparse", "negative_blank", "zero_prefix"))
def test_tree_preparation_against_independently_normalized_signed_vectors(bits, kind):
    size = 1 << bits
    if kind == "dense":
        integers = np.array([(-1 if index % 3 else 1)*(index % 7+1) for index in range(size)])
    else:
        integers = np.zeros(size, dtype=int)
        if kind == "sparse":
            integers[0], integers[size//2+1], integers[-1] = -3, 4, -12
        elif kind == "negative_blank":
            integers[0] = -1
        else:
            integers[size//2], integers[-2], integers[-1] = -2, -3, 6
    length = math.sqrt(sum(int(value)**2 for value in integers))
    target = np.array([int(value)/length for value in integers])
    word = c.compile_real_state(target)
    blank = np.zeros(size)
    blank[0] = 1
    actual = c.execute(word, bits, blank)
    # Vector equality includes global and relative phases, including negative
    # leaves. Comparing only probabilities or density matrices would miss them.
    np.testing.assert_allclose(actual, target, atol=4e-12, rtol=0)
    assert len(word) == 3*((bits-1)*(1 << bits)+1)
    inverse = c.execute(c.inverse(word), bits, actual)
    np.testing.assert_allclose(inverse, blank, atol=8e-12, rtol=0)


def test_tree_word_is_unitary_on_all_inputs_and_passes_unchanged_roundoff_guard():
    target, _ = field.vacuum_vector(16)
    word = c.compile_real_state(target)
    executed = c.execute(word, 8)
    restored, diagnostics = c.stabilize_executed_unitary(executed)
    assert diagnostics["raw_unitarity_defect"] <= 1e-10
    assert diagnostics["projection_change"] <= 5.1e-11
    np.testing.assert_allclose(executed[:, 0], target, atol=4e-12, rtol=0)
    np.testing.assert_allclose(restored[:, 0], target, atol=4e-12, rtol=0)
    np.testing.assert_allclose(restored.conj().T@restored, np.eye(256), atol=2e-13, rtol=0)


def test_tree_schedule_charges_zero_angles_and_preserves_complete_source_phase():
    blank = np.zeros(128)
    blank[0] = 1
    dense = np.arange(1, 129, dtype=float)
    dense /= np.linalg.norm(dense)
    empty_word, full_word = c.compile_real_state(blank), c.compile_real_state(dense)
    assert [(row["kind"], row["targets"]) for row in empty_word] == [(row["kind"], row["targets"]) for row in full_word]
    assert len(empty_word) == len(full_word) == 2307
    assert sum(row["angle"] == 0 for row in empty_word) == 127
    identity = c.execute(empty_word, 7)
    np.testing.assert_allclose(identity, np.eye(128), atol=4e-12, rtol=0)


def test_four_qubit_basis_word_retains_every_column_and_inverse_phase():
    # Exact finite Walsh values, with an independently supplied column sign.
    basis = np.array([[(-1)**((row & column).bit_count())/4
                       for column in range(16)] for row in range(16)])
    basis[:, 5] *= -1
    word = c.compile_real_orthogonal(basis)
    actual = c.execute(word, 4)
    np.testing.assert_allclose(actual, basis, atol=8e-12, rtol=0)
    np.testing.assert_allclose(c.execute(c.inverse(word), 4)@actual, np.eye(16), atol=2e-11, rtol=0)


def test_expanded_capacity_keeps_a_declared_upper_limit():
    beyond = np.zeros(8192)
    beyond[0] = 1
    with pytest.raises(ValueError):
        c.compile_real_state(beyond)
    with pytest.raises(ValueError):
        c.compile_real_orthogonal(np.eye(32))
