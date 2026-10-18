"""Controls for action positivity, phase freedom and the Hamiltonian distinction."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

import numpy as np
import pytest

from . import construction as c
from ..record_laws import family as f


@pytest.mark.parametrize("phase", tuple(c.PHASES.values()))
def test_actual_source_selective_maps_and_binary_kraus_equality(phase):
    # Selective equality on every diagonal matrix unit proves equality of all
    # finite adaptive classical histories; unread marginal equality is weaker.
    assert c.source_selective_equality(phase) == 156
    assert np.array_equal(c.binary_channel(c.phase_table(phase)),
                          c.group_channel(c.native_binary_groups(phase), 4))


def test_binary_classification_includes_noncoin_and_nontrivial_phase_laws():
    tables = [f.extract_table(build()) for build in f.INSTRUMENTS.values()]
    tables += [table for name, table in f.generic_tables().items() if name != "outside"]
    tables += [f.table(c.I, 0, 0, 0, 0), f.table(0, 0, c.I*F(1, 4), 0, 0)]
    for table in tables:
        cert, symmetric, skew = c.local_certificate(table)
        assert cert["symmetric_action_nonnegative"]
        assert np.array_equal(symmetric, c.adjoint(symmetric))
        assert np.array_equal(skew, -c.adjoint(skew))
        assert cert["self_adjoint_Liouville_action"] == all(table[key].is_real() for key in ("k", "g", "h"))
    assert not c.local_certificate(tables[-1])[0]["skew_is_Hamiltonian_commutator"]


@pytest.mark.parametrize("phase", tuple(c.PHASES.values()))
def test_exact_symmetric_form_for_complex_operator(phase):
    table = c.phase_table(phase)
    channel = c.binary_channel(table)
    _, symmetric, _ = c.local_certificate(table)
    vector = np.array([c.S(F(i-3, 7), 0, F(2*i+1, 11)) for i in range(16)], dtype=object)
    actual = sum(value.conj()*result for value, result in zip(vector, symmetric@vector, strict=True))
    difference = vector[5]-vector[10]
    expected = F(1, 2)*difference.conj()*difference
    for a, b in product(range(4), repeat=2):
        if a != b:
            index = 4*a+b
            expected += (1-channel[index, index].real())*vector[index].conj()*vector[index]
    assert actual == expected
    assert actual.sign() >= 0


@pytest.mark.parametrize("phase", tuple(c.PHASES.values()))
def test_hamiltonian_decomposition_with_overlapping_seams(phase):
    result = c.network_certificate(phase)
    assert result["matrix_units_checked"] == 64
    assert result["exact_commutator_equality"]
    assert result["two_site_mixed_difference"] == str(c.imaginary(phase))
    assert result["native_carrier_seams"] == [list(c.native.E), list(c.native.F)]


def test_skew_part_of_general_channel_is_not_always_a_hamiltonian():
    # Build a physical Kraus instrument independently of its table: E followed
    # by a phase -i on full-seam output. Its skew part violates the Leibniz rule.
    original = f.INSTRUMENTS["E"]()
    changed = {label: {x: (y, amplitude*(-c.I if f.key(f.STATES[y]) == (1, 1) else 1))
                       for x, (y, amplitude) in branch.items()}
               for label, branch in original.items()}
    table = f.extract_table(changed)
    assert table == f.table(c.I, 0, 0, 0, 0)
    assert f.completeness(changed) and f.positive_semidefinite(table)
    cert, _, skew = c.local_certificate(table)
    assert not cert["skew_is_Hamiltonian_commutator"]
    assert skew[1, 1]-skew[3, 3]-skew[13, 13] == -c.I


def test_purity_loss_is_independent_of_record_phase():
    for phase in c.PHASES.values():
        output = c.binary_channel(c.phase_table(phase))[:, 5].reshape(4, 4)
        assert sum((output@output)[i, i] for i in range(4)) == c.S(F(1, 2))
        assert output[1, 1] == output[2, 2] == c.S(F(1, 2))


def test_same_fixed_space_does_not_fix_action_or_interference():
    tables = [c.phase_table(phase) for phase in c.PHASES.values()]
    assert all(c.local_certificate(table)[0]["fixed_operator_dimension"] == 5 for table in tables)
    # The same independently specified fixed space: the even 2x2 block and a
    # scalar on the odd block. Check every spanning vector for each law.
    fixed = [(0,), (15,), (3,), (12,), (5, 10)]
    for table in tables:
        channel = c.binary_channel(table)
        for support in fixed:
            vector = np.array([c.S(index in support) for index in range(16)], dtype=object)
            assert np.array_equal(channel@vector, vector)
    zero, quarter, half, _ = tables
    assert 1-zero["g"].real() == F(1, 2)
    assert 1-quarter["g"].real() == 1
    assert 1-half["g"].real() == F(3, 2)
    probability = lambda table: f.plus_probability(
        f.apply_table(table, f.prepared("phase_pi_over_2", f.oi.Z)), f.oi.Z)
    assert probability(zero) == F(3, 8)
    assert probability(quarter) == F(1, 8)


@pytest.mark.parametrize("phase", (0, F(1, 2), 2))
def test_nonunit_phases_rejected(phase):
    with pytest.raises(ValueError):
        c.phase_instrument(phase)


def test_nonpositive_table_rejected_before_action_construction():
    with pytest.raises(ValueError):
        c.binary_channel(f.generic_tables()["outside"])


def test_retained_receipt_replays_exactly():
    retained = json.loads(Path(__file__).with_name("results.json").read_text())
    assert c.run() == retained
