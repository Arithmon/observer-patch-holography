"""Adversarial checks for complete clocks, exact field memory and references."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as sp

from . import construction as c
from ..federated_quantum_sector import construction as q
from ..primitive_observer import runtime


@pytest.mark.parametrize('law', ('C', 'E'))
def test_clock_records_are_original_native_currents_for_every_basis_state(law):
    vertices, edges = c.geometry()
    expected = [sp.zeros(4) for _ in c.MARKS]
    for i, vertex in enumerate(vertices):
        expected = [sp.zeros(4) for _ in c.MARKS]
        for edge, coin in product(edges, (0, 1)):
            final, current, _ = q.native_image(1 << vertex, edge, coin)
            destination = vertices.index(final.bit_length()-1)
            mark = tuple(int(destination == v)-int(i == v) for v in (0, 2))
            expected[c.MARKS.index(mark)][destination, destination] += sp.Rational(1, 6)
            assert current == ((vertex == edge[0])-(vertices[destination] == edge[0]))
        initial = sp.zeros(4); initial[i, i] = 1
        assert tuple(c.unvector(m*c.vector(initial)) for m in c.source_maps(law)) == tuple(expected)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_every_quantum_matrix_unit_and_hidden_source_record_group(law):
    instrument = c.source_instrument(law)
    for i, j in product(range(4), repeat=2):
        unit = sp.zeros(4); unit[i, j] = 1
        actual = runtime._apply(instrument, c._dense(unit), 1)
        for output, exact in zip(actual, c.source_maps(law), strict=True):
            assert np.linalg.norm(output-c._dense(c.unvector(exact*c.vector(unit)))) < 1e-14
    # A no-event path includes genuine remote occupation transfer.
    initial = sp.zeros(4); initial[1, 1] = 1
    assert c.unvector(c.source_maps(law)[0]*c.vector(initial))[3, 3] == sp.Rational(1, 6)


@pytest.mark.parametrize('law', ('C', 'E'))
def test_complete_first_event_resolvent_and_every_finite_timeout(law):
    maps = c.source_maps(law)
    trace = sp.Matrix([[int(i//4 == i % 4) for i in range(16)]])
    for cap in range(9):
        completed = sum((sum(maps[1:], sp.zeros(16))*maps[0]**n for n in range(cap)), sp.zeros(16))
        assert trace*(completed+maps[0]**cap) == trace
        assert trace*maps[0]**cap == sp.Rational(5, 6)**cap*trace
    report = c.event_certificate(law)
    assert report['stationary_probability_consecutive_union_ticks_use_different_clocks'] == '1/6'
    assert report['largest_no_local_tick_probability_over_two_following_union_events'] == ['2/3']*2


def test_clock_waiting_normalization_cannot_erase_remote_motion():
    # Local event hazard is constant, but R/(1-q) is still the wrong event map.
    q0, mark = c.source_maps('C')[:2]
    assert q0 != sp.Rational(5, 6)*sp.eye(16)
    assert c.event_maps('C')[1] != 6*c.source_maps('C')[2]
    # Both candidates give the same completed tick instrument, not merely means.
    assert c.event_maps('C') == c.event_maps('E')


@pytest.mark.parametrize('law', ('C', 'E'))
def test_entangled_reference_event_assemblage_is_preserved(law):
    # Four-dimensional work entangled with a four-dimensional untouched reference.
    initial = sp.zeros(16)
    for i, j in product(range(4), repeat=2):
        initial[5*i, 5*j] = sp.Rational(1, 4)
    destinations = (1, 0, 3, 2)
    total_reference = sp.zeros(4)
    for matrix, destination in zip(c.event_maps(law), destinations, strict=True):
        actual = c._apply_exact(matrix, initial, 4)
        reference = sp.Matrix(4, 4, lambda r, s: sum(actual[4*i+r, 4*i+s] for i in range(4)))
        expected = sp.zeros(16)
        expected[4*destination:4*destination+4, 4*destination:4*destination+4] = reference
        assert actual == expected
        assert reference != sp.trace(reference)*sp.eye(4)/4
        total_reference += reference
    assert total_reference == sp.eye(4)/4


def test_equal_coarse_state_does_not_determine_clock_instrument():
    result = c.memory_witnesses()
    assert result['first_event_law_total_variation'] == '1'
    assert result['exact_word_probability'] == '5/72'
    assert result['memoryless_word_probability'] == '1/16'


def test_full_quantum_memory_identity_for_complex_units_and_non_native_cp_intervention():
    # A nonunitary amplitude-damping service, plus complex phases, tests more
    # than native invariant diagonal sectors or the production examples.
    k0 = sp.kronecker_product(sp.diag(1, sp.Rational(1, 2)), sp.eye(2))
    k1 = sp.kronecker_product(sp.Matrix([[0, 1], [0, 0]]), sp.eye(2))
    damping = c.superoperator(((F(1), k0.tolist()), (F(3, 4), k1.tolist())))
    phase = c.superoperator(((F(1), sp.diag(1, sp.I, -1, -sp.I).tolist()),))
    maps = (c.source_maps('C')[0], phase, damping, c.source_maps('E')[2])
    for i, j in product(range(4), repeat=2):
        initial = sp.zeros(4); initial[i, j] = 1+sp.I
        predicted = c.memory_forecast(initial, maps)
        fine = c.vector(initial)
        assert predicted[0] == c.reduced(initial)
        for step, matrix in enumerate(maps, start=1):
            fine = matrix*fine
            assert predicted[step] == c.reduced(c.unvector(fine))


@pytest.mark.parametrize('law,plus', (('C', '35/54'), ('E', '14/27')))
def test_actual_full_runtime_retains_no_event_quantum_information(law, plus):
    result = c.execute_experiment(law)
    assert result['complete_record_count'] == 250
    assert result['complete_field_plus_probability'] == plus
    assert result['first_event_timeout_probability'] == '125/216'
    assert result['second_event_timeout_probability'] == '25/27'
    assert result['full_state_and_checkpoint_residual_diagnostic'] < 1e-12
    # No-event quantum survival differs while the complete clock law agrees.
    no_event_plus = next(row['probability'] for row in result['complete_record_law']
                         if row['attempt_marks_and_readout'] == [0, 0, 0, 0])
    assert no_event_plus == ('7/16' if law == 'C' else '133/432')


@pytest.mark.parametrize('law', ('C', 'E'))
def test_live_record_control_changes_state_without_reset_and_all_references_survive(law):
    controlled = c.execute_experiment(law, feedback=True, cap=2)
    uncontrolled = c.execute_experiment(law, feedback=False, cap=2)
    reference = c.execute_experiment(law, reference=True, cap=2)
    assert controlled['complete_field_plus_probability'] != uncontrolled['complete_field_plus_probability']
    assert reference['complete_field_plus_probability'] == '1/2'
    assert reference['reference_dimension'] == 2
    assert reference['complete_record_count'] == 50
    assert reference['exact_memory_matches_every_reduced_terminal_matrix']


@pytest.mark.parametrize('bad', (True, 0, -1, 1.0, '2', 9))
def test_invalid_attempt_caps_are_not_silently_coerced(bad):
    with pytest.raises(ValueError):
        c.experiment('C', cap=bad)


def test_exactness_shape_and_controller_contracts():
    with pytest.raises(ValueError, match='cap at most3'):
        c.experiment('C', cap=4)
    for law in ('c', True, None):
        with pytest.raises(ValueError):
            c.source_maps(law)
    for mode in (True, 'inferred', None):
        with pytest.raises(ValueError):
            c.field_unitary(mode)
    for feedback in (1, 0, 'yes'):
        with pytest.raises(ValueError):
            c.experiment('C', feedback=feedback)
    for initial in (sp.zeros(3), sp.eye(4)*0.5, sp.eye(4)*sp.oo):
        with pytest.raises(ValueError):
            c.memory_forecast(initial, ())
    with pytest.raises(ValueError):
        c.memory_forecast(sp.eye(4)/4, (sp.eye(15),))
    with pytest.raises(ValueError):
        c.memory_forecast(sp.eye(4)/4, (sp.eye(16),)*9)
    with pytest.raises(ValueError):
        c.superoperator(())
