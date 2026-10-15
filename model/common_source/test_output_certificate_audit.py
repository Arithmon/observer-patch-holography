"""Independent controls for the exact comparison of retained numerical outputs.

Fixtures start from explicit Fock preparations and one Ramsey read. They do
not reuse the comparator's embedding, inner products, or source dispatcher.
The low-level comparison trusts a locally computed theorem enclosure; the
canonical producer owns replay and external proof provenance.
"""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import exact_output as exact
from . import protocol, smooth_history as smooth


FAST = smooth.Resolution(cutoff=4, preparation=2, degree=16,
                         steps_per_unit=2, bits=96)


@pytest.fixture(scope='module')
def vacuum():
    program = protocol.Program(dimension=2)
    enclosure = smooth.enclose(program, dimension=2, resolution=FAST)
    vector = np.zeros(16, complex)
    vector[0] = 1
    return program, enclosure, {(): vector}


@pytest.fixture(scope='module')
def ramsey():
    program = protocol.Program(operations=(protocol.Operation('smooth_read'),))
    enclosure = smooth.enclose(program, dimension=2, resolution=FAST)
    angle = 1/math.sqrt(8)
    vectors = {}
    for bit in (0, 1):
        vector = np.zeros(16, complex)
        phase = (1+1j*(2*bit-1))/2
        vector[0] = phase*math.cos(angle)
        vector[4] = phase*((-1)**bit)*math.sin(angle)
        vectors[(bit,)] = vector
    return program, enclosure, vectors


def _density(vectors):
    return {record: np.outer(vector, vector.conj())
            for record, vector in vectors.items()}


def _compare(program, enclosure, vectors, density=False):
    if density:
        return exact.compare_densities(program, _density(vectors), enclosure)
    return exact.compare_vectors(program, vectors, enclosure)


@pytest.mark.parametrize('density', [False, True])
def test_raw_mass_is_retained_and_distance_is_not_capped_at_one(vacuum, density):
    program, enclosure, vectors = vacuum
    result = _compare(program, enclosure, {(): 2*vectors[()]}, density)
    assert result['source_mass_real'] == 4
    assert result['source_mass_imag'] == 0
    # 1/2 ||4|0><0| - |0><0|||_1 = 3/2, not a normalized-state distance.
    assert result['distance_lower'] <= F(3, 2) <= result['distance_upper']
    assert result['distance_upper'] >= F(3, 2)
    assert result['distance_lower'] > 1


@pytest.mark.parametrize('density', [False, True])
def test_zero_output_is_not_silently_replaced_or_normalized(vacuum, density):
    program, enclosure, vectors = vacuum
    result = _compare(program, enclosure, {(): 0*vectors[()]}, density)
    assert result['source_mass_real'] == 0
    assert result['distance_lower'] <= F(1, 2) <= result['distance_upper']
    assert result['distance_lower'] > F(49, 100)


@pytest.mark.parametrize('density', [False, True])
def test_independent_ramsey_output_and_branch_normalization_control(ramsey, density):
    program, enclosure, vectors = ramsey
    correct = _compare(program, enclosure, vectors, density)
    assert correct['distance_upper'] < F(1, 10**10)
    wrong = {record: vector/math.sqrt(float(np.vdot(vector, vector).real))
             for record, vector in vectors.items()}
    changed = _compare(program, enclosure, wrong, density)
    assert changed['distance_lower'] > F(49, 100)
    assert changed['source_mass_real'] > F(199, 100)


@pytest.mark.parametrize('phase', [1, -1, 1j, -1j])
@pytest.mark.parametrize('density', [False, True])
def test_complete_record_phase_is_harmless(ramsey, phase, density):
    program, enclosure, vectors = ramsey
    changed = {record: phase*vector if record == (0,) else vector.copy()
               for record, vector in vectors.items()}
    result = _compare(program, enclosure, changed, density)
    assert result['distance_upper'] < F(1, 10**10)


@pytest.mark.parametrize('mutation', ['scalar_phase', 'swapped_records'])
@pytest.mark.parametrize('density', [False, True])
def test_unchanged_record_probabilities_do_not_hide_wrong_quantum_output(ramsey, mutation, density):
    program, enclosure, vectors = ramsey
    if mutation == 'scalar_phase':
        wrong = {record: vector.copy() for record, vector in vectors.items()}
        for vector in wrong.values():
            vector[4] *= -1
    else:
        wrong = {(0,): vectors[(1,)].copy(), (1,): vectors[(0,)].copy()}
    assert all(np.vdot(wrong[key], wrong[key]) == np.vdot(value, value)
               for key, value in vectors.items())
    result = _compare(program, enclosure, wrong, density)
    # The general matrix comparator uses a Frobenius lower bound, while
    # vector inputs admit the tighter exact rank-two trace-norm formula.
    assert result['distance_lower'] > (F(2, 5) if density else F(3, 5))


@pytest.mark.parametrize('density', [False, True])
def test_factorial_metric_and_gauge_fast_index_are_both_required(density):
    program = protocol.Program(dimension=4)
    enclosure = smooth.enclose(program, dimension=4, resolution=FAST)
    # (|00> - sqrt(2)|20>/36 + 2|11>/36 - sqrt(2)|02>/36)
    # with a blank four-dimensional gauge factor and norm squared 163/162.
    vector = np.zeros(64, complex)
    vector[0] = 1
    vector[8] = vector[32] = -math.sqrt(2)/36
    vector[20] = F(1, 18)
    vector /= math.sqrt(F(163, 162))
    correct = _compare(program, enclosure, {(): vector}, density)
    assert correct['distance_upper'] < F(1, 10**10)
    missing_factorial = vector.copy()
    missing_factorial[8] /= math.sqrt(2)
    missing_factorial[32] /= math.sqrt(2)
    wrong_gauge_order = vector.reshape(4, 4, 4).transpose(2, 0, 1).ravel()
    for wrong in (missing_factorial, wrong_gauge_order):
        result = _compare(program, enclosure, {(): wrong}, density)
        assert result['distance_lower'] > F(1, 100)


@pytest.mark.parametrize('density', [False, True])
def test_even_zero_outcomes_cannot_be_omitted(density):
    program = protocol.Program(operations=(protocol.Operation('gauge_read'),))
    enclosure = smooth.enclose(program, dimension=2, resolution=FAST)
    vector = np.zeros(16, complex)
    vector[0] = 1
    with pytest.raises((ValueError, ArithmeticError)):
        _compare(program, enclosure, {(0,): vector}, density)


@pytest.mark.parametrize('mutation', ['missing', 'extra', 'boolean_record'])
@pytest.mark.parametrize('density', [False, True])
def test_record_menu_is_strict(ramsey, mutation, density):
    program, enclosure, vectors = ramsey
    changed = deepcopy(vectors)
    if mutation == 'missing':
        changed.pop((1,))
    elif mutation == 'extra':
        changed[(2,)] = np.zeros(16, complex)
    else:
        changed = {(False,): changed[(0,)], (1,): changed[(1,)]}
    with pytest.raises((ValueError, ArithmeticError)):
        _compare(program, enclosure, changed, density)


@pytest.mark.parametrize('bad', ['nan', 'inf', 'bool', 'object', 'masked', 'shape'])
@pytest.mark.parametrize('density', [False, True])
def test_original_arrays_are_validated_before_rational_conversion(vacuum, bad, density):
    program, enclosure, vectors = vacuum
    value = _density(vectors)[()] if density else vectors[()].copy()
    if bad in ('nan', 'inf'):
        value.flat[0] = float(bad)
    elif bad == 'bool':
        value = value.astype(bool)
    elif bad == 'object':
        value = value.astype(object)
    elif bad == 'masked':
        value = np.ma.array(value, mask=np.zeros(value.shape, bool))
    else:
        value = value.reshape(-1) if density else value.reshape(16, 1)
    method = exact.compare_densities if density else exact.compare_vectors
    with pytest.raises((ValueError, ArithmeticError, TypeError)):
        method(program, {(): value}, enclosure)


@pytest.mark.parametrize('field', ['normalization', 'vector_error', 'preparation_error', 'branches'])
def test_unresealed_theorem_mutation_is_not_accepted(vacuum, field):
    program, enclosure, vectors = vacuum
    changed = deepcopy(enclosure)
    if field == 'branches':
        changed.branches[()].coefficients[(0, 0, 0)] = (2, 0)
    else:
        changed = replace(changed, **{field: getattr(changed, field)+F(1, 3)})
    with pytest.raises((ValueError, ArithmeticError)):
        exact.compare_vectors(program, vectors, changed)


def test_matching_shapes_do_not_erase_the_declared_current(vacuum):
    program, enclosure, vectors = vacuum
    with pytest.raises((ValueError, ArithmeticError)):
        exact.compare_vectors(replace(program, native_current=1), vectors, enclosure)


def test_exact_mass_keeps_an_amplitude_whose_square_underflows_binary64(vacuum):
    program, enclosure, vectors = vacuum
    tiny = float.fromhex('0x1p-600')
    vector = vectors[()].copy()
    vector[1] = tiny
    assert tiny*tiny == 0
    result = exact.compare_vectors(program, {(): vector}, enclosure)
    assert result['source_mass_real'] == 1+F(tiny)**2
    assert result['source_mass_real'] > 1
    assert result['distance_lower'] > 0
    assert result['distance_upper'] < F(1, 10**100)


def test_comparison_does_not_mutate_inputs_or_consume_randomness(ramsey):
    program, enclosure, vectors = ramsey
    before = {record: vector.tobytes() for record, vector in vectors.items()}
    state = np.random.get_state()
    first = exact.compare_vectors(program, vectors, enclosure)
    second = exact.compare_vectors(program, vectors, enclosure)
    assert first == second
    assert before == {record: vector.tobytes() for record, vector in vectors.items()}
    following = np.random.get_state()
    assert state[0] == following[0] and np.array_equal(state[1], following[1])
    assert state[2:] == following[2:]
