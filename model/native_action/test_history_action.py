"""Complete native histories, including adaptive access and a quantum reference."""
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest

from . import construction as native_action
from . import history_action as h
from ..primitive_observer import runtime


@pytest.mark.parametrize('law', ('C', 'E'))
def test_exact_completion_is_unitary_and_keeps_the_source_isometry(law):
    _, w, d = h.completion_blocks(law)
    assert np.array_equal(w.T@w, 2*np.eye(4, dtype=object))
    assert np.array_equal(d, d.T)
    assert np.array_equal(d@d, d)
    assert not np.any(d@w)
    # These four exact block identities prove U^2=I and U J=V for all inputs.
    assert np.array_equal(d@d+F(1, 2)*(w@w.T), np.eye(len(w), dtype=object))


@pytest.mark.parametrize('law', ('C', 'E'))
def test_all_input_matrix_units_preserve_adaptive_complete_histories(law):
    for left, right in product(range(8), repeat=2):
        assert h.history_unit(law, left, right) == h.direct_history_unit(law, left, right)
    matrices = [matrix for _, matrix in h.controlled_history_numerators(law)]
    assert np.array_equal(sum((m.T@m for m in matrices)), 4*np.eye(8, dtype=object))


@pytest.mark.parametrize('law', ('C', 'E'))
def test_actual_runtime_history_matches_dilation_with_entangled_reference(law):
    def instrument(seam):
        groups = native_action.full_groups(law, seam, ports=3)
        return native_action._instrument(groups, visible=True)

    identity = runtime.Instrument(kraus=((np.eye(8),),))
    phase = runtime.Instrument(kraus=((np.diag([(-1.)**(x & 1) for x in range(8)]),),))
    first, second, alternate = instrument((0, 1)), instrument((1, 2)), instrument((0, 2))
    psi = np.zeros(16, complex)
    psi[2], psi[5] = 1/np.sqrt(2), 1j/np.sqrt(2)
    rho = np.outer(psi, psi.conj())
    reduced = np.trace(rho.reshape(8, 2, 8, 2), axis1=1, axis2=3)
    program = runtime.Program('native_action_reference_history', 8, reduced, (
        runtime.Event('first_repair', instruments=(first,)),
        runtime.Event('local_phase_feedback', control=0, instruments=(identity, phase, identity)),
        runtime.Event('current_selected_seam', control=0, instruments=(second, second, alternate)),
    ))
    actual = runtime.execute(program, initial=rho)
    expected = {word: np.zeros((16, 16), complex) for word in product(range(3), repeat=2)}
    for labels, numerator in h.controlled_history_numerators(law):
        matrix = np.kron(np.asarray(numerator, dtype=float)/2, np.eye(2))
        word = tuple(label[1]+1 for label in labels)
        expected[word] += matrix@rho@matrix.conj().T
    assert len(actual.branches) == len(expected) == 9
    assert abs(sum(np.trace(value) for value in expected.values())-1) < 1e-14
    for branch in actual.branches:
        assert np.linalg.norm(branch.density-expected[branch.record], ord='nuc') < 1e-14
    prefix = runtime.execute(program, initial=rho, stop_after=1)
    resumed = runtime.execute(program, checkpoint=runtime.checkpoint_from_json(
        program, runtime.checkpoint_to_json(prefix)))
    assert runtime.complete_distance(actual, resumed) < 1e-14


@pytest.mark.parametrize('law', ('C', 'E'))
def test_same_native_endpoint_does_not_fix_intermediate_completion_readouts(law):
    _, w, d = h.completion_blocks(law)
    w = np.asarray(w, dtype=float)/np.sqrt(2)
    u = np.block([[np.asarray(d, dtype=float), w], [w.T, np.zeros((4, 4))]])
    p_plus, p_minus = (np.eye(len(u))+u)/2, (np.eye(len(u))-u)/2
    initial = np.zeros(len(u))
    initial[-3] = 1
    endpoints, intermediate = [], []
    for branch in (0, 2):
        # Both A_m=pi P_-+2pi m P_+ are positive; the displayed spectral
        # formula evaluates their exponential without a matrix-log routine.
        def evolve(parameter):
            return (np.exp(-1j*np.pi*parameter)*p_minus
                    +np.exp(-2j*np.pi*branch*parameter)*p_plus)@initial
        endpoints.append(evolve(1))
        intermediate.append(np.linalg.norm(evolve(F(1, 3))[-4:])**2)
    assert np.linalg.norm(endpoints[0]-endpoints[1]) < 1e-14
    assert np.linalg.norm(endpoints[0]-u@initial) < 1e-14
    assert abs(intermediate[0]-.75) < 1e-14
    assert abs(intermediate[1]) < 1e-14
    # This is a read of an auxiliary completion flag between native events.
    # No permitted native-boundary protocol can distinguish the two actions.


@pytest.mark.parametrize('args', (('X', (0, 1), 2), ('C', (1, 0), 2),
    ('C', (0, 2), 2), ('E', (0, 1), True), ('C', (0, 0), 2), ('C', [0, 1], 2)))
def test_completion_rejects_invalid_native_domains(args):
    with pytest.raises(ValueError):
        h.kraus_numerators(*args)


def test_receipt_scopes_representation_and_its_exact_evidence():
    result = h.run()
    assert not result['scope']['federation_required']
    assert not result['scope']['autonomous_source_hamiltonian_derived']
    assert all(row['all_history_maps_equal'] for row in result['native_instrument_history_actions'])
