"""Noncommutative replay, native cover and hostile certificate controls."""
from copy import deepcopy

import pytest

from . import presentation as p
from ..clock_defect_energy import defect_history as history


@pytest.fixture(scope='module')
def packet():
    return p.certificate()


def test_original_two_complex_reduces_to_order_two_with_odd_exchange(packet):
    result = p.verify(packet)
    assert result['initial_generators'] == 235
    assert result['native_relators'] == 495
    assert result['verified_Tietze_eliminations'] == 234
    assert result['identity_substitutions'] == 167
    assert result['single_letter_substitutions'] == 67
    assert packet['remaining_generators'] == [1]
    assert packet['remaining_relators'] == [[-1, -1]]
    assert result['exchange_is_nontrivial_generator']
    assert result['higher_cells_or_continuous_surface_classification_used'] is False
    assert result['physical_statistics_selected'] is False


def test_native_relators_reconstructed_independently_in_edge_order(packet):
    pairs, edges, cells, tree, chords, _ = history.geometry()
    variables = {edges[edge]: i+1 for i, edge in enumerate(chords)}
    actual = []
    for _, vertices in cells:
        loop = [*vertices, vertices[0]]
        word = []
        for i in range(len(loop)-1):
            a, b = loop[i:i+2]
            pair = (min(a, b), max(a, b))
            if pair in variables:
                word.append(variables[pair]*(1 if a < b else -1))
        actual.append(word)
    assert actual == packet['native_presentation']['relators']


def test_replay_does_not_call_the_elimination_search(packet, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('producer search called by verifier')
    monkeypatch.setattr(p, 'reduce_presentation', forbidden)
    assert p.verify(packet)['fundamental_group'] == 'Z/2'


def test_free_group_order_is_not_replaced_by_abelianized_exponent_sums():
    commutator = (1, 2, -1, -2)
    assert sum(commutator) == 0
    assert p.free_reduce(commutator) == commutator
    result = p.reduce_presentation(2, [commutator])
    assert result['eliminations'] == []
    assert result['remaining_generators'] == [1, 2]
    assert len(result['remaining_relators'][0]) == 4
    assert p.free_reduce((1, 2, -2, -1)) == ()


def test_unused_free_generator_is_retained_and_alternative_valid_elimination_replays():
    relators = [[2, -1, -1]]
    # b=a^2 defines b but imposes no relation on a.
    result = p.reduce_presentation(2, relators)
    assert result['remaining_generators'] == [1]
    assert result['remaining_relators'] == []
    move = {'generator': 2, 'relator': [-2, 1, 1], 'replacement': [1, 1]}
    assert p.replay_eliminations(2, relators, [move]) == ((1,), ())
    # Both unused generators survive an empty relator inventory.
    assert p.replay_eliminations(2, [], []) == ((1, 2), ())


@pytest.mark.parametrize('damage', ['missing_move', 'duplicate_move', 'wrong_relator',
                                     'wrong_replacement', 'wrong_final_group', 'free_generator',
                                     'changed_cover', 'changed_native_cell', 'boolean'])
def test_mutated_proof_cannot_be_admitted(packet, damage):
    bad = deepcopy(packet)
    if damage == 'missing_move':
        bad['eliminations'].pop()
    elif damage == 'duplicate_move':
        bad['eliminations'].insert(1, deepcopy(bad['eliminations'][0]))
    elif damage == 'wrong_relator':
        bad['eliminations'][0]['relator'] = [1, 1]
    elif damage == 'wrong_replacement':
        move = next(m for m in bad['eliminations'] if m['replacement'])
        move['replacement'][0] *= -1
    elif damage == 'wrong_final_group':
        bad['remaining_relators'] = [[-1]]
    elif damage == 'free_generator':
        bad['remaining_generators'].append(2)
    elif damage == 'changed_cover':
        bad['native_presentation']['cover_parities'][0] ^= 1
    elif damage == 'changed_native_cell':
        bad['native_presentation']['relators'][0] += [1, 2, -1, -2]
    elif damage == 'boolean':
        bad['eliminations'][0]['generator'] = True
    with pytest.raises(ValueError):
        p.verify(bad)


def test_order_corruption_with_same_abelianized_replacement_is_rejected():
    # The native greedy trace happens to use only single-letter replacements.
    # This explicit free-group control tests noncommutative order directly.
    relation = [-3, 1, 2]
    move = {'generator': 3, 'relator': relation, 'replacement': [2, 1]}
    with pytest.raises(ValueError, match='noncommutative'):
        p.replay_eliminations(3, [relation], [move])


def test_omitting_face_fillings_changes_the_native_presentation(packet):
    _, _, cells, _, _, _ = history.geometry()
    original = packet['native_presentation']['relators']
    squares = [word for word, (kind, _) in zip(original, cells) if kind == 'disjoint_square']
    result = p.reduce_presentation(235, squares)
    remaining, relations = p.replay_eliminations(235, squares, result['eliminations'])
    assert len(remaining) >= 20
    assert not (len(remaining) == 1 and len(relations) == 1)
    bad = deepcopy(packet)
    bad['native_presentation']['relators'] = squares
    with pytest.raises(ValueError, match='original native'):
        p.verify(bad)


@pytest.mark.parametrize('word', [[0], [True], [1.0], [236], [1]*4097])
def test_malformed_or_oversized_words_fail_before_reduction(word):
    with pytest.raises(ValueError):
        p.reduce_presentation(235, [word])
