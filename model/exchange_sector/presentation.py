"""Exact noncommutative presentation proof for the native exchange complex.

A spanning-tree presentation of the ORIGINAL two-complex is reduced by
elementary Tietze eliminations to <g | g^2>.  No abelianized rank, higher
product cells or continuous sphere classification enters this certificate.
"""
from __future__ import annotations

from collections import Counter
import json

from ..clock_defect_energy import defect_history as history


MAX_GENERATORS = 235
MAX_RELATORS = 512
MAX_WORD = 4096
SCHEMA = 'oph.native_exchange_presentation.v1'


def _word(value, generators):
    if (type(value) not in (tuple, list) or len(value) > MAX_WORD
            or any(type(x) is not int or not 1 <= abs(x) <= generators for x in value)):
        raise ValueError('a bounded word of nonzero signed integer generators required')
    return tuple(value)


def free_reduce(word):
    """Cancel adjacent inverse letters, without commuting any letters."""
    stack = []
    for letter in word:
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


def inverse(word):
    return tuple(-x for x in reversed(word))


def _relator(word):
    word = free_reduce(word)
    while len(word) > 1 and word[0] == -word[-1]:
        word = word[1:-1]
    if not word:
        return ()
    # Cyclic conjugates and inverses have the same normal closure.
    return min(w[i:]+w[:i] for w in (word, inverse(word)) for i in range(len(w)))


def _relations(relators, generators):
    if (type(generators) is not int or not 1 <= generators <= MAX_GENERATORS
            or type(relators) not in (tuple, list) or len(relators) > MAX_RELATORS):
        raise ValueError('bounded finite generator and relator inventories required')
    return {_relator(_word(w, generators)) for w in relators}-{()}


def _substitute(word, generator, replacement):
    backwards = inverse(replacement)
    expanded = []
    for letter in word:
        expanded.extend(replacement if letter == generator else
                        backwards if letter == -generator else (letter,))
        if len(expanded) > MAX_WORD:
            raise ValueError('finite presentation word budget exceeded')
    return free_reduce(expanded)


def _replacement(relator, generator):
    places = [i for i, letter in enumerate(relator) if abs(letter) == generator]
    if len(places) != 1:
        raise ValueError('an elimination requires exactly one occurrence of its generator')
    place = places[0]
    rest = relator[place+1:]+relator[:place]
    # Rotate to g^s rest=1; hence g=rest^-1 for s=1 and g=rest for s=-1.
    return inverse(rest) if relator[place] > 0 else rest


def reduce_presentation(generators, relators):
    """Generate a bounded exact elimination trace, without claiming its group."""
    current = _relations(relators, generators)
    remaining = set(range(1, generators+1))
    moves = []
    while True:
        choices = []
        for relation in current:
            counts = Counter(map(abs, relation))
            choices.extend((len(relation), relation, i)
                           for i, letter in enumerate(relation) if counts[abs(letter)] == 1)
        if not choices:
            break
        _, relation, place = min(choices)
        generator = abs(relation[place])
        replacement = _replacement(relation, generator)
        moves.append({'generator': generator, 'relator': list(relation),
                      'replacement': list(replacement)})
        current = {_relator(_substitute(w, generator, replacement))
                   for w in current-{relation}}-{()}
        remaining.remove(generator)
    return {'eliminations': moves, 'remaining_generators': sorted(remaining),
            'remaining_relators': [list(w) for w in sorted(current)]}


def replay_eliminations(generators, relators, moves):
    """Verify each elementary group-presentation transformation independently.

    The verifier does not rerun the producer's greedy elimination selection.
    Any valid elimination order is admitted. Empty, duplicate, inverse and
    conjugate relator removals preserve their normal closure in a free group.
    """
    current = _relations(relators, generators)
    remaining = set(range(1, generators+1))
    if type(moves) is not list or len(moves) > generators:
        raise ValueError('a bounded list of elementary eliminations required')
    for move in moves:
        if type(move) is not dict or set(move) != {'generator', 'relator', 'replacement'}:
            raise ValueError('complete elementary elimination record required')
        generator = move['generator']
        if type(generator) is not int or generator not in remaining:
            raise ValueError('an uneliminated integer generator required')
        relation = _word(move['relator'], generators)
        replacement = _word(move['replacement'], generators)
        if relation not in current:
            raise ValueError('the eliminated relator is absent from the current presentation')
        if sum(abs(x) == generator for x in relation) != 1:
            raise ValueError('the generator must occur exactly once')
        if any(abs(x) == generator or abs(x) not in remaining for x in replacement):
            raise ValueError('replacement uses an eliminated or self-referential generator')
        # Independent equation check: substitution into the chosen relation
        # must freely reduce to the empty word, not merely abelianize to zero.
        if _substitute(relation, generator, replacement):
            raise ValueError('replacement does not solve the noncommutative relator')
        current = {_relator(_substitute(w, generator, replacement))
                   for w in current-{relation}}-{()}
        remaining.remove(generator)
    return tuple(sorted(remaining)), tuple(sorted(current))


def native_presentation():
    """Construct the tree presentation and cover parity from native cells."""
    pairs, edges, cells, tree, chords, _ = history.geometry()
    if (len(pairs), len(edges), len(cells), len(tree), len(chords)) != (66, 300, 495, 65, 235):
        raise ArithmeticError('the declared native two-complex census changed')
    if (len(set(tree)) != 65 or set(tree) & set(chords)
            or set(tree) | set(chords) != set(range(300))):
        raise ArithmeticError('tree and chord inventories do not partition native edges')
    adjacency = {i: [] for i in range(66)}
    for index in tree:
        a, b = edges[index]
        adjacency[a].append(b)
        adjacency[b].append(a)
    parity, queue = {0: 0}, [0]

    def edge_parity(a, b):
        source, target = pairs[a], pairs[b]
        removed, added = set(source)-set(target), set(target)-set(source)
        if len(removed) != 1 or len(added) != 1:
            raise ArithmeticError('a declared movement is not a one-particle step')
        old, new = next(iter(removed)), next(iter(added))
        return int(tuple(new if x == old else x for x in source) != target)

    for a in queue:
        for b in adjacency[a]:
            if b not in parity:
                parity[b] = parity[a] ^ edge_parity(a, b)
                queue.append(b)
    if len(parity) != 66:
        raise ArithmeticError('the declared tree is not connected')
    positions = {edge: i for i, edge in enumerate(edges)}
    variables = {edge: i+1 for i, edge in enumerate(chords)}

    def path_word(path):
        word = []
        for a, b in zip(path, path[1:]):
            index = positions[tuple(sorted((a, b)))]
            if index in variables:
                word.append(variables[index] if a < b else -variables[index])
        return tuple(word)

    relators = tuple(path_word(v+v[:1]) for _, v in cells)
    cover = tuple(parity[a] ^ edge_parity(a, b) ^ parity[b] for a, b in (edges[i] for i in chords))
    if any(sum(cover[abs(x)-1] for x in w) % 2 for w in relators):
        raise ArithmeticError('a native deformation cell does not lift closed')
    _, path = history.exchange_loop()
    if path[0] != path[-1]:
        raise ArithmeticError('the exchange witness is not closed')
    exchange = path_word(path)
    if sum(cover[abs(x)-1] for x in exchange) % 2 != 1:
        raise ArithmeticError('the native exchange witness lost its odd ordered lift')
    return {'generators': len(chords), 'relators': [list(w) for w in relators],
            'tree_edges': list(tree), 'chord_edges': list(chords),
            'cover_parities': list(cover), 'exchange_word': list(exchange)}


def certificate():
    initial = native_presentation()
    reduced = reduce_presentation(initial['generators'], initial['relators'])
    result = {'schema': SCHEMA, 'native_presentation': initial, **reduced}
    verify(result)
    return result


def verify(packet):
    """Replay a complete native certificate against the actual source cells."""
    fields = {'schema', 'native_presentation', 'eliminations',
              'remaining_generators', 'remaining_relators'}
    if type(packet) is not dict or set(packet) != fields or packet['schema'] != SCHEMA:
        raise ValueError('complete native presentation certificate required')
    initial = native_presentation()
    supplied = packet['native_presentation']
    # JSON equality distinguishes true from1, unlike Python numeric equality.
    if json.dumps(supplied, sort_keys=True) != json.dumps(initial, sort_keys=True):
        raise ValueError('certificate does not describe the original native two-complex')
    remaining, relators = replay_eliminations(initial['generators'], initial['relators'],
                                            packet['eliminations'])
    if json.dumps(packet['remaining_generators']) != json.dumps(list(remaining)):
        raise ValueError('remaining generator inventory disagrees with replay')
    if json.dumps(packet['remaining_relators']) != json.dumps([list(w) for w in relators]):
        raise ValueError('remaining relators disagree with replay')
    if len(remaining) != 1 or relators != ((-remaining[0], -remaining[0]),):
        raise ValueError('replay did not produce exactly the presentation <g | g^2>')
    generator = remaining[0]
    word = tuple(initial['exchange_word'])
    for move in packet['eliminations']:
        word = _substitute(word, move['generator'], tuple(move['replacement']))
    if (any(abs(x) != generator for x in word) or len(word) % 2 != 1
            or initial['cover_parities'][generator-1] != 1):
        raise ValueError('remaining generator is not the odd native exchange class')
    return {'initial_generators': initial['generators'], 'native_relators': len(initial['relators']),
            'verified_Tietze_eliminations': len(packet['eliminations']),
            'identity_substitutions': sum(not move['replacement'] for move in packet['eliminations']),
            'single_letter_substitutions': sum(len(move['replacement']) == 1 for move in packet['eliminations']),
            'fundamental_group': 'Z/2', 'remaining_generator': generator,
            'exchange_is_nontrivial_generator': True,
            'noncommutative_free_word_verification': True,
            'higher_cells_or_continuous_surface_classification_used': False,
            'physical_statistics_selected': False}


if __name__ == '__main__':
    packet = certificate()
    print(json.dumps({'certificate': packet, 'verification': verify(packet)}, indent=2))
