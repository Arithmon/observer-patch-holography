"""Producer for shortest two-record path histories, using Dyck words."""
from fractions import Fraction as F
import hashlib
from math import comb

from . import codec


def words(depth):
    def walk(up, down, prefix):
        if up == down == depth-1:
            yield tuple(prefix+[depth-1])
        if up < depth-1:
            yield from walk(up+1, down, prefix+[up+1])
        if down < up:
            yield from walk(up, down+1, prefix+[down])
    return sorted(walk(0, 0, []))


def make():
    controls = []
    for depth in range(1, 11):
        stream = hashlib.sha256()
        for word in words(depth):
            state = [[F(i == j) for j in range(2)] for i in range(depth+1)]
            first = state[-1][:] if depth == 1 else None
            first_slot = 0 if depth == 1 else None
            for k, a in enumerate(word, 1):
                mean = [(x+y)/2 for x, y in zip(state[a], state[a+1])]
                state[a] = state[a+1] = mean
                if a == depth-1 and first is None:
                    first, first_slot = state[-1][:], k
            stream.update(codec.canonical([word, first_slot, list(map(str, first)), list(map(str, state[-1]))]))
        controls.append({"depth": depth, "means": 2*depth-1, "samples": 2*depth,
                         "successful_words": comb(2*depth-2, depth-1)//depth,
                         "word_response_sha256": stream.hexdigest(),
                         "first_response": ["0", str(F(1, 2**(depth-1)))],
                         "final_response": [str(F(1, 2**depth)), str(F(depth+1, 2**(depth+1)))],
                         "first_record_noise_gain": str(F((depth+5)*2**depth, 4)),
                         "second_record_noise_gain": str(2**(depth-1))})
    return controls
