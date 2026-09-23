"""Exact zero-error central-channel criterion, including misleading controls.

This criterion concerns a supplied operational channel. An algebraic seam
inclusion or an invertible matrix of probabilities is not such a decoder.
"""
from fractions import Fraction as Q


def classify(rows):
    if type(rows) is not list or not rows or type(rows[0]) is not list or not rows[0]:
        raise ValueError("channel shape")
    width = len(rows[0])
    parsed = []
    for row in rows:
        if type(row) is not list or len(row) != width:
            raise ValueError("channel row shape")
        values = []
        for entry in row:
            if type(entry) is not str or len(entry) >= 128:
                raise ValueError("channel rational type/size")
            value = Q(entry)
            if str(value) != entry or value < 0:
                raise ValueError("channel probability")
            values.append(value)
        if sum(values) != 1:
            raise ValueError("channel normalization")
        parsed.append(values)
    decoder = []
    collisions = []
    for output in range(width):
        possible = [source for source, row in enumerate(parsed) if row[output] > 0]
        decoder.append(possible[0] if len(possible) == 1 else None)
        if len(possible) > 1:
            collisions.append([output, possible[0], possible[1]])
    if not collisions:
        # Independently multiply channel and deterministic decoder. Every
        # conditional output law must return the intervened source with mass 1.
        for source, row in enumerate(parsed):
            for guess in range(len(parsed)):
                mass = sum(row[y] for y, d in enumerate(decoder) if d == guess)
                if mass != int(source == guess):
                    raise ValueError("decoder reconstruction")
    return {"inputs": len(parsed), "outputs": width,
            "decoder": None if collisions else decoder,
            "collisions": collisions}


def controls():
    count = 12
    identity = [[str(int(x == y)) for y in range(count)] for x in range(count)]
    noisy = [[str(Q(11, 12)*int(x == y)+Q(1, 144))
              for y in range(count)] for x in range(count)]
    cases = {"identity": identity,
             "permutation": [[str(int(y == (x+1) % count)) for y in range(count)]
                             for x in range(count)],
             "scalar": [["1"] for _ in range(count)],
             "full_rank_noisy": noisy,
             "redundant_outputs": [[str(Q(1, 2) if y//2 == x else Q(0))
                                    for y in range(2*count)] for x in range(count)]}
    result = {name: classify(rows) for name, rows in cases.items()}
    if any((result[name]["decoder"] is not None) != expected for name, expected in
           {"identity": True, "permutation": True, "scalar": False,
            "full_rank_noisy": False, "redundant_outputs": True}.items()):
        raise ValueError("channel separation controls")
    # P = 11/12 I + 1/144 J has determinant (11/12)^11, yet no
    # single-shot exact decoder: every output occurs for every source input.
    inverse = [[Q(12, 11)*int(x == y)-Q(1, 132) for y in range(count)]
               for x in range(count)]
    for x in range(count):
        for z in range(count):
            if sum(Q(noisy[x][y])*inverse[y][z] for y in range(count)) != int(x == z):
                raise ValueError("noisy full-rank control")
    result["full_rank_noisy"]["signed_inverse_is_stochastic"] = all(
        value >= 0 for row in inverse for value in row)
    return result
