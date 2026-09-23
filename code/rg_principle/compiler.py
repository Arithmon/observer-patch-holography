"""An intentionally simple, resource-explicit finite Boolean compiler."""

def bits(value, width):
    return tuple((value >> j) & 1 for j in range(width))


def lower_cnot(program):
    """Eliminate CNOT with one clean constant-control bit, including setup/cleanup."""
    if not any(len(gate) == 2 for gate in program["gates"]):
        return program
    control = program["n"] + program["m"] + program["work"]
    gates = [[control]]
    gates.extend([control, *gate] if len(gate) == 2 else list(gate) for gate in program["gates"])
    gates.append([control])
    return {**program, "work": program["work"] + 1, "gates": gates}


def compile_table(n, m, table):
    """Retain inputs, XOR f(input) into outputs, and clean all work bits."""
    if type(n) is not int or type(m) is not int or not 0 <= n <= 8 or not 1 <= m <= 8:
        raise ValueError("bounded reference compiler requires 0<=n<=8, 1<=m<=8")
    if (type(table) is not list or len(table) != 2**n or
            any(type(row) is not list or len(row) != m or
                any(type(b) is not int or b not in (0, 1) for b in row) for row in table)):
        raise ValueError("complete binary truth table required")
    work = max(n - 2, 0)
    gates = []
    for value, output in enumerate(table):
        negate = [j for j, b in enumerate(bits(value, n)) if b == 0]
        gates.extend([j] for j in negate)
        for j, bit in enumerate(output):
            if not bit:
                continue
            target = n + j
            if n < 3:
                gates.append(list(range(n)) + [target])
            else:
                ladder = [[0, 1, n + m]]
                for k in range(2, n - 1):
                    ladder.append([n + m + k - 2, k, n + m + k - 1])
                gates.extend(ladder)
                gates.append([n + m + work - 1, n - 1, target])
                gates.extend(reversed(ladder))
        gates.extend([j] for j in reversed(negate))
    return lower_cnot({"n": n, "m": m, "work": work, "table": table, "gates": gates})


def cases():
    for code in range(256):
        yield f"boolean3_{code:03d}", compile_table(3, 1, [[(code >> i) & 1] for i in range(8)])
    for n in (0, 1, 2, 4, 5):
        table = [[sum(bits(i, n)) % 2, int(all(bits(i, n))), int(i == 0)]
                 for i in range(2**n)]
        yield f"multiple_{n}", compile_table(n, 3, table)
    # Transport diagnostics expose distinct designated intermediate versions.
    yield "retained_copy", lower_cnot({"n": 1, "m": 2, "work": 0,
                                      "table": [[0, 0], [1, 1]], "gates": [[0, 1], [1, 2]]})
    yield "clean_conjunction", compile_table(4, 1, [[int(i == 15)] for i in range(16)])
