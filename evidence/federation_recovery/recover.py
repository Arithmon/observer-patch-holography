"""Recover a bounded observed federation from baseline event rows only.

``recover(rows)`` performs no I/O and accepts no support, input header, source
specification or receipt.  A row contains eight integers: opcode, two input
registers, output register, owner, two consumed writer IDs, and twice the
scalar value.  The first seven are unsigned 64-bit integers.  The last may
be signed int64 or its unsigned two's-complement representation, as emitted
by the lossless tape decoder.  Event IDs are consecutive row positions.

The operation alphabet and the baseline preparation/control grammar below
are supplied, not discovered laws.  All reported population, routing and
order data are reconstructed from actual reads and writes within that
grammar.  In particular, an unexecuted seam cannot be recovered.  Replaying
the observed program does not identify the complete topology or a general
population/routing algorithm.  The bounded interface deliberately accepts
the baseline cohort only, not intervention tapes.
"""

from collections import Counter


ABSENT = (1 << 64) - 1
MAX_EVENTS = 100_000
MAX_REGISTERS = 30_000
MAX_LOGICAL_NODES = 512
OP_NAMES = (
    "ZERO", "SCRATCH", "PUBLISH", "EXPORT", "RESET", "MEAN", "CAPTURE",
    "START", "ADD", "COMMIT", "ADDRESS",
)
ZERO, SCRATCH, PUBLISH, EXPORT, RESET, MEAN, CAPTURE, START, ADD, COMMIT, ADDRESS = range(11)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _bits(mask):
    result = []
    while mask:
        low = mask & -mask
        result.append(low.bit_length() - 1)
        mask ^= low
    return result


def _rows(source):
    result = []
    for row in source:
        require(len(result) < MAX_EVENTS, "bounded recovery event limit")
        require(isinstance(row, (tuple, list)) and len(row) == 8, "event width/type")
        require(all(type(x) is int for x in row), "event fields must be exact integers")
        require(all(0 <= x <= ABSENT for x in row[:7]), "unsigned event metadata")
        require(-(1 << 63) <= row[7] <= ABSENT, "scalar outside signed/encoded int64")
        value = row[7] - (1 << 64) if row[7] >= (1 << 63) else row[7]
        require(row[0] < len(OP_NAMES), "unknown primitive")
        result.append(tuple(row[:7]) + (value,))
    require(bool(result), "empty event history")
    return result


def recover(source):
    """Validate baseline rows and return JSON-compatible observed structure.

    Input custody, level attribution and comparison with a hidden generating
    specification belong outside this function.  Rejection is ``ValueError``;
    a successful return never certifies absence of unobserved seams.
    """
    rows = _rows(source)
    # Cell fields: value, latest writer, owner, immutable, logical ancestry,
    # unique transported logical version (or None), and actual owner path.
    cells = []
    operations = Counter()
    reads = writes = maximum = 0
    cursor = 0

    def execute(kind, expected, *, fresh=False, immutable=False, origin=None, path=()):
        nonlocal cursor, reads, writes, maximum
        require(cursor < len(rows), "truncated operation word")
        op, a, b, out, owner, wa, wb, value = rows[cursor]
        require(op == kind, "unexpected primitive in baseline grammar")
        require(value == expected, "incorrect primitive arithmetic")
        mask = 0
        for key, writer in ((a, wa), (b, wb)):
            if key == ABSENT:
                require(writer == ABSENT, "writer supplied without a read")
                continue
            require(key < len(cells), "read of absent register")
            cell = cells[key]
            require(writer == cell[1] and writer < cursor, "stale/incorrect consumed writer")
            require(op == MEAN or cell[2] == owner, "nonlocal register read")
            mask |= cell[4]
            reads += 1
        require((a != ABSENT or wa == ABSENT) and (b != ABSENT or wb == ABSENT), "absent read")
        if op in (PUBLISH, COMMIT):
            require(origin is not None and origin < MAX_LOGICAL_NODES, "logical-node limit")
            mask |= 1 << origin
        targets = (a, b) if op == MEAN else (out,)
        require(op != MEAN or (a != b and out == a and owner == ABSENT), "pair-mean write interface")
        if op != MEAN:
            require(owner < len(carrier_ports), "unknown owner")
        for target in targets:
            require(target <= len(cells), "nonconsecutive register allocation")
            who = cells[target][2] if op == MEAN else owner
            if target == len(cells):
                require(fresh and len(cells) < MAX_REGISTERS, "unexpected/too many fresh registers")
                cells.append([value, cursor, who, immutable, mask, origin, tuple(path)])
            else:
                require(not fresh and not cells[target][3], "overwrite of protected/fresh register")
                require(cells[target][2] == who, "nonlocal register write")
                cells[target] = [value, cursor, who, False, mask, origin, tuple(path)]
            writes += 1
        operations[OP_NAMES[op]] += 1
        maximum = max(maximum, abs(value))
        cursor += 1

    def peek(kind=None):
        require(cursor < len(rows), "truncated baseline history")
        row = rows[cursor]
        if kind is not None:
            require(row[0] == kind, "unexpected primitive in baseline grammar")
        return row

    def no_reads(row):
        require(row[1] == row[2] == row[5] == row[6] == ABSENT, "preparation has inputs")

    def one_read(row):
        require(row[1] != ABSENT and row[2] == row[6] == ABSENT, "one-input interface")
        require(row[1] < len(cells), "read of absent register")

    # Infer the owner/port census from the opening SCRATCH block.  The block
    # order and twelve ports per owner are part of the accepted grammar.
    opening = 0
    while opening < len(rows) and rows[opening][0] == SCRATCH:
        opening += 1
    require(opening > 0 and opening % 12 == 0, "twelve-port preparation blocks")
    carrier_ports = [list(range(i, i + 12)) for i in range(0, opening, 12)]
    for index in range(opening):
        row = peek(SCRATCH)
        no_reads(row)
        require(row[3] == index and row[4] == index // 12, "prepared port owner/order")
        execute(SCRATCH, row[7], fresh=True)

    zero_registers = []
    for owner in range(len(carrier_ports)):
        row = peek(ZERO)
        no_reads(row)
        require(row[3] == len(cells) and row[4] == owner, "one zero per owner")
        zero_registers.append(row[3])
        execute(ZERO, 0, fresh=True, immutable=True)

    population = []
    while cursor < len(rows) and rows[cursor][0] == SCRATCH:
        row = peek(SCRATCH)
        no_reads(row)
        require(row[3] == len(cells) and row[4] < len(carrier_ports), "population accumulator")
        population.append({"site": len(population), "owner": row[4], "accumulator": row[3],
                           "address_registers": [], "address_coefficients": []})
        execute(SCRATCH, 0, fresh=True)
    n = len(population)
    require(n > 0 and len({p["owner"] for p in population}) == n, "nonempty injective population")
    require(n <= MAX_LOGICAL_NODES, "logical-node limit")
    for site in population:
        for _ in range(6):
            row = peek(ADDRESS)
            no_reads(row)
            require(row[3] == len(cells) and row[4] == site["owner"] and row[7] % 2 == 0,
                    "six integral prepared address coefficients per site")
            site["address_registers"].append(row[3])
            site["address_coefficients"].append(row[7] // 2)
            execute(ADDRESS, row[7], fresh=True, immutable=True)

    # The physical role of these six coefficients is supplied.  This checks
    # the observed baseline split without consulting a golden-coordinate
    # generator or asserting a population law at any unobserved regulator.
    positive, negative = (0, 1, 4, 5, 8, 9), (3, 2, 7, 6, 11, 10)
    occupied = {p["owner"]: p for p in population}
    for owner, ports in enumerate(carrier_ports):
        expected = [0] * 12
        if owner in occupied:
            for axis, value in enumerate(occupied[owner]["address_coefficients"]):
                expected[positive[axis]] = 2 * max(0, value)
                expected[negative[axis]] = 2 * max(0, -value)
        require([cells[p][0] for p in ports] == expected, "baseline antipodal preparation")

    versions = []

    def record_version(site, layer, eid, register):
        index = len(versions)
        cell = cells[register]
        require(cell[5] == index, "logical-version identity")
        versions.append({"index": index, "site": site, "layer": layer, "event": eid,
                         "register": register, "scaled_value": cell[0], "ancestors": _bits(cell[4])})

    for s, site in enumerate(population):
        row = peek(PUBLISH)
        no_reads(row)
        require(row[3] == len(cells) and row[4] == site["owner"], "baseline initial publication")
        site["initial_payload_register"] = row[3]
        eid = cursor
        execute(PUBLISH, 2 * (s + 1), fresh=True, immutable=True, origin=s, path=(site["owner"],))
        record_version(s, 0, eid, row[3])

    observed = set()
    port_partners = {}
    carrier_pair_seams = {}
    hops = []
    logical_reads = []
    edges = set()
    rounds = 0
    reference_menu = None
    while cursor < len(rows):
        layer = rounds + 1
        require((layer + 1) * n <= MAX_LOGICAL_NODES, "logical-node limit")
        for site in population:
            row = peek(START)
            one_read(row)
            require(row[1] == zero_registers[site["owner"]] and row[3] == site["accumulator"]
                    and row[4] == site["owner"], "layer start interface/order")
            execute(START, 2)
        consumed = [set() for _ in population]
        while peek()[0] != COMMIT:
            row = peek()
            if row[0] == EXPORT:
                one_read(row)
                archive, source_port = row[1], row[3]
                require(source_port < opening and cells[archive][3] and cells[archive][5] is not None,
                        "export of protected logical payload to a port")
                origin = cells[archive][5]
                require((layer - 1) * n <= origin < layer * n, "wrong-layer transported payload")
                require(row[4] == cells[source_port][2] == cells[archive][2], "export owner")
                path = cells[archive][6]
                require(path and path[-1] == row[4], "transport path owner")
                start = cursor
                execute(EXPORT, cells[archive][0], origin=origin, path=path)
                row = peek(RESET)
                one_read(row)
                target_port, target_owner = row[3], row[4]
                require(target_port < opening and target_port != source_port and target_owner < len(carrier_ports),
                        "receiver port")
                require(target_owner == cells[target_port][2] and target_owner != cells[source_port][2]
                        and row[1] == zero_registers[target_owner], "receiver reset locality")
                require(target_owner not in path, "cyclic transported path")
                execute(RESET, 0)
                row = peek(MEAN)
                require(row[1] == source_port and row[2] == target_port and row[3] == source_port
                        and row[4] == ABSENT, "mean endpoints")
                total = cells[source_port][0] + cells[target_port][0]
                require(total % 2 == 0 and cells[target_port][0] == 0, "exact pair mean with reset receiver")
                seam = tuple(sorted((source_port, target_port)))
                carrier_pair = tuple(sorted((cells[source_port][2], cells[target_port][2])))
                require(port_partners.get(source_port, target_port) == target_port
                        and port_partners.get(target_port, source_port) == source_port,
                        "observed port has inconsistent gluing partner")
                require(carrier_pair_seams.get(carrier_pair, seam) == seam,
                        "multiple observed seams between one carrier pair")
                port_partners[source_port] = target_port
                port_partners[target_port] = source_port
                carrier_pair_seams[carrier_pair] = seam
                execute(MEAN, total // 2, origin=origin, path=path)
                observed.add(seam)
                row = peek(CAPTURE)
                one_read(row)
                capture = row[3]
                require(row[1] == target_port and capture == len(cells) and row[4] == target_owner,
                        "factor-two protected capture")
                execute(CAPTURE, 2 * cells[target_port][0], fresh=True, immutable=True,
                        origin=origin, path=path + (target_owner,))
                for endpoint in (source_port, target_port):
                    row = peek(RESET)
                    one_read(row)
                    owner = cells[endpoint][2]
                    require(row[1] == zero_registers[owner] and row[3] == endpoint and row[4] == owner,
                            "hop cleanup interface/order")
                    execute(RESET, 0)
                hops.append({"events": list(range(start, cursor)), "source_version": origin,
                             "source_register": archive, "source_port": source_port,
                             "target_port": target_port, "capture_register": capture})
            elif row[0] == ADD:
                a, b, out, owner = row[1:5]
                require(a < len(cells) and b < len(cells) and a == out, "accumulation interface")
                target = next((p["site"] for p in population if p["accumulator"] == out), None)
                require(target is not None and owner == population[target]["owner"], "population target")
                origin = cells[b][5]
                require(cells[b][3] and origin is not None and (layer - 1) * n <= origin < layer * n,
                        "read of latest protected logical payload")
                require(origin not in consumed[target], "duplicate logical read")
                consumed[target].add(origin)
                logical_reads.append({"event": cursor, "source": origin, "target": layer * n + target,
                                      "via_register": b, "path": list(cells[b][6])})
                edges.add((origin, layer * n + target))
                execute(ADD, cells[a][0] + cells[b][0])
            else:
                raise ValueError("unexpected primitive outside a complete transport word")
        menu = [[source - (layer - 1) * n for source in sorted(items)] for items in consumed]
        require(all(menu), "each population site must have an observed read")
        if reference_menu is None:
            reference_menu = menu
        else:
            require(menu == reference_menu, "observed read menu changes between baseline layers")
        for s, site in enumerate(population):
            row = peek(COMMIT)
            one_read(row)
            require(row[1] == site["accumulator"] and row[3] == len(cells) and row[4] == site["owner"],
                    "complete baseline layer commit/order")
            eid = cursor
            execute(COMMIT, cells[row[1]][0], fresh=True, immutable=True,
                    origin=layer * n + s, path=(site["owner"],))
            record_version(s, layer, eid, row[3])
        rounds += 1
    require(rounds > 0, "no completed logical layer")

    # Recover the strict logical order from consumed-writer ancestry.  Check
    # independently inside this decoder that the reconstructed read edges
    # generate precisely that order, excluding serial/control ancestry.
    ancestors = []
    heights = []
    incoming = [[] for _ in versions]
    for a, b in edges:
        incoming[b].append(a)
    for version in versions:
        index = version["index"]
        mask = 1 << index
        for parent in incoming[index]:
            require(parent < index, "logical edge is not forward")
            mask |= ancestors[parent]
        require(_bits(mask) == version["ancestors"], "read-edge and actual-writer order disagree")
        ancestors.append(mask)
        heights.append(1 + max((heights[parent] for parent in incoming[index]), default=0))
    descendants = [0] * len(versions)
    for target, mask in enumerate(ancestors):
        for source_index in _bits(mask & ~(1 << target)):
            descendants[source_index] |= 1 << target
    histogram = Counter()
    strict = 0
    for target, mask in enumerate(ancestors):
        for source_index in _bits(mask & ~(1 << target)):
            strict += 1
            interior = (descendants[source_index] & mask & ~(1 << target)).bit_count()
            histogram[str(interior)] += 1

    return {
        "schema": "oph.federation_recovery.observed.v1",
        "grammar": {"cohort": "baseline", "value_encoding": "signed_int64_half_units",
                    "ports_per_owner": 12, "address_coefficients_per_site": 6,
                    "primitive_semantics": "supplied finite register grammar"},
        "scope": {"event_only": True, "observed_routing_recovered": True,
                  "population_records_recovered": True, "logical_poset_recovered": True,
                  "full_port_gluing_recovered": False, "general_population_law_recovered": False,
                  "general_routing_algorithm_recovered": False, "physical_identification": False},
        "program": [list(row[:5]) + [row[7] if row[0] in (ZERO, SCRATCH, PUBLISH, ADDRESS) else None]
                    for row in rows],
        "census": {"events": len(rows), "registers": len(cells), "scalar_reads": reads,
                   "scalar_writes": writes, "operation_counts": {name: operations[name] for name in OP_NAMES},
                   "max_abs_scaled_scalar": maximum},
        "carrier_ports": carrier_ports,
        "population": population,
        "observed_gluing": [list(edge) for edge in sorted(observed)],
        "transport_hops": hops,
        "logical_versions": versions,
        "logical_reads": logical_reads,
        "logical_edges": [list(edge) for edge in sorted(edges)],
        "poset": {"nodes": len(versions), "strict_relations": strict,
                  "reflexive_relations": strict + len(versions), "layer_sizes": [n] * (rounds + 1),
                  "max_chain_length": max(heights),
                  "open_interval_cardinality_histogram": dict(sorted(histogram.items(), key=lambda item: int(item[0])))},
        "invariants": {"owners": len(carrier_ports), "ports": opening, "population_sites": n,
                       "rounds": rounds, "observed_glued_seams": len(observed),
                       "logical_nodes": len(versions), "logical_reads": len(logical_reads),
                       "logical_strict_relations": strict, "logical_height": max(heights)},
    }


def instantiate(recovered):
    """Execute the recovered finite program, computing writers and results.

    Each instruction is ``[op, input_a, input_b, output, owner, literal]``.
    A literal is present only for initial preparation; every other result
    is computed from the current cells.  The output uses uint64 value bits,
    permitting exact byte comparison with decoded source tapes.  This is
    instantiation of an observed finite instruction sequence under supplied
    semantics, not inference of a source generator or full gluing graph.
    """
    require(isinstance(recovered, dict) and recovered.get("schema") == "oph.federation_recovery.observed.v1",
            "recovered program schema")
    program = recovered.get("program")
    require(isinstance(program, list) and 0 < len(program) <= MAX_EVENTS, "bounded finite program")
    cells = []  # signed value, writer, owner, immutable
    result = []
    for eid, instruction in enumerate(program):
        require(isinstance(instruction, list) and len(instruction) == 6, "instruction width/type")
        op, a, b, out, owner, literal = instruction
        require(all(type(x) is int and 0 <= x <= ABSENT for x in instruction[:5]), "instruction metadata")
        require(op < len(OP_NAMES), "unknown instruction")
        initial = op in (ZERO, SCRATCH, PUBLISH, ADDRESS)
        if initial:
            require(type(literal) is int and -(1 << 63) <= literal < (1 << 63), "initial literal")
            require(a == b == ABSENT and out == len(cells) and owner != ABSENT, "initial instruction")
            value = literal
            require(op != ZERO or value == 0, "zero initialization")
        else:
            require(literal is None, "computed instruction contains an output literal")
            require(a < len(cells), "instruction input absent")
            if op in (MEAN, ADD):
                require(b < len(cells), "second instruction input absent")
            else:
                require(b == ABSENT, "extra instruction input")
            if op in (EXPORT, COMMIT):
                value = cells[a][0]
            elif op == RESET:
                require(cells[a][0] == 0 and cells[a][3], "reset source is not protected zero")
                value = 0
            elif op == MEAN:
                total = cells[a][0] + cells[b][0]
                require(total % 2 == 0, "half-unit mean domain")
                value = total // 2
            elif op == CAPTURE:
                value = 2 * cells[a][0]
            elif op == START:
                require(cells[a][0] == 0 and cells[a][3], "start source is not protected zero")
                value = 2
            elif op == ADD:
                value = cells[a][0] + cells[b][0]
            else:
                raise ValueError("unknown arithmetic instruction")
        require(-(1 << 63) <= value < (1 << 63), "instantiated scalar overflow")
        for key in (a, b):
            if key != ABSENT:
                require(key < len(cells), "instruction input absent")
                require(op == MEAN or cells[key][2] == owner, "nonlocal instruction read")
        wa = ABSENT if a == ABSENT else cells[a][1]
        wb = ABSENT if b == ABSENT else cells[b][1]
        result.append((op, a, b, out, owner, wa, wb, value & ABSENT))
        targets = (a, b) if op == MEAN else (out,)
        if op == MEAN:
            require(a != b and out == a and owner == ABSENT, "pair-mean instruction interface")
        fresh = initial or op in (CAPTURE, COMMIT)
        immutable = op in (ZERO, PUBLISH, ADDRESS, CAPTURE, COMMIT)
        for target in targets:
            require(target <= len(cells), "nonconsecutive instruction output")
            who = cells[target][2] if op == MEAN else owner
            if target == len(cells):
                require(fresh and len(cells) < MAX_REGISTERS, "unexpected/too many instruction registers")
                cells.append([value, eid, who, immutable])
            else:
                require(not fresh and not cells[target][3] and cells[target][2] == who,
                        "protected/nonlocal instruction overwrite")
                cells[target] = [value, eid, who, False]
    return result
