"""Replay immutable logical versions in reusable owner-local scalar slots.

The original event tape remains the source of every value and consumed writer.
No event is omitted. A separate last-use pass rejects an overwrite before its
final admitted read. This is an offline allocation certificate, not a bounded
memory claim for this Python verifier or for its retained input/output logs.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import struct

ABSENT = 2**64 - 1
ZERO, SCRATCH, PUBLISH, EXPORT, RESET, MEAN, CAPTURE, START, ADD, COMMIT, ADDRESS = range(11)


def require(ok, message):
    if not ok:
        raise ValueError(message)


@dataclass(frozen=True)
class Cell:
    virtual: int
    value: int
    writer: int
    owner: int


def last_uses(rows):
    result = {}
    for eid, row in enumerate(rows):
        for key in row[1:3]:
            if key != ABSENT:
                result[key] = eid
    return result


def replay(rows, carriers, sites, *, slot_override=None):
    """Return the complete mapped tape and its checked allocation census.

    `slot_override` is solely an adversarial test seam. The default layout is
    the 23-slot coloring in Geometry.SourceRoutingStorage. Dead bindings are
    discarded; future reads of a discarded version must fail.
    """
    require(0 < sites <= carriers, "carrier/site bounds")
    rows = list(rows)
    last = last_uses(rows)
    cells, locations = {}, {}
    generations, axes = {}, {}
    mapped, virtual_next, recycled, peak = bytearray(), 0, 0, 0
    occupied, captures, commits = set(), 0, 0
    reads = writes = 0

    def read(key, writer, eid, op, actor):
        if key == ABSENT:
            require(writer == ABSENT, "writer on absent read")
            return ABSENT, None
        require(key in locations, "read of retired or absent version")
        physical = locations[key]
        cell = cells[physical]
        require(cell.virtual == key and cell.writer == writer < eid, "stale writer")
        require(op == MEAN or cell.owner == actor, "remote local read")
        return physical, cell.value

    for eid, row in enumerate(rows):
        require(len(row) == 8, "event width")
        op, a, b, out, actor, wa, wb, encoded = row
        require(0 <= op <= ADDRESS, "opcode")
        value = encoded if encoded < 2**63 else encoded - 2**64
        pa, va = read(a, wa, eid, op, actor)
        pb, vb = read(b, wb, eid, op, actor)
        reads += (a != ABSENT) + (b != ABSENT)
        if op == MEAN:
            require(a != b and va is not None and vb is not None, "mean operands")
            require((va + vb) % 2 == 0 and value == (va + vb)//2, "mean value")
            require(out == a and actor == ABSENT, "mean interface")
            targets = [(a, pa, cells[pa].owner), (b, pb, cells[pb].owner)]
            mapped_out = pa
        else:
            require(0 <= actor < carriers, "owner")
            if op in (ZERO, SCRATCH, ADDRESS, PUBLISH):
                require(a == b == ABSENT, "preparation reads")
            elif op == EXPORT:
                require(va is not None and b == ABSENT and value == va, "export value")
            elif op == RESET:
                require(va == 0 and b == ABSENT and value == 0, "reset value")
            elif op == CAPTURE:
                require(va is not None and b == ABSENT and value == 2*va, "capture value")
            elif op == START:
                require(va == 0 and b == ABSENT and value == 2, "layer start")
            elif op == ADD:
                require(va is not None and vb is not None and value == va + vb, "add value")
            elif op == COMMIT:
                # A local commit intervention is inherited unchanged. Its
                # permitted value is independently checked by the source oracle.
                require(va is not None and b == ABSENT, "commit interface")

            if out == virtual_next:
                require(op in (ZERO, SCRATCH, ADDRESS, PUBLISH, CAPTURE, COMMIT), "allocation opcode")
                if op == ZERO:
                    require(value == 0, "prepared zero")
                    local = 12
                elif op == SCRATCH:
                    local = out % 12 if out < 12*carriers else 19
                    require(out >= 12*carriers or out//12 == actor, "port owner")
                elif op == ADDRESS:
                    axis = axes.get(actor, 0)
                    require(axis < 6, "address width")
                    axes[actor] = axis + 1
                    local = 13 + axis
                elif op == PUBLISH:
                    require(actor not in generations, "duplicate initial payload")
                    generations[actor] = 0
                    occupied.add(actor)
                    local = 20
                elif op == CAPTURE:
                    local = 22
                    captures += 1
                else:
                    require(actor in generations, "commit without initial payload")
                    generations[actor] += 1
                    local = 20 + generations[actor] % 2
                    commits += 1
                physical = 23*actor + local
                if slot_override is not None:
                    physical = slot_override(eid, op, actor, local, physical)
                require(0 <= physical < 23*carriers and physical//23 == actor, "physical owner/slot")
                if physical in cells:
                    old = cells[physical].virtual
                    require(last.get(old, -1) < eid, "overwrite of a live version")
                    require(op in (CAPTURE, COMMIT), "static allocation alias")
                    del locations[old]
                    recycled += 1
                locations[out] = physical
                virtual_next += 1
            else:
                require(out in locations, "write to retired/absent register")
                physical = locations[out]
                require(op in (EXPORT, RESET, START, ADD), "protected overwrite")
                require(physical % 23 < 12 or physical % 23 == 19, "mutable destination")
            targets = [(out, physical, actor)]
            mapped_out = physical
        for virtual, physical, who in targets:
            require(physical//23 == who, "write ownership")
            require(physical not in cells or cells[physical].owner == who, "owner drift")
            cells[physical] = Cell(virtual, value, eid, who)
            writes += 1
        mapped.extend(struct.pack("<8Q", op, pa, pb, mapped_out, actor, wa, wb, encoded))
        peak = max(peak, len(cells))
        require(len(locations) == len(cells), "binding inventory")

    require(len(occupied) == sites, "occupied host inventory")
    require(all(axes.get(owner) == 6 for owner in occupied), "missing addresses")
    require(peak <= 14*carriers + 9*sites, "working store bound")
    return bytes(mapped), {
        "events": len(rows), "reads": reads, "writes": writes,
        "logical_allocations": virtual_next, "physical_slots_used": peak,
        "slot_bound": 14*carriers + 9*sites,
        "per_carrier_slot_bound": 23,
        "recycled_allocations": recycled, "captures": captures, "commits": commits,
        "mapped_decoded_sha256": hashlib.sha256(mapped).hexdigest(),
    }


def check_refinement(source_rows, data, carriers):
    """Compare semantic fields independently of the allocation algorithm.

    A physical-store execution alone may implement a different program. Exact
    opcodes, values, owners and consumed writers bind it to the checked source
    word, hence preserve its consumed-writer graph, including absent operands.
    """
    require(len(data) == 64*len(source_rows), "target/source event count")
    for source, target in zip(source_rows, struct.iter_unpack("<8Q", data), strict=True):
        require(all(source[i] == target[i] for i in (0,4,5,6,7)),
                "target/source semantic field drift")
        require(all((source[i] == ABSENT) == (target[i] == ABSENT) for i in (1,2,3)),
                "target/source operand presence drift")
        for field in (1,2,3):
            address = source[field]
            if address < 12*carriers:
                expected = 23*(address//12) + address%12
                require(target[field] == expected, "target/source port identity drift")
    return len(source_rows)


def replay_physical(data, carriers):
    """Independent target-store replay: no logical-register map or last-use data.

    Every consumed writer and scalar arithmetic operation is reconstructed
    directly on the mapped local store. Source preparation and commit-bias
    validation belong to the separately replayed original program.
    """
    require(carriers > 0, "target carrier count")
    require(len(data) % 64 == 0, "partial mapped event")
    cells = {}
    for eid, row in enumerate(struct.iter_unpack("<8Q", data)):
        op, a, b, out, owner, wa, wb, encoded = row
        require(0 <= op <= ADDRESS, "target opcode")
        value = encoded if encoded < 2**63 else encoded - 2**64
        operands = []
        for address, writer in ((a, wa), (b, wb)):
            if address == ABSENT:
                require(writer == ABSENT, "target absent writer")
                operands.append(None)
            else:
                require(0 <= address < 23*carriers, "target read address bound")
                require(address in cells and cells[address][1] == writer < eid, "target stale writer")
                require(op == MEAN or address//23 == owner, "target remote read")
                operands.append(cells[address][0])
        x, y = operands
        if op == MEAN:
            require(a != b and a//23 != b//23 and a%23 < 12 and b%23 < 12,
                    "target mean ports")
            require(out == a and owner == ABSENT, "target mean interface")
        else:
            require(0 <= out < 23*carriers and 0 <= owner < carriers and out//23 == owner,
                    "target local destination")
            allowed = {ZERO: (12,), SCRATCH: (*range(12),19), ADDRESS: tuple(range(13,19)),
                       PUBLISH: (20,), CAPTURE: (22,), COMMIT: (20,21),
                       EXPORT: tuple(range(12)), RESET: tuple(range(12)), START: (19,), ADD: (19,)}
            require(out%23 in allowed[op], "target destination kind")
            if op in (ZERO, SCRATCH, ADDRESS, PUBLISH):
                require(a == b == ABSENT and out not in cells, "target preparation interface")
            elif op == ADD:
                require(a == out and b != ABSENT and b%23 in (20,21,22), "target add interface")
            else:
                require(a != ABSENT and b == ABSENT, "target unary interface")
            if op in (EXPORT, RESET, START, ADD):
                require(out in cells, "target unprepared mutable destination")
            if op in (RESET, START):
                require(a%23 == 12, "target zero operand")
            elif op == EXPORT:
                require(a%23 in (20,21,22), "target export operand")
            elif op == CAPTURE:
                require(a%23 < 12, "target capture operand")
            elif op == COMMIT:
                require(a%23 == 19, "target commit operand")
        if op == ZERO:
            require(value == 0, "target prepared zero")
        if op == EXPORT:
            require(x is not None and value == x, "target export")
        elif op == RESET:
            require(x == value == 0, "target reset")
        elif op == START:
            require(x == 0 and value == 2, "target start")
        elif op == CAPTURE:
            require(x is not None and value == x*2, "target capture")
        elif op == ADD:
            require(x is not None and y is not None and value == x+y, "target addition")
        elif op == MEAN:
            require(x is not None and y is not None and x+y == 2*value, "target mean")
        targets = (a, b) if op == MEAN else (out,)
        for address in targets:
            require(op == MEAN or address//23 == owner, "target remote write")
            cells[address] = (value, eid)
    return {"events": len(data)//64, "physical_slots_used": len(cells)}
