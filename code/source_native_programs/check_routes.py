"""Independent complete route-witness validation; no search implementation import."""
import hashlib
import json

from . import codec


def require(condition, message):
    if not condition:
        raise ValueError(message)


def edges(level):
    require(type(level) is int and level in (3,4,5), "support level")
    support = codec.load(codec.ROOT/codec.SUPPORTS[level])
    schema = "oph.source_routing.w12_support.v1" if level == 3 else "oph.source_read_routing.support.v1"
    require(support["schema"] == schema and
            support["level"] == level and support["carriers"] == 20*4**level,
            "captured support identity")
    result = {tuple(sorted((12*c+a,12*c+b))) for c in range(support["carriers"])
              for a,b in support["intra_carrier_seams"]}
    result.update(tuple(sorted((12*c+a,12*d+b))) for c,a,d,b in support["glued_pairs"])
    return result


def check_route(row, sites, level, carrier, direction, allowed, protected=None):
    require(type(row) is dict and set(row) == {"carrier","direction","pairs"}, "route fields")
    codec.equal([row["carrier"],row["direction"]],[carrier,direction], "route inventory/order")
    pairs = row["pairs"]
    require(type(pairs) is list and len(pairs) >= 3, "paired route too short")
    require(all(type(p) is list and len(p) == 2 and
                all(type(v) is int and 0 <= v < 12*20*4**level for v in p) for p in pairs),
            "invalid physical port")
    flat = [v for p in pairs for v in p]
    require(len(set(flat)) == len(flat), "route alias or intersecting rails")
    require(all(tuple(sorted((a[r],b[r]))) in allowed
                for a,b in zip(pairs,pairs[1:]) for r in (0,1)), "unsupported route edge")
    record = {12*carrier,12*carrier+1}
    if protected is None:
        protected = {12*c+r for c in range(1,2*sites+1) for r in (0,1)}
    if direction == "read":
        require(pairs[0] == [12*carrier+5,12*carrier+9] and set(pairs[-1]) == {3,5},
                "read endpoints")
        forbidden = protected | {1,2,4,6,7,9}
        require(all(tuple(sorted((12*carrier+r,pairs[0][r]))) in allowed for r in (0,1)),
                "unsupported source export")
    else:
        require(pairs[0] == [4,9] and set(pairs[-1]) == record, "write endpoints")
        forbidden = (protected-record) | {1,2,3,5,6,7}
        require(tuple(sorted(record)) in allowed, "destination has no reset rung")
        require(all(tuple(sorted((s,t))) in allowed for s,t in zip((2,1),pairs[0])),
                "unsupported accumulator export")
    require(not set(flat)&forbidden, "route clobbers protected live record/core")
    require(tuple(sorted(pairs[0])) in allowed, "bus root has no reset rung")
    return len(pairs)


def check(rows, sites, level, retain=False):
    require(type(sites) is int and sites >= 1 and 2*sites+1 <= 20*4**level,
            "two-bank capacity")
    allowed = edges(level)
    protected = {12*c+r for c in range(1,2*sites+1) for r in (0,1)}
    # Eight distinct core ports, both operand polarities, and retirement rectangle.
    core = [(4,9),(2,4),(1,9),(6,2),(7,1),
            (3,4),(5,9),(3,9),(5,4)]
    require(all(tuple(sorted(e)) in allowed for e in core), "unsupported arithmetic core")
    rows = iter(rows)
    end = object()  # JSON null is a record, never an end-of-stream marker.
    digest = hashlib.sha256()
    total = maximum = count = reversed_ends = 0
    kept = {}
    for carrier in range(1,2*sites+1):
        for direction in ("read","write"):
            row = next(rows,end)
            length = check_route(row,sites,level,carrier,direction,allowed,protected)
            digest.update((codec.compact(row)+"\n").encode("ascii"))
            count += 1
            total += length
            maximum = max(maximum,length)
            reversed_ends += row["pairs"][-1][0] > row["pairs"][-1][1]
            if retain:
                kept[carrier,direction] = row["pairs"]
    require(next(rows,end) is end, "extra route witness")
    summary = {"sites":sites,"level":level,"routes":count,"paired_cells":total,
               "max_paired_cells":maximum,"reversed_endpoints":reversed_ends,
               "routes_sha256":digest.hexdigest(),
               "support_sha256":hashlib.sha256((codec.ROOT/codec.SUPPORTS[level]).read_bytes()).hexdigest(),
               "scope":"complete_two_bank_topology_witness_not_native_execution"}
    return (summary,kept) if retain else summary


def read_rows(path):
    # Reuse only the strict JSON decoder's parser policy, never route semantics.
    from source_read_acceptance.codec import unique, reject
    with path.open(encoding="ascii",newline="") as stream:
        for line in stream:
            row = json.loads(line,object_pairs_hook=unique,parse_float=reject,parse_constant=reject)
            require(line == codec.compact(row)+"\n", "noncanonical route line")
            yield row
