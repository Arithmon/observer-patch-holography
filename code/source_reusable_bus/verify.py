"""Independent captured-edge, integer-tape and exact linear-map verification."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
if __package__:
    from . import codec
else:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_reusable_bus import codec


def require(condition, label):
    if not condition:
        raise ValueError(label)


def operation(k):
    cycle, phase = divmod(k, 568)
    if phase < 8:
        transfer, rail = divmod(phase, 2)
        source = (0, 1, 0)[cycle] if transfer == 0 else transfer+1
        target = transfer+2
        return 2*source+rail, 2*target+rail
    stage = (phase-8) % 7
    if stage == 0:
        return 4, 5
    hop, rail = divmod(stage-1, 2)
    return 4+2*hop+rail, 6+2*hop+rail


def support_edges():
    support = codec.load(codec.ROOT/codec.SUPPORT)
    edges = set()
    for carrier in range(support["carriers"]):
        for u, v in support["intra_carrier_seams"]:
            edges.add(tuple(sorted((12*carrier+u, 12*carrier+v))))
    for left, lp, right, rp in support["glued_pairs"]:
        edges.add(tuple(sorted((12*left+lp, 12*right+rp))))
    return edges, 12*support["carriers"]


def layout(ports):
    # Each rail path is listed independently of the producer's flat register map.
    positive = (1, 14, 23, 45)
    negative = (5, 40, 38, 39)
    expected = [0, 9, 7, 4]
    for p, m in zip(positive, negative):
        expected.extend((p, m))
    codec.equal(ports, expected, "declared rail embedding")
    require(len(set(ports)) == 12, "rail aliases")
    edges, vertices = support_edges()
    require(all(0 <= p < vertices for p in ports), "port outside support")
    require(all(p//12 == 0 for p in ports[:6]), "source carrier")
    require(all(p//12 == 3 for p in ports[10:]), "receiver carrier")
    for k in range(1704):
        u, v = operation(k)
        require(tuple(sorted((ports[u], ports[v]))) in edges, "non-native seam")
    return {"source_carrier": 0, "receiver_carrier": 3,
            "participating_carriers": sorted({p//12 for p in ports}),
            "positive_path": list(positive), "negative_path": list(negative),
            "registered_port_count": vertices}


def averaged_rows(rows, u, v):
    row = tuple((a+b)/2 for a, b in zip(rows[u], rows[v]))
    rows[u] = rows[v] = row


@lru_cache(maxsize=1)
def ideal_maps():
    # Columns are arbitrary initial source amplitudes, not selected payload values.
    rows = [(F(0), F(0)) for _ in range(12)]
    rows[0], rows[1] = (F(1), F(0)), (F(-1), F(0))
    rows[2], rows[3] = (F(0), F(1)), (F(0), F(-1))
    reads, scrub_ends = [], []
    for k in range(1704):
        averaged_rows(rows, *operation(k))
        if k % 568 == 7:
            require(all(rows[2*i] == tuple(-v for v in rows[2*i+1]) for i in range(6)),
                    "ideal paired encoding")
            reads.append(tuple((p-m)/2 for p, m in zip(rows[10], rows[11])))
        if k % 568 == 567:
            scrub_ends.append(tuple(rows))
    return tuple(reads), tuple(scrub_ends)


def contraction_certificate():
    # Six arbitrary amplitudes exercise the full scrub operator, including archives.
    rows = [tuple(F((1 if rail == 0 else -1) if j == i else 0)
                  for j in range(6)) for i in range(6) for rail in range(2)]
    for k in range(8, 15):
        averaged_rows(rows, *operation(k))
    expected = [(1,0,0,0,0,0), (0,1,0,0,0,0), (0,0,0,F(1,2),0,0),
                (0,0,0,F(1,4),F(1,2),0),
                (0,0,0,F(1,8),F(1,4),F(1,2)),
                (0,0,0,F(1,8),F(1,4),F(1,2))]
    require([rows[2*i] for i in range(6)] == expected, "scrub matrix differs from theorem")
    require(all(rows[2*i] == tuple(-v for v in rows[2*i+1]) for i in range(6)),
            "scrub destroys paired encoding")
    require(all(all(w >= 0 for w in row) and sum(row) <= F(7,8)
                for row in expected[2:]), "scrub max-norm contraction")
    return [[str(v) for v in row] for row in expected]


def decode_receiver(positive, negative, use, grid):
    scaled = (positive-negative)*2**(use+3)
    cuts = (-2*grid, 0, 2*grid)
    if scaled in cuts:
        return None
    return ("-3/2", "-1/2", "1/2", "3/2")[sum(scaled > cut for cut in cuts)]


def read_bounds(cycle, coefficients):
    record, use = ((0,1), (1,1), (0,2))[cycle]
    gain = F(1, 2**(use+3))
    epsilon = F(3,2)*F(7,8)**80
    ideal_error = F(3,2)*sum(abs(weight-(gain if j == record else 0))
                           for j, weight in enumerate(coefficients))
    require(ideal_error <= epsilon, "universal ideal transfer error")
    means = 568*cycle+8
    implementation = F(1,2**18) + means*(F(1,2**28)+F(1,2**21)) + F(1,2**16)
    total = epsilon+implementation
    require(total < gain/2, "insufficient decoding margin")
    return {"record": record, "use": use, "means_before_read": means,
            "ideal_coefficients": list(map(str, coefficients)),
            "ideal_uniform_error": str(ideal_error), "cleanup_bound": str(epsilon),
            "initial_plus_rounding_disturbance_readout_bound": str(implementation),
            "total_error_bound": str(total), "half_spacing": str(gain/2),
            "margin": str(gain/2-total),
            "ideal_other_record_influence": coefficients[1-record] != 0}


def check_case(case, payloads, maps, scrub_ends):
    q = 2**20
    values = [2*q]*12
    for record, amplitude in enumerate(payloads):
        values[2*record] += int(amplitude*q)
        values[2*record+1] -= int(amplitude*q)
    codec.equal(case["payloads"], list(map(str, payloads)), "payload enumeration")
    codec.equal(case["initial_units"], values, "preparation or payload-dependent blank")
    require(len(case["tape"]) == 1704, "missing/extra mean or cleanup")
    writers = [-1-i for i in range(12)]
    ancestors = [1 << i for i in range(12)]
    reads, audit = [], []
    max_rounding = F(0)
    for k, event in enumerate(case["tape"]):
        require(isinstance(event, list) and len(event) == 7 and
                all(type(v) is int for v in event), "invalid tape cell")
        u, v = operation(k)
        exact = F(values[u]+values[v], 2)
        floor = exact.numerator // exact.denominator
        output = min((floor, floor+1), key=lambda n: (abs(F(n)-exact), n % 2))
        expected = [u,v,values[u],values[v],writers[u],writers[v],output]
        codec.equal(event, expected, "native mean, rounding, value or consumed writer")
        local_error = abs(F(output)-exact)/q
        require(local_error <= F(1,2*q), "rounding allowance")
        max_rounding = max(max_rounding, local_error)
        require(0 <= output <= 4*q, "raw-load storage range")
        ancestry = ancestors[u] | ancestors[v] | (1 << (12+k))
        values[u] = values[v] = output
        writers[u] = writers[v] = k
        ancestors[u] = ancestors[v] = ancestry
        if k % 568 == 7:
            cycle = k//568
            record, use = ((0,1), (1,1), (0,2))[cycle]
            ideal = sum(w*a for w, a in zip(maps[cycle], payloads))
            observed = F(values[10]-values[11], 2*q)
            require(abs(observed-ideal) <= F(k+1,2*q), "global rounded error bound")
            decoded = decode_receiver(values[10], values[11], use, q)
            require(decoded == str(payloads[record]), "receiver decoded wrong payload")
            reads.append({"after_mean": k, "record": record, "use": use,
                          "local_units": values[10:12], "writers": writers[10:12],
                          "decoded": decoded})
            lineage = ancestors[10] | ancestors[11]
            initial_ports = [i for i in range(12) if lineage & (1 << i)]
            if cycle > 0:
                require(all(i in initial_ports for i in range(4)),
                        "cleanup erased raw cross-record ancestry")
            audit.append({"decoded": decoded, "observed_contrast": str(observed),
                          "observed_rounding_error": str(abs(observed-ideal)),
                          "initial_ancestor_ports": initial_ports,
                          "mean_ancestor_count": (lineage >> 12).bit_count()})
        if k % 568 == 567:
            epsilon = F(3,2)*F(7,8)**80
            for i in range(4,12):
                ideal = sum(w*a for w,a in zip(scrub_ends[k//568][i], payloads))
                require(abs(ideal) <= epsilon, "ideal cleanup residual")
                require(abs(F(values[i],q)-2-ideal) <= F(k+1,2*q),
                        "rounded cleanup allowance")
    codec.equal(case["reads"], reads, "receiver read interface")
    codec.equal(case["final_units"], values, "retained final loads")
    codec.equal(case["final_writers"], writers, "retained final writers")
    resources = {"scalar_registers": 12, "preparation_writes": 12, "means": 1704,
                 "mean_reads": 3408, "cross_carrier_means": 729,
                 "receiver_scalar_samples": 6, "writes": 3420,
                 "scalar_storage_bits_bound": 276, "mean_adder_bits_bound": 24,
                 "receiver_arithmetic_bits_bound": 29,
                 "expanded_local_schedule_bits": 13632, "global_port_map_bits": 168,
                 "program_counter_bits": 11, "version_use_counter_bits": 4,
                 "read_schedule_bits": 42, "receiver_address_bits": 8,
                 "request_record_bits": 3}
    codec.equal(case["resources"], resources, "resource accounting")
    require(set(case) == {"payloads", "initial_units", "tape", "reads", "final_units",
                          "final_writers", "resources"}, "case fields")
    return {"payloads": list(map(str,payloads)), "reads": audit,
            "max_observed_rounding_error": str(max_rounding), "resources": resources}


def verify(data=None):
    if data is None:
        data = codec.load(codec.HERE/"controls.json")
    expected_fields = {"schema", "source_sha256", "global_ports", "precision", "baseline",
                       "cleanup_sweeps", "requests", "tape_columns", "assumed_error_bounds",
                       "executions"}
    require(set(data) == expected_fields, "control fields")
    codec.equal(data["schema"], "oph.source_reusable_bus.controls.v1", "schema")
    codec.equal(data["source_sha256"], codec.pins(), "source pins")
    for key, value in {"precision":20, "baseline":"2", "cleanup_sweeps":80,
                       "requests":[0,1,0],
                       "tape_columns":["u","v","input_u","input_v","writer_u","writer_v","output"],
                       "assumed_error_bounds":{"initial":"1/262144", "extra_per_mean":"1/268435456",
                                               "terminal_per_rail":"1/65536"}}.items():
        codec.equal(data[key], value, f"declared {key}")
    support = layout(data["global_ports"])
    matrix = contraction_certificate()
    maps, scrub_ends = ideal_maps()
    bounds = [read_bounds(i, row) for i, row in enumerate(maps)]
    require(bounds[2]["ideal_other_record_influence"], "false raw version independence")
    for state in scrub_ends:
        require(all(F(3,2)*sum(abs(v) for v in row) <= F(3,2)*F(7,8)**80
                    for row in state[4:]), "uniform cleanup coefficient bound")
    payloads = list(product((F(-3,2),F(-1,2),F(1,2),F(3,2)), repeat=2))
    require(len(data["executions"]) == len(payloads), "missing payload history")
    cases = [check_case(case, pair, maps, scrub_ends)
             for case, pair in zip(data["executions"],payloads)]
    return {"schema":"oph.source_reusable_bus.verified.v1", "source_sha256":codec.pins(),
            "support":support, "scrub_matrix":matrix, "uniform_read_bounds":bounds,
            "histories":cases, "retained_means":1704*16,
            "scope":"Finite declared four-level records; raw influence and supplied control remain."}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = verify()
    if args.write_receipt:
        (codec.HERE/"receipt.json").write_bytes(codec.canonical(result))
    else:
        codec.equal(codec.load(codec.HERE/"receipt.json"), result, "retained receipt")
    print(f"Verified {len(result['histories'])} histories, {result['retained_means']} means, "
          "captured seams, all writer custody, and positive universal decoding margins.")
