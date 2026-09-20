"""Write complete rounded-mean tapes; the separate verifier certifies them."""
from fractions import Fraction as F
from itertools import product
if __package__:
    from . import codec
else:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from source_reusable_bus import codec

LEVELS = tuple(F(k, 2) for k in (-3, -1, 1, 3))
PORTS = [0, 9, 7, 4, 1, 5, 14, 40, 23, 38, 45, 39]


def copy_edges(s, t):
    return [(2*s, 2*t), (2*s+1, 2*t+1)]


def schedule():
    result, read_positions = [], []
    for record in (0, 1, 0):
        result.extend(copy_edges(record, 2))
        for cell in (2, 3, 4):
            result.extend(copy_edges(cell, cell+1))
        read_positions.append(len(result)-1)
        for _ in range(80):
            result.append((4, 5))
            for cell in (2, 3, 4):
                result.extend(copy_edges(cell, cell+1))
    return result, read_positions


def round_mean(u, v):
    quotient, remainder = divmod(u+v, 2)
    return quotient + (remainder if quotient % 2 else 0)


def decode_local(positive, negative, use, grid):
    observation = F(positive-negative, 2*grid)
    scores = [abs(observation-level/2**(use+3)) for level in LEVELS]
    winners = [i for i, score in enumerate(scores) if score == min(scores)]
    return str(LEVELS[winners[0]]) if len(winners) == 1 else None


def execute(payloads, precision=20):
    grid = 2**precision
    baseline = 2*grid
    state = [baseline]*12
    for record, payload in enumerate(payloads):
        delta = int(payload*grid)
        if F(delta, grid) != payload:
            raise ValueError("payload is not represented on this grid")
        state[2*record] += delta
        state[2*record+1] -= delta
    initial = state.copy()
    writers = [-1-i for i in range(12)]
    tape, reads = [], []
    word, positions = schedule()
    uses = [0, 0]
    for k, (u, v) in enumerate(word):
        output = round_mean(state[u], state[v])
        tape.append([u, v, state[u], state[v], writers[u], writers[v], output])
        state[u] = state[v] = output
        writers[u] = writers[v] = k
        if k in positions:
            record = (0, 1, 0)[positions.index(k)]
            uses[record] += 1
            reads.append({"after_mean": k, "record": record, "use": uses[record],
                          "local_units": state[10:12], "writers": writers[10:12],
                          "decoded": decode_local(*state[10:12], uses[record], grid)})
    return {"payloads": list(map(str, payloads)), "initial_units": initial,
            "tape": tape, "reads": reads, "final_units": state,
            "final_writers": writers,
            "resources": {"scalar_registers": 12, "preparation_writes": 12,
                          "means": len(tape), "mean_reads": 2*len(tape),
                          "cross_carrier_means": sum(PORTS[u]//12 != PORTS[v]//12 for u,v in word),
                          "receiver_scalar_samples": 6, "writes": 12+2*len(tape),
                          "scalar_storage_bits_bound": 12*(precision+3),
                          "mean_adder_bits_bound": precision+4,
                          "receiver_arithmetic_bits_bound": precision+9,
                          "expanded_local_schedule_bits": 8*len(tape),
                          "global_port_map_bits": 12*14,
                          "program_counter_bits": len(tape).bit_length(),
                          "version_use_counter_bits": 4,
                          "read_schedule_bits": 42, "receiver_address_bits": 8,
                          "request_record_bits": 3}}


def build():
    return {"schema": "oph.source_reusable_bus.controls.v1",
            "source_sha256": codec.pins(), "global_ports": PORTS,
            "precision": 20, "baseline": "2", "cleanup_sweeps": 80,
            "requests": [0, 1, 0],
            "tape_columns": ["u", "v", "input_u", "input_v", "writer_u", "writer_v", "output"],
            "assumed_error_bounds": {"initial": "1/262144", "extra_per_mean": "1/268435456",
                                     "terminal_per_rail": "1/65536"},
            "executions": [execute(pair) for pair in product(LEVELS, repeat=2)]}


if __name__ == "__main__":
    (codec.HERE/"controls.json").write_bytes(codec.canonical(build()))
    print("Wrote complete rounded tapes; independent verification is required.")
