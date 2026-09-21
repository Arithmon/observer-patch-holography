"""Independent full-word census, scalar intervention replay and policy validation.

Does not import the producer, its elimination routine or its dynamic planner.
The census uses integer-scaled scalar basis runs and determinant tests. Online
transition masses are reconstructed by summing all accepted full words.
"""
import argparse
from fractions import Fraction as F
import hashlib
from itertools import product
from math import comb
from pathlib import Path

from source_temporal_acceptance.verify import adjacency, check_certificate, rational, require
from . import codec


def determined(pivot, full, targets):
    if full:
        return True
    if pivot is None:
        return all(not any(t) for t in targets)
    return all(pivot[0]*t[1] == pivot[1]*t[0] for t in targets)


def classify(size, word, targets):
    # Independent scalar preparations; at step k each entry is scaled by 2**k.
    columns = [[int(i == j) for i in range(size)] for j in range(2)]
    initial = [column[-1] for column in columns]
    pivot, full = initial if any(initial) else None, False
    first = 0 if determined(pivot, full, targets) else None
    for k, a in enumerate(word, 1):
        row = []
        for column in columns:
            mean_scaled = column[a]+column[a+1]
            for i in range(size):
                column[i] *= 2
            column[a] = column[a+1] = mean_scaled
            row.append(column[-1])
        if any(row):
            if pivot is None:
                pivot = row
            elif pivot[0]*row[1] != pivot[1]*row[0]:
                full = True
        if first is None and determined(pivot, full, targets):
            first = k
    return first


def census(spec, prefixes):
    h, weights = spec["horizon"], spec["weights"]
    first_counts, roots = [0]*(h+1), [0]*len(weights)
    prefix_mass = {p: 0 for p in prefixes}
    stream = hashlib.sha256()
    accepted_count = 0
    for word in product(range(len(weights)), repeat=h):
        first = classify(len(spec["ports"]), word, spec["targets"])
        weight = 1
        for a in word:
            weight *= weights[a]
        stream.update(codec.canonical([word, first, weight]))
        if first is None:
            continue
        accepted_count += 1
        first_counts[first] += weight
        roots[word[0]] += weight
        for k in range(h+1):
            prefix = word[:k]
            if prefix in prefix_mass:
                prefix_mass[prefix] += weight
    return first_counts, roots, prefix_mass, accepted_count, stream.hexdigest()


def scalar_trace(spec, word, payload):
    ports = spec["ports"]
    values = {p: F(3) for p in ports}
    values[ports[0]] += payload[0]
    values[ports[1]] += payload[1]
    writers = {p: -p-1 for p in ports}
    samples, events = ["3"], []
    for k, edge in enumerate(word):
        a, b = ports[edge:edge+2]
        after = (values[a]+values[b])/2
        events.append([a, b, writers[a], writers[b], str(values[a]), str(values[b]), str(after)])
        values[a] = values[b] = after
        writers[a] = writers[b] = k
        samples.append(str(values[ports[-1]]))
    return {"payload": list(payload), "samples": samples, "events": events}


def check_history(spec, item, selector, prefixes, accepted_mass):
    require(set(item) == {"word", "decisions", "responses", "certificates", "payloads"}, "history fields")
    word, h = item["word"], spec["horizon"]
    require(type(word) is list and len(word) == h and all(
        type(a) is int and 0 <= a < len(spec["weights"]) for a in word), "native word")
    require(type(item["decisions"]) is list and len(item["decisions"]) == h, "decision count")
    prefix, factor, probability = (), 1, F(1)
    for k, a in enumerate(word):
        masses = [prefixes[prefix+(b,)]//factor for b in range(len(spec["weights"]))]
        require(all(prefixes[prefix+(b,)] % factor == 0 for b in range(len(spec["weights"]))),
                "nonintegral continuation mass")
        total = sum(masses)
        require(total > 0, "dead prefix")
        ticket = (0, total//2, total-1)[selector]
        chosen = next(b for b in range(len(masses)) if ticket < sum(masses[:b+1]))
        codec.equal(item["decisions"][k], {"masses": masses, "ticket": ticket, "edge": chosen},
                    "continuation mass or ticket")
        require(a == chosen, "policy choice")
        probability *= F(masses[a], total)
        prefix += (a,)
        factor *= spec["weights"][a]
    require(probability == F(factor, accepted_mass), "telescoping policy law")
    expected_payloads = [scalar_trace(spec, word, p) for p in product((-1, 0, 1), repeat=2)]
    codec.equal(item["payloads"], expected_payloads, "scalar samples, values, writers or interventions")
    columns = [scalar_trace(spec, word, p)["samples"] for p in ((1, 0), (0, 1))]
    rows = [[F(column[k])-3 for column in columns] for k in range(h+1)]
    codec.equal(item["responses"], [list(map(str, row)) for row in rows], "native responses")
    require(type(item["certificates"]) is list and len(item["certificates"]) == h+1, "prefix certificates")
    first = None
    for k, certs in enumerate(item["certificates"]):
        require(type(certs) is list and len(certs) == len(spec["targets"]), "record grammar")
        accepted = [check_certificate(rows[:k+1], target, cert)
                    for target, cert in zip(spec["targets"], certs)]
        if all(accepted) and first is None:
            first = k
        for target, cert, ok in zip(spec["targets"], certs, accepted):
            if not ok:
                continue
            coefficients = [rational(x) for x in cert["coefficients"]]
            for payload in item["payloads"]:
                decoded = sum(c*(F(s)-3) for c, s in zip(coefficients, payload["samples"][:k+1]))
                require(decoded == sum(a*b for a, b in zip(target, payload["payload"])), "decoded intervention")
    require(first == classify(len(spec["ports"]), word, spec["targets"]), "first publication")
    require(first is not None, "unpublished selected word")


def check_case(item, spec):
    require(type(item) is dict, "case object")
    codec.equal(item.get("spec"), spec, "fixed source experiment")
    histories = item.get("histories")
    require(type(histories) is list and len(histories) in (0, 3), "history controls")
    prefixes = {()}
    for row in histories:
        require(type(row) is dict and type(row.get("word")) is list, "history shape")
        word = row["word"]
        require(len(word) == spec["horizon"] and all(type(a) is int and 0 <= a < len(spec["weights"])
                                                    for a in word), "word shape")
        for k in range(len(word)):
            prefixes.update(tuple(word[:k])+ (a,) for a in range(len(spec["weights"])))
    first, roots, masses, count, stream = census(spec, prefixes)
    total = sum(first)
    expected = {"spec": spec, "full_word_count": len(spec["weights"])**spec["horizon"],
                "reference_mass": sum(spec["weights"])**spec["horizon"],
                "accepted_mass": total, "first_edge_masses": roots, "first_publication_masses": first,
                "first_publication_expectation": str(F(sum(k*v for k, v in enumerate(first)), total)) if total else None,
                "selected_mean_count": spec["horizon"] if total else None,
                "selected_mean_reads": 2*spec["horizon"] if total else None,
                "selected_mean_writes": 2*spec["horizon"] if total else None,
                "selected_sample_count": spec["horizon"]+1 if total else None,
                "histories": histories}
    codec.equal(item, expected, "full-word census, reference, publication or work")
    require(len(histories) == (3 if total else 0), "missing or infeasible selected histories")
    for selector, history in enumerate(histories):
        check_history(spec, history, selector, masses, total)
    return {"name": spec["name"], "horizon": spec["horizon"], "accepted_word_count": count,
            "accepted_reference_probability": str(F(total, expected["reference_mass"])),
            "full_census_sha256": stream}


def pipeline_extensions(depth):
    """Independent enumeration of a two-wave precedence poset."""
    requirements = {}
    for i in range(1, depth):
        requirements[("first", i)] = {("first", i-1)} if i > 1 else set()
    for i in range(depth):
        parents = {("second", i-1)} if i else set()
        if i < depth-1:
            parents.add(("first", i+1))
        requirements[("second", i)] = parents

    def extend(done, word):
        if len(done) == len(requirements):
            yield tuple(word)
        for event, parents in requirements.items():
            if event not in done and parents <= done:
                yield from extend(done | {event}, word+[event[1]])
    return sorted(extend(set(), []))


def check_pipelines(items):
    require(type(items) is list and len(items) == 10, "pipeline depth family")
    receipts = []
    for depth, item in enumerate(items, 1):
        words = pipeline_extensions(depth)
        require(len(words) == len(set(words)), "duplicate event histories")
        count = comb(2*depth-2, depth-1)//depth
        require(len(words) == count, "Catalan extension census")
        stream = hashlib.sha256()
        for word in words:
            columns = [[int(i == j) for i in range(depth+1)] for j in range(2)]
            first = [F(0), F(1)] if depth == 1 else None
            first_slot = 0 if depth == 1 else None
            for k, edge in enumerate(word, 1):
                for values in columns:
                    joined = values[edge]+values[edge+1]
                    for i in range(depth+1): values[i] *= 2
                    values[edge] = values[edge+1] = joined
                if first is None and edge == depth-1:
                    first = [F(values[-1], 2**k) for values in columns]
                    first_slot = k
            last = [F(values[-1], 2**len(word)) for values in columns]
            require(first == [0, F(1, 2**(depth-1))], "first native wave")
            require(last == [F(1, 2**depth), F(depth+1, 2**(depth+1))], "second native wave")
            require(word.count(0) == 1 and all(word.count(a) == 2 for a in range(1, depth)),
                    "native crossing work")
            stream.update(codec.canonical([word, first_slot, list(map(str, first)), list(map(str, last))]))
        expected = {"depth": depth, "means": 2*depth-1, "samples": 2*depth,
                    "successful_words": count, "word_response_sha256": stream.hexdigest(),
                    "first_response": ["0", str(F(1, 2**(depth-1)))],
                    "final_response": [str(F(1, 2**depth)), str(F(depth+1, 2**(depth+1)))],
                    "first_record_noise_gain": str(F((depth+5)*2**depth, 4)),
                    "second_record_noise_gain": str(2**(depth-1))}
        codec.equal(item, expected, "native pipeline theorem control")
        # At depths <=5 this scans every minimal-length word, not just the poset.
        exhaustive = None
        if depth <= 5:
            exhaustive = [w for w in product(range(depth), repeat=2*depth-1)
                          if classify(depth+1, w, [[1, 0], [0, 1]]) is not None]
            require(exhaustive == words, "optimal-language completeness")
        receipts.append({"depth": depth, "positive_native_words": len(words),
                         "all_words_checked": depth**(2*depth-1) if exhaustive is not None else None,
                         "necessity_scope": "exhaustive" if exhaustive is not None else "analytic theorem"})
    return receipts


def verify(packet):
    require(type(packet) is dict and set(packet) == {"schema", "scope", "pins", "cases", "pipelines"}, "packet schema")
    codec.equal(packet["schema"], 1, "schema version")
    codec.equal(packet["pins"], codec.pins(), "source pins")
    codec.equal(packet["scope"], codec.SCOPE, "scientific scope")
    specs = codec.cases()
    require(type(packet["cases"]) is list and len(packet["cases"]) == len(specs), "complete experiment family")
    edges = adjacency()
    for spec in specs:
        if spec["name"] != "chain":
            require(all(frozenset((a, b)) in edges for a, b in zip(spec["ports"], spec["ports"][1:])),
                    "captured seam provenance")
    results = [check_case(item, spec) for item, spec in zip(packet["cases"], specs)]
    return {"schema": 1, "controls_sha256": codec.digest(packet), "cases": results,
            "pipelines": check_pipelines(packet["pipelines"]),
            "scope": codec.SCOPE, "source_pins_sha256": codec.digest(packet["pins"])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls", type=Path, default=codec.HERE/"controls.json")
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    receipt = verify(codec.load(args.controls))
    target = codec.HERE/"receipt.json"
    if args.write_receipt:
        target.write_bytes(codec.canonical(receipt))
    else:
        codec.equal(codec.load(target), receipt, "verification receipt")
    print(codec.digest(receipt))


if __name__ == "__main__":
    main()
