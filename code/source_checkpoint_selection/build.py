"""Producer: memoized native response-state planning and exact online choices."""
import argparse
from fractions import Fraction as F
from functools import cache
from itertools import product

from source_temporal_acceptance.build import certificate
from . import codec, pipeline


def basis(rows):
    pivots = {}
    for row in rows:
        v = list(map(F, row))
        for j, r in sorted(pivots.items()):
            scale = v[j]
            v = [a-scale*b for a, b in zip(v, r)]
        j = next((j for j, value in enumerate(v) if value), None)
        if j is not None:
            pivots[j] = tuple(x/v[j] for x in v)
    if len(pivots) == 2:
        return ((F(1), F(0)), (F(0), F(1)))
    return tuple(pivots.values())


class Planner:
    def __init__(self, spec):
        ports, weights, targets = spec["ports"], spec["weights"], spec["targets"]
        if (type(spec["horizon"]) is not int or spec["horizon"] < 0
                or type(ports) not in (list, tuple) or len(ports) < 2
                or not all(type(p) is int and p >= 0 for p in ports)
                or len(set(ports)) != len(ports)
                or type(weights) not in (list, tuple) or len(weights) != len(ports)-1
                or not all(type(w) is int and w > 0 for w in weights)
                or type(targets) not in (list, tuple)
                or not all(type(t) in (list, tuple) and len(t) == 2
                           and all(type(x) in (int, F) for x in t) for t in targets)):
            raise ValueError("invalid native checkpoint specification")
        self.spec = spec
        self.size = len(spec["ports"])
        self.weights = tuple(spec["weights"])
        self.initial = tuple((F(i == 0), F(i == 1)) for i in range(self.size))
        self.targets = basis(spec["targets"])

        @cache
        def mass(left, state, observed):
            if type(left) is not int or left < 0:
                raise ValueError("invalid remaining horizon")
            if self.complete(observed):
                return sum(self.weights)**left
            if not left:
                return 0
            return sum(w*mass(left-1, *self.advance(state, observed, a))
                       for a, w in enumerate(self.weights))
        self.mass = mass

    def complete(self, observed):
        return basis(observed+self.targets) == observed

    def advance(self, state, observed, edge):
        updated = list(state)
        mean = tuple((a+b)/2 for a, b in zip(state[edge], state[edge+1]))
        updated[edge] = updated[edge+1] = mean
        return tuple(updated), basis(observed+(updated[-1],))

    def select(self, selector):
        """Each selector chooses a ticket in the exact current integer CDF."""
        state, observed = self.initial, basis((self.initial[-1],))
        if selector not in ("first", "middle", "last"):
            raise ValueError("invalid ticket selector")
        if self.mass(self.spec["horizon"], state, observed) == 0:
            raise ValueError("infeasible checkpoint")
        word, decisions = [], []
        for left in range(self.spec["horizon"], 0, -1):
            masses = [w*self.mass(left-1, *self.advance(state, observed, a))
                      for a, w in enumerate(self.weights)]
            total = sum(masses)
            if total == 0:
                raise ValueError("infeasible checkpoint")
            ticket = {"first": 0, "middle": total//2, "last": total-1}[selector]
            cumulative = 0
            for a, value in enumerate(masses):
                cumulative += value
                if ticket < cumulative:
                    break
            decisions.append({"masses": masses, "ticket": ticket, "edge": a})
            word.append(a)
            state, observed = self.advance(state, observed, a)
        return word, decisions

    def first_distribution(self):
        totals = [0]*(self.spec["horizon"]+1)

        def visit(left, state, observed, multiplicity):
            if self.complete(observed):
                totals[self.spec["horizon"]-left] += multiplicity*sum(self.weights)**left
            elif left:
                for a, weight in enumerate(self.weights):
                    visit(left-1, *self.advance(state, observed, a), multiplicity*weight)
        visit(self.spec["horizon"], self.initial, basis((self.initial[-1],)), 1)
        return totals


def history(spec, word, decisions):
    ports = spec["ports"]
    state = [[F(i == j) for j in range(2)] for i in range(len(ports))]
    responses = [state[-1][:]]
    for a in word:
        mean = [(x+y)/2 for x, y in zip(state[a], state[a+1])]
        state[a] = state[a+1] = mean
        responses.append(state[-1][:])
    payloads = []
    for payload in product((-1, 0, 1), repeat=2):
        values = [F(3)]*len(ports)
        values[0] += payload[0]
        values[1] += payload[1]
        writers = [-p-1 for p in ports]
        samples, events = [values[-1]], []
        for k, a in enumerate(word):
            value = (values[a]+values[a+1])/2
            events.append([ports[a], ports[a+1], writers[a], writers[a+1],
                           str(values[a]), str(values[a+1]), str(value)])
            values[a] = values[a+1] = value
            writers[a] = writers[a+1] = k
            samples.append(values[-1])
        payloads.append({"payload": list(payload), "samples": list(map(str, samples)),
                         "events": events})
    return {"word": word, "decisions": decisions,
            "responses": [list(map(str, row)) for row in responses],
            "certificates": [[certificate(responses[:k+1], t) for t in spec["targets"]]
                             for k in range(len(responses))], "payloads": payloads}


def case(spec):
    planner = Planner(spec)
    root_basis = basis((planner.initial[-1],))
    total = planner.mass(spec["horizon"], planner.initial, root_basis)
    roots = ([w*planner.mass(spec["horizon"]-1, *planner.advance(planner.initial, root_basis, a))
              for a, w in enumerate(planner.weights)] if spec["horizon"] else [0]*len(planner.weights))
    first = planner.first_distribution()
    histories = [history(spec, *planner.select(s)) for s in ("first", "middle", "last")] if total else []
    return {"spec": spec, "full_word_count": len(planner.weights)**spec["horizon"],
            "reference_mass": sum(planner.weights)**spec["horizon"],
            "accepted_mass": total, "first_edge_masses": roots, "first_publication_masses": first,
            "first_publication_expectation": str(F(sum(k*v for k, v in enumerate(first)), total)) if total else None,
            "selected_mean_count": spec["horizon"] if total else None,
            "selected_mean_reads": 2*spec["horizon"] if total else None,
            "selected_mean_writes": 2*spec["horizon"] if total else None,
            "selected_sample_count": spec["horizon"]+1 if total else None,
            "histories": histories}


def build():
    return {"schema": 1, "scope": codec.SCOPE, "pins": codec.pins(),
            "cases": [case(spec) for spec in codec.cases()], "pipelines": pipeline.make()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=codec.HERE/"controls.json")
    args = parser.parse_args()
    packet = build()
    from pathlib import Path
    Path(args.output).write_bytes(codec.canonical(packet))
    print(codec.digest(packet))


if __name__ == "__main__":
    main()
