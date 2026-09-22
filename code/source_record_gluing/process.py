"""Execute finite clock-budget experiments in RG's generated process theory.

There is no input radius or read menu. A unit flight has speed one; flight
duration is the square root of the emitted rational duration_squared.
The caller chooses only a spatial sampling cutoff and an observation clock.
"""
from fractions import Fraction as Q
from itertools import product


def execute(q, layers=2, intervention=None):
    sites = list(product(range(q), repeat=3))
    delta_squared = Q(2, 2*q+1)
    state = {(-1, i): i+1+int(intervention == i) for i in range(len(sites))}
    for i in range(len(sites)):
        yield ["prepare", i, list(sites[i]), state[-1, i]]
    # The process comparison is between one primitive flight and the clock,
    # not between a requested radius and a coordinate-labelled target graph.
    flights = []
    incoming = {i: [] for i in range(len(sites))}
    for receiver, y in enumerate(sites):
        for source, x in enumerate(sites):
            squared = Q(sum((a-b)**2 for a, b in zip(x, y)), q*q)
            if squared <= delta_squared:
                incoming[receiver].append(source)
                if source != receiver:
                    flights.append((squared, source, receiver))
    flights.sort()
    for layer in range(layers):
        copied = {}
        for squared, source, receiver in flights:
            value = state[layer-1, source]
            copied[source, receiver] = value
            yield ["fork", layer, source, receiver, value]
        for squared, source, receiver in flights:
            yield ["flight", layer, source, receiver, str(squared), copied[source, receiver]]
        for squared, source, receiver in flights:
            # The wait is sqrt(delta_squared)-sqrt(squared), not their
            # difference. Its nonnegativity is checked without floating point.
            yield ["wait", layer, source, receiver, str(delta_squared), str(squared)]
        for receiver in range(len(sites)):
            value = 0
            for source in incoming[receiver]:
                sample = state[layer-1, source] if source == receiver else copied[source, receiver]
                before = value
                value += sample
                yield ["accumulate", layer, source, receiver, sample, before, value]
            state[layer, receiver] = value
            yield ["commit", layer, receiver, value]
        yield ["checkpoint", layer, [state[layer, i] for i in range(len(sites))]]
