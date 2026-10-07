"""Independent completeness replay for finite public record sections.

This module imports neither the record solver nor the capacity evaluator.
It enumerates the Cartesian product when small. For larger diagrams with
injective readouts it reconstructs each component from a root atom, then
checks every interface, including loops and edges outside the spanning tree.
"""
from collections.abc import Mapping, Sequence
from itertools import product
from math import prod


class ReplayBudgetExceeded(ValueError):
    """The independent procedure cannot certify completeness at this budget."""


def _strings(values):
    if (not isinstance(values, Sequence) or isinstance(values, (str, bytes))
            or not values or any(type(x) is not str or not x for x in values)
            or len(set(values)) != len(values)):
        raise ValueError("invalid atom alphabet")
    return tuple(values)


def _read_diagram(observers, interfaces):
    if (not isinstance(observers, Mapping) or not observers
            or any(type(x) is not str or not x for x in observers)):
        raise ValueError("invalid observers")
    domains = {x: _strings(values) for x, values in observers.items()}
    if not isinstance(interfaces, Sequence) or isinstance(interfaces, (str, bytes)):
        raise ValueError("invalid interface declarations")
    edges = []
    for seam in interfaces:
        if not isinstance(seam, Mapping):
            raise ValueError("invalid interface")
        left, right = seam["left_observer"], seam["right_observer"]
        if type(left) is not str or type(right) is not str or left not in domains or right not in domains:
            raise ValueError("unknown endpoint")
        a, b = seam["left_readout"], seam["right_readout"]
        if (not isinstance(a, Mapping) or not isinstance(b, Mapping)
                or set(a) != set(domains[left]) or set(b) != set(domains[right])):
            raise ValueError("readouts must be total")
        image = list(a.values())+list(b.values())
        if any(type(y) is not str or not y for y in image):
            raise ValueError("invalid interface atom")
        if "interface_atoms" in seam and not set(image) <= set(_strings(seam["interface_atoms"])):
            raise ValueError("invalid readout codomain")
        edges.append((left, right, dict(a), dict(b)))
    return domains, edges


def _expected(domains, edges, limit):
    ids = sorted(domains)

    def consistent(section):
        return all(a[section[left]] == b[section[right]] for left, right, a, b in edges)

    if prod(len(domains[x]) for x in ids) <= limit:
        return {values for values in product(*(domains[x] for x in ids))
                if consistent(dict(zip(ids, values)))}

    adjacency = {x: [] for x in ids}
    for left, right, a, b in edges:
        if left == right:
            continue
        if len(set(a.values())) != len(a) or len(set(b.values())) != len(b):
            raise ReplayBudgetExceeded("noninjective diagram exceeds exhaustive replay budget")
        adjacency[left].append((right, a, {y: x for x, y in b.items()}))
        adjacency[right].append((left, b, {y: x for x, y in a.items()}))
    roots, seen = [], set()
    for root in ids:
        if root in seen:
            continue
        roots.append(root)
        seen.add(root)
        queue = [root]
        for current in queue:
            for neighbor, _, _ in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
    if prod(len(domains[x]) for x in roots) > limit:
        raise ReplayBudgetExceeded("root assignments exceed independent replay budget")
    expected = set()
    for seeds in product(*(domains[x] for x in roots)):
        assignment = dict(zip(roots, seeds))
        queue = list(roots)
        possible = True
        for current in queue:
            for neighbor, read_current, inverse_neighbor in adjacency[current]:
                value = inverse_neighbor.get(read_current[assignment[current]])
                if value is None or (neighbor in assignment and assignment[neighbor] != value):
                    possible = False
                    break
                if neighbor not in assignment:
                    assignment[neighbor] = value
                    queue.append(neighbor)
            if not possible:
                break
        if possible and consistent(assignment):
            expected.add(tuple(assignment[x] for x in ids))
    return expected


def verify_public_record_sections(observers, interfaces, sections, *, max_candidates=100_000):
    """Verify exactly the complete record set, not merely valid listed records.

    Malformed inputs, duplicate/extra/omitted sections return False. An
    unsupported enumeration size raises ReplayBudgetExceeded; it does not
    declare the supplied section list invalid or certify an empty set.
    """
    if type(max_candidates) is not int or max_candidates < 1:
        raise ValueError("max_candidates must be a positive integer")
    try:
        domains, edges = _read_diagram(observers, interfaces)
        ids = sorted(domains)
        if not isinstance(sections, Sequence) or isinstance(sections, (str, bytes)):
            return False
        actual = set()
        for section in sections:
            if (not isinstance(section, Mapping) or set(section) != set(ids)
                    or any(type(section[x]) is not str or section[x] not in domains[x] for x in ids)):
                return False
            values = tuple(section[x] for x in ids)
            if values in actual:
                return False
            actual.add(values)
        return actual == _expected(domains, edges, max_candidates)
    except ReplayBudgetExceeded:
        raise
    except (KeyError, TypeError, ValueError):
        return False
