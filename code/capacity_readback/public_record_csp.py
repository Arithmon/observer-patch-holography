#!/usr/bin/env python3
"""Exact finite record gluing, with no label coercion or recursive search.

Self-interfaces are equalizers on one atom, not relations between two
independent copies of that atom. Arc consistency only prunes impossible
values; branching still checks globally obstructed cycles.
"""
from __future__ import annotations

from collections import deque
import re
from typing import Any, Mapping, Sequence


def _label(value: Any) -> bool:
    return type(value) is str and bool(value)


def _atoms(values: Any, name: str) -> tuple[str, ...]:
    if (not isinstance(values, Sequence) or isinstance(values, (str, bytes))
            or not values or any(not _label(x) for x in values)
            or len(set(values)) != len(values)):
        raise ValueError(f"{name} must be a nonempty sequence of distinct string atoms")
    return tuple(values)


def _escape(value: str) -> str:
    return value.replace("%", "%25").replace("|", "%7C").replace("=", "%3D")


def encode_section_id(section: Mapping[str, str]) -> str:
    """Injective identifier, preserving legacy bytes for delimiter-free labels."""
    if (not isinstance(section, Mapping) or not section
            or any(not _label(x) or not _label(y) for x, y in section.items())):
        raise ValueError("a section needs nonempty string observer and atom labels")
    return "|".join(f"{_escape(x)}={_escape(section[x])}" for x in sorted(section))


def decode_section_id(identifier: str) -> dict[str, str]:
    """Decode only the canonical spelling; reject aliases and duplicate keys."""
    if not _label(identifier):
        raise ValueError("invalid section identifier")
    escapes = {"25": "%", "7C": "|", "3D": "="}
    section: dict[str, str] = {}
    for entry in identifier.split("|"):
        if entry.count("=") != 1:
            raise ValueError("invalid section identifier")
        left, right = (re.sub(r"%(25|7C|3D)", lambda m: escapes[m[1]], part)
                       for part in entry.split("="))
        if left in section:
            raise ValueError("duplicate observer in section identifier")
        section[left] = right
    if encode_section_id(section) != identifier:
        raise ValueError("section identifier is not canonical")
    return section


def _diagram(observers, interfaces):
    if (not isinstance(observers, Mapping) or not observers
            or any(not _label(x) for x in observers)):
        raise ValueError("observers must have nonempty string labels")
    domains = {x: _atoms(observers[x], f"atoms of {x}") for x in sorted(observers)}
    if not isinstance(interfaces, Sequence) or isinstance(interfaces, (str, bytes)):
        raise ValueError("interfaces must be a sequence of readout declarations")
    constraints = []
    for interface in interfaces:
        if not isinstance(interface, Mapping):
            raise ValueError("each interface must be a readout declaration")
        left, right = interface.get("left_observer"), interface.get("right_observer")
        if not _label(left) or not _label(right) or left not in domains or right not in domains:
            raise ValueError("interface references an unknown observer")
        maps = []
        codomain = (set(_atoms(interface["interface_atoms"], "interface atoms"))
                    if "interface_atoms" in interface else None)
        for endpoint, field in [(left, "left_readout"), (right, "right_readout")]:
            readout = interface.get(field)
            if not isinstance(readout, Mapping) or set(readout) != set(domains[endpoint]):
                raise ValueError("atom readout maps must be total on endpoint atoms")
            if any(not _label(y) for y in readout.values()):
                raise ValueError("readout values must be nonempty string atom labels")
            if codomain is not None and not set(readout.values()) <= codomain:
                raise ValueError("readout does not land in the declared interface atoms")
            maps.append(dict(readout))
        constraints.append((left, right, *maps))
    return domains, constraints


def _propagate(domains, arcs, incoming):
    """Delete only atoms without support across a declared directed interface."""
    queue = deque(range(len(arcs)))
    queued = set(queue)
    while queue:
        index = queue.popleft()
        queued.remove(index)
        source, target, read_source, read_target = arcs[index]
        supported = {read_target[x] for x in domains[target]}
        reduced = tuple(x for x in domains[source] if read_source[x] in supported)
        if not reduced:
            return False
        if reduced != domains[source]:
            domains[source] = reduced
            for affected in incoming[source]:
                if affected not in queued:
                    queue.append(affected)
                    queued.add(affected)
    return True


def public_global_sections_csp(
    observers: Mapping[str, Sequence[str]],
    interfaces: Sequence[Mapping[str, Any]],
    *,
    max_search_nodes: int = 100_000,
    max_sections: int = 100_000,
) -> list[dict[str, str]]:
    """Enumerate every compatible section or refuse an exhausted resource bound.

    A self-interface first restricts its single domain to equal readouts.
    Binary interfaces propagate supported domains to a fixed point. Surviving
    multivalued domains are split by an explicit depth-first stack, so arc
    consistency is never confused with global satisfiability. The output is
    a complete list, never a truncated prefix; its order is not contractual.
    """
    if any(type(n) is not int or n < 1 for n in (max_search_nodes, max_sections)):
        raise ValueError("search and section limits must be positive integers")
    domains, constraints = _diagram(observers, interfaces)
    incoming = {x: [] for x in domains}
    neighbors = {x: [] for x in domains}
    arcs = []
    for left, right, read_left, read_right in constraints:
        neighbors[left].append(right)
        neighbors[right].append(left)
        if left == right:
            domains[left] = tuple(x for x in domains[left] if read_left[x] == read_right[x])
        else:
            incoming[right].append(len(arcs))
            arcs.append((left, right, read_left, read_right))
            incoming[left].append(len(arcs))
            arcs.append((right, left, read_right, read_left))
    if any(not values for values in domains.values()):
        return []
    sections = []
    stack = [domains]
    allocated = 1
    while stack:
        current = stack.pop()
        if not _propagate(current, arcs, incoming):
            continue
        choices = [x for x in current if len(current[x]) > 1]
        if not choices:
            if len(sections) == max_sections:
                raise ValueError("section limit exhausted; no complete record set certified")
            sections.append({x: values[0] for x, values in current.items()})
            continue
        observer = max(choices, key=lambda x: (
            sum(len(current[y]) == 1 for y in neighbors[x]), len(neighbors[x]), x))
        for atom in reversed(current[observer]):
            allocated += 1
            if allocated > max_search_nodes:
                raise ValueError("search budget exhausted; no complete record set certified")
            branch = current.copy()
            branch[observer] = (atom,)
            stack.append(branch)
    return sections
