"""Root replay checked against direct equations on unequal local alphabets."""
from itertools import permutations, product

import pytest

from public_record_csp import public_global_sections_csp
from verify_public_record_sections import verify_public_record_sections


@pytest.mark.parametrize("extra", ["none", "cycle", "parallel", "self"])
def test_all_injective_chain_maps_with_disconnected_component(extra):
    observers = {"a": ["a0", "a1"], "b": ["b0", "b1", "b2"],
                 "c": ["c0", "c1"], "d": ["d0", "d1"]}
    ids = tuple(observers)
    images = ("x", "y", "z")

    def edge(left, right, a, b):
        return {"left_observer": left, "right_observer": right,
                "left_readout": dict(zip(observers[left], a)),
                "right_readout": dict(zip(observers[right], b)),
                "interface_atoms": [*images, "unused"]}

    all_maps = [tuple(permutations(images, size)) for size in (2, 3, 3, 2)]
    empty = nonempty = 0
    for a, b, c, d in product(*all_maps):
        seams = [edge("a", "b", a, b), edge("b", "c", c, d)]
        if extra == "cycle":
            seams.append(edge("c", "a", images[:2], images[:2]))
        elif extra == "parallel":
            seams.append(edge("a", "b", images[:2], images))
        elif extra == "self":
            seams.append(edge("b", "b", images, ("x", "z", "y")))

        # Enumerate all 24 assignments independently of both implementations.
        expected = []
        for values in product(*(observers[x] for x in ids)):
            section = dict(zip(ids, values))
            if all(s["left_readout"][section[s["left_observer"]]] ==
                   s["right_readout"][section[s["right_observer"]]] for s in seams):
                expected.append(section)
        result = public_global_sections_csp(observers, seams)
        assert {tuple(s[x] for x in ids) for s in result} == {
            tuple(s[x] for x in ids) for s in expected}
        assert len(result) == len(expected)

        # Four root tuples (a,d) force the alternate replay: Cartesian size is
        # 24. Roots need not have the same alphabet size as their neighbors,
        # and their injective readouts need not have identical images.
        assert verify_public_record_sections(observers, seams, expected, max_candidates=4)
        if expected:
            nonempty += 1
            assert not verify_public_record_sections(
                observers, seams, expected[:-1], max_candidates=4)
        else:
            empty += 1
        # A well-formed extra local assignment still must satisfy every seam.
        invalid = next((dict(zip(ids, values))
                        for values in product(*(observers[x] for x in ids))
                        if dict(zip(ids, values)) not in expected), None)
        assert invalid is not None
        assert not verify_public_record_sections(
            observers, seams, [*expected, invalid], max_candidates=4)
    assert nonempty > 0
    if extra != "none":
        assert empty > 0
