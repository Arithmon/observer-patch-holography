"""Portable numerical replay must preserve exact scientific claims."""

import pytest

from .audit import compare_replay


def test_roundoff_is_tolerated_but_claims_and_custody_are_exact():
    retained = {
        "affine_class_excluded_exact": True,
        "bound_exact": "1/100",
        "source_bindings": {"model/source.py": "a" * 64},
        "diagnostics": [0.027005827417535756, 14.8126],
    }
    current = {
        **retained,
        "diagnostics": [0.02700582741753576, 14.81260000001],
    }
    compare_replay(retained, current)

    for key, value in (("affine_class_excluded_exact", False),
                       ("bound_exact", "1/101"),
                       ("source_bindings", {"model/source.py": "b" * 64}),
                       ("diagnostics", [0.02700582741753576, 15.0])):
        with pytest.raises(ValueError):
            compare_replay(retained, {**current, key: value})

    with pytest.raises(ValueError):
        compare_replay(retained, {key: value for key, value in current.items()
                                  if key != "source_bindings"})
