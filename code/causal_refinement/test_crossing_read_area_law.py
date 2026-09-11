"""Exact crossing counts, mutation guards and receipt checks for the seam-count area law."""
from copy import deepcopy
from fractions import Fraction
from math import pi, sqrt
from pathlib import Path
import re
import types

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
SMALL = (5, 8)


def module(name):
    path = HERE / (name + ".py")
    obj = types.ModuleType(name)
    obj.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), obj.__dict__)
    return obj


producer = module("crossing_read_area_law")
verifier = module("verify_crossing_read_area_law")


@pytest.fixture(scope="module")
def families():
    return producer.evidence_families()


@pytest.fixture(scope="module")
def fresh(families):
    return {q: producer.level_report(q, families[q]) for q in SMALL}


@pytest.fixture(scope="module")
def graphs():
    out = {}
    for q in SMALL:
        values, A, B, sites = producer.population(q)
        indptr, indices = producer.read_relation(q, values, A, B, sites)
        out[q] = (values, A, B, sites, indptr, indices)
    return out


@pytest.fixture(scope="module")
def stored():
    return verifier.load()


def midplane(values, sites):
    return producer.site_sides(values, sites, 1, Fraction(1, 2))


# --------------------------------------------------------------------------
# Receipt against fresh rebuilds
# --------------------------------------------------------------------------


def test_stored_receipt_matches_fresh_small_levels(stored, fresh):
    rows = {r["q"]: r for r in stored["levels"]}
    for q in SMALL:
        assert producer.canonical(rows[q]) == producer.canonical(fresh[q])


def test_stored_receipt_verifies_on_small_levels(stored):
    result = verifier.verify(stored, SMALL)
    assert [r["q"] for r in result["levels_rebuilt"]] == list(SMALL)
    assert result["levels_in_receipt"][:2] == list(SMALL)


def test_independent_rebuild_equals_producer(fresh, families):
    for q in SMALL:
        assert verifier.canonical(verifier.rebuild_level(q, families[q])) == producer.canonical(fresh[q])


def test_relation_equals_evidence_package(fresh):
    for q in SMALL:
        assert fresh[q]["equals_evidence_receipt_digest"] is True
        assert fresh[q]["equals_evidence_receipt_edges"] is True


def test_scaled_relation_reproduces_evidence_builder(graphs):
    for q in SMALL:
        values, A, B, sites, indptr, indices = graphs[q]
        indptr2, indices2 = producer.scaled_read_relation(q, values, A, B, sites, Fraction(1))
        assert np.array_equal(indptr, indptr2) and np.array_equal(indices, indices2)


# --------------------------------------------------------------------------
# Mutation guards on the count
# --------------------------------------------------------------------------


def test_halving_the_radius_changes_the_count(graphs):
    values, A, B, sites, indptr, indices = graphs[8]
    side = midplane(values, sites)
    full = producer.crossing_count(indptr, indices, side)
    indptr_h, indices_h = producer.scaled_read_relation(8, values, A, B, sites, Fraction(1, 4))
    half = producer.crossing_count(indptr_h, indices_h, side)
    assert 0 < half < full
    assert producer.pair_count(indptr_h) < producer.pair_count(indptr)


def test_moving_the_plane_changes_the_count(graphs):
    values, A, B, sites, indptr, indices = graphs[8]
    counts = [producer.crossing_count(indptr, indices, producer.site_sides(values, sites, 1, pos))
              for pos in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 4))]
    assert len(set(counts)) == 3
    values, A, B, sites, indptr, indices = graphs[5]
    assert producer.crossing_count(indptr, indices, producer.site_sides(values, sites, 1, Fraction(1, 5))) != \
        producer.crossing_count(indptr, indices, midplane(values, sites))


def test_flipping_the_plane(graphs):
    values, A, B, sites, indptr, indices = graphs[8]
    quarter = producer.crossing_count(indptr, indices, producer.site_sides(values, sites, 1, Fraction(1, 4)))
    mirrored = producer.crossing_count(indptr, indices, producer.site_sides(values, sites, 1, Fraction(3, 4)))
    assert quarter != mirrored
    side = midplane(values, sites)
    assert producer.crossing_count(indptr, indices, -side) == producer.crossing_count(indptr, indices, side)


def test_self_reads_never_cross_but_are_pairs(graphs):
    for q in SMALL:
        values, A, B, sites, indptr, indices = graphs[q]
        side = midplane(values, sites)
        base = producer.crossing_count(indptr, indices, side)
        assert producer.crossing_count(indptr, indices, side, include_self=True) == base
        assert producer.pair_count(indptr, include_self=True) == producer.pair_count(indptr) + q ** 3
        rows = np.repeat(np.arange(q ** 3), np.diff(indptr))
        keep = indices != rows
        pruned_ptr = np.zeros(q ** 3 + 1, dtype=np.int64)
        pruned_ptr[1:] = np.cumsum(np.bincount(rows[keep], minlength=q ** 3))
        assert producer.crossing_count(pruned_ptr, indices[keep], side) == base
        assert int(pruned_ptr[-1]) == 2 * producer.pair_count(indptr)


def test_ordered_count_is_twice_unordered(graphs):
    for q in SMALL:
        values, A, B, sites, indptr, indices = graphs[q]
        side = midplane(values, sites)
        assert producer.crossing_count(indptr, indices, side, ordered=True) == \
            2 * producer.crossing_count(indptr, indices, side)


def test_no_site_on_any_plane(graphs):
    for q in SMALL:
        values, A, B, sites, indptr, indices = graphs[q]
        for axis, position in producer.PLANES + ((1, Fraction(1, 3)),):
            assert not np.any(producer.site_sides(values, sites, axis, position) == 0)


def test_midplane_symmetry_across_axes(fresh):
    for q in SMALL:
        counts = {c["plane"]: c["count"] for c in fresh[q]["crossings"]}
        assert counts["x1 = L/2"] == counts["x2 = L/2"] == counts["x3 = L/2"]
        assert fresh[q]["midplane_counts_agree_across_axes"] is True


# --------------------------------------------------------------------------
# Analytic identities
# --------------------------------------------------------------------------


def gauss(f, lo, hi, nodes=24):
    x, w = np.polynomial.legendre.leggauss(nodes)
    u = 0.5 * (hi - lo) * x + 0.5 * (hi + lo)
    return 0.5 * (hi - lo) * float(np.dot(w, f(u)))


@pytest.mark.parametrize("a", [0.5, 1.0, 2.5])
def test_cap_volume_integral_identity(a):
    cap = lambda u: pi * (a - u) ** 2 * (2 * a + u) / 3
    assert abs(gauss(cap, 0.0, a) - pi * a ** 4 / 4) <= 1e-12 * a ** 4
    depth = np.linspace(0.0, a, 7)
    assert np.allclose(cap(depth), pi * (2 * a ** 3 - 3 * a ** 2 * depth + depth ** 3) / 3)
    for n in (0.3, 2.0):
        assert abs(n ** 2 * gauss(cap, 0.0, a) - pi / 4 * n ** 2 * a ** 4) <= 1e-12 * a ** 4


@pytest.mark.parametrize("q", [5, 8, 13, 21, 34, 55])
def test_golden_identity_and_dictionary(q):
    L = 2.0 / sqrt((1.0 + sqrt(5.0)) / 2.0 + 2.0)
    n, a2 = q ** 3 / L ** 3, L ** 2 / q
    assert abs(n ** 2 * a2 ** 2 * L ** 2 - q ** 4) <= 1e-9 * q ** 4
    assert abs(L ** 2 * pi / 4 * n ** 2 * a2 ** 2 - producer.infinite_plane_prediction(q)) <= 1e-9 * q ** 4
    ell = producer.planck_length_over_L(q) * L
    assert abs(pi / 4 * n ** 2 * a2 ** 2 - 1.0 / (4 * ell ** 2)) <= 1e-9 * n ** 2 * a2 ** 2
    assert abs(ell ** 2 - L ** 2 / (pi * q ** 4)) <= 1e-15


@pytest.mark.parametrize("q", [5, 8, 13, 21, 34, 55])
def test_box_corrected_matches_closed_form(q):
    box = producer.box_corrected_prediction(q, Fraction(1, 2))
    assert abs(box - producer.box_corrected_closed_form(q)) <= 1e-12 * box
    assert box < producer.infinite_plane_prediction(q)
    quarter = producer.box_corrected_prediction(q, Fraction(1, 4))
    if producer.cap_truncated(q, Fraction(1, 4)):
        assert quarter < box
    else:
        assert abs(quarter - box) <= 1e-12 * box
    assert abs(producer.box_corrected_prediction(q, Fraction(3, 4)) - quarter) <= 1e-12 * box


def test_box_correction_closer_than_infinite_plane(stored):
    for row in stored["levels"]:
        mid = row["crossings"][0]
        assert mid["plane"] == "x1 = L/2"
        assert abs(mid["ratio_to_box_corrected"] - 1) < abs(mid["ratio_to_infinite_plane"] - 1)
        assert row["mean_neighbours_excluding_self"] < row["continuum_ball_neighbours"]
        assert abs(row["mean_neighbours_excluding_self"] / row["box_corrected_ball_neighbours"] - 1) < \
            abs(row["mean_neighbours_excluding_self"] / row["continuum_ball_neighbours"] - 1)


# --------------------------------------------------------------------------
# Verifier discipline
# --------------------------------------------------------------------------


@pytest.mark.parametrize("mutation", ["count", "sites_on_plane", "digest_flag", "edges", "float_format",
                                      "ratio", "extra", "scope", "pin", "evidence", "quadrature",
                                      "plane", "level_order", "self_reads"])
def test_receipt_forgeries_rejected(stored, mutation):
    item = deepcopy(stored)
    row = item["levels"][0]
    if mutation == "count":
        row["crossings"][0]["count"] += 1
    elif mutation == "sites_on_plane":
        row["crossings"][1]["sites_on_plane"] = 1
    elif mutation == "digest_flag":
        row["equals_evidence_receipt_digest"] = False
    elif mutation == "edges":
        row["undirected_edges"] -= 1
    elif mutation == "float_format":
        row["infinite_plane_prediction"] = pi * 5 ** 4 / 4
    elif mutation == "ratio":
        row["crossings"][0]["ratio_to_box_corrected"] = 1.0
    elif mutation == "extra":
        row["undocumented"] = True
    elif mutation == "scope":
        item["scope"]["planck_length_physically_identified"] = True
    elif mutation == "pin":
        item["source_pins"]["Lean/Geometry/CrossingReadAreaLaw.lean"] = "0" * 64
    elif mutation == "evidence":
        item["evidence_receipt_sha256"] = "0" * 64
    elif mutation == "quadrature":
        item["quadrature"]["nodes_per_piece"] = 32
    elif mutation == "plane":
        item["planes"][0]["position_over_L"] = "1/3"
    elif mutation == "level_order":
        item["levels"] = item["levels"][::-1]
    elif mutation == "self_reads":
        row["self_reads"] = 0
    with pytest.raises(ValueError):
        verifier.verify(item, SMALL)


@pytest.mark.parametrize("raw", ['{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":-Infinity}'])
def test_strict_json_parser(tmp_path, raw):
    path = tmp_path / "forged.json"
    path.write_bytes(raw.encode("ascii"))
    with pytest.raises(ValueError):
        verifier.load(path)


def test_verifier_is_independent():
    source = (HERE / "verify_crossing_read_area_law.py").read_text()
    assert re.search(r"^\s*(import|from)\s+(crossing_read_area_law|build_causal_poset)\b", source, re.M) is None
    assert "sys.path" not in source and "exec(" not in source and "importlib" not in source
    assert not hasattr(verifier, "bcp") and not hasattr(verifier, "level_report")


def test_exact_sign_rules_agree():
    rng = np.random.default_rng(20260911)
    c = rng.integers(-10 ** 6, 10 ** 6, size=20000)
    b = rng.integers(-10 ** 6, 10 ** 6, size=20000)
    c[:5] = [0, 0, 3, -3, 9]
    b[:5] = [0, 2, 0, 0, -4]
    exact = verifier.sqrt5_sign(c, b)
    floats = np.sign(c + b * sqrt(5.0))
    decided = np.abs(c + b * sqrt(5.0)) > 1e-3
    assert np.array_equal(exact[decided], floats[decided])
    assert np.array_equal(exact[:5], [0, 1, 1, -1, 1])
    a = rng.integers(-10 ** 5, 10 ** 5, size=5000)
    bb = rng.integers(-10 ** 5, 10 ** 5, size=5000)
    assert np.array_equal(verifier.phi_sign(a, bb), producer.bcp.phi_sign_array(a, bb))
