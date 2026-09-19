from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import experiment as producer
import readouts
import verify


def test_exact_dyadic_encoding_and_overflow():
    values = np.array([0, 1, 2**63, 2**127 + 7, 2**192 - 1], dtype=object)
    assert np.array_equal(values, producer.unlimbs(producer.limbs(values)))
    with pytest.raises(AssertionError):
        producer.limbs(np.array([2**192], dtype=object))
    # The nearest-integer update used by the abandoned partial run is excluded.
    state = np.array([0, 1], dtype=object)
    numerator = state[0] + state[1]
    assert numerator == 1  # both outputs are 1 / 2, neither 0 nor 1


@pytest.mark.parametrize("level", [3, 4, 5])
def test_complete_matching_schedules_and_readback(level):
    g = dict(np.load(producer.HERE / f"geometry/geometry_l{level}.npz"))
    c = len(g["faces"])
    for wiring, expected in (("w12", 36*c-30), ("w3", 63*c//2), ("isolated", 30*c)):
        schedule = producer.phases(g, wiring)
        assert sum(len(e) for _, e in schedule) == expected
        for _, edges in schedule:
            assert len(np.unique(edges)) == 2 * len(edges)
    lap, projector, frame = producer.carrier(g["seams"])
    np.testing.assert_allclose(lap @ projector, (5-np.sqrt(5))*projector, atol=1e-13)
    np.testing.assert_allclose(frame@frame.T, 4*projector, atol=1e-13)
    assert np.linalg.matrix_rank(frame) == 3


def test_refinement_is_four_children_and_state_preserving():
    for level in (4, 5):
        g = dict(np.load(producer.HERE / f"geometry/geometry_l{level}.npz"))
        coarse = dict(np.load(producer.HERE / f"geometry/geometry_l{level-1}.npz"))
        parents, children, weights = g["parent"], g["children"], g["expectation_weights"]
        np.testing.assert_array_equal(parents[children], np.arange(len(children))[:, None].repeat(4, axis=1))
        np.testing.assert_allclose(weights[children].sum(axis=1), 1, atol=2e-15)
        np.testing.assert_allclose(g["areas"][children].sum(axis=1), coarse["areas"], atol=2e-15)
        np.testing.assert_allclose(weights, g["areas"] / coarse["areas"][parents], atol=1e-13)


@pytest.fixture
def small_trace(tmp_path, monkeypatch):
    """Small executor fixture; the full geodesic census is checked separately.

    Use two initial carriers and the same four-child copy convention. This
    exercises every serialized phase and law without copying the 60 MB archive
    for each deliberate corruption.
    """
    original = producer.HERE
    seams = np.load(original / "geometry/geometry_l3.npz")["seams"]
    (tmp_path / "geometry").mkdir()
    (tmp_path / "SPECIFICATION.md").write_bytes((original / "SPECIFICATION.md").read_bytes())
    (tmp_path / "geometry/geometry.json").write_text("{}\n")
    geometries = {}
    for level, count in ((3, 2), (4, 8), (5, 32)):
        glued = np.array([[a, p, a+1, p] for a in range(0, count, 2) for p in range(12)])
        g = {"faces": np.zeros((count, 3), dtype=int), "seams": seams, "w12": glued}
        if level > 3:
            g["parent"] = np.arange(count)//4
        np.savez_compressed(tmp_path / f"geometry/geometry_l{level}.npz", **g)
        geometries[level] = g
    monkeypatch.setattr(producer, "HERE", tmp_path)
    monkeypatch.setattr(verify, "HERE", tmp_path)
    trace = tmp_path / "trace"
    producer.build_trace(trace)
    return trace, geometries


def reseal(trace, manifest):
    """Adversary updates all hashes: semantic rejection must still hold."""
    import hashlib
    chain = hashlib.sha256(producer.canonical(manifest["binding"])).hexdigest()
    for item in manifest["chunks"]:
        item["previous"] = chain
        item["sha256"] = producer.sha(trace / item["file"])
        item.pop("chain", None)
        chain = hashlib.sha256(producer.canonical(item)).hexdigest()
        item["chain"] = chain
    manifest["final_chain"] = chain
    producer.save_json(trace / "trace.json", manifest)


def test_small_trace_independent_replay(small_trace):
    directory, g = small_trace
    replay = verify.verify_trace(directory, g)
    assert replay[-1]["exact_total_and_descent_checks"]
    # The two writes of an action read old versions, not one another.
    first = replay[0]["chunks"][1]["first_event"]
    assert np.array_equal(replay[1][first], replay[1][first+1])
    assert first not in replay[1][first+1]


@pytest.mark.parametrize("mutation, message", [
    ("mean", "noncanonical seam mean"),
    ("stale", "stale or unauthenticated read"),
    ("omit", "missing or changed seam attempts"),
    ("copy", "refinement componentwise copy"),
    ("join_writer", "refinement current parent writers"),
    ("position", "writer rank-three readback"),
    ("law", "trace contract binding"),
    ("missing_phase", "complete scheduled phase count"),
])
def test_rehashed_semantic_corruptions_fail(small_trace, mutation, message):
    import json
    trace, g = small_trace
    manifest = json.loads((trace / "trace.json").read_text())
    item = next(x for x in manifest["chunks"] if x["kind"] == ("copy" if mutation in ("copy", "join_writer") else "intra"))
    data = dict(np.load(trace / item["file"]))
    if mutation in ("mean", "copy"):
        data["values"][0, 0] ^= np.uint64(1)
    elif mutation in ("stale", "join_writer"):
        data["parents"].flat[0] = -1
    elif mutation == "omit":
        data["endpoints"] = data["endpoints"][1:]
    elif mutation == "position":
        data["positions"][0, 0] += .01
    elif mutation == "law":
        manifest["binding"]["law"] = "integer_nearest_agreement"
    elif mutation == "missing_phase":
        del manifest["chunks"][1]
    np.savez_compressed(trace / item["file"], **data)
    reseal(trace, manifest)
    with pytest.raises(ValueError, match=message):
        verify.verify_trace(trace, g)


def test_broken_hash_fails(small_trace):
    trace, g = small_trace
    path = trace / "phase_001.npz"
    with path.open("ab") as stream:
        stream.write(b"tampered")
    with pytest.raises(ValueError, match="trace chunk hash"):
        verify.verify_trace(trace, g)


def test_interval_algorithms_against_exhaustive_oracle():
    rng = np.random.default_rng(776)
    for size in range(2, 25):
        parents = [np.asarray([], dtype=int)]
        for v in range(1, size):
            parents.append(np.unique([v-1, int(rng.integers(v))]))
        for bottom in (0, size//2):
            members, row = readouts.interval(parents, bottom, size-1)
            verified, count, height = verify.count_interval(parents, bottom, size-1)
            reach = np.eye(size, dtype=bool)
            for v in range(size):
                for p in parents[v]:
                    reach[:, v] |= reach[:, p]
            expected = np.flatnonzero(reach[bottom] & reach[:, size-1])
            oracle = sum(bool(reach[a, b]) for a in expected for b in expected if a != b)
            assert np.array_equal(members, expected) and np.array_equal(verified, expected)
            assert count == row["strict_pairs"] == oracle
            verify.check_order_row(row, members, count, height)
            row["strict_pairs"] += 1
            with pytest.raises(ValueError, match="interval counts"):
                verify.check_order_row(row, members, count, height)


def test_dimension_references_and_undefined_cases():
    for fraction, dimension in ((1, 1), (.5, 2), (8/35, 3), (.1, 4)):
        assert readouts.dimension(fraction) == pytest.approx(dimension, abs=1e-12)
    assert readouts.dimension(0) is None
    assert readouts.dimension(-1) is None
    assert readouts.dimension(1.1) is None
