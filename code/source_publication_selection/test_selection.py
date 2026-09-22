"""Semantic forgery rejection and independent native/geometry controls."""
from copy import deepcopy
from decimal import Decimal, localcontext
from itertools import product
from functools import lru_cache
import os
import shutil
import subprocess
import sys

import pytest

from . import build, codec, verify


@pytest.fixture(scope="module")
def packet():
    return codec.load(codec.HERE/"controls.json")


def test_committed_receipt(packet):
    codec.equal(verify.verify(packet), codec.load(codec.HERE/"receipt.json"), "receipt")


def test_producer_independence():
    script = '''
import importlib.abc
import sys
class RejectProducer(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in {"source_publication_selection.build", "source_publication_selection.exterior", "source_publication_selection.locality", "source_checkpoint_selection.build",
                        "source_checkpoint_selection.pipeline", "source_temporal_acceptance.build"}:
            raise AssertionError("producer import: " + fullname)
sys.meta_path.insert(0, RejectProducer())
from source_publication_selection import verify
verify.main()
'''
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("path,value", (
    (("schema",), True), (("schema",), 1.0),
    (("locality", "intervals", 0, "events"), True),
    (("locality", "intervals", 11, "ordering_fraction"), "1/10"),
    (("locality", "intervals", 11, "incomparable_pairs"), 0),
    (("locality", "stencils", 0, "target_reachable"), True),
    (("locality", "stencils", 0, "minimum_steps"), 3),
    (("locality", "stencils", 1, "minimum_steps"), 0),
    (("locality", "stencils", 1, "stencil_size_including_self"), 7),
    (("locality", "stencils", 5, "euclidean_strict_timelike"), False),
    (("locality", "stencils"), []),
    (("locality", "intervals"), []),
    (("scope", "M1_derived"), True),
    (("scope", "meaning_population_radius_selected"), True),
    (("scope", "all_axioms_instantiated"), True),
    (("native", 0, "potential"), "1"),
    (("native", 0, "selected_first"), ["0", "1/2", "1/2"]),
    (("native", 0, "equal_prior_menu_posterior"), ["1", "0"]),
    (("native", 0, "aggregate_potential"), "1/2"),
    (("native", 0, "weights"), [True, 1, 1]),
    (("native", 0, "levels", 0, "words"), 0),
    (("native", 0, "levels", 0, "published"), 1),
    (("native", 0, "levels", 4, "mass"), "1/10"),
    (("native", 0, "levels", 8, "published"), 674),
    (("native", 0, "levels", 8, "census_sha256"), "0"*64),
    (("native", 0, "levels", 5, "cylinder_total"), "1/2"),
    (("native", 0, "levels", 5, "selected_published"), "1"),
    (("native", 0, "levels", 5, "conditional_first"), ["0", "2/3", "1/3"]),
    (("geometry", 0, "ordered_reads_including_self"), [7919, 2784, 1153]),
    (("geometry", 0, "ordered_reads_including_self"), [7919, 1153, 2785]),
    (("geometry", 0, "sites"), 124),
    (("geometry", 0, "q"), 5.0),
    (("geometry", 0, "population"), "grid"),
    (("geometry", 0, "one_dimensional_distance_types"), 1),
    (("geometry", 4, "all_ordered_pairs"), 9261),
    (("paths", 0, "sources"), [0, 1]),
    (("paths", 0, "levels", 8, "minimum_fixed"), "0"),
    (("paths", 0, "levels", 8, "census_sha256"), "0"*64),
    (("paths", 0, "levels", 8, "words"), 1),
    (("paths", 4, "levels", 6, "rank_lost"), 0),
))
def test_resealed_semantic_forgery(packet, path, value):
    forged = deepcopy(packet)
    target = forged
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    forged["sha256"] = codec.digest({k: v for k, v in forged.items() if k != "sha256"})
    with pytest.raises(ValueError):
        verify.verify(forged)


@pytest.mark.parametrize("field", ("native", "paths", "geometry", "locality", "scope", "source_pins", "sha256"))
def test_omitted_field(packet, field):
    forged = deepcopy(packet)
    del forged[field]
    with pytest.raises(ValueError):
        verify.verify(forged)


@pytest.mark.parametrize("field", ("native", "paths", "geometry"))
def test_omitted_case(packet, field):
    forged = deepcopy(packet)
    forged[field].pop()
    forged["sha256"] = codec.digest({k: v for k, v in forged.items() if k != "sha256"})
    with pytest.raises(ValueError):
        verify.verify(forged)


def test_stale_pins_and_unsealed_edit(packet):
    forged = deepcopy(packet)
    first = next(iter(forged["source_pins"]))
    forged["source_pins"][first] = "0"*64
    forged["sha256"] = codec.digest({k: v for k, v in forged.items() if k != "sha256"})
    with pytest.raises(ValueError, match="source pins"):
        verify.verify(forged)
    forged = deepcopy(packet)
    forged["native"][0]["potential"] = "1"
    with pytest.raises(ValueError, match="digest"):
        verify.verify(forged)


@pytest.mark.parametrize("raw", ('{"x":0,"x":1}', '{"x":0.5}', '{"x":NaN}', '{"x":Infinity}'))
def test_json_rejection(tmp_path, raw):
    path = tmp_path/"bad.json"
    path.write_text(raw, encoding="ascii")
    with pytest.raises(ValueError):
        codec.load(path)


@pytest.mark.parametrize("word", ([True], [-1], [3], [1.0], "120", None))
def test_native_word_validation(word):
    with pytest.raises(ValueError):
        verify.experiment(word)


@pytest.mark.parametrize("weights,horizon", (
    ((0, 1, 1), 3), ((-1, 1, 1), 3), ((True, 1, 1), 3), ((1.0, 1, 1), 3),
    ((1, 1), 3), ((1, 1, 1), -1), ((1, 1, 1), True), ((1, 1, 1), 13)))
def test_native_specification_validation(weights, horizon):
    with pytest.raises(ValueError):
        build.native(weights, horizon)


@pytest.mark.parametrize("q,pop", ((True, "golden"), (5.0, "grid"), (1, "golden"), (5, "torus")))
def test_population_validation(q, pop):
    with pytest.raises(ValueError):
        build.counts(q, pop)


def scalar_leaves(value, path=()):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from scalar_leaves(child, (*path, key))
    elif isinstance(value, list):
        for key, child in enumerate(value):
            yield from scalar_leaves(child, (*path, key))
    else:
        yield path, value


def test_every_scalar_field_resealed(packet, monkeypatch):
    # Cache only the independently recomputed expected experiments. Every
    # forged packet still traverses the actual schema, custody and semantic
    # comparisons. No expected values come from the producer or the packet.
    from . import check_exterior, check_locality
    for module, name in ((codec, "pins"), (verify, "native"), (verify, "geometry"),
                         (check_exterior, "case"), (check_locality, "interval"),
                         (check_locality, "stencil")):
        monkeypatch.setattr(module, name, lru_cache(maxsize=None)(getattr(module, name)))
    verify.verify(packet)
    count = 0
    for path, value in scalar_leaves(packet):
        if path == ("sha256",):
            continue  # Resealing would overwrite this edit; tested separately.
        forged = deepcopy(packet)
        parent = forged
        for key in path[:-1]:
            parent = parent[key]
        if type(value) is bool:
            bad = not value
        elif type(value) is int:
            bad = value+1
        elif type(value) is str:
            bad = value+"!"
        elif value is None:
            bad = "not-null"
        else:
            raise AssertionError((path, value))
        parent[path[-1]] = bad
        forged["sha256"] = codec.digest({k: v for k, v in forged.items() if k != "sha256"})
        with pytest.raises(ValueError):
            verify.verify(forged)
        count += 1
    assert count > 1000
    verify.verify(packet)  # Rejections cannot come from poisoning the verifier.
    print(f"individually rejected {count} resealed scalar-field corruptions")


def test_real_cli_rejects_artifacts_and_source_edits(tmp_path, packet):
    clone = tmp_path/"repository"
    paths = set(codec.pins()) | {"code/source_publication_selection/controls.json",
                                "code/source_publication_selection/receipt.json"}
    for name in paths:
        destination = clone/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(codec.ROOT/name, destination)
    package = clone/"code/source_publication_selection"
    env = dict(os.environ, PYTHONPATH=str(clone/"code"))

    def run():
        return subprocess.run([sys.executable, "-m", "source_publication_selection.verify"],
                              cwd=clone, env=env, capture_output=True, text=True, timeout=90)

    good = run()
    assert good.returncode == 0, good.stderr
    assert '"M1_derived":false' in good.stdout
    changes = []
    for raw in (b'{}', b'{"x":0,"x":1}', b'{"x":0.5}', b'{"x":NaN}'):
        changes.append((package/"controls.json", raw, None))
    forged = deepcopy(packet)
    forged["native"][0]["potential"] = "1"
    forged["sha256"] = codec.digest({k: v for k, v in forged.items() if k != "sha256"})
    changes.append((package/"controls.json", codec.canonical(forged), "native semantics"))
    receipt = codec.load(package/"receipt.json")
    receipt["M1_derived"] = True
    changes.append((package/"receipt.json", codec.canonical(receipt), "receipt"))
    source = clone/"Lean/Geometry/SourcePublicationPathInvariant.lean"
    changes.append((source, source.read_bytes()+b'\n-- deliberate source tampering\n', "source pins"))
    for path, bad, diagnostic in changes:
        original = path.read_bytes()
        try:
            path.write_bytes(bad)
            result = run()
            assert result.returncode != 0, (path, result.stdout)
            assert '"verified_controls_sha256"' not in result.stdout
            if diagnostic:
                assert diagnostic in result.stderr
        finally:
            path.write_bytes(original)
    restored = run()
    assert restored.returncode == 0, restored.stderr


@pytest.mark.parametrize("q,population", tuple(product((2, 3, 5), ("golden", "grid"))))
def test_direct_point_pair_census(q, population):
    # A third algorithm, enumerating actual 3D sites and every ordered pair.
    # High-precision decimal comparisons are only a test cross-check; the
    # committed verifier certifies boundaries by rational enclosures.
    with localcontext() as ctx:
        ctx.prec = 80
        D = Decimal
        phi = (1+D(5).sqrt())/2
        axis = [(phi*b)%1 if population == "golden" else D(b)/q for b in range(q)]
        sites = list(product(axis, repeat=3))
        thresholds = [1/D(q).sqrt(), 1/D(q), 1/(D(q)*D(q).sqrt())]
        counts = [0, 0, 0]
        for x, y in product(sites, repeat=2):
            squared = sum((a-b)**2 for a, b in zip(x, y))
            for i, cutoff in enumerate(thresholds):
                # Rational grid equality is handled exactly, including the
                # possible last decimal place in a repeating coordinate.
                if population == "grid" and i == 1:
                    numerator = sum((int((a*q).to_integral_value())-int((b*q).to_integral_value()))**2
                                    for a, b in zip(x, y))
                    inside = numerator <= q
                else:
                    assert abs(squared-cutoff) > D("1e-60")
                    inside = squared <= cutoff
                counts[i] += inside
    assert counts == build.counts(q, population)["ordered_reads_including_self"]
    assert counts == verify.geometry(q, population)["ordered_reads_including_self"]
