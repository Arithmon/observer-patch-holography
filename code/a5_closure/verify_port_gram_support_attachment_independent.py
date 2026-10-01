"""Independent reconstruction of seed attachments and support-chain data.

This module deliberately does not import the certificate producer and never
reads its candidate list. Its permutation search is driven by partial triangle
incidence checks; chain pushforwards are recomputed from the source faces.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
FACES_PATH = ROOT / "Lean/ObserverPatchHolography/CoreAxioms.lean"
GEOMETRY_PATH = ROOT / "code/source_selection_model/geometry.py"
RECEIPT_PATH = ROOT / "code/m1_source_realization/receipt.json"
DEGREE_VERIFIER = ROOT / "code/a5_closure/verify_galois_port_frame_independent.py"


def need(value, label):
    if not value:
        raise ValueError(label)


def cyclic(face):
    a, b, c = face
    return min((a, b, c), (b, c, a), (c, a, b))


def chain(faces):
    out = Counter()
    for face in faces:
        fwd = (face, face[1:] + face[:1], face[2:] + face[:2])
        rev = ((face[0], face[2], face[1]),
               (face[2], face[1], face[0]),
               (face[1], face[0], face[2]))
        canonical = min(fwd + rev)
        out[canonical] += 1 if canonical in fwd else -1
    return Counter({key: value for key, value in out.items() if value})


def read_local():
    text = FACES_PATH.read_text()
    block = text.split("def orientedFaces :", 1)[1].split("/--", 1)[0]
    return [tuple(map(int, x)) for x in re.findall(r"\((\d+), (\d+), (\d+)\)", block)]


def read_geometry_seed():
    text = GEOMETRY_PATH.read_text()
    block = text.split("def seed():", 1)[1].split("\n\n\ndef subdivide", 1)[0]
    return [tuple(map(int, x)) for x in re.findall(r"\((\d+), (\d+), (\d+)\)", block)]


def automorphisms(faces):
    # Each mapped oriented face must still be a source face whenever its
    # three vertices have been assigned. This is independent of the producer's
    # distance-table refinement strategy.
    face_set = {cyclic(f) for f in faces}
    incident = [set() for _ in range(12)]
    adjacent = [set() for _ in range(12)]
    for f in faces:
        for v in f:
            incident[v].add(f)
            adjacent[v].update(set(f) - {v})
    found = []

    def extend(mapping, used):
        if len(mapping) == 12:
            p = tuple(mapping[i] for i in range(12))
            if {cyclic(tuple(p[v] for v in f)) for f in faces} == face_set:
                found.append(p)
            return
        # Choose the vertex constrained by the most already mapped neighbors.
        v = max((x for x in range(12) if x not in mapping),
                key=lambda x: (sum(n in mapping for n in adjacent[x]), len(incident[x]), -x))
        for w in range(12):
            if w in used or len(adjacent[v]) != len(adjacent[w]):
                continue
            if any((u in adjacent[v]) != (image in adjacent[w])
                   for u, image in mapping.items()):
                continue
            trial = {**mapping, v: w}
            valid = True
            for f in incident[v]:
                if all(x in trial for x in f):
                    if cyclic(tuple(trial[x] for x in f)) not in face_set:
                        valid = False
                        break
            if valid:
                extend(trial, used | {w})

    extend({}, set())
    return sorted(set(found))


def verify_support(tower):
    need(len(tower) == 4, "tower row count")
    levels = []
    for level, row in enumerate(tower):
        s = row["support"]
        faces, section = [tuple(f) for f in s["faces"]], row["section"]
        need(s["level"] == level and len(faces) == 20 * 4**level, "support census")
        need(len(section) == s["vertices"] and len(set(section)) == len(section), "section")
        # The section cycle is pushed by b to the entire oriented support cycle.
        section_chain = [tuple(section[v] for v in f) for f in faces]
        vertex_projection = [tuple(row["carriers"][c][0] for c in f) for f in section_chain]
        need(chain(vertex_projection) == chain(faces), "b_* z != [S]")
        if level:
            prev = tower[level - 1]
            support_parent = s["coarsen"]
            nerve_parent = row["coarsen"]
            need(len(support_parent) == s["vertices"], "support parent map")
            need(all(nerve_parent[section[v]] == prev["section"][support_parent[v]]
                     for v in range(s["vertices"])), "section naturality")
            pushed = [tuple(nerve_parent[c] for c in f) for f in section_chain]
            pushed = [f for f in pushed if len(set(f)) == 3]
            prev_chain = [tuple(prev["section"][v] for v in f)
                          for f in prev["support"]["faces"]]
            need(chain(pushed) == chain(prev_chain), "section cycle refinement square")
        levels.append({"level": level, "vertices": s["vertices"],
                       "faces": len(faces), "carriers": len(row["carriers"]),
                       "bridge_degree": 1})
    return levels


def verify():
    local = read_local()
    seed = read_geometry_seed()
    packet = json.loads(RECEIPT_PATH.read_text())
    tower = packet["evidence"]["topology"]
    support = [tuple(f) for f in tower[0]["support"]["faces"]]
    need(local == seed == support, "literal ordered seed equality")
    need(len(local) == 20 and len({frozenset(f) for f in local}) == 20, "local faces")
    group = automorphisms(local)
    need(len(group) == 60, "Aut+ census")
    section = tower[0]["section"]
    z0 = [tuple(section[v] for v in f) for f in support]
    target = chain(z0)
    for perm in group:
        h = [section[perm[v]] for v in range(12)]
        image = [tuple(h[v] for v in f) for f in local]
        need(chain(image) == target, "exact integer-chain image")
    # Composition is a free transitive action on the 60 bijective seed maps.
    sample = group[0]
    images = {tuple(g[sample[v]] for v in range(12)) for g in group}
    need(len(images) == 60, "free transitive group action")
    support_levels = verify_support(tower)

    # Degree data comes from the separately maintained independent Galois
    # verifier, pinned here by source bytes. It does not supply attachments.
    verifier_sha = hashlib.sha256(DEGREE_VERIFIER.read_bytes()).hexdigest()
    expected_sha = "0b805033eb3e7c27c812a2d0919504796f73980d5d6d75f23dd68925045fbbea"
    need(verifier_sha == expected_sha, "pinned degree verifier changed")
    run = subprocess.run([sys.executable, "-B", str(DEGREE_VERIFIER)], cwd=ROOT,
                         check=True, capture_output=True, text=True)
    degree_receipt = json.loads(run.stdout)
    degrees = {name: degree_receipt["verified"][name]["degree"] for name in ("plus", "minus")}
    need(degrees == {"plus": 1, "minus": 7}, "independent local degrees")
    return {"verified": True, "seed_relation": "literal oriented equality",
            "candidate_count": 60, "exact_chain_pushforwards": 60,
            "stabilizer_size": 1, "orbit_count": 1, "support_levels": support_levels,
            "degree_receipt_sha256": verifier_sha, "local_degrees": degrees,
            "producer_candidate_list_consumed": False}


if __name__ == "__main__":
    print(json.dumps(verify(), sort_keys=True))
