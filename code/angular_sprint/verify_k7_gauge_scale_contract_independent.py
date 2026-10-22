#!/usr/bin/env python3
"""Independent branch enumeration and symbolic RG transport verification.

Never imports the scale-contract producer, its axial covariance audit, or its
clock reduction. The frozen semantic digest covers source interpretations;
the numerical fields excluded from that digest are reconstructed below.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import re
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "code/angular_sprint/runtime/k7_gauge_scale_contract.json"
PINS = {
    "model/primitive_observer/transfer.py": "1d87c7a73f6949b3b8507ecd2953454b13178c50",
    "model/common_source/slow_geometry.py": "63685aaa951622dce316afcbd5037d3f305f5d58",
    "model/observer/native.py": "2fd5415ed5ca3f19ba44e1514068609222121a33",
    "evidence/native_geometric_source_20260925/spec.json": "1f1424332973665c7ae1333f3ea1f29717382602",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py": "0cbce7722440a80ac7b8c1914862135fddcb4022",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py": "fee7f6228031a6d9f4e9e05976e4d49e9d39ade1",
    "Lean/Screen/RGRepresentationFrontier.lean": "c35eb3bc1f8e5f8bd4633bdea07911049bb985d0",
    "model/native_circulation/construction.py": "ddd84d7390cef65d389a390c2f1318b0e2e34e3b",
    "model/native_clock/clock.py": "569fa77e414b11f87c891774fd19f4642529b798",
    "code/angular_sprint/external/k7_gauge_scale/ERRATUM_v3.5.md": "0770f07ebaf48df4d7ce026e89da4dca3dbd0563",
}
SEMANTIC_SHA256 = "008ec0d9786ae814ce522bbe6db04ad0f4d3005c9d9af32502f4a654a43704ba"
DYNAMIC = {
    "native_signature": ("branches", "covariance_diagonal", "counter_diffusion",
                         "one_attempt_law", "unit_extensions"),
    "conditional_plane_transport": ("k", "b", "cofactor", "beta_annihilator",
        "k_annihilator",
        "mixed_component_transport_delta_D", "wrong_branch_transport_delta_D",
        "threshold_basis_delta_D"),
}


def semantic_digest(doc):
    frozen = deepcopy(doc)
    del frozen["source_pins"]
    for section, keys in DYNAMIC.items():
        for key in keys:
            del frozen[section][key]
    return hashlib.sha256(json.dumps(frozen, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def pairs(rows):
    doc = {}
    for key, value in rows:
        if key in doc:
            raise ValueError("duplicate JSON key")
        doc[key] = value
    return doc


def load(path):
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def reconstruct_native():
    from model.observer import native
    import ast
    from itertools import combinations
    closure = ROOT / "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py"
    tree = ast.parse(closure.read_text())
    node = next(statement.value for statement in tree.body
                if isinstance(statement, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == "ORIENTED_BASE_FACES"
                        for target in statement.targets))
    faces = ast.literal_eval(node)
    seams = tuple(sorted({tuple(sorted((face[i],face[(i+1)%3])))
                          for face in faces for i in range(3)}))
    # Reconstruct the Q(sqrt(5)) chart in coordinates, independently of the
    # producer's sorted graph traversal and rational current bookkeeping.
    phi = sp.Rational(1,2)+sp.sqrt(5)/2
    minus_phi = -phi
    chart = ((-1,phi,0),(1,phi,0),(-1,minus_phi,0),(1,minus_phi,0),
             (0,-1,phi),(0,1,phi),(0,-1,minus_phi),(0,1,minus_phi),
             (phi,0,-1),(phi,0,1),(minus_phi,0,-1),(minus_phi,0,1))
    vertices = [sp.Matrix(row) for row in chart]
    radius2 = sp.expand(vertices[0].dot(vertices[0]))
    weights = [(vertices[a].cross(vertices[b])/radius2).applyfunc(sp.simplify)
               for a,b in seams]
    # Rebuild the 66 configurations rather than consuming the producer kernel.
    states = [tuple(int(i == a)+int(i == b) for i in range(12))
              for a,b in combinations(range(12),2)]
    covariance = sp.zeros(3)
    counts = [0,0]
    branches = 0
    for state in states:
        drift = sp.zeros(3,1)
        for seam, weight in zip(seams, weights, strict=True):
            for coin in (0,1):
                out, record = native.native(state,seam,coin)
                if out not in states or record[0] != coin:
                    raise ValueError("native binary carrier changed")
                current = record[1]
                increment = weight*current
                drift += increment
                covariance += increment*increment.T
                counts[int(bool(current))] += 1
                branches += 1
        if drift.applyfunc(sp.simplify) != sp.zeros(3,1):
            raise ValueError("native counter has drift")
    covariance = (covariance/branches).applyfunc(sp.simplify)
    variance = covariance[0,0]
    if covariance != variance*sp.eye(3):
        raise ValueError("native counter is not isotropic")
    return {"branches": branches, "covariance_diagonal": str(variance),
            "counter_diffusion": str(variance/2),
            "one_attempt_law": [str(Q(n,branches)) for n in counts],
            "unit_extensions": [
                {"length_per_counter_unit": str(length), "time_per_attempt": str(time),
                 "diffusion_in_assigned_units": str(variance*length**2/(2*time))}
                for length,time in ((1,1),(2,1),(1,2))]}


def reconstruct_plane():
    text = (ROOT/"Lean/Screen/RGRepresentationFrontier.lean").read_text()
    def read(name):
        body = re.search(r"^def "+name+r"\b[^\n]*:=\n(.*?)(?=\n\n)",
                         text, re.M | re.S).group(1)
        # The custody gate fixes the trusted grammar; symbolic rational parsing
        # is a different route from the producer's Fraction AST evaluator.
        return sp.sympify(" ".join(body.split()).replace("^", "**"),
                          locals={"nG": sp.Integer(3), "nH": sp.Integer(1)},
                          rational=True)
    k = sp.Matrix([read(name) for name in
                   ("u1WeylIndexPerCopy","su2WeylIndexPerCopy","su3WeylIndexPerCopy")])
    b = sp.Matrix([read(name) for name in ("bY","b2","b3")])
    x = sp.Matrix(sp.symbols("x0:3")); t = sp.Symbol("t")
    determinant = lambda y: sp.Matrix.hstack(y,k,b).det()
    if sp.expand(determinant(x-t*b)-determinant(x)) != 0:
        raise ValueError("conditional RG invariance fails")
    cofactor = k.cross(b)
    mixed = sp.Matrix([k[0]-b[0],k[1],k[2]])
    wrong = sp.Matrix([11,1,-3])
    if determinant(mixed) == 0 or determinant(wrong) == 0:
        raise ValueError("transport controls lost discrimination")
    return {"k": list(map(str,k)), "b": list(map(str,b)),
            "cofactor": list(map(str,cofactor)), "beta_annihilator": str(cofactor.dot(b)),
            "k_annihilator": str(cofactor.dot(k)),
            "mixed_component_transport_delta_D": str(determinant(mixed)),
            "wrong_branch_transport_delta_D": str(-determinant(wrong)),
            "threshold_basis_delta_D": list(map(str,cofactor))}


def verify(doc):
    pins = {}
    for path, expected in PINS.items():
        raw = (ROOT/path).read_bytes()
        actual = hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
        if actual != expected:
            raise ValueError(f"source pin changed: {path}")
        pins[path] = {"git_blob": expected, "sha256": hashlib.sha256(raw).hexdigest()}
    if doc["source_pins"] != pins:
        raise ValueError("source custody mismatch")
    if semantic_digest(doc) != SEMANTIC_SHA256:
        raise ValueError("reviewed source semantics changed")
    # JSON encoding distinguishes false from 0 and prevents lossy coercions.
    for section, expected in (("native_signature",reconstruct_native()),
                              ("conditional_plane_transport",reconstruct_plane())):
        actual = {key: doc[section][key] for key in DYNAMIC[section]}
        if json.dumps(actual,sort_keys=True) != json.dumps(expected,sort_keys=True):
            raise ValueError(f"independent reconstruction differs: {section}")
    return True


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receipt", type=Path, default=OUT)
    args = parser.parse_args()
    verify(load(args.receipt))
    print("PASS independent native enumeration and symbolic conditional RG transport")


if __name__ == "__main__":
    main()
