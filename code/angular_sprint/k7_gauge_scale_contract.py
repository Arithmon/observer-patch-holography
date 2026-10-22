#!/usr/bin/env python3
"""Exact native-scale typing and conditional RG-plane transport certificate.

The unit extensions concern the declared native counter signature only.
The plane identity assumes the frozen matter branch and common subtraction
prescription; it selects neither a physical kinetic form nor an RG map.
"""
from __future__ import annotations

import argparse
import ast
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
BASELINES = {
    "oph": "879ed37ad2c1933ecb6dd8f54271362ee3f4d41c",
    "oph_integration": "7b0997b4a8a425679d9e2946f82c86cd92289059",
    "k7": "26d18fdd2fce61cc0231ec3485f8f5fe7efd5fc6",
    "k7_lean": "70dc589bd7c3f77db54a2acb6e67986a7ef69468",
}


def blob(raw):
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def custody(root=ROOT):
    rows = {}
    for path, expected in PINS.items():
        raw = (root / path).read_bytes()
        if blob(raw) != expected:
            raise ValueError(f"source pin changed: {path}")
        rows[path] = {"git_blob": expected, "sha256": hashlib.sha256(raw).hexdigest()}
    return rows


def rational_expression(text, **names):
    """Read only rational arithmetic from the pinned Lean definitions."""
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return Q(node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return Q(names[node.id])
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return a+b
            if isinstance(node.op, ast.Sub): return a-b
            if isinstance(node.op, ast.Mult): return a*b
            if isinstance(node.op, ast.Div): return a/b
            if isinstance(node.op, ast.Pow) and b.denominator == 1: return a**int(b)
        raise ValueError("unsupported source arithmetic")
    return visit(ast.parse(" ".join(text.split()).replace("^", "**"), mode="eval").body)


def oph_columns():
    text = (ROOT / "Lean/Screen/RGRepresentationFrontier.lean").read_text()
    def definition(name):
        match = re.search(r"^def " + name + r"\b[^\n]*:=\n(.*?)(?=\n\n)", text, re.M | re.S)
        if match is None: raise ValueError(f"missing definition {name}")
        return rational_expression(match.group(1), nG=3, nH=1)
    return (tuple(definition(n) for n in ("u1WeylIndexPerCopy", "su2WeylIndexPerCopy", "su3WeylIndexPerCopy")),
            tuple(definition(n) for n in ("bY", "b2", "b3")))


def determinant(x, k, b):
    return sum(x[i]*(k[(i+1)%3]*b[(i+2)%3]-k[(i+2)%3]*b[(i+1)%3]) for i in range(3))


def reconstruct_native_signature():
    """Enumerate only the pinned primitive, face graph and exact frame."""
    from itertools import combinations
    from model.observer import native

    closure = ROOT / "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py"
    tree = ast.parse(closure.read_text())
    faces_node = next(node.value for node in tree.body
                      if isinstance(node, ast.Assign)
                      and any(isinstance(target, ast.Name) and target.id == "ORIENTED_BASE_FACES"
                              for target in node.targets))
    faces = ast.literal_eval(faces_node)
    seams = tuple(sorted({tuple(sorted((face[i], face[(i+1)%3])))
                          for face in faces for i in range(3)}))
    if len(seams) != 30:
        raise ValueError("pinned face graph no longer has thirty seams")

    # Q(sqrt(5)) frame, expressed as exact pairs (rational, sqrt(5) coefficient).
    phi = sp.Rational(1, 2) + sp.sqrt(5)/2
    neg_phi = -phi
    frame = ((-1,phi,0),(1,phi,0),(-1,neg_phi,0),(1,neg_phi,0),
             (0,-1,phi),(0,1,phi),(0,-1,neg_phi),(0,1,neg_phi),
             (phi,0,-1),(phi,0,1),(neg_phi,0,-1),(neg_phi,0,1))
    radius2 = sp.simplify(sum(x*x for x in frame[0]))
    triangle = next(t for t in combinations(range(12),3)
                    if all(tuple(sorted(e)) in seams for e in combinations(t,2)))
    oriented = tuple(zip(triangle, triangle[1:]+triangle[:1]))
    weights = {}
    for seam in seams:
        a,b = seam
        weights[seam] = sp.Matrix([sp.simplify(x/radius2) for x in
            sp.Matrix(frame[a]).cross(sp.Matrix(frame[b]))])
    states = [tuple(int(i in (a,b)) for i in range(12))
              for a in range(12) for b in range(a+1,12)]
    covariance = sp.zeros(3)
    drift_by_state = [sp.zeros(3,1) for _ in states]
    counts = [0,0]
    branches = 0
    for state_index,state in enumerate(states):
        for seam in seams:
            weight = weights[seam]
            for coin in (0,1):
                out, record = native.native(state,seam,coin)
                current = record[1]
                if out not in states or record[0] != coin:
                    raise ValueError("native binary carrier changed")
                increment = weight*current
                drift_by_state[state_index] += increment
                covariance += increment*increment.T
                counts[int(bool(current))] += 1
                branches += 1
    if any(row.applyfunc(sp.simplify) != sp.zeros(3,1) for row in drift_by_state):
        raise ValueError("native counter has drift")
    covariance = (covariance/branches).applyfunc(sp.simplify)
    variance = covariance[0,0]
    if covariance != variance*sp.eye(3) or variance != sp.Rational(4,99):
        raise ValueError("native counter covariance changed")
    # Clock record law is the independent branch count for this same primitive.
    return {"branches": branches, "covariance_diagonal": str(variance),
            "counter_diffusion": str(variance/2),
            "one_attempt_law": [str(Q(n,branches)) for n in counts],
            "triangle": list(triangle)}


def build():
    pins = custody()
    native = reconstruct_native_signature()
    variance = Q(native["covariance_diagonal"])
    diffusion = variance/2
    k, b = oph_columns()
    cofactor = tuple(determinant(tuple(Q(i == j) for i in range(3)), k, b) for j in range(3))
    if determinant(k,k,b) != 0 or determinant(b,k,b) != 0:
        raise ValueError("transport identity failed")
    mixed = tuple(k[i] - b[i]*(i == 0) for i in range(3))
    wrong_branch = (Q(11), Q(1), Q(-3))  # MSSM in Y normalization, control only.
    branch_drift = -determinant(wrong_branch,k,b)
    if determinant(mixed,k,b) == 0 or branch_drift == 0:
        raise ValueError("hostile transport controls lost discrimination")
    return {
        "schema": "arithmon.oph.k7_gauge_scale_contract.v1",
        "baselines": BASELINES, "source_pins": pins,
        "native_signature": {
            "scope": "uniform native attempts; binary total-two current counter and its recovered axial readout",
            "branches": native["branches"], "covariance_diagonal": str(variance),
            "counter_diffusion": str(diffusion), "generator": "D Delta_R3",
            "one_attempt_law": native["one_attempt_law"],
            "unit_extensions": [
                {"length_per_counter_unit": str(l), "time_per_attempt": str(t),
                 "diffusion_in_assigned_units": str(diffusion*l*l/t)}
                for l,t in ((Q(1),Q(1)),(Q(2),Q(1)),(Q(1),Q(2)))],
            "arbitrary_positive_unit_extensions_preserve_native_records": True,
            "physical_continuum_length_identified": False,
            "rg_subtraction_scale_identified": False,
            "laboratory_energy_identified": False,
            "countermodel_proves_universal_physical_nonselection": False,
        },
        "conditional_plane_transport": {
            "completion": {"nG": 3, "nH": 1}, "normalization": "Y",
            "k": list(map(str,k)), "b": list(map(str,b)), "cofactor": list(map(str,cofactor)),
            "identity": "det(x-t*b,k,b)=det(x,k,b) for every x and common scalar t",
            "proof": "det(b,k,b)=0; determinant linear in first column",
            "beta_annihilator": str(sum(c*d for c,d in zip(cofactor,b))),
            "k_annihilator": str(sum(c*d for c,d in zip(cofactor,k))),
            "translation_kernel_dimension": 2,
            "zero_test_depends_only_on_transverse_component": True,
            "common_absolute_mu_needed_for_plane_zero_test_under_hypotheses": False,
            "requires_common_scheme_branch_and_matching": True,
            "matter_kinetic_branch_physically_selected": False,
            "mixed_component_transport_delta_D": str(determinant(mixed,k,b)),
            "wrong_branch_transport_delta_D": str(branch_drift),
            "threshold_basis_delta_D": list(map(str,cofactor)),
            "arbitrary_threshold_or_scheme_shift_preserves_plane": False,
        },
        "k7_source_contract": {
            "erratum_path": "publications/ERRATUM_v3.5.md", "erratum_section": "8.2",
            "alpha_em_scale": "q^2=0", "weak_and_strong_scale": "M_Z",
            "three_couplings_common_scale_claimed": False,
            "three_thirteenths_gut_boundary_and_dependent_rge_chain": "withdrawn",
            "source_clean_transport_to_common_subtraction_prescription": "not established",
        },
        "interface_status": "NO_ADMISSIBLE_K7_GAUGE_VECTOR_ESTABLISHED_BY_PINNED_SOURCES",
        "narrowed_missing_datum": "source-clean common subtraction prescription and compatible running/threshold transport of the heterogeneous K7 inputs; absolute common mu cancels only on the fixed one-loop branch",
        "trust_boundary": {"public_data_used_as_selection_input": False,
                           "target_driven_bridge_choices": False,
                           "k7_coupling_vector_constructed": False,
                           "physical_rg_map_derived": False},
    }


def canonical(doc):
    return (json.dumps(doc,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode()


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--check",action="store_true")
    args=parser.parse_args(); data=canonical(build())
    if args.check:
        if OUT.read_bytes()!=data: raise SystemExit("scale contract receipt differs from live reconstruction")
        print("PASS exact native scale / conditional RG-plane transport")
    else:
        OUT.write_bytes(data); print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__": main()
