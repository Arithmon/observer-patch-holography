"""Native complement repair selects a finite bulk/edge entropy completion.

On the binary twelve-port total-two source, choose a nonempty cap R with
at least three complementary ports B and a connected induced B graph.
Only seams inside B are attempted, uniformly and independently, with fair
placement coins. All attempt records are discarded. The same C/E primitives
then converge to

  E_R(X) = direct_sum_j Tr_B(P_j X P_j) tensor I_(b_j)/b_j,
  a_j = binomial(|R|,j), b_j = binomial(|B|,2-j).

This map preserves the actual number-preserving cap algebra
direct_sum_j M_(a_j) tensor I_(b_j). It is an orthogonal trace-preserving
conditional expectation, not an assigned reference or a cap-state reset.
For a fixed retained cap marginal, E_R(rho) is the unique entropy maximum:

  S(E_R(rho))-S(rho) = D(rho || E_R(rho)) >= 0.

The equality follows since log(E_R(rho)) is in the retained algebra on
its support; positivity puts rho on that support. Equality in relative
entropy forces rho=E_R(rho). Its entropy is H(p)+sum p_j S(rho_R,j)
plus sum p_j log(b_j). Thus the dimensions and central weights in the
finite bulk/edge identity follow from this native protocol. They describe
complement multiplicity, with no physical-area or boost identification.

Proof of convergence: apply the existing native signed-difference reduction
to B. Zero outside difference gives a connected fixed-number exclusion
chain and preserves all cap matrix units within one charge. Every other
outside difference is nonconstant: a constant +/-1 would require a full
B configuration, forbidden by |B|>2. Connectedness gives C killing at
least 1/(2M), and E killing at least 1/M. Every local channel is a positive
self-adjoint Hilbert--Schmidt contraction. Therefore, using the inherited
graph certificate gamma_graph and gamma=min(gamma_graph,1/(2M)),

  ||Phi_L^n-E_R||_(2->2) <= (1-gamma)^n,
  dist((Phi_L^n tensor id)(rho),(E_R tensor id)(rho))
      <= min(1,33*(1-gamma)^n), L=C,E.

The reference-uniform bound follows by purifying a source input with
Schmidt rank at most 66, exactly as in QUANTUM_MIXING.md. Any later common
complete instrument preserves it, retaining all failures and outcomes.
The parameter n counts discarded native attempts, not physical time.
Allowing both empty and full B invalidates this completion theorem:
their cross-charge coherence can be fixed by every local repair.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import comb, lcm
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from .construction import graph_certificate, preparation_budget
from ..observer import native, quantum_mixing

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "cut_results.json"
CAP = (0, 1, 5, 7)
SOURCE_PATHS = (
    "model/native_equilibration/cut_completion.py",
    "model/native_equilibration/construction.py",
    "model/observer/native.py", "model/observer/quantum_mixing.py",
    "model/observer/QUANTUM_MIXING.md",
    "model/common_source/slow_geometry.py",
)


@dataclass(frozen=True)
class Sector:
    charge: int
    cap_states: tuple
    complement_states: tuple
    indices: tuple


def _cap(cap):
    if (type(cap) is not tuple or not 0 < len(cap) < 12
            or any(type(x) is not int or not 0 <= x < 12 for x in cap)
            or tuple(sorted(set(cap))) != cap):
        raise ValueError("a nonempty proper increasing tuple of native port integers is required")
    return cap


def cut_sectors(cap=CAP):
    """Tensor factors enumerated directly from the native binary source."""
    return _cut_sectors(_cap(cap))


@lru_cache(maxsize=32)
def _cut_sectors(cap):
    complement = tuple(i for i in range(12) if i not in cap)
    pairs = { (tuple(x[i] for i in cap), tuple(x[i] for i in complement)): k
              for k, x in enumerate(native.states(binary=True)) }
    result = []
    for charge in range(max(0, 2-len(complement)), min(len(cap), 2)+1):
        left = tuple(sorted({a for a, b in pairs if sum(a) == charge}))
        right = tuple(sorted({b for a, b in pairs if sum(a) == charge}))
        indices = tuple(tuple(pairs[a, b] for b in right) for a in left)
        if len(left) != comb(len(cap), charge) or len(right) != comb(len(complement), 2-charge):
            raise ArithmeticError("source sector factorization failed")
        result.append(Sector(charge, left, right, indices))
    return tuple(result)


def cut_certificate(cap=CAP, *, precision_bits=10):
    """Exact native graph/gap and reference-uniform finite-attempt budget."""
    cap = _cap(cap)
    complement = tuple(i for i in range(12) if i not in cap)
    if len(complement) < 3:
        raise ValueError("at least three complement ports exclude empty/full coherence")
    edges = tuple((complement.index(a), complement.index(b))
                  for a, b in quantum_mixing.seams() if a in complement and b in complement)
    graph = graph_certificate(len(complement), edges)
    gamma = min(F(graph["quantum_gap_lower_per_attempt"]), F(1, 2*len(edges)))
    budget = preparation_budget(66, gamma, precision_bits)
    epsilon = F(budget["trace_distance_upper"])
    # Audenaert's bound in nats is epsilon*log(65)+h(epsilon), increasing
    # up to epsilon=65/66. Here epsilon=33/2^(7+precision_bits)<1/2,
    # log(1/epsilon)<(precision_bits+2)log2, log65<7log2, and
    # -(1-epsilon)log(1-epsilon)<=epsilon. The exact inequality log2<7/10
    # follows already from the degree-four positive Taylor sum for exp(7/10).
    entropy_upper = epsilon*(1+F(7, 10)*(precision_bits+9))
    return {"cap": list(cap), "complement": list(complement),
            "internal_native_seams": len(edges), "source_dimension": 66,
            "fixed_algebra_dimension": sum(len(s.cap_states)**2 for s in cut_sectors(cap)),
            "sector_dimensions": [{"cap_charge": s.charge, "bulk_dimension": len(s.cap_states),
                                   "complement_dimension": len(s.complement_states)}
                                  for s in cut_sectors(cap)],
            "graph_certificate": graph, "gap_lower_per_attempt": str(gamma),
            "preparation": budget,
            "source_entropy_deficit_upper_nats": str(entropy_upper)}


def completion_counts(cap=CAP):
    """Integer matrix-unit representation of the independently defined E_R.

    It remains a valid conditional expectation outside the convergence
    domain; cut_certificate separately rejects unproved native selection.
    """
    sectors = cut_sectors(cap)
    scale = lcm(*(len(s.complement_states) for s in sectors))
    rows, cols, data = [], [], []
    for sector in sectors:
        a, b = len(sector.cap_states), len(sector.complement_states)
        ix = sector.indices
        for i, j, u, v in product(range(a), range(a), range(b), range(b)):
            rows.append(66*ix[i][v]+ix[j][v])
            cols.append(66*ix[i][u]+ix[j][u])
            data.append(scale//b)
    return coo_matrix((data, (rows, cols)), shape=(4356, 4356), dtype=np.int64).tocsr(), scale


def completed_palm(cap=CAP):
    """Exact E_R(sigma) in native binary-state order, from native seam degrees.

    This preserves the source Palm cap marginal and gives a faithful
    reference whose modular flow leaves the cap algebra invariant. Finite
    preparation and a compiled modular word need composed error bounds.
    """
    basis = native.states(binary=True)
    degree = tuple(sum(x[a] != x[b] for a, b in quantum_mixing.seams()) for x in basis)
    total = sum(degree)
    result = [F(0)]*66
    for sector in cut_sectors(cap):
        for row in sector.indices:
            value = F(sum(degree[i] for i in row), total*len(row))
            for i in row:
                result[i] = value
    if sum(result) != 1 or min(result) <= 0:
        raise ArithmeticError("native Palm completion failed normalization or faithfulness")
    return tuple(result)


def native_cut_counts(cap=CAP, law="C"):
    """Actual source Kraus records, including cross-charge coherences."""
    cap = _cap(cap)
    if law not in ("C", "E"):
        raise ValueError("candidate C or E required")
    basis = native.states(binary=True)
    positions = {x: i for i, x in enumerate(basis)}
    edges = tuple((a, b) for a, b in quantum_mixing.seams() if a not in cap and b not in cap)
    if not edges:
        raise ValueError("the complement has no native repair seam")
    rows, cols, data = [], [], []
    for seam, coin in product(edges, (0, 1)):
        images = []
        for x in basis:
            y, (_, current) = native.native(x, seam, coin)
            images.append((positions[y], current if law == "C" else x[seam[0]]-x[seam[1]]))
        for i, (a, left_record) in enumerate(images):
            for j, (b, right_record) in enumerate(images):
                if left_record == right_record:
                    rows.append(66*a+b)
                    cols.append(66*i+j)
                    data.append(1)
    return coo_matrix((data, (rows, cols)), shape=(4356, 4356), dtype=np.int64).tocsr(), 2*len(edges)


def _zero(matrix):
    matrix.eliminate_zeros()
    return matrix.nnz == 0


def exact_audit(cap=CAP):
    """All 4,356 matrix units: idempotence, invariance and cap preservation."""
    e, es = completion_counts(cap)
    if not _zero(e@e-es*e) or not _zero(e-e.T):
        raise ArithmeticError("the tensor completion is not an orthogonal projection")
    trace = np.zeros(4356, dtype=np.int64)
    trace[np.arange(66)*67] = 1
    if not np.array_equal(trace@e, es*trace):
        raise ArithmeticError("the completion changed total trace")
    laws = {}
    for law in ("C", "E"):
        p, ps = native_cut_counts(cap, law)
        if (not _zero(p-p.T) or not _zero(p@e-ps*e)
                or not _zero(e@p-ps*e) or not np.array_equal(trace@p, ps*trace)):
            raise ArithmeticError("native cut repair failed the fixed-algebra identity")
        # One integral digest fixes all actual Kraus entries, without BLAS output.
        triples = sorted(zip(*p.nonzero(), strict=True))
        digest = hashlib.sha256(json.dumps([[int(i), int(j), int(p[i, j])]
                                            for i, j in triples], separators=(",", ":")).encode()).hexdigest()
        laws[law] = {"native_channel_integer_sha256": digest, "denominator": ps,
                     "all_matrix_units": 4356, "completion_fixed_and_cap_preserved": True}
    return {"completion_denominator": es, "completion_rank": int(e.diagonal().sum()//es),
            "exact_idempotence_and_trace_preservation": True, "laws": laws}


def run():
    completed = completed_palm()
    degree = tuple(sum(x[a] != x[b] for a, b in quantum_mixing.seams())
                   for x in native.states(binary=True))
    palm_change = sum(abs(p-F(d, sum(degree))) for p, d in zip(completed, degree, strict=True))/2
    return {"schema": "oph.native_equilibration.cut_completion.v1",
            "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
            "native_cut": cut_certificate(), "all_matrix_unit_audit": exact_audit(),
            "completed_palm_reference": {
                "minimum_eigenvalue": str(min(completed)),
                "half_trace_distance_from_original_palm": str(palm_change),
                "cap_marginal_preserved": True,
                "bare_cap_algebra_modular_invariant": True,
                "per_microstate_eigenvalues_by_cap_state": [
                    {"cap_charge": s.charge, "complement_multiplicity": len(s.complement_states),
                     "eigenvalues": [str(completed[row[0]]) for row in s.indices]}
                    for s in cut_sectors()]},
            "derived_entropy_contract": {
                "completion": "E_R(rho)=direct_sum_j Tr_B(P_j rho P_j) tensor I_(b_j)/b_j",
                "unique_maximum": "S(E_R(rho))-S(rho)=D(rho||E_R(rho)) on the fixed cap marginal fibre",
                "bulk": "H(p)+sum_j p_j S(rho_R,j)",
                "edge": "sum_j p_j log b_j; b_j=binomial(12-|R|,2-j)",
                "modular_split": "K_completion=K_bulk+Z, Z_j=log(b_j) I, on the state's support",
                "stationarity": "dS_sigma[X]=0 if sigma=E_R(sigma) is faithful and E_R(X)=0",
                "complete_future_instrument": "same half-trace bound, arbitrary untouched reference, every outcome retained"},
            "scope": {"native_completion_selected": True, "cap_and_complement_only_scheduler_supplied": True,
                      "coherent_cap_observables_declared_accessible": False,
                      "input_cap_marginal_or_sector_probabilities_selected": False,
                      "retained_preparation_tape_covered": False,
                      "physical_edge_area_volume_stress_or_BW_boost_identified": False,
                      "generalized_entropy_stationarity_under_geometry_variations_derived": False,
                      "Palm_reference_preserved_by_this_tracial_completion": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, sort_keys=True)+"\n"
    if args.check:
        if args.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit("native cut-completion receipt drift")
    else:
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    print("native cut entropy completion: all checks passed")


if __name__ == "__main__":
    main()
