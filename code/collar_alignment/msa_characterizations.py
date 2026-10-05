#!/usr/bin/env python3
"""Numerical verification of the Markov-split alignment (MSA) characterizations.

Implements the finite alignment objects introduced for paper-audit issue 001
(GitHub #543), historically labelled `prop:msachar` and `cor:msareduction`.
Those labels are absent from the reorganized active spacetime paper. The
finite Gibbs connection and numerical scope are in GIBBS_SECTOR_AUDIT.md.

Collar model: H = oplus_alpha  H_A (x) H_{bL^alpha} (x) H_{bR^alpha} (x) H_D.
A state is represented blockwise as a list of (weight, block density matrix,
(dA, dL, dR, dD)) triples, which enforces [rho, P_alpha] = 0 by construction;
`embed_blocks` produces the full direct-sum matrix when needed.

Checked characterizations, for faithful blockwise states:

1. EC-aligned normal form  rho = oplus p_a rho_{A bL}^a (x) rho_{bR D}^a
2. modular splitting       log rho in M_L + M_R (blockwise K = h_L (x) 1 + 1 (x) h_R)
3. Takesaki criterion      [log rho, x (x) 1] in B(H_{A bL}) (x) 1  for all x
4. entropic criterion      I(A bL : bR D) = 0 on every block

together with the implication (any of 1-4)  =>  I(A:D|B) = 0, and the
Bell-pair counterexample showing that I(A:D|B) = 0 alone implies none of 1-4.
"""

from __future__ import annotations

import sys
from math import prod
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from quantum_information import (
    conditional_mutual_information as _conditional_mutual_information,
    density_matrix, dimensions, direct_sum_state, faithful_density_matrix, finite_real_scalar,
    gibbs_sectors,
    faithful_log as logm_psd, mutual_information, partial_trace, probabilities,
    one_sided_projection as _one_sided_projection, von_neumann_entropy,
)

LN2 = float(np.log(2.0))


# ---------------------------------------------------------------------------
# linear-algebra helpers
# ---------------------------------------------------------------------------

def dagger(m: np.ndarray) -> np.ndarray:
    return m.conj().T


def conditional_mutual_information(rho, dims, part_a, part_b, part_d):
    """Preserve the collar's historical D-subsystem keyword."""
    return _conditional_mutual_information(rho,dims,part_a,part_b,part_d)


def one_sided_projection(k, d_left, d_right):
    return _one_sided_projection(k,d_left,d_right)


# ---------------------------------------------------------------------------
# blockwise state representation
# ---------------------------------------------------------------------------

Block = tuple[float, np.ndarray, tuple[int, int, int, int]]


def validated_blocks(blocks):
    """Reject empty, malformed or unnormalized block families before scoring."""
    if not isinstance(blocks, (list, tuple)) or not blocks:
        raise ValueError("nonempty collar block family required")
    if any(not isinstance(b, (list, tuple)) or len(b) != 3 for b in blocks):
        raise ValueError("a collar block needs weight, state and four dimensions")
    probabilities([b[0] for b in blocks])
    out = []
    for p, rho, dims in blocks:
        dims = dimensions(dims)
        a = density_matrix(rho)
        if len(dims) != 4 or prod(dims) != len(a):
            raise ValueError("four collar dimensions must match each block")
        out.append((p, a, dims))
    return out


def block_dims(block: Block) -> list[int]:
    return list(block[2])


def embed_blocks(blocks: list[Block]) -> tuple[np.ndarray, list[np.ndarray]]:
    """Direct-sum embedding; returns (rho_full, center projectors P_alpha)."""
    blocks = validated_blocks(blocks)
    rho = direct_sum_state([b[0] for b in blocks], [b[1] for b in blocks])
    projectors, offset = [], 0
    for _, state, _ in blocks:
        proj = np.zeros_like(rho)
        proj[offset:offset+len(state), offset:offset+len(state)] = np.eye(len(state))
        projectors.append(proj)
        offset += len(state)
    return rho, projectors


def collar_cmi(blocks: list[Block]) -> float:
    """I(A:D|B) of the blockwise state; entropies add over blocks."""
    total = 0.0
    for p, rho_a, dims in validated_blocks(blocks):
        if p <= 0.0:
            continue
        d = list(dims)
        total += p * conditional_mutual_information(rho_a, d, [0], [1, 2], [3])
    # Classical block label contributes equally to S(AB), S(BD), S(B), S(ABD)
    # and cancels in the CMI combination, so the blockwise sum is exact.
    return total


# ---------------------------------------------------------------------------
# characterization 4: entropic criterion (no faithfulness needed)
# ---------------------------------------------------------------------------

def entropic_alignment_defect(blocks: list[Block]) -> float:
    """max_alpha I(A bL^alpha : bR^alpha D); zero iff EC-aligned (item 4)."""
    worst = 0.0
    for p, rho_a, dims in validated_blocks(blocks):
        if p <= 0.0:
            continue
        d = list(dims)
        worst = max(worst, mutual_information(rho_a, d, [0, 1], [2, 3]))
    return worst


def is_ec_aligned(blocks: list[Block], tol: float = 1e-9) -> bool:
    tol = finite_real_scalar(tol, "alignment tolerance")
    if tol <= 0:
        raise ValueError("alignment tolerance must be finite and positive")
    return entropic_alignment_defect(blocks) < tol


# ---------------------------------------------------------------------------
# characterization 2: modular splitting  log rho in M_L + M_R
# ---------------------------------------------------------------------------

def modular_splitting_defect(blocks: list[Block]) -> float:
    """max_alpha || log rho^alpha - Proj_{M_L + M_R}(log rho^alpha) ||_inf.

    Zero iff log rho in M_L + M_R blockwise (item 2). Requires faithful blocks.
    """
    worst = 0.0
    for p, rho_a, dims in validated_blocks(blocks):
        if p <= 0.0:
            continue
        d_a, d_l, d_r, d_d = dims
        k = logm_psd(rho_a)
        proj = one_sided_projection(k, d_a * d_l, d_r * d_d)
        worst = max(worst, float(np.linalg.norm(k - proj, ord=2)))
    return worst


# ---------------------------------------------------------------------------
# characterization 3: Takesaki commutator criterion
# ---------------------------------------------------------------------------

def takesaki_defect(blocks: list[Block]) -> float:
    """Complete matrix-unit test of the Takesaki commutator criterion.

    Every matrix unit E_ij of the left factor is tested. By complex
    linearity, vanishing on this basis is equivalent to vanishing for all
    left observables. No sample count or random seed can weaken the check.
    The returned maximum is a finite-basis diagnostic, not a proof about
    exact zeros inferred from floating-point tolerances.
    """
    worst = 0.0
    for p, rho_a, dims in validated_blocks(blocks):
        if p == 0:
            continue
        d_a, d_l, d_r, d_d = dims
        dl_tot, dr_tot = d_a*d_l, d_r*d_d
        k = logm_psd(rho_a)
        for i in range(dl_tot):
            for j in range(dl_tot):
                x = np.zeros((dl_tot, dl_tot))
                x[i, j] = 1
                x_full = np.kron(x, np.eye(dr_tot))
                comm = k @ x_full - x_full @ k
                comm4 = comm.reshape(dl_tot, dr_tot, dl_tot, dr_tot)
                comm_l = np.trace(comm4, axis1=1, axis2=3)/dr_tot
                resid = comm - np.kron(comm_l, np.eye(dr_tot))
                worst = max(worst, float(np.linalg.norm(resid, ord=2)))
    return worst


# ---------------------------------------------------------------------------
# state constructors
# ---------------------------------------------------------------------------

def random_density(dim: int, rng: np.random.Generator, floor: float = 0.05) -> np.ndarray:
    """Random faithful density matrix with spectral floor."""
    g = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
    rho = g @ dagger(g)
    rho = rho / np.trace(rho).real
    rho = (1.0 - floor) * rho + floor * np.eye(dim) / dim
    return (rho + dagger(rho)) / 2.0


def random_aligned_blocks(rng: np.random.Generator,
                          spec: list[tuple[int, int, int, int]] | None = None
                          ) -> list[Block]:
    """Random EC-aligned exact Markov state: blockwise product across the cut."""
    if spec is None:
        spec = [(2, 2, 2, 2), (2, 2, 2, 2)]
    weights = rng.dirichlet(np.ones(len(spec)))
    blocks: list[Block] = []
    for w, dims in zip(weights, spec):
        d_a, d_l, d_r, d_d = dims
        left = random_density(d_a * d_l, rng)
        right = random_density(d_r * d_d, rng)
        blocks.append((float(w), np.kron(left, right), dims))
    return blocks


def random_generic_blocks(rng: np.random.Generator,
                          spec: list[tuple[int, int, int, int]] | None = None
                          ) -> list[Block]:
    """Random faithful blockwise state with generic cross-cut correlation."""
    if spec is None:
        spec = [(2, 2, 2, 2)]
    weights = rng.dirichlet(np.ones(len(spec)))
    return [(float(w), random_density(int(np.prod(d)), rng), d)
            for w, d in zip(weights, spec)]


def bell_state(dim: int = 2) -> np.ndarray:
    vec = np.eye(dim).reshape(-1) / np.sqrt(dim)
    return np.outer(vec, vec.conj())


def bell_counterexample() -> list[Block]:
    """rho = Phi(A, bR) (x) Phi(bL, D): exactly Markov, maximally misaligned.

    Single center block, qubits ordered (A, bL, bR, D).
    """
    phi = bell_state(2)
    rho = np.kron(phi, phi)  # order (A, bR, bL, D)
    rho = rho.reshape([2] * 8)
    perm = [0, 2, 1, 3]      # -> (A, bL, bR, D)
    rho = rho.transpose(perm + [p + 4 for p in perm]).reshape(16, 16)
    return [(1.0, rho, (2, 2, 2, 2))]


def gibbs_blocks(hamiltonians: list[tuple[np.ndarray, tuple[int, int, int, int]]],
                 central_energies: list[float], beta: float = 1.0) -> list[Block]:
    """Faithful direct-sum Gibbs state, retaining relative sector partition functions.

    Scalar origins cannot change alignment within a sector. See
    GIBBS_SECTOR_AUDIT.md for the energy/log-state identity and precision scope.
    """
    if not isinstance(hamiltonians, (list, tuple)) or not hamiltonians:
        raise ValueError("one central energy is required per nonempty sector")
    operators, shapes = [], []
    for entry in hamiltonians:
        if not isinstance(entry, (list, tuple)) or len(entry) != 2:
            raise ValueError("a Hamiltonian and four subsystem dimensions are required")
        operator, dims = entry
        dims = dimensions(dims)
        if len(dims) != 4 or np.shape(operator) != (prod(dims),)*2:
            raise ValueError("sector Hamiltonian and subsystem dimensions must match")
        operators.append(operator)
        shapes.append(dims)
    return [(p, state, dims) for (p, state), dims in
            zip(gibbs_sectors(operators, central_energies, beta), shapes)]


def one_sided_hamiltonian(rng: np.random.Generator,
                          dims: tuple[int, int, int, int]) -> np.ndarray:
    """H = h_L (x) 1 + 1 (x) h_R: every cross-cut coupling routed through the center."""
    d_a, d_l, d_r, d_d = dims
    h_l = rng.normal(size=(d_a * d_l, d_a * d_l))
    h_l = (h_l + h_l.T) / 2.0
    h_r = rng.normal(size=(d_r * d_d, d_r * d_d))
    h_r = (h_r + h_r.T) / 2.0
    return np.kron(h_l, np.eye(d_r * d_d)) + np.kron(np.eye(d_a * d_l), h_r)


def cross_cut_coupling(rng: np.random.Generator,
                       dims: tuple[int, int, int, int]) -> np.ndarray:
    """A traceless non-central bL--bR interaction term crossing the cut."""
    d_a, d_l, d_r, d_d = dims
    g_l = rng.normal(size=(d_l, d_l))
    g_l = (g_l + g_l.T) / 2.0
    g_l -= np.trace(g_l) / d_l * np.eye(d_l)
    g_r = rng.normal(size=(d_r, d_r))
    g_r = (g_r + g_r.T) / 2.0
    g_r -= np.trace(g_r) / d_r * np.eye(d_r)
    return np.kron(np.kron(np.eye(d_a), g_l), np.kron(g_r, np.eye(d_d)))
