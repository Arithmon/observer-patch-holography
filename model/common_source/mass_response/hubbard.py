"""Small periodic Hubbard controls for a separately declared coherent law.

The source preparation supplies coefficients; assigning t=clock/2 and
U=interaction to the spinful Hubbard chain is an explicit architecture choice.
These sparse numerical ground-state calculations are diagnostics, distinct
from both the executed four-mode ring and the thermodynamic gap certificate.
Neither a finite-volume addition energy nor an eigensolver residual identifies
an observed particle mass.
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations
import json
import math

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags, eye, kron
from scipy.sparse.linalg import eigsh

from . import model


MAX_SECTOR_DIMENSION = 4900


def _length(value):
    if type(value) is not int or value not in (4, 6, 8):
        raise ValueError("this finite diagnostic supports periodic lengths 4, 6, 8")
    return value


def _coupling(value, name, *, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float)):
        raise ValueError(f"finite real {name} required")
    value = float(value)
    if not math.isfinite(value) or value < 0 or (positive and value == 0):
        raise ValueError(f"finite {'positive' if positive else 'nonnegative'} {name} required")
    return value


def spin_basis(length, count):
    """Occupation bit strings for one fixed spin population, in ascending order."""
    _length(length)
    if type(count) is not int or not 0 <= count <= length:
        raise ValueError("integer spin population between zero and length required")
    return _spin_basis(length, count)


@lru_cache(maxsize=27)
def _spin_basis(length, count):
    return tuple(sorted(sum(1 << j for j in occupied)
                        for occupied in combinations(range(length), count)))


@lru_cache(maxsize=27)
def _unit_hopping(length, count):
    basis = spin_basis(length, count)
    index = {bits: j for j, bits in enumerate(basis)}
    rows, columns, values = [], [], []
    for column, bits in enumerate(basis):
        for a in range(length):
            b = (a + 1) % length
            for origin, target in ((a, b), (b, a)):
                if not (bits >> origin) & 1 or (bits >> target) & 1:
                    continue
                low, high = sorted((origin, target))
                between = ((1 << high) - 1) ^ ((1 << (low + 1)) - 1)
                parity = (bits & between).bit_count() % 2
                rows.append(index[bits ^ (1 << origin) ^ (1 << target)])
                columns.append(column)
                values.append(1.0 if parity else -1.0)
    return coo_matrix((values, (rows, columns)),
                      shape=(len(basis), len(basis)), dtype=float).tocsr()


def sector_hamiltonian(length, n_up, n_down, *, hopping=None, repulsion=None):
    """Return H=-t sum(c†_j c_{j+1}+h.c.)+U sum(n_up n_down).

    All up-spin modes precede all down-spin modes. The tensor product basis
    orders up occupations first and down occupations second. The wrapped bond
    carries its fermionic occupation sign. No chemical-potential centering is
    included; it cancels from the second particle-number difference below.
    """
    up, down = spin_basis(length, n_up), spin_basis(length, n_down)
    source = model.parameters()
    t = _coupling(source.clock / 2 if hopping is None else hopping,
                  "hopping", positive=True)
    u = _coupling(source.interaction if repulsion is None else repulsion, "repulsion")
    dimension = len(up) * len(down)
    if dimension > MAX_SECTOR_DIMENSION:
        raise ValueError("sector exceeds the small-volume diagnostic budget")
    hopping_matrix = t * (
        kron(_unit_hopping(length, n_up), eye(len(down), format="csr"), format="csr")
        + kron(eye(len(up), format="csr"), _unit_hopping(length, n_down), format="csr"))
    doubles = np.fromiter(((a & b).bit_count() for a in up for b in down),
                          dtype=float, count=dimension)
    result = csr_matrix(hopping_matrix + diags(u * doubles, format="csr"))
    result.eliminate_zeros()
    return result


def sector_ground_state(length, n_up, n_down, *, hopping=None, repulsion=None):
    """Numerical lowest eigenpair diagnostic, with a reproducible start vector.

    The residual measures eigenpair consistency; it is not an interval proof of
    eigenvalue ordering. The vector is deliberately omitted from JSON results.
    """
    h = sector_hamiltonian(length, n_up, n_down,
                           hopping=hopping, repulsion=repulsion)
    if h.shape[0] == 1:
        energy, vector = float(h[0, 0]), np.ones(1)
    else:
        start = np.random.default_rng(118 + 31 * length + 7 * n_up + n_down).normal(
            size=h.shape[0])
        energies, vectors = eigsh(h, k=1, which="SA", v0=start,
                                  tol=1e-12, maxiter=20000)
        energy, vector = float(energies[0]), vectors[:, 0]
    return {
        "n_up": n_up, "n_down": n_down,
        "dimension": int(h.shape[0]), "nonzero_entries": int(h.nnz),
        "energy": energy,
        "eigenpair_residual_norm": float(np.linalg.norm(h @ vector - energy * vector)),
        "eigenvector_norm_error": abs(float(np.linalg.norm(vector)) - 1.0),
        "method": ("one-dimensional sector, direct scalar evaluation" if h.shape[0] == 1
                   else "scipy.sparse.linalg.eigsh, smallest algebraic eigenvalue"),
        "deterministic_start": True,
        "eigenvalue_ordering_certified": False,
    }


def free_charge_gap(length, *, hopping=None):
    """Exact free-chain shell formula, evaluated with ordinary floating point.

    At half filling the periodic chain has zero addition-energy difference for
    L divisible by four, and 4 t sin(pi/L) for L=4m+2. The latter is a finite
    shell effect; both sequences approach zero in the thermodynamic limit.
    """
    _length(length)
    t = _coupling(model.parameters().clock / 2 if hopping is None else hopping,
                  "hopping", positive=True)
    return 0.0 if length % 4 == 0 else 4 * t * math.sin(math.pi / length)


def _charge_row(length, t, u):
    count = length // 2
    sectors = {
        "half_filled": sector_ground_state(length, count, count, hopping=t, repulsion=u),
        "one_added": sector_ground_state(length, count + 1, count, hopping=t, repulsion=u),
        "one_removed": sector_ground_state(length, count - 1, count, hopping=t, repulsion=u),
    }
    e0, eplus, eminus = (sectors[key]["energy"]
                         for key in ("half_filled", "one_added", "one_removed"))
    return {
        "E0": e0, "Eplus": eplus, "Eminus": eminus,
        "charge_gap": eplus + eminus - 2 * e0,
        "particle_hole_energy_residual": eplus - eminus - u,
        "sectors": sectors,
    }


def finite_volume_controls(lengths=(4, 6, 8), *, verify_free=(4, 6)):
    """Return small-volume numerical controls, never a size extrapolation.

    Inputs are source-computed native-reference coefficients. Choosing this
    spinful nearest-neighbor law and its half-filled ground-state sector remains
    explicit. The thermodynamic certificate concerns the infinite version of
    this law; it does not certify these floating-point spectra.
    """
    lengths, verify_free = tuple(lengths), tuple(verify_free)
    if not lengths or len(set(lengths)) != len(lengths):
        raise ValueError("nonempty distinct finite lengths required")
    for length in (*lengths, *verify_free):
        _length(length)
    source = model.parameters()
    t, u = source.clock / 2, source.interaction
    rows = []
    for length in lengths:
        row = _charge_row(length, t, u)
        row.update({
            "length": length, "total_population": length,
            "charge_gap_over_clock": row["charge_gap"] / source.clock,
            "free_charge_gap_analytic": free_charge_gap(length, hopping=t),
            "free_shell_formula": "0" if length % 4 == 0 else "4*t*sin(pi/L)",
        })
        row["free_numerical_check"] = (
            _charge_row(length, t, 0.0) if length in verify_free else None)
        rows.append(row)
    return {
        "source": source.record(),
        "architecture": {
            "law": "periodic one-dimensional spinful repulsive Hubbard chain",
            "hopping": t, "repulsion": u,
            "coefficient_map": {"t": "native clock/2", "U": "native interaction"},
            "half_filled_ground_state_is_contract": True,
            "native_reference_is_hubbard_vacuum": False,
            "distinct_from_executed_four_mode_ring": True,
        },
        "rows": rows,
        "interpretation": {
            "numerical_diagnostic": True, "certified_spectrum": False,
            "thermodynamic_extrapolation_performed": False,
            "finite_addition_gap_is_physical_mass": False,
            "isolated_electron_pole_established": False,
            "three_spatial_dimensions_established": False,
            "standard_model_species_identified": False,
            "free_finite_gap_is_shell_effect": True,
            "time_parameter": "chosen coherent Hamiltonian parameter; source clock frequency reported",
        },
    }


if __name__ == "__main__":
    print(json.dumps(finite_volume_controls(), indent=2, sort_keys=True))
