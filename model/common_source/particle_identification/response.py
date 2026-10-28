"""Complete finite CAR and density responses of the executed ring law.

The Lehmann identities are exact finite-dimensional identities. Their arrays
here are floating evaluations, not outward numerical enclosures. No measured
mass, continuum particle pole, Lorentz spin or Standard Model identification
is supplied or inferred. ``transfer_read`` does execute a new neutral source
instrument on the caller's continuing work/reference state.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.linalg import eigh

from ..mass_response import history, model
from ...observer.source_scalar import compiler


VISIBILITY_TOLERANCE = 1e-12
DEGENERACY_TOLERANCE = 1e-12


def _real(value, label):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"finite real {label} required")
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"finite real {label} required") from exc
    if not math.isfinite(out) or np.ndim(value) != 0 or isinstance(value, (str, bytes)):
        raise ValueError(f"finite real {label} required")
    return out


def _freeze(value):
    raw = np.asarray(value, complex)
    # A bytes-backed array cannot be made writable by its consumer.
    return np.frombuffer(raw.tobytes(), dtype=complex).reshape(raw.shape)


def _matrix_record(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


@dataclass(frozen=True)
class Pole:
    energy_cost: float
    frequency: float
    charge_change: int
    target_number: int
    multiplicity: int
    residue: np.ndarray

    def record(self):
        trace = float(np.trace(self.residue).real)
        return {"energy_cost": self.energy_cost, "frequency": self.frequency,
                "charge_change": self.charge_change, "target_number": self.target_number,
                "energy_eigenspace_multiplicity": self.multiplicity,
                "residue": _matrix_record(self.residue),
                "local_weights": np.diag(self.residue).real.tolist(),
                "total_operator_weight": trace,
                "numerically_visible": trace > VISIBILITY_TOLERANCE}


@dataclass(frozen=True)
class FiniteResponse:
    interaction: bool
    chemical_potential: float
    hamiltonian: np.ndarray
    ground_energy: float
    ground: np.ndarray
    density_means: np.ndarray
    addition: tuple[Pole, ...]
    removal: tuple[Pole, ...]
    density: tuple[Pole, ...]


def _sector(hamiltonian, number):
    indices = [i for i in range(16) if i.bit_count() == number]
    energies, block = eigh(hamiltonian[np.ix_(indices, indices)])
    vectors = np.zeros((16, len(indices)), complex)
    vectors[indices] = block
    return energies, vectors


def _poles(energies, vectors, channels, ground_energy, charge, *, removal=False):
    amplitudes = vectors.conj().T @ channels
    groups = []
    for index, energy in enumerate(energies):
        if not groups or abs(energy - energies[groups[-1][0]]) > DEGENERACY_TOLERANCE:
            groups.append([])
        groups[-1].append(index)
    result = []
    for group in groups:
        rows = amplitudes[group]
        # Addition: <g|c_i|m><m|c_j*|g>.
        # Removal:  <g|c_j*|m><m|c_i|g>.
        residue = rows.T @ rows.conj() if removal else rows.conj().T @ rows
        energy = float(np.mean(energies[group]) - ground_energy)
        result.append(Pole(energy, -energy if removal else energy, charge,
                           2 + charge, len(group), _freeze(residue)))
    return tuple(result)


def finite_response(*, interaction=True, chemical_potential=0.0):
    """All N=2-ground CAR addition/removal and centered-density residues.

    ``H(mu)=H(0)+mu*N`` follows the existing model's sign convention.
    The ground is the N=2 ground, even if a supplied mu makes another charge
    sector globally preferable. Costs are then signed energy differences.
    Equal-energy eigenspaces are summed before exporting residues, removing
    the arbitrary eigenbasis within the ring's degenerate neutral subspace.
    No pole or small residue is discarded.
    """
    if type(interaction) is not bool:
        raise ValueError("Boolean interaction control required")
    mu = _real(chemical_potential, "chemical potential")
    h = model.ring_hamiltonian(interaction=interaction, chemical_potential=mu)
    energy, ground = model.sector_ground(h, 2)
    ground = np.asarray(ground, complex)
    ladders = model.annihilators()
    numbers = tuple(a.conj().T @ a for a in ladders)
    means = np.array([np.vdot(ground, n @ ground).real for n in numbers])
    create = np.column_stack([a.conj().T @ ground for a in ladders])
    destroy = np.column_stack([a @ ground for a in ladders])
    centered = np.column_stack([(n - mean * np.eye(16)) @ ground
                               for n, mean in zip(numbers, means, strict=True)])
    e3, v3 = _sector(h, 3)
    e1, v1 = _sector(h, 1)
    e2, v2 = _sector(h, 2)
    return FiniteResponse(interaction, mu, _freeze(h), energy, _freeze(ground),
                          _freeze(means),
                          _poles(e3, v3, create, energy, 1),
                          _poles(e1, v1, destroy, energy, -1, removal=True),
                          _poles(e2, v2, centered, energy, 0))


def charged_correlation(response, time, *, kind="addition"):
    """C+(t)=<c_i(t)c_j*> or C-(t)=<c_j* c_i(t)> in the N=2 ground.

    Addition frequencies are +Delta+ and removal frequencies -Delta-.
    The parameter is the same programmed model/modular parameter used by
    the continuing ring; it is not host time or calibrated laboratory time.
    """
    if type(response) is not FiniteResponse or kind not in ("addition", "removal"):
        raise ValueError("a finite response and addition/removal channel required")
    time = _real(time, "model parameter")
    poles = response.addition if kind == "addition" else response.removal
    return sum((np.exp(-1j * pole.frequency * time) * pole.residue for pole in poles),
               np.zeros((4, 4), complex))


def density_correlation(response, time):
    """Connected <delta n_i(t) delta n_j> with only charge-zero targets."""
    if type(response) is not FiniteResponse:
        raise ValueError("a finite response required")
    time = _real(time, "model parameter")
    return sum((np.exp(-1j * pole.frequency * time) * pole.residue for pole in response.density),
               np.zeros((4, 4), complex))


def green_function(response, frequency):
    """Finite retarded CAR resolvent, evaluated strictly above the real axis."""
    if type(response) is not FiniteResponse or isinstance(frequency, (bool, np.bool_)):
        raise ValueError("a finite response and upper-half-plane frequency required")
    try:
        z = complex(frequency)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("an upper-half-plane frequency required") from exc
    if not math.isfinite(z.real) or not math.isfinite(z.imag) or z.imag <= 0:
        raise ValueError("an upper-half-plane frequency required")
    return sum((pole.residue / (z - pole.frequency)
                for pole in (*response.addition, *response.removal)), np.zeros((4, 4), complex))


def _sites(site_a, site_b, quadrature):
    if (type(site_a) is not int or type(site_b) is not int
            or not 0 <= site_a < 4 or not 0 <= site_b < 4
            or type(quadrature) is not str or quadrature not in ("X", "Y", "Z")):
        raise ValueError("two cell-local sites 0..3 and an X/Y/Z quadrature required")
    return site_a, 4 + site_b


def transfer_operator(site_a=0, site_b=0, *, quadrature="X"):
    """Even two-cell transfer X=T+T*, Y=-i(T-T*), or Z=n_A-n_B.

    A occupies the low four bits; T=C_A* C_B raises A and lowers B.
    C_A=I_B tensor c_A, C_B=c_B tensor (-1)^N_A. The parity factor is
    essential: distinct cells are fermion modes, not commuting odd qubits.
    Every returned observable commutes with total charge and fermion parity.
    This is global occupation-U(1) invariance. If that symmetry were gauged,
    a separated transfer would need a gauge link/dressing; none is supplied.
    """
    _sites(site_a, site_b, quadrature)
    a, b = model.annihilators()[site_a], model.annihilators()[site_b]
    if quadrature == "Z":
        return np.kron(np.eye(16), a.conj().T @ a) - np.kron(b.conj().T @ b, np.eye(16))
    parity = np.diag((-1.) ** model.number_diagonal())
    transfer = np.kron(b, a.conj().T @ parity)
    return transfer + transfer.conj().T if quadrature == "X" else -1j * (transfer - transfer.conj().T)


def transfer_kraus(site_a=0, site_b=0, *, quadrature="X"):
    """Independent finite-matrix positive Kraus maps for the transfer read.

    Q^3=Q, hence sqrt((I+sQ)/2)=I/sqrt(2)+sQ/2+(1/2-1/sqrt(2))*Q^2.
    These target matrices are not consulted by the actual source executor.
    """
    q = transfer_operator(site_a, site_b, quadrature=quadrature)
    common = np.eye(256)/math.sqrt(2)+(.5-1/math.sqrt(2))*(q@q)
    return common+q/2, common-q/2


def transfer_basis_word(site_a=0, site_b=0, *, quadrature="X"):
    """B with B* (n_A-n_B) B=the requested transfer quadrature.

    For i<j, X=(X_i Z...Z X_j+Y_i Z...Z Y_j)/2 and
    Y=(X_i Z...Z Y_j-Y_i Z...Z X_j)/2. Native pulses implement exp(i a P).
    Thus B_X=exp(+i*pi*Y/4) and B_Y=exp(-i*pi*X/4). The paired Pauli
    factors commute; their complete product is number conserving.
    """
    i, j = _sites(site_a, site_b, quadrature)
    if quadrature == "Z":
        return []
    requests = (("X", "Y", math.pi/8), ("Y", "X", -math.pi/8)) if quadrature == "X" else (
        ("X", "X", -math.pi/8), ("Y", "Y", -math.pi/8))
    word = []
    for left, right, angle in requests:
        label = ["I"] * 8
        label[i], label[j] = left, right
        label[i+1:j] = ["Z"] * (j-i-1)
        word += compiler.pauli_rotation(label, angle)
    return word


def transfer_read(state, site_a=0, site_b=0, *, quadrature="X"):
    """Execute E0=(I+quadrature)/2 and E1=(I-quadrature)/2 on retained work.

    The first eight bits are the same two fermion cells. Higher work bits
    (including the existing clocks) and every reference column are retained.
    The positive Kraus outputs are B* sqrt(E_D,b) B times the input; no
    selection, reset, re-preparation or probability normalization occurs.
    Q has eigenvalues -1,0,+1, so this is a binary POVM, not a projector
    read of Q's zero eigenspace. A fresh neutral pointer is a declared service.
    """
    _sites(site_a, site_b, quadrature)
    if np.ma.isMaskedArray(state):
        raise ValueError("masked retained states are inadmissible")
    raw = np.asarray(state)
    if (raw.ndim not in (1, 2) or raw.shape[0] < 256
            or raw.shape[0] > (1 << 15) or raw.shape[0] & (raw.shape[0]-1)):
        raise ValueError("eight fermion bits and at most seven spectator bits required")
    qubits = raw.shape[0].bit_length()-1
    basis = transfer_basis_word(site_a, site_b, quadrature=quadrature)
    rotated, _ = history.execute_word(basis, qubits, state)
    diagonal = history.diagonal_read(rotated, [.5, 1., 0., .5], targets=(site_a, 4+site_b))
    inverse = compiler.inverse(basis)
    outputs = tuple(history.execute_word(inverse, qubits, branch)[0] for branch in diagonal.outputs)
    word = [*basis, *diagonal.word, *inverse]
    resources = {**compiler.resources(word), "fresh_pointer_blanks": 1,
                 "committed_classical_bits": 1, "work_preparations": 0,
                 "work_qubits": qubits, "reference_columns": outputs[0].shape[1],
                 "endpoint_total_fermion_charge_preserved": True,
                 "charge_domain": "completed logical instrument, not every intermediate processor pulse",
                 "physical_service_error_certified": False}
    return history.Readout(outputs, tuple(word), resources)


def transfer_visibility(response, *, site_a=0, site_b=0):
    """CAR visibility of the SAME (2,2)<->(3,1) ground-arm interferometer.

    This is a matrix-element/Lehmann certificate. It does not reinterpret a
    final transfer expectation, after other instruments, as a pole fit.
    """
    if type(response) is not FiniteResponse:
        raise ValueError("a finite response required")
    _sites(site_a, site_b, "X")
    e1, g1 = model.sector_ground(response.hamiltonian, 1)
    e3, g3 = model.sector_ground(response.hamiltonian, 3)
    g2 = response.ground
    reference = np.kron(g2, g2)
    charged = np.kron(g1, g3)
    parity = np.diag((-1.) ** model.number_diagonal())
    ladders = model.annihilators()
    addition = np.array([np.vdot(g3, c.conj().T @ g2) for c in ladders])
    removal = np.array([np.vdot(g1, c @ g2) for c in ladders])
    element = addition[site_a] * removal[site_b]
    # Check the factorization against the globally graded, even operator.
    t = np.kron(ladders[site_b], ladders[site_a].conj().T @ parity)
    direct = np.vdot(charged, t @ reference)
    pair = float(e1+e3-2*response.ground_energy)
    addition_cost = float(e3-response.ground_energy)
    removal_cost = float(e1-response.ground_energy)
    symmetric_gap = (response.chemical_potential == 0.0
                     and math.isclose(addition_cost, removal_cost, rel_tol=0., abs_tol=1e-12))
    return {"reference_arm": "N_A=2,N_B=2", "charged_arm": "N_A=3,N_B=1",
            "total_charge": 4, "selected_sites": [site_a, site_b],
            "transfer_direction": "T=C_A* C_B raises A and lowers B",
            "global_occupation_U1_invariant": True,
            "local_electromagnetic_gauge_invariance_established": False,
            "gauge_link_or_dressing_supplied": False,
            "apparatus_pair_has_derived_spacetime_locality": False,
            "transfer_matrix_element": [float(element.real), float(element.imag)],
            "transfer_weight": float(abs(element)**2),
            "factorization_residual": float(abs(direct-element)),
            "all_site_pair_weights": (abs(addition[:, None]*removal[None, :])**2).tolist(),
            "pair_frequency": pair, "pair_frequency_over_clock": pair/model.parameters().clock,
            "particle_hole_gap_identification_admitted": symmetric_gap,
            "single_gap_over_clock_under_particle_hole_symmetry": (
                pair/(2*model.parameters().clock) if symmetric_gap else None),
            "addition_ground_cost": addition_cost,
            "removal_ground_cost": removal_cost,
            "X_expectation_on_unread_two_arm_state": "Re(m*exp(+i*pair_frequency*s))",
            "Y_expectation_on_unread_two_arm_state": "Im(m*exp(+i*pair_frequency*s))",
            "m_definition": "<g3_A,g1_B|T|g2_A,g2_B>",
            "coherent_odd_charge_preparation_required": False,
            "final_history_expectation_is_independent_pole_measurement": False}


def charged_sector_symmetry(response):
    """Simple charged sectors admit no nontrivial degenerate color triplet.

    The exact native ratios prove b>0 and a-b>0. For N=1 or N=3 the
    centered quartic term vanishes, leaving the four distinct energies
    (-a-b,-a+b,a-b,a+b)+mu*N. Every charge-preserving symmetry commuting
    with this H is diagonal in that sector's energy basis. Its unitary
    commutant is U(1)^4, irrespective of the basis used to name the modes.
    This is a finite elementary-pole statement, not an exclusion of gauge
    descriptions, composites or other continuing/refined architectures.
    """
    if type(response) is not FiniteResponse:
        raise ValueError("a finite response required")
    p = model.parameters()
    w0, w1, w2 = p.weights
    inter_ratio, imbalance_ratio = w2/w1, w1*w1/(w0*w2)
    if not inter_ratio > 1 or not imbalance_ratio > 1:
        raise ArithmeticError("the exact simple-spectrum branch changed")
    omega, imbalance = p.intra+p.inter, p.intra-p.inter
    bare = (-omega, -imbalance, imbalance, omega)
    ranks = {name: [int(np.count_nonzero(np.linalg.eigvalsh(pole.residue) > VISIBILITY_TOLERANCE))
                    for pole in poles]
             for name, poles in (("addition", response.addition), ("removal", response.removal))}
    return {"sector_numbers": [1, 3], "unshifted_energies": list(bare),
            "exact_energy_forms": ["-(a+b)", "-(a-b)", "a-b", "a+b"],
            "strict_order_proof": {"exp(b)": str(inter_ratio), "exp(a-b)": str(imbalance_ratio),
                                   "both_exact_ratios_exceed_one": True},
            "energy_eigenspace_multiplicities": [1, 1, 1, 1],
            "minimum_level_separation": min(bare[i+1]-bare[i] for i in range(3)),
            "charged_residue_ranks_at_visibility_tolerance": ranks,
            "complex_commutant_dimension_per_sector": 4,
            "unitary_commutant_per_sector": "U(1)^4",
            "nontrivial_energy_and_charge_preserving_color_triplet_by_relabeling": False,
            "scope": "the single ring's N=1 or N=3 elementary charged sectors only",
            "composite_or_refined_constructions_excluded": False}


def report(*, interaction=True, chemical_potential=0.0):
    """JSON-ready complete finite response certificate and its limits."""
    result = finite_response(interaction=interaction, chemical_potential=chemical_potential)
    h, g = result.hamiltonian, result.ground
    ladders = model.annihilators()
    integrated = sum((p.residue for p in (*result.addition, *result.removal)), np.zeros((4, 4), complex))
    first = sum((p.frequency*p.residue for p in (*result.addition, *result.removal)), np.zeros((4, 4), complex))
    moment = np.array([[np.vdot(g, ((a@h-h@a)@b.conj().T+b.conj().T@(a@h-h@a)) @ g)
                        for b in ladders] for a in ladders])
    centered = [a.conj().T@a-mean*np.eye(16) for a, mean in zip(ladders, result.density_means, strict=True)]
    covariance = np.array([[np.vdot(g, a@b@g) for b in centered] for a in centered])
    density_sum = sum((p.residue for p in result.density), np.zeros((4, 4), complex))
    visible = lambda poles: [p.energy_cost for p in poles if np.trace(p.residue).real > VISIBILITY_TOLERANCE]
    all_poles = (*result.addition, *result.removal, *result.density)
    n = np.diag(model.number_diagonal())
    return {"schema": "oph.common_source.particle_identification.finite_response.v1",
            "law": "same four-mode spinless interacting SSH ring as mass_response",
            "source_parameters": model.parameters().record(), "interaction_enabled": interaction,
            "chemical_potential": result.chemical_potential,
            "chemical_potential_convention": "H(mu)=H(0)+mu*N",
            "reference_number": 2, "reference_energy": result.ground_energy,
            "reference_is_global_ground": bool(np.linalg.eigvalsh(h).min() >= result.ground_energy-1e-12),
            "addition": [p.record() for p in result.addition],
            "removal": [p.record() for p in result.removal],
            "connected_density": [p.record() for p in result.density],
            "density_means": result.density_means.real.tolist(),
            "lowest_CAR_visible_addition_cost": min(visible(result.addition)),
            "lowest_CAR_visible_removal_cost": min(visible(result.removal)),
            "lowest_visible_neutral_density_cost": min(visible(result.density)),
            "sum_rules": {"CAR_integral_identity_residual": float(np.max(abs(integrated-np.eye(4)))),
                          "CAR_first_moment_residual": float(np.max(abs(first-moment))),
                          "density_covariance_residual": float(np.max(abs(density_sum-covariance))),
                          "total_number_density_response_residual": float(max(abs(p.residue.sum()) for p in result.density)),
                          "minimum_residue_eigenvalue": float(min(np.linalg.eigvalsh(p.residue).min() for p in all_poles)),
                          "charge_commutator_residual": float(np.max(abs(n@h-h@n)))},
            "neutral_transfer_apparatus": transfer_visibility(result),
            "finite_charged_sector_symmetry": charged_sector_symmetry(result),
            "visibility_tolerance": VISIBILITY_TOLERANCE,
            "degenerate_energy_grouping_tolerance": DEGENERACY_TOLERANCE,
            "all_poles_and_dark_residues_retained": True,
            "spectral_identity_is_exact_finite_dimensional": True,
            "numerical_arrays_are_outward_certified": False,
            "source_transfer_instrument_available": True,
            "CAR_statistics": "encoded fermionic canonical anticommutation relations",
            "charge_meaning": "integer occupation relative to N=2, not identified electromagnetic charge",
            "spatial_spin_representation": "spinless finite ring; no Lorentz spin representation derived",
            "infinite_volume_stable_particle_established": False,
            "Lorentzian_interacting_particle_pole_established": False,
            "Standard_Model_species_identified": False,
            "empirical_mass_inputs": []}
