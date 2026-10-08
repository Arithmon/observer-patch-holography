#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.linalg import eigh

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from quantum_information.gibbs import _numeric
from quantum_information.recovery import (
    matrix_inv_sqrt_psd, matrix_sqrt_psd, petz_recovery, state_fidelity,
    trace_distance, validated_state, UnresolvedPetzSupport,
)
from quantum_information.positive_rounding import _is_exact_psd, round_positive_gram
from quantum_information.states import (
    ATOL, conditional_mutual_information as _cmi, partial_trace,
    von_neumann_entropy as _entropy,
)

# Numerical recovery diagnostics do not require the optional circuit/cloud SDKs.
if TYPE_CHECKING:
    from qiskit import QuantumCircuit


PAULI_MATRICES = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def qiskit_density_to_q0_order(rho: np.ndarray, num_qubits: int) -> np.ndarray:
    dims = [2] * num_qubits
    tensor = rho.reshape(*(dims + dims))
    ket_axes = list(range(num_qubits))
    bra_axes = list(range(num_qubits, 2 * num_qubits))
    perm = list(reversed(ket_axes)) + list(reversed(bra_axes))
    return np.transpose(tensor, perm).reshape(2**num_qubits, 2**num_qubits)


def circuit_density_q0_order(circuit: QuantumCircuit, *, return_diagnostics=False):
    """Construct an exactly PSD matrix from the numerical circuit amplitudes.

    The Gram-rounding bound covers these supplied amplitudes, not simulator
    error relative to an ideal circuit. No tomography projection is applied.
    """
    from qiskit.quantum_info import Statevector

    _qubit_count(circuit.num_qubits)
    if type(return_diagnostics) is not bool:
        raise ValueError("return_diagnostics must be Boolean")
    amplitudes = _numeric(Statevector.from_instruction(circuit).data, "circuit amplitudes")
    if amplitudes.shape != (2**circuit.num_qubits,):
        raise ValueError("statevector dimensions must match the circuit")
    # Reverse tensor axes before forming the Gram matrix. Constructing an
    # ordinary floating outer product first can lose exact positivity.
    q0 = amplitudes.reshape([2]*circuit.num_qubits).transpose().reshape(-1, 1)
    rho, bound = round_positive_gram(q0)
    rho = validated_state(rho)
    diagnostic = {
        "state_construction": "outward_rounded_statevector_gram",
        "gram_rounding_trace_bound": bound,
        "rounding_bound_scope": "supplied_numerical_statevector",
        "trace_normalization_tolerance": ATOL,
        "circuit_simulation_error_certified": False,
    }
    return (rho, diagnostic) if return_diagnostics else rho


def partial_trace_q0_order(rho: np.ndarray, keep: list[int], num_qubits: int) -> np.ndarray:
    _qubit_count(num_qubits)
    return partial_trace(validated_state(rho), [2]*num_qubits, keep)


def von_neumann_entropy(rho: np.ndarray, base: float = 2.0) -> float:
    base = _numeric(base, "entropy base", real=True)
    if base.ndim != 0 or base <= 1:
        raise ValueError("entropy base must be a finite real number greater than one")
    return _entropy(validated_state(rho))/math.log(float(base))


def conditional_mutual_information(rho: np.ndarray) -> float:
    return _cmi(validated_state(rho), [2, 2, 2], [0], [1], [2])/math.log(2)


def project_to_physical_density_matrix(rho: np.ndarray, *, return_rounding_bound=False):
    """Explicit tomography estimator: normalized positive spectral part.

    This is not validation, a channel, or the nearest trace-one PSD matrix.
    Only tomography calls it; recovery and distances never repair inputs.
    Outward Gram rounding certifies PSD of the returned entries when ordinary
    reconstruction loses it. Trace normalization retains the shared ATOL
    convention; the rounding bound is not a statistical error certificate.
    """
    rho = _numeric(rho, "tomographic estimate")
    if type(return_rounding_bound) is not bool:
        raise ValueError("return_rounding_bound must be Boolean")
    if (rho.ndim != 2 or not len(rho) or rho.shape[0] != rho.shape[1]
            or np.linalg.norm(rho-rho.conj().T) > ATOL
            or abs(np.trace(rho)-1) > ATOL):
        raise ValueError("a Hermitian trace-one tomographic estimate is required")
    herm = (rho + rho.conj().T) / 2.0
    evals, evecs = eigh(herm)
    evals = np.clip(np.real_if_close(evals), 0.0, None)
    total = float(np.sum(evals))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("tomographic positive part exceeds numerical range")
    weights = evals/total
    state = validated_state((evecs*weights) @ evecs.conj().T)
    rounding_bound = None  # No Gram-rounding certificate was needed/emitted.
    if not _is_exact_psd(state):
        state, rounding_bound = round_positive_gram(evecs*np.sqrt(weights))
        state = validated_state(state)
    return (state, rounding_bound) if return_rounding_bound else state


def pauli_expectation(rho: np.ndarray, pauli_string_q0: str) -> float:
    rho = validated_state(rho)
    if (not isinstance(pauli_string_q0, str) or not pauli_string_q0
            or any(p not in "IXYZ" for p in pauli_string_q0)
            or rho.shape != (2**len(pauli_string_q0),)*2):
        raise ValueError("Pauli label must match the state space")
    op = PAULI_MATRICES[pauli_string_q0[0]]
    for char in pauli_string_q0[1:]:
        op = np.kron(op, PAULI_MATRICES[char])
    return float(np.real_if_close(np.trace(rho @ op)))


def low_weight_observable_mismatch(rho: np.ndarray, sigma: np.ndarray) -> float:
    labels = []
    paulis = "XYZ"
    for weight in (1, 2):
        for positions in itertools.combinations(range(3), weight):
            for ops in itertools.product(paulis, repeat=weight):
                label = ["I", "I", "I"]
                for pos, op in zip(positions, ops):
                    label[pos] = op
                labels.append("".join(label))
    diffs = [abs(pauli_expectation(rho, label) - pauli_expectation(sigma, label)) for label in labels]
    return float(np.mean(diffs))


def fawzi_renner_fidelity_lower_bound(cmi_bits: float) -> float:
    """Lower bound on optimal *squared* recovery fidelity.

    Fawzi--Renner (https://arxiv.org/abs/1410.0664, Eq. 3) uses root
    fidelity >= 2**(-I/2). Squaring matches ``state_fidelity`` here.
    The theorem guarantees a B -> BC recovery channel; it does not certify
    this benchmark's particular unrotated Petz map at nonzero CMI.
    """
    cmi_bits = _numeric(cmi_bits, "CMI in bits", real=True)
    if cmi_bits.ndim != 0 or cmi_bits < -ATOL:
        raise ValueError("CMI must be nonnegative, apart from declared spectral roundoff")
    return float(2 ** (-max(float(cmi_bits), 0.0)))


def basis_rotation(circuit: QuantumCircuit, qubit: int, basis: str) -> None:
    if basis == "X":
        circuit.h(qubit)
    elif basis == "Y":
        circuit.sdg(qubit)
        circuit.h(qubit)
    elif basis == "Z":
        return
    else:
        raise ValueError(f"Unsupported basis {basis}")


def measurement_bases(num_qubits: int) -> list[str]:
    _qubit_count(num_qubits)
    return ["".join(chars) for chars in itertools.product("XYZ", repeat=num_qubits)]


def add_measurement_basis(circuit: QuantumCircuit, basis_q0: str) -> QuantumCircuit:
    from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

    qreg = QuantumRegister(circuit.num_qubits, "q")
    creg = ClassicalRegister(circuit.num_qubits, "c")
    measured = QuantumCircuit(qreg, creg, name=f"{circuit.name}__{basis_q0}")
    measured.compose(circuit, qubits=qreg, inplace=True)
    for qubit, basis in enumerate(basis_q0):
        basis_rotation(measured, qubit, basis)
    measured.measure(qreg, creg)
    return measured


def bitstring_to_q0_order(bitstring: str) -> str:
    return bitstring[::-1]


def _qubit_count(value):
    if type(value) is not int or value <= 0:
        raise ValueError("a positive integer qubit count is required")


def _count_totals(counts, pauli_q0):
    if (not isinstance(pauli_q0, str) or not pauli_q0
            or any(p not in "IXYZ" for p in pauli_q0)
            or not isinstance(counts, dict) or not counts):
        raise ValueError("a Pauli label and nonempty counts are required")
    acc = total = 0
    for bitstring, count in counts.items():
        if (not isinstance(bitstring, str) or len(bitstring) != len(pauli_q0)
                or any(bit not in "01" for bit in bitstring)
                or type(count) is not int or count < 0):
            raise ValueError("matching binary outcomes and nonnegative integer counts required")
        total += count
        bits_q0 = bitstring_to_q0_order(bitstring)
        eigenvalue = 1
        for bit, char in zip(bits_q0, pauli_q0):
            if char != "I":
                eigenvalue *= 1 if bit == "0" else -1
        acc += eigenvalue * count
    if total <= 0:
        raise ValueError("each measurement setting requires positive shots")
    return acc, total


def expectation_from_counts(counts: dict[str, int], pauli_q0: str) -> float:
    acc, total = _count_totals(counts, pauli_q0)
    return acc/total


def reconstruct_density_matrix(
    counts_by_basis: dict[str, dict[str, int]],
    num_qubits: int,
    *, return_diagnostics: bool = False,
):
    """Complete local-Pauli tomography with pooled shots and explicit repair."""
    required = measurement_bases(num_qubits)
    if not isinstance(counts_by_basis, dict) or set(counts_by_basis) != set(required):
        raise ValueError("complete local-Pauli measurement settings are required")
    if type(return_diagnostics) is not bool:
        raise ValueError("return_diagnostics must be Boolean")
    shots = {basis: _count_totals(counts, "I"*num_qubits)[1]
             for basis, counts in counts_by_basis.items()}
    expectations: dict[str, float] = {"I" * num_qubits: 1.0}
    for pauli_q0 in itertools.product("IXYZ", repeat=num_qubits):
        label = "".join(pauli_q0)
        if label == "I" * num_qubits:
            continue
        support = [
            basis
            for basis in counts_by_basis
            if all(p == "I" or p == b for p, b in zip(label, basis))
        ]
        totals = [_count_totals(counts_by_basis[basis], label) for basis in support]
        expectations[label] = sum(t[0] for t in totals)/sum(t[1] for t in totals)

    rho = np.zeros((2**num_qubits, 2**num_qubits), dtype=complex)
    for label_q0, value in expectations.items():
        op = PAULI_MATRICES[label_q0[0]]
        for char in label_q0[1:]:
            op = np.kron(op, PAULI_MATRICES[char])
        rho += value * op
    rho /= 2**num_qubits
    state, rounding_bound = project_to_physical_density_matrix(rho, return_rounding_bound=True)
    diagnostics = {
        "estimator": "pooled_pauli_linear_inversion_then_normalized_positive_part_with_psd_rounding",
        "gram_rounding_trace_bound": rounding_bound,
        "trace_normalization_tolerance": ATOL,
        "complete_settings": len(required),
        "shots_by_basis": shots,
        "raw_min_eigenvalue": float(np.linalg.eigvalsh(rho)[0]),
        "negative_spectral_mass": float(-sum(v for v in np.linalg.eigvalsh(rho) if v < 0)),
        "state_correction_frobenius": float(np.linalg.norm(state-rho)),
        "statistical_error_certified": False,
    }
    return (state, expectations, diagnostics) if return_diagnostics else (state, expectations)


def build_structured_family(theta: float) -> QuantumCircuit:
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(3, name=f"structured_theta_{theta:.2f}")
    qc.h(1)
    qc.cx(1, 0)
    if abs(theta) > 1e-12:
        qc.cry(theta, 1, 2)
    return qc


def build_ghz() -> QuantumCircuit:
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(3, name="ghz_control")
    qc.h(0)
    qc.cx(0, 1)
    qc.cx(1, 2)
    return qc


def build_random_control(seed: int, depth: int) -> QuantumCircuit:
    from qiskit import QuantumCircuit
    from qiskit.circuit.random import random_circuit

    random_qc = random_circuit(3, depth=depth, max_operands=2, measure=False, seed=seed)
    qc = QuantumCircuit(3, name=f"random_seed_{seed}")
    qc.compose(random_qc, inplace=True)
    return qc


def choose_random_control(depth: int, seeds: list[int]) -> tuple[QuantumCircuit, dict[str, float]]:
    candidates = []
    for seed in seeds:
        circ = build_random_control(seed, depth)
        rho, construction = circuit_density_q0_order(circ, return_diagnostics=True)
        candidates.append(
            {
                "seed": seed,
                "circuit": circ,
                "exact_cmi_bits": conditional_mutual_information(rho),
                "circuit_state": construction,
            }
        )
    chosen = max(candidates, key=lambda item: item["exact_cmi_bits"])
    return chosen["circuit"], {
        "seed": chosen["seed"],
        "exact_cmi_bits": chosen["exact_cmi_bits"],
        "circuit_state": chosen["circuit_state"],
        "candidate_summary": [
            {"seed": item["seed"], "exact_cmi_bits": item["exact_cmi_bits"],
             "circuit_state": item["circuit_state"]} for item in candidates
        ],
    }


def state_catalog(random_depth: int, random_seeds: list[int]) -> tuple[list[QuantumCircuit], dict]:
    random_circ, random_meta = choose_random_control(random_depth, random_seeds)
    circuits = [
        build_structured_family(0.0),
        build_structured_family(0.6),
        build_structured_family(1.0),
        build_ghz(),
        random_circ,
    ]
    return circuits, {"random_control_selection": random_meta}


def analyze_state(rho: np.ndarray) -> dict:
    """Retain valid information even when Petz support cannot be resolved."""
    cmi_bits = conditional_mutual_information(rho)
    result = {
        "cmi_bits": cmi_bits,
        "recovery_map": "unrotated_petz_with_reference_state_kernel_completion",
        "input_state_policy": "validated_without_normalization",
        "fidelity_convention": "squared_uhlmann",
        "fawzi_renner_fidelity_lower_bound": fawzi_renner_fidelity_lower_bound(cmi_bits),
        "fawzi_renner_bound_scope": "optimal_recovery_over_B_to_BC_channels",
    }
    try:
        recovered = petz_recovery(rho)
    except UnresolvedPetzSupport as error:
        result.update(petz_status="unresolved_support", petz_unavailable_reason=str(error),
                      petz_fidelity=None, petz_trace_distance=None, petz_observable_mismatch=None)
    else:
        result.update(petz_status="available", petz_unavailable_reason=None,
                      petz_fidelity=state_fidelity(rho, recovered),
                      petz_trace_distance=trace_distance(rho, recovered),
                      petz_observable_mismatch=low_weight_observable_mismatch(rho, recovered))
    return result


def run_sampler(
    circuits: list[QuantumCircuit],
    mode: str,
    shots: int,
    transpile_seed: int,
    credentials_file: Path,
    backend_name: str | None,
) -> tuple[dict, str | None]:
    from qiskit import transpile
    from qiskit_aer import AerSimulator
    from qiskit_ibm_runtime import SamplerV2

    from ibm_runtime_common import get_service

    if mode == "local":
        backend = AerSimulator()
        service = None
    else:
        service = get_service(credentials_file)
        if backend_name:
            backend = service.backend(backend_name)
        else:
            backend = service.least_busy(operational=True, simulator=False, min_num_qubits=3)
            backend_name = backend.name

    isa_circuits = transpile(
        circuits,
        backend=backend,
        optimization_level=1,
        seed_transpiler=transpile_seed,
    )
    sampler = SamplerV2(mode=backend)
    job = sampler.run(isa_circuits, shots=shots)
    result = job.result()
    counts_by_name = {}
    for circuit, pub_result in zip(circuits, result):
        keys = list(pub_result.data.keys())
        if len(keys) != 1:
            raise RuntimeError(f"Unexpected classical data keys: {keys}")
        bit_array = getattr(pub_result.data, keys[0])
        counts_by_name[circuit.name] = bit_array.get_counts()
    metadata = {
        "backend_name": getattr(backend, "name", str(backend)),
        "job_id": None if mode == "local" else job.job_id(),
        "mode": mode,
        "shots": shots,
        "transpile_seed": transpile_seed,
    }
    if service is not None:
        metadata["active_instance"] = service.active_instance()
    return {"counts_by_name": counts_by_name, "run_metadata": metadata}, backend_name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Stage 1A / 1B Markov fingerprint and recovery-map benchmark."
    )
    parser.add_argument(
        "--mode",
        choices=["local", "hardware"],
        default="local",
        help="Whether to run on local Aer or IBM hardware.",
    )
    parser.add_argument(
        "--local-testing",
        action="store_true",
        help="Alias for --mode local.",
    )
    parser.add_argument("--shots", type=int, default=512)
    parser.add_argument("--transpile-seed", type=int, default=7)
    parser.add_argument("--random-depth", type=int, default=3)
    parser.add_argument(
        "--random-seeds",
        type=int,
        nargs="+",
        default=list(range(10)),
    )
    parser.add_argument(
        "--credentials-file",
        type=Path,
        default=Path("IBM_cloud.txt"),
    )
    parser.add_argument("--backend", type=str, default=None)
    parser.add_argument(
        "--outdir",
        type=Path,
        required=True,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    from ibm_runtime_common import ensure_dir, write_json

    mode = "local" if args.local_testing else args.mode
    outdir = ensure_dir(args.outdir)
    # A failed rerun must not overwrite counts beside an older success report.
    for name in ("acquired_counts.json", "summary.json", "summary_pretty.txt"):
        if (outdir / name).exists():
            raise FileExistsError(f"Stage 1 output already exists at {outdir / name}; choose a fresh --outdir")

    circuits, catalog_meta = state_catalog(args.random_depth, args.random_seeds)
    bases = measurement_bases(3)
    measured = []
    measured_index = {}
    for circuit in circuits:
        measured_index[circuit.name] = {}
        for basis in bases:
            full = add_measurement_basis(circuit, basis)
            measured_index[circuit.name][basis] = full.name
            measured.append(full)

    sampler_output, resolved_backend = run_sampler(
        circuits=measured,
        mode=mode,
        shots=args.shots,
        transpile_seed=args.transpile_seed,
        credentials_file=args.credentials_file,
        backend_name=args.backend,
    )

    run_context = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "experiment": "stage1_markov_fingerprint",
        "mode": mode,
        "backend": resolved_backend,
        "shots": args.shots,
        "transpile_seed": args.transpile_seed,
        "random_depth": args.random_depth,
        "random_seeds": args.random_seeds,
        "catalog": catalog_meta,
        "run_metadata": sampler_output["run_metadata"],
    }
    # Persist the returned evidence before regrouping counts or evaluating
    # either reference or sampled states. Unexpected analysis errors still
    # raise, but cannot erase the acquired data needed to replay the failure.
    write_json(outdir / "acquired_counts.json", {
        **run_context,
        "measured_index": measured_index,
        "counts_by_name": sampler_output["counts_by_name"],
    })

    exact_analysis = {}
    for circuit in circuits:
        rho, construction = circuit_density_q0_order(circuit, return_diagnostics=True)
        exact_analysis[circuit.name] = analyze_state(rho)
        exact_analysis[circuit.name]["circuit_state"] = construction

    counts_by_state = {}
    flat_counts = sampler_output["counts_by_name"]
    for state_name, mapping in measured_index.items():
        counts_by_state[state_name] = {
            basis: flat_counts[circuit_name] for basis, circuit_name in mapping.items()
        }

    reconstructed_analysis = {}
    for circuit in circuits:
        rho_recon, expectations, tomography = reconstruct_density_matrix(
            counts_by_state[circuit.name], 3, return_diagnostics=True)
        reconstructed_analysis[circuit.name] = analyze_state(rho_recon)
        reconstructed_analysis[circuit.name]["tomography"] = tomography
        reconstructed_analysis[circuit.name]["selected_expectations"] = {
            label: expectations[label]
            for label in ["ZZI", "IZZ", "ZIZ", "XXX", "YYY", "ZZZ"]
            if label in expectations
        }

    structured_fidelities = [reconstructed_analysis[f"structured_theta_{theta}"]["petz_fidelity"]
                            for theta in ("0.00", "0.60", "1.00")]
    summary = {
        **run_context,
        "exact_analysis": exact_analysis,
        "reconstructed_analysis": reconstructed_analysis,
        "tomography_counts_by_state": counts_by_state,
        "fingerprint_checks": {
            "structured_theta_0.00_lt_random_control": reconstructed_analysis["structured_theta_0.00"][
                "cmi_bits"
            ]
            < reconstructed_analysis[f"random_seed_{catalog_meta['random_control_selection']['seed']}"][
                "cmi_bits"
            ],
            "structured_theta_0.00_lt_ghz": reconstructed_analysis["structured_theta_0.00"]["cmi_bits"]
            < reconstructed_analysis["ghz_control"]["cmi_bits"],
            "recovery_improves_as_cmi_drops": (
                None if any(fidelity is None for fidelity in structured_fidelities)
                else structured_fidelities[0] >= structured_fidelities[1] >= structured_fidelities[2]
            ),
        },
    }

    write_json(outdir / "summary.json", summary)
    (outdir / "summary_pretty.txt").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
