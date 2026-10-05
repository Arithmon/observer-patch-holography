#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from scipy.linalg import eigh

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


def circuit_density_q0_order(circuit: QuantumCircuit) -> np.ndarray:
    from qiskit.quantum_info import DensityMatrix, Statevector

    rho = DensityMatrix(Statevector.from_instruction(circuit)).data
    return qiskit_density_to_q0_order(rho, circuit.num_qubits)


def partial_trace_q0_order(rho: np.ndarray, keep: list[int], num_qubits: int) -> np.ndarray:
    dims = [2] * num_qubits
    keep = sorted(keep)
    trace_out = [idx for idx in range(num_qubits) if idx not in keep]
    tensor = rho.reshape(*(dims + dims))
    current_n = num_qubits
    for idx in sorted(trace_out, reverse=True):
        tensor = np.trace(tensor, axis1=idx, axis2=idx + current_n)
        current_n -= 1
    final_dim = 2 ** len(keep)
    return tensor.reshape(final_dim, final_dim)


def von_neumann_entropy(rho: np.ndarray, base: float = 2.0) -> float:
    evals = np.linalg.eigvalsh((rho + rho.conj().T) / 2.0)
    evals = np.clip(np.real_if_close(evals), 0.0, None)
    total = float(np.sum(evals))
    if total <= 0:
        return 0.0
    evals = evals / total
    nonzero = evals[evals > 1e-12]
    if len(nonzero) == 0:
        return 0.0
    return float(-np.sum(nonzero * np.log(nonzero) / np.log(base)))


def conditional_mutual_information(rho: np.ndarray) -> float:
    s_ab = von_neumann_entropy(partial_trace_q0_order(rho, [0, 1], 3))
    s_bc = von_neumann_entropy(partial_trace_q0_order(rho, [1, 2], 3))
    s_b = von_neumann_entropy(partial_trace_q0_order(rho, [1], 3))
    s_abc = von_neumann_entropy(rho)
    return float(s_ab + s_bc - s_b - s_abc)


def project_to_physical_density_matrix(rho: np.ndarray) -> np.ndarray:
    herm = (rho + rho.conj().T) / 2.0
    evals, evecs = eigh(herm)
    evals = np.clip(np.real_if_close(evals), 0.0, None)
    total = float(np.sum(evals))
    if total <= 0:
        return np.eye(rho.shape[0], dtype=complex) / rho.shape[0]
    return (evecs @ np.diag(evals / total) @ evecs.conj().T).astype(complex)


def matrix_sqrt_psd(rho: np.ndarray) -> np.ndarray:
    evals, evecs = eigh((rho + rho.conj().T) / 2.0)
    evals = np.clip(np.real_if_close(evals), 0.0, None)
    return evecs @ np.diag(np.sqrt(evals)) @ evecs.conj().T


def matrix_inv_sqrt_psd(rho: np.ndarray, cutoff: float = 1e-10) -> np.ndarray:
    evals, evecs = eigh((rho + rho.conj().T) / 2.0)
    inv_sqrt = np.array([1.0 / np.sqrt(v) if v > cutoff else 0.0 for v in evals], dtype=float)
    return evecs @ np.diag(inv_sqrt) @ evecs.conj().T


def state_fidelity(rho: np.ndarray, sigma: np.ndarray) -> float:
    """Squared Uhlmann fidelity, ``||sqrt(rho) sqrt(sigma)||_1**2``."""
    sqrt_rho = matrix_sqrt_psd(project_to_physical_density_matrix(rho))
    inner = sqrt_rho @ project_to_physical_density_matrix(sigma) @ sqrt_rho
    evals = np.linalg.eigvalsh((inner + inner.conj().T) / 2.0)
    evals = np.clip(np.real_if_close(evals), 0.0, None)
    fidelity = np.sum(np.sqrt(evals))
    return float(np.real_if_close(fidelity * fidelity))


def trace_distance(rho: np.ndarray, sigma: np.ndarray) -> float:
    delta = (rho - sigma + (rho - sigma).conj().T) / 2.0
    evals = np.linalg.eigvalsh(delta)
    return float(0.5 * np.sum(np.abs(np.real_if_close(evals))))


def pauli_expectation(rho: np.ndarray, pauli_string_q0: str) -> float:
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


def petz_recovery(rho_abc: np.ndarray) -> np.ndarray:
    rho_ab = partial_trace_q0_order(rho_abc, [0, 1], 3)
    rho_bc = partial_trace_q0_order(rho_abc, [1, 2], 3)
    rho_b = partial_trace_q0_order(rho_abc, [1], 3)
    sqrt_bc = matrix_sqrt_psd(rho_bc)
    inv_sqrt_b = matrix_inv_sqrt_psd(rho_b)
    whitened_ab = np.kron(np.eye(2), inv_sqrt_b) @ rho_ab @ np.kron(np.eye(2), inv_sqrt_b)
    lifted = np.kron(whitened_ab, np.eye(2))
    embed_bc = np.kron(np.eye(2), sqrt_bc)
    recovered = embed_bc @ lifted @ embed_bc
    return project_to_physical_density_matrix(recovered)


def fawzi_renner_fidelity_lower_bound(cmi_bits: float) -> float:
    """Lower bound on optimal *squared* recovery fidelity.

    Fawzi--Renner (https://arxiv.org/abs/1410.0664, Eq. 3) uses root
    fidelity >= 2**(-I/2). Squaring matches ``state_fidelity`` here.
    The theorem guarantees a B -> BC recovery channel; it does not certify
    this benchmark's particular unrotated Petz map at nonzero CMI.
    """
    return float(2 ** (-max(cmi_bits, 0.0)))


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


def expectation_from_counts(counts: dict[str, int], pauli_q0: str) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    acc = 0.0
    for bitstring, count in counts.items():
        bits_q0 = bitstring_to_q0_order(bitstring)
        eigenvalue = 1.0
        for bit, char in zip(bits_q0, pauli_q0):
            if char != "I":
                eigenvalue *= 1.0 if bit == "0" else -1.0
        acc += eigenvalue * count
    return float(acc / total)


def reconstruct_density_matrix(
    counts_by_basis: dict[str, dict[str, int]],
    num_qubits: int,
) -> tuple[np.ndarray, dict[str, float]]:
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
        values = [expectation_from_counts(counts_by_basis[basis], label) for basis in support]
        expectations[label] = float(np.mean(values)) if values else 0.0

    rho = np.zeros((2**num_qubits, 2**num_qubits), dtype=complex)
    for label_q0, value in expectations.items():
        op = PAULI_MATRICES[label_q0[0]]
        for char in label_q0[1:]:
            op = np.kron(op, PAULI_MATRICES[char])
        rho += value * op
    rho /= 2**num_qubits
    return project_to_physical_density_matrix(rho), expectations


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
        rho = circuit_density_q0_order(circ)
        candidates.append(
            {
                "seed": seed,
                "circuit": circ,
                "exact_cmi_bits": conditional_mutual_information(rho),
            }
        )
    chosen = max(candidates, key=lambda item: item["exact_cmi_bits"])
    return chosen["circuit"], {
        "seed": chosen["seed"],
        "exact_cmi_bits": chosen["exact_cmi_bits"],
        "candidate_summary": [
            {"seed": item["seed"], "exact_cmi_bits": item["exact_cmi_bits"]} for item in candidates
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
    recovered = petz_recovery(rho)
    cmi_bits = conditional_mutual_information(rho)
    fidelity = state_fidelity(rho, recovered)
    return {
        "cmi_bits": cmi_bits,
        "fidelity_convention": "squared_uhlmann",
        "petz_fidelity": fidelity,
        "petz_trace_distance": trace_distance(rho, recovered),
        "petz_observable_mismatch": low_weight_observable_mismatch(rho, recovered),
        "fawzi_renner_fidelity_lower_bound": fawzi_renner_fidelity_lower_bound(cmi_bits),
        "fawzi_renner_bound_scope": "optimal_recovery_over_B_to_BC_channels",
    }


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

    exact_analysis = {}
    for circuit in circuits:
        rho = circuit_density_q0_order(circuit)
        exact_analysis[circuit.name] = analyze_state(rho)

    sampler_output, resolved_backend = run_sampler(
        circuits=measured,
        mode=mode,
        shots=args.shots,
        transpile_seed=args.transpile_seed,
        credentials_file=args.credentials_file,
        backend_name=args.backend,
    )

    counts_by_state = {}
    flat_counts = sampler_output["counts_by_name"]
    for state_name, mapping in measured_index.items():
        counts_by_state[state_name] = {
            basis: flat_counts[circuit_name] for basis, circuit_name in mapping.items()
        }

    reconstructed_analysis = {}
    for circuit in circuits:
        rho_recon, expectations = reconstruct_density_matrix(counts_by_state[circuit.name], 3)
        reconstructed_analysis[circuit.name] = analyze_state(rho_recon)
        reconstructed_analysis[circuit.name]["selected_expectations"] = {
            label: expectations[label]
            for label in ["ZZI", "IZZ", "ZIZ", "XXX", "YYY", "ZZZ"]
            if label in expectations
        }

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "experiment": "stage1_markov_fingerprint",
        "mode": mode,
        "backend": resolved_backend,
        "shots": args.shots,
        "transpile_seed": args.transpile_seed,
        "random_depth": args.random_depth,
        "catalog": catalog_meta,
        "run_metadata": sampler_output["run_metadata"],
        "exact_analysis": exact_analysis,
        "reconstructed_analysis": reconstructed_analysis,
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
                reconstructed_analysis["structured_theta_0.00"]["petz_fidelity"]
                >= reconstructed_analysis["structured_theta_0.60"]["petz_fidelity"]
                >= reconstructed_analysis["structured_theta_1.00"]["petz_fidelity"]
            ),
        },
    }

    write_json(outdir / "summary.json", summary)
    (outdir / "summary_pretty.txt").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
