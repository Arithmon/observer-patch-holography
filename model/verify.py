"""Verify artifact bindings and replay the runnable model receipts."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from .replay_empirical import EMPIRICAL_REPLAYS

ROOT = Path(__file__).resolve().parents[1]
PACKETS = ("model/observer", "model/observer/p_response", "model/local_tick",
           "model/observer/scalar_response", "model/observer/source_scalar",
           "model/predictions", "model/observer/source_discrimination",
           "model/common_source", "model/spatial_interaction", "model/relativistic_scalar",
           "model/scalar_observations", "model/observer_reconstruction",
           "model/interacting_reconstruction", "model/coupled_reconstruction",
           "model/anharmonic_observations", "model/anharmonic_transfer",
           "model/audit_20261015", "model/occupation_control",
           "model/identification_design", "model/primitive_observer",
           "model/native_quantum_history", "model/native_clock",
           "model/native_fluctuations", "model/empirical_diagnostics",
           "model/native_relational_readout", "model/native_circulation",
           "model/native_interface", "model/measurement_bounds")
REPLAYS = (
    ("model.native_circulation", None, "model/native_circulation/results.json"),
    ("model.native_interface", None, "model/native_interface/results.json"),
    ("model.measurement_bounds", None, "model/measurement_bounds/results.json"),
    ("model.native_relational_readout", None, "model/native_relational_readout/results.json"),
    ("model.native_quantum_history", None, "model/native_quantum_history/results.json"),
    ("model.native_clock", None, "model/native_clock/results.json"),
    ("model.native_fluctuations", None, "model/native_fluctuations/results.json"),
    ("model.empirical_diagnostics", None, "model/empirical_diagnostics/results.json"),
    ("model.primitive_observer.transfer", None, "model/primitive_observer/transfer_results.json"),
    ("model.common_source.open_observer", None, "model/common_source/open_observer_results.json"),
    ("model.occupation_control", None, "model/occupation_control/results.json"),
    ("model.identification_design", None, "model/identification_design/results.json"),
    ("model.primitive_observer", None, "model/primitive_observer/results.json"),
    ("model.anharmonic_transfer", None, "model/anharmonic_transfer/results.json"),
    ("model.anharmonic_observations", None, "model/anharmonic_observations/results.json"),
    ("model.coupled_reconstruction", None, "model/coupled_reconstruction/results.json"),
    ("model.interacting_reconstruction", None, "model/interacting_reconstruction/results.json"),
    ("model.observer_reconstruction", None, "model/observer_reconstruction/results.json"),
    ("model.scalar_observations", None, "model/scalar_observations/results.json"),
    ("model.scalar_observations.occupation", None, "model/scalar_observations/occupation_results.json"),
    ("model.scalar_observations.experimental_readout", None,
     "model/scalar_observations/experimental_readout/results.json"),
    ("model.scalar_observations.experimental_readout.executable", None,
     "model/scalar_observations/experimental_readout/executable_results.json"),
    ("model.relativistic_scalar", None, "model/relativistic_scalar/results.json"),
    ("model.spatial_interaction", None, "model/spatial_interaction/results.json"),
    ("model.observer.source_selection", None, "model/observer/results/source_selection.json"),
    ("model.observer.coherent_continuation", None, "model/observer/results/coherent_continuation.json"),
    ("model.observer", "native", "model/observer/results/native.json"),
    ("model.observer.quantum_mixing_receipt", None, "model/observer/results/quantum_mixing.json"),
    ("model.observer", "capacity", "model/observer/results/capacity.json"),
    ("model.observer", "run", "model/observer/results/run.json"),
    ("model.observer.p_response", None, "model/observer/p_response/results.json"),
    ("model.local_tick", None, "model/local_tick/results.json"),
    ("model.observer.scalar_response", None, "model/observer/scalar_response/results.json"),
    ("model.observer.scalar_response", "--continue-field", "model/observer/scalar_response/continuation_results.json"),
    ("model.observer.scalar_response.local_evolution", None,
     "model/observer/scalar_response/local_results.json"),
    ("model.observer.source_scalar", None, "model/observer/source_scalar/results.json"),
    ("model.predictions", None, "model/predictions/results.json"),
    ("model.observer.source_discrimination", None,
     "model/observer/source_discrimination/results.json"),
    ("model.common_source", None, "model/common_source/results.json"),
    ("model.common_source.composition", None, "model/common_source/composition_results.json"),
    ("model.common_source.smooth_history", None, "model/common_source/smooth_history_results.json"),
    ("model.common_source.quantum_native", None, "model/common_source/quantum_native_results.json"),
    ("model.common_source.encoded_observer", None, "model/common_source/encoded_observer_results.json"),
    ("model.common_source.infinite_history", None, "model/common_source/infinite_results.json"),
    ("model.common_source.sharp_output", None, "model/common_source/sharp_output_results.json"),
    ("model.common_source.bridge", None, "model/common_source/bridge_results.json"),
    ("model.common_source.smooth_bridge", None, "model/common_source/smooth_results.json"),
    ("model.common_source.thermal_bridge", None, "model/common_source/thermal_results.json"),
    ("model.common_source.causality", None, "model/common_source/causality_results.json"),
    ("model.common_source.native_fields", None, "model/common_source/native_fields_results.json"),
    ("model.observer.record_recovery_receipt", None, "model/observer/results/record_recovery.json"),
)


def replay_groups(replays=None):
    """Partition the declaration without splitting the empirical cache worker."""
    replays = REPLAYS if replays is None else tuple(replays)
    if len(set(replays)) != len(replays):
        raise ValueError("Replay declarations must be unique")
    groups, index, empirical_groups = [], 0, 0
    while index < len(replays):
        if replays[index:index+len(EMPIRICAL_REPLAYS)] == EMPIRICAL_REPLAYS:
            groups.append(EMPIRICAL_REPLAYS)
            empirical_groups += 1
            index += len(EMPIRICAL_REPLAYS)
        else:
            if replays[index] in EMPIRICAL_REPLAYS:
                raise ValueError("Empirical replays must remain contiguous and in their canonical order")
            groups.append((replays[index],))
            index += 1
    if empirical_groups != 1:
        raise ValueError("Exactly one complete empirical replay group is required")
    return tuple(groups)


def selected_replay_groups(shard_index=None, shard_count=None, *, replays=None):
    """Keep the default order; isolate expensive shared-cache work when sharded.

    Shard zero owns the empirical group. Other groups retain declaration order
    and are distributed round-robin over the remaining shards. This prevents
    another costly certificate from joining the dominant empirical worker.
    Every new declared replay is included automatically.
    """
    groups = replay_groups(replays)
    if shard_index is None and shard_count is None:
        return groups
    if (type(shard_index) is not int or type(shard_count) is not int
            or not 1 <= shard_count <= len(groups)
            or not 0 <= shard_index < shard_count):
        raise ValueError(f"Supply both shard arguments: 1 <= count <= {len(groups)}, 0 <= index < count")
    if shard_count == 1:
        return groups
    if shard_index == 0:
        return (EMPIRICAL_REPLAYS,)
    remaining = tuple(group for group in groups if group != EMPIRICAL_REPLAYS)
    return tuple(group for index, group in enumerate(remaining)
                 if index % (shard_count-1) == shard_index-1)


def verify_artifacts():
    covered = set()
    for name in PACKETS:
        folder = ROOT / name
        manifest = json.loads((folder / "manifest.json").read_text())
        files = manifest["files"]
        if not files:
            raise ValueError(f"Empty artifact manifest: {name}")
        for relative, expected in files.items():
            target = (folder / relative).resolve()
            if not target.is_relative_to(folder.resolve()):
                raise ValueError(f"Artifact outside its package: {name}/{relative}")
            actual = hashlib.sha256(target.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f"Artifact digest mismatch: {name}/{relative}")
            covered.add(target)
    for name in PACKETS:
        for target in (ROOT / name).rglob("*"):
            if (target.is_file() and target.suffix in {".py", ".md", ".json"}
                    and target.name != "manifest.json"
                    and target.resolve() not in covered):
                raise ValueError(f"Unbound artifact: {target.relative_to(ROOT)}")
    return len(covered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hashes-only", action="store_true",
                        help="check artifact integrity without recalculation")
    parser.add_argument("--shard-index", type=int,
                        help="zero-based replay shard; requires --shard-count")
    parser.add_argument("--shard-count", type=int,
                        help="number of replay shards; every shard checks all artifact hashes")
    args = parser.parse_args()
    try:
        groups = selected_replay_groups(args.shard_index, args.shard_count)
    except ValueError as error:
        parser.error(str(error))
    count = verify_artifacts()
    print(f"Verified {count} model artifacts", flush=True)
    if args.hashes_only:
        return
    environment = dict(os.environ)
    environment.setdefault("OPENBLAS_NUM_THREADS", "1")
    environment.setdefault("OMP_NUM_THREADS", "1")
    selected_count = sum(map(len, groups))
    if args.shard_count is not None:
        print(f"Replaying shard {args.shard_index+1}/{args.shard_count}: "
              f"{selected_count} receipts", flush=True)
    for group in groups:
        if group == EMPIRICAL_REPLAYS:
            # Keep every canonical check, sharing only the expensive source
            # compilation caches in this one child. Failures still propagate.
            subprocess.run([sys.executable, "-m", "model.replay_empirical"],
                           cwd=ROOT, env=environment, check=True)
            continue
        (module, command, receipt), = group
        argv = [sys.executable, "-m", module]
        if command is not None:
            argv.append(command)
        argv.extend(["--check", receipt])
        print(f"Replaying {receipt}", flush=True)
        subprocess.run(argv, cwd=ROOT, env=environment, check=True)
    if args.shard_count is None:
        print(f"All {selected_count} model receipts reproduced", flush=True)
    else:
        print(f"All {selected_count} receipts in shard {args.shard_index+1}/{args.shard_count} reproduced",
              flush=True)


if __name__ == "__main__":
    main()
