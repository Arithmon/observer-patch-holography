"""Execute the pinned simulator and compare independently derived read semantics."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

from . import simulator_verifier as check
from .verify import HERE, verify


def canonicalize_source_aliases(packet, source, simulator):
    """Remove only proven duplicate-case names for the same pinned file.

    On case-insensitive macOS volumes, PosixPath treats import names such as
    PatchNet and patchnet as distinct even though both open the same file.
    Every canonical pin must exist; separate files and unknown imports fail.
    """
    check.equal(packet["sha256"], hashlib.sha256(check.canonical(
        {k: v for k, v in packet.items() if k != "sha256"})).hexdigest(), "live packet digest")
    actual, expected = packet["source"]["files"], source["files"]
    check.equal({name: actual[name] for name in expected if name in actual},
                expected, "live canonical source closure")
    aliases = set(actual)-set(expected)
    for alias in aliases:
        matches = [name for name in expected if name.casefold() == alias.casefold()]
        if len(matches) != 1:
            raise ValueError("unrecognized live source alias: "+alias)
        name = matches[0]
        canonical_path, alias_path = simulator/name, simulator/alias
        if (not alias_path.samefile(canonical_path)
                or alias_path.read_bytes() != canonical_path.read_bytes()
                or hashlib.sha256(alias_path.read_bytes()).hexdigest() != expected[name]
                or actual[alias] != expected[name]):
            raise ValueError("live source alias is not the same pinned file: "+alias)
    check.equal({**packet["source"], "files": expected}, source, "live source identity")
    if not aliases:
        return packet
    normalized = {**packet, "source": {**packet["source"], "files": dict(expected)}}
    normalized["sha256"] = hashlib.sha256(check.canonical(
        {k: v for k, v in normalized.items() if k != "sha256"})).hexdigest()
    return normalized


def live(simulator):
    simulator = Path(simulator).resolve()
    source = check.strict_load(HERE/"source_snapshot.json")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=simulator).decode().strip()
    check.equal(revision, source["revision"], "live simulator revision")
    for name, digest in source["files"].items():
        raw = (simulator/name).read_bytes()
        committed = subprocess.check_output(["git", "show", revision+":"+name], cwd=simulator)
        if raw != committed or hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError("live simulator source differs from pinned commit: "+name)
    with tempfile.TemporaryDirectory(prefix="oph-operation-reads-") as directory:
        output = Path(directory)/"capture.json"
        run = subprocess.run([sys.executable, "-m", "oph_fpe.bulk.primitive_source_reads", str(output)],
                             cwd=simulator, capture_output=True, text=True, timeout=180)
        if run.returncode:
            raise ValueError("native capture failed: "+run.stderr)
        packet = check.strict_load(output)
    check.verify(packet, source_root=simulator)
    packet = canonicalize_source_aliases(packet, source, simulator)
    result = verify(packet)["native_derivation"]
    frozen_packet = check.strict_load(HERE/"capture.json")
    frozen = verify(frozen_packet)["native_derivation"]
    # Numerical expm/normalization can differ at the last bits across libraries.
    # Both raw executions must independently pass their complete IEEE replay
    # and rational exponential enclosure. Exact source structure and derived
    # read statements must agree; hash equality of cross-platform floats is
    # neither required nor substituted for those checks.
    for actual, expected in zip(packet["cases"], frozen_packet["cases"]):
        for field in ("carriers", "cycles", "support_level", "seams", "order"):
            check.equal(actual[field], expected[field], "live source structure: "+field)
    for field in ("cases", "quantum", "M1_derived", "history_extension_is_spatial_refinement"):
        check.equal(result[field], frozen[field], "live derived semantics: "+field)
    print("Live native capture, source custody, complete replay and derived read semantics verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simulator", type=Path)
    live(parser.parse_args().simulator)
