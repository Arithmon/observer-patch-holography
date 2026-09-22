"""Execute the pinned simulator and compare independently derived read semantics."""
import argparse
from fractions import Fraction
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile

from . import simulator_verifier as check
from .verify import HERE, verify


def compare_replays(packet, frozen_packet, result, frozen):
    """Bind verified replay semantics to the retained numerical preparation.

    The absolute 2^-40 comparison budget is a reproducibility tolerance for
    two binary64 executions, not a physical error bound or an exact identity.
    Both packets must separately pass the independent replay verifier first.
    """
    def close(values, expected, label):
        if len(values) != len(expected):
            raise ValueError("live numeric shape: "+label)
        for a,b in zip(values,expected):
            if abs(Fraction(check.hex_number(a))-Fraction(check.hex_number(b))) > Fraction(1,2**40):
                raise ValueError("live numeric reproduction: "+label)

    check.equal(len(packet['cases']),len(frozen_packet['cases']),"live case count")
    for actual, expected in zip(packet["cases"], frozen_packet["cases"]):
        for field in ("carriers", "cycles", "support_level", "seams", "order"):
            check.equal(actual[field], expected[field], "live source structure: "+field)
        for field in ('initial_hex','final_hex'):
            close(actual[field],expected[field],field)
        check.equal(len(actual['observer_records']),len(expected['observer_records']),"live record count")
        for observed, retained in zip(actual['observer_records'],expected['observer_records']):
            check.equal(observed[:3],retained[:3],"live record identity")
            close(observed[4],retained[4],'observer full port state')
    for field in ('laplacian','step','phase_ports','uniform_snapshot_feedback'):
        check.equal(packet['quantum'][field],frozen_packet['quantum'][field],"live quantum source: "+field)
    for field in ("cases", "quantum", "instruments", "M1_derived", "history_extension_is_spatial_refinement"):
        check.equal(result[field], frozen[field], "live derived semantics: "+field)


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
    result = verify(packet)["native_derivation"]
    check.verify(packet, source_root=simulator)
    frozen_packet = check.strict_load(HERE/"capture.json")
    frozen = verify(frozen_packet)["native_derivation"]
    # Numerical expm/normalization can differ at the last bits across libraries.
    # Both raw executions must independently pass their complete IEEE replay
    # and rational exponential enclosure. Exact source structure and derived
    # read statements must agree, and every ledger/snapshot value must reproduce
    # within the stated comparison budget. Summary equality alone is not enough.
    compare_replays(packet, frozen_packet, result, frozen)
    print("Live native capture, source custody, complete replay and derived read semantics verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("simulator", type=Path)
    live(parser.parse_args().simulator)
