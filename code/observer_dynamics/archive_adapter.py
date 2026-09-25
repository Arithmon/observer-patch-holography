"""Execute preserved producers against archived, hash-checked source bytes.

The historical producers requested seven source blobs with ``git show``. The
archive contains exactly those blobs, so no simulator checkout or git history is
required. Only that specific read command is intercepted; other commands retain
normal subprocess behavior. No scientific implementation or receipt is altered.
"""
from contextlib import contextmanager
from functools import wraps
import hashlib
import importlib.machinery
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from unittest.mock import patch

RER = Path(__file__).resolve().parents[2]
ARCHIVE = RER / "evidence/observer_dynamics_20260925"
SIM = ARCHIVE / "sim-analysis"
VENDOR = ARCHIVE / "oph-physics-sim"
COMMIT = "14d16994b12adb2320b994044bc75db02378b8b1"


def portable_metadata(module):
    """Keep archived default producers' relative provenance paths in POSIX form.

    The retained JSON uses repository-relative slash paths. Historical producers
    used str(Path.relative_to(...)), whose separators depend on the host. Adapt
    only these generated metadata fields; comparisons, numerical values, input
    hashes and immutable producer bytes remain unchanged.
    """
    filename = getattr(module, "__file__", None)
    if not filename:
        return
    try:
        relative = Path(filename).resolve().relative_to(SIM.resolve()).as_posix()
    except ValueError:
        return
    functions = {
        "codex/observer_cmb/observer_spectrum.py": "make_output",
        "codex/observer_cmb/measured_comparison.py": "load_measurements",
        "codex/boltzmann_data/audit_inputs.py": "verify_downloads",
        "codex/boltzmann/summarize.py": "build",
    }
    name = functions.get(relative)
    if name is None:
        return
    original = getattr(module, name)
    if getattr(original, "__observer_posix_metadata__", False):
        return

    @wraps(original)
    def portable(*args, **kwargs):
        result = original(*args, **kwargs)
        if name == "make_output":
            paths = result["inputs_relative_to_codex"]
            result["inputs_relative_to_codex"] = {key: value.replace("\\", "/") for key, value in paths.items()}
        elif name == "load_measurements":
            for row in result[1]:
                row["path"] = row["path"].replace("\\", "/")
        elif name == "verify_downloads":
            for row in result:
                row["path_relative_to_codex"] = row["path_relative_to_codex"].replace("\\", "/")
        elif name == "build":
            paths = result["observational_input_hashes"]
            result["observational_input_hashes"] = {key.replace("\\", "/"): value for key, value in paths.items()}
        return result

    portable.__observer_posix_metadata__ = True
    setattr(module, name, portable)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_manifest(directory=ARCHIVE):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text())
    total = 0
    for relative, row in manifest["files"].items():
        path = directory / relative
        if not path.resolve().is_relative_to(directory.resolve()):
            raise ValueError(f"Evidence path escapes archive: {relative}")
        if path.stat().st_size != row["bytes"] or sha(path) != row["sha256"]:
            raise ValueError(f"Evidence digest mismatch: {relative}")
        total += row["bytes"]
    return {"files": len(manifest["files"]), "bytes": total,
            "manifest_sha256": sha(directory / "manifest.json")}


@contextmanager
def frozen_sources():
    manifest = json.loads((ARCHIVE / "manifest.json").read_text())
    original = subprocess.check_output
    original_exec = importlib.machinery.SourceFileLoader.exec_module

    def exec_archived(loader, module):
        original_exec(loader, module)
        portable_metadata(module)

    def read(args, *positional, **keywords):
        if isinstance(args, (list, tuple)) and len(args) == 5 and args[:2] == ["git", "-C"] and args[3] == "show":
            repository = Path(args[2]).resolve()
            if repository == VENDOR.resolve():
                commit, relative = args[4].split(":", 1)
                if len(commit) < 7 or not COMMIT.startswith(commit):
                    raise ValueError("Requested simulator revision is not archived")
                row = manifest["files"].get("oph-physics-sim/" + relative)
                if row is None or row.get("origin_commit") != COMMIT:
                    raise ValueError(f"Unarchived simulator source: {relative}")
                data = (VENDOR / relative).read_bytes()
                if hashlib.sha256(data).hexdigest() != row["sha256"]:
                    raise ValueError(f"Frozen source digest mismatch: {relative}")
                return data.decode(keywords.get("encoding") or "utf-8") if keywords.get("text") or keywords.get("encoding") else data
        return original(args, *positional, **keywords)

    sys.path.insert(0, str(VENDOR))
    try:
        with patch.object(subprocess, "check_output", read), patch.object(
            importlib.machinery.SourceFileLoader, "exec_module", exec_archived
        ):
            yield
    finally:
        sys.path.remove(str(VENDOR))


def load(relative):
    path = SIM / relative
    name = "observer_archive_" + relative.replace("/", "_").replace(".", "_")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    sys.path.insert(0, str(path.parent))
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(str(path.parent))
    return module
