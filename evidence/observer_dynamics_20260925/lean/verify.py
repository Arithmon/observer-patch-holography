#!/usr/bin/env python3
"""Build the canonical modules and independently print every new theorem's axioms."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LEAN = ROOT / 'Lean'
MODULES = {
    'ObserverPatchHolography.EinsteinBranch.EdgeCenterTiltFreedom':
        ('ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltFreedom.lean',
         'OPH.EinsteinBranch.EdgeCenterTiltFreedom'),
    'NativeRepairFourCycle':
        ('Thermodynamics/NativeRepairFourCycle.lean',
         'OPH.Thermodynamics.NativeRepairFourCycle'),
    'TemporalFilterMonotonicity':
        ('Thermodynamics/TemporalFilterMonotonicity.lean',
         'OPH.Thermodynamics.TemporalFilterMonotonicity'),
    'Geometry.SourceNetClockSeparation':
        ('Geometry/SourceNetClockSeparation.lean', 'OPH.SourceNetClockSeparation'),
    'Geometry.FiniteLayerQuadrature':
        ('Geometry/FiniteLayerQuadrature.lean', 'OPH.FiniteLayerQuadrature'),
    'Geometry.GoldenSourceFourier':
        ('Geometry/GoldenSourceFourier.lean', 'OPH.GoldenSourceFourier'),
}
ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}
DEPENDENCIES = (
    'Geometry/GoldenSourceAssignment.lean',
    'Geometry/SourcePopulationQuadrature.lean',
    'Geometry/SourceNetLayeredOrder.lean',
    'Geometry/SourceNetCausalCone.lean',
    'ObserverPatchHolography/EinsteinBranch/EdgeCenterTiltCocycle.lean',
)
CONFIGURATIONS = ('lakefile.lean', 'lean-toolchain', 'lake-manifest.json')


def compiler_identity(report: str) -> dict[str, str]:
    """Architecture is provenance; version and source commit identify this compiler."""
    match = re.fullmatch(r'Lean \(version ([^,]+), [^,]+, commit ([0-9a-f]+), [^)]+\)', report)
    if not match:
        raise ValueError('Unrecognized Lean compiler version report')
    return {'version': match.group(1), 'commit': match.group(2)}


def check_receipt(expected: dict, actual: dict) -> None:
    """Compare every proof pin while permitting an equivalent compiler platform."""
    for receipt in (expected, actual):
        if receipt['compiler_identity'] != compiler_identity(receipt['lean_version']):
            raise ValueError('Compiler identity disagrees with its provenance report')
    if expected['compiler_identity'] != actual['compiler_identity']:
        raise ValueError('Lean version or source commit differs from the frozen receipt')
    # Preserve the original complete platform report; --check never rewrites it.
    expected_comparison = {key: value for key, value in expected.items() if key != 'lean_version'}
    actual_comparison = {key: value for key, value in actual.items() if key != 'lean_version'}
    if expected_comparison != actual_comparison:
        raise ValueError('Compiler receipt differs from fresh verification')



def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str]) -> str:
    result = subprocess.run(command, cwd=LEAN, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(result.stdout)
    return result.stdout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Check the frozen receipt without rewriting it.')
    parser.add_argument('--full-build', action='store_true', help='Build every default library first.')
    args = parser.parse_args()
    version = run(['lake', 'env', 'lean', '--version']).strip()
    identity = compiler_identity(version)
    dependency_rows = [{'source': 'Lean/' + relative, 'sha256': sha(LEAN / relative)}
                       for relative in DEPENDENCIES]
    configuration_pins = {'Lean/' + relative: sha(LEAN / relative)
                          for relative in CONFIGURATIONS}
    source_rows = []
    names = []
    for module, (relative, namespace) in MODULES.items():
        source = LEAN / relative
        declared = re.findall(r'^(?:theorem|lemma)\s+(\w+)', source.read_text(), re.M)
        names.extend(namespace + '.' + name for name in declared)
        source_rows.append({'module': module, 'source': 'Lean/' + relative,
                            'sha256': sha(source), 'declaration_count': len(declared)})
    run(['lake', 'build'] + ([] if args.full_build else list(MODULES)))
    with tempfile.TemporaryDirectory(prefix='oph-observer-dynamics-lean-') as temp:
        audit = Path(temp) / 'AxiomAudit.lean'
        audit.write_text(''.join('import ' + module + '\n' for module in MODULES)
                         + '\n' + ''.join('#print axioms ' + name + '\n' for name in names))
        output = run(['lake', 'env', 'lean', str(audit)])
    dependencies = {}
    for name, values in re.findall(r"'([^']+)' depends on axioms:\s*\[(.*?)\]", output, re.S):
        actual = {value.strip() for value in values.split(',') if value.strip()}
        if not actual <= ALLOWED:
            raise SystemExit('Unexpected proof axioms for ' + name + ': ' + str(actual - ALLOWED))
        dependencies[name] = sorted(actual)
    for name in re.findall(r"'([^']+)' does not depend on any axioms", output):
        dependencies[name] = []
    if set(dependencies) != set(names):
        raise SystemExit('Missing or additional declaration axiom report')
    for row in source_rows + dependency_rows:
        if row['sha256'] != sha(ROOT / row['source']):
            raise SystemExit('Proof source changed during verification')
    for relative, pinned in configuration_pins.items():
        if sha(ROOT / relative) != pinned:
            raise SystemExit('Lean configuration changed during verification')
    log = '\n'.join(name + ': ' + ', '.join(dependencies[name]) for name in names) + '\n'
    receipt = {
        'schema': 'oph.observer-dynamics.canonical-lean.v1',
        'lean_version': version,
        'compiler_identity': identity,
        'build_scope': 'six_named_canonical_modules',
        'local_dependency_sources': dependency_rows,
        'project_configuration_sha256': configuration_pins,
        'build_exit_code': 0,
        'sources': source_rows,
        'declaration_count': len(names),
        'allowed_axioms': sorted(ALLOWED),
        'declaration_axioms': dependencies,
        'verifier_sha256': sha(Path(__file__)),
        'axiom_log_sha256': hashlib.sha256(log.encode()).hexdigest(),
        'scope': ('Finite rational repair witness, survival-cocycle parameter freedom, pointwise '
                  'filter inequality, layered-clock separation, exact discrete-layer parity sums, '
                  'finite FLRW polynomial mismatch, and finite golden-orbit Fourier identities. '
                  'No stochastic CLT, geometric FLRW integral, fixed-mode Fourier asymptotic, '
                  'or cosmological source identification is formalized here.'),
    }
    if args.check:
        try:
            check_receipt(json.loads((HERE / 'compiler_receipt.json').read_text()), receipt)
        except ValueError as error:
            raise SystemExit(str(error)) from error
        if (HERE / 'axiom_audit.log').read_text() != log:
            raise SystemExit('Axiom log differs from fresh verification')
    else:
        (HERE / 'compiler_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        (HERE / 'axiom_audit.log').write_text(log)
    print(f'{len(names)} declarations compiled and audited; only ordinary Lean foundations.')
    if args.full_build:
        print('All default Lake libraries also built successfully.')


if __name__ == '__main__':
    main()
