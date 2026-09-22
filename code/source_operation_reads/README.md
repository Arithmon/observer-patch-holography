# Source operations to reads

The [derivation](DERIVATION.md) classifies the registered simulator's complete
pair-repair algebra, reconstructs actual numerical observer records, and
proves its read limits. It includes normalized source witnesses, the
nonclosure of intensity-only recurrent propagation, and a native feedback
equivariance control. Its axiom table states what the executable supplies.

The compact capture is produced by [simulator PR #16](https://github.com/muellerberndt/oph-physics-sim/pull/16).
`source_snapshot.json` pins its commit and transitive local Python sources.
`simulator_verifier.py` is a byte-identical copy of that contribution's
independent standard-library verifier. It imports no producer or simulator.
The receipt additionally binds this derivation and the transitive local Lean
proof sources. Raw event tapes are reconstructed rather than archived.

With `PYTHONPATH=code`, from the OPH repository root:

```sh
python -m source_operation_reads.verify
python -m pytest -q code/source_operation_reads
python -m source_operation_reads.check_live /path/to/pinned/simulator
```

The live check requires the exact simulator revision and source bytes,
executes the native driver afresh, and independently checks its primitive and
observer output. It compares exact structure and derived read semantics
across platforms. Last-bit numerical differences must pass the rational
unitary enclosure and complete IEEE ledger replay; hashes of floating-point
expm results are not claimed to be universal across numerical libraries.

The formal target is `Geometry.SourceOperationReadsAxiomAudit`. Matching
composition, all-word dependencies, decoder obstruction and scalar/phase
arithmetic are kernel checked. Normalized rank, finite monoid size, matrix
Taylor bounds and their executable identification are analytic and independently
tested. No test census is substituted for the all-word proof.

Global seeds and custody hashes are not included in the local numeric read
channel. Exact real inverses are distinguished from thresholded, rounded
execution. The chosen driver is distinct from the shared-vertex W12 family
and the group-valued repair engine. It neither instantiates complete A1--A3
nor derives M1's population, radius or canonical temporal instrument.

To reproduce the retained source revision before the companion is merged:

```sh
git clone https://github.com/MarioPoneder/oph-physics-sim.git
cd oph-physics-sim
git checkout cfe291881fe09accdd6e80600db335b6505adb64
python -m pip install -e ".[test]"
```

The capture is about 72 KB and contains replay inputs, not expanded event
tapes. The complete seven-case attempt denominator is retained. CI checks
both the offline proof artifact and a fresh native execution on Linux and
Windows. The canonical support criterion used here is already present in
`source_read_selection`; this contribution does not depend on an unmerged
stronger theorem.
