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
Every initial, final and recorded ledger value must also reproduce within
absolute 2^-40. This is an explicit cross-platform comparison tolerance,
not a physical precision claim. Matching only ranks and graph summaries
does not establish reproduction of the retained numerical preparation.

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
git checkout 84d46b97b07965406c608e1465bec18b2a308372
python -m pip install -e ".[test]"
```

The capture is 118 KB and contains replay inputs, not expanded event
tapes. The complete seven-case attempt denominator is retained. CI checks
both the offline proof artifact and a fresh native execution on Linux and
Windows. The canonical support criterion used here is already present in
`source_read_selection`; this contribution does not depend on an unmerged
stronger theorem.

The [instrument continuation](INSTRUMENT_DERIVATION.md) goes beyond a missing
receipt: it proves affine obstructions for the native no-op rule and phase
lift, constructs two CPTP extensions of the ideal scalar repairs, and
executes a subsequent unitary read that distinguishes them. It derives a
recurrence influence cone and a 5N+1-dimensional normalized erasure kernel.
Influence is not promoted to an exact read. A separate algebra-generation
argument rules out identifying the coherent probes with central records;
an exact lumpability test rules out carrier-only Markov state for the
measured extension. These controls do not select a new canonical model.

Two constructive completions prevent that erasure with explicitly added
records: a minimal 5N+1-coordinate linear difference archive reverses the
deterministic ideal population means, and retained random-unitary flags
allow recovery of the quantum input. The latter matches the mean only in
expectation; its flags alone contain no input information. Neither archive
is claimed to be present in the native source or to supply a local M1 read.
