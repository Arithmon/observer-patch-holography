# Finite federation log recovery

This package tests how much a bounded self-reading federation can infer from
its public event records. Carriers have local scalar state, twelve ports,
protected readback records and declared feedback operations. The decoder
receives authenticated event rows only. It reconstructs the observed
population and routes, checks consumed writers and arithmetic, and computes
the finite logical ancestry relation. It cannot read the generating support
or population specification.

`recover.py` supplies `recover(rows)` and `instantiate(recovered)`. Recovery
accepts the retained baseline grammar, with explicit finite resource limits.
The recovered program stores instruction operands and preparation literals;
instantiation recomputes every derived value and writer. Two iterations
reproduce the observed instruction program. This is a finite trace round
trip, not a fixed point of complete federation recovery.

`custody.py` authenticates and decodes the retained compressed/differential
q3 histories outside the blind interface. It also compiles the declared C++
producer and executes both full gluing inputs under all four retained
variants. `topology.py` compares the two graph completions. `verify.py`
independently reconstructs the event and graph invariants without importing
the recovery or topology implementation. `build.py` generates `receipt.json`.
The receipt pins source and retained archive bytes. The paper's scientific
statement is in `paper/tex_fragments/FEDERATION_LOG_RECOVERY.tex`.

The named identification failure is the complete gluing's triangle count.
The counterexample class consists of connected finite twelve-port gluing
inputs to the declared feedback program. The alternative graph is not
asserted to be a canonical icospherical refinement. Opcode meanings,
preparation, half-unit arithmetic, input-log completeness and level-three
provenance are supplied. No general population law, source-selected repair
law, physical clock, quantum instrument or cosmic uniqueness follows.

From the repository root, with the pinned Python dependencies and a C++17
compiler available:

```sh
python3 evidence/federation_recovery/build.py --check
python3 evidence/federation_recovery/verify.py
python3 -m pytest -q evidence/federation_recovery
```

To regenerate the certificate after an intentional source change, run
`python3 evidence/federation_recovery/build.py` and independently verify it.
Native runs use temporary directories and retain only digests; no production
q13/q21 tape is generated or downloaded. Existing routing histories and
contributor implementation files are inputs and are not rewritten.
