# Protected and live source records

`build_population.py` executes the frozen `specification.json`: one initialized
twelve-port observer-like carrier at each q=13 golden record, one active
pair-mean repair per nonconstant carrier, local read/write records and an
authenticated operation history. It compares live readback positions with
separately retained initialized address records. The source population and
protected write footprint are declared; no routing or physical selection is
inferred. The operation index is not physical time.

Run from the research repository root with Python and NumPy:

```sh
python3 -B code/source_population/build_population.py
python3 -B code/source_population/verify_population.py
python3 -m pytest -q code/source_population/test_population.py
```

The verifier imports the hash-checked source-archive reconstruction to obtain
the original population and carrier incidence. It imports no producer. It
replays actual rational writes and independently encloses every pair-distance
comparison using a rational bracket for sqrt(5). The source proof is in
`paper/tex_fragments/SOURCE_RECORD_PROTECTION.tex`; the linear invariance
characterization is kernel-checked in `Lean/Geometry/SourceRecordProtection.lean`.
The receipt is a scoped finite experiment, not a certificate that every
population law or physical spacetime is impossible.
