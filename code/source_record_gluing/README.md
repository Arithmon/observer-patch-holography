# Reference execution of the proposed RG source law

The [sufficiency theorem](../source_selection_model/RECORD_GLUING.md) derives
a round record-access cone and its causal/count limit from one explicit
additional source law. The law connects the source response's intrinsic
three-dimensional ideal to scalable, retained classical record transport.
M1's finite golden construction is a regulator of that limit; it is not
asserted to be a uniquely selected microscopic population.

This package executes the law's **reference process theory**, not a native
implementation in `oph-physics-sim`. A physical source must still satisfy the
proposed law. The existing simulator and the scalar-seam countermodel must
not be silently relabelled as implementations.

For q=2,3,5 the program prepares q^3 independent records and executes two
layers. Primitive flight duration, compared with the observation clock,
determines each read; no radius or parent menu is input. Every preparation,
immutable fork, flight, wait, local addition, commit and checkpoint enters
the stream. The audit log linearizes a parallel process: its line number is
not elapsed physical time. All arithmetic and squared times are exact.

The baseline, every single-site initial +1 intervention and every first-layer
commit +1 intervention give 323 executions and 5,344,106 events. An
intermediate intervention occurs before that version's sole immutable commit.
`capture.json` retains their hashes and operation counts in a compact packet;
the complete streams are reproducible, not committed as
millions of lines. The independent `check_process.py` imports no producer.
It reconstructs addresses, legal flights, clock budgets, immutable writers,
consumed values, local arithmetic, whole-layer completeness and resource
counts. The live verifier also checks intervention differences against path
multiplicities extracted from the consumed-writer graph. These linear
reference experiments do not prove arbitrary microphysical source behavior.
The reference's additive diagnostic is distinct from the theorem's lossless
tuple payload, which RG's finite classical maps construct analytically.
Initial-input equality alone cannot certify ancestry: two distinct stored
versions can have identical dependence on every initial input. The commit
interventions test those additional independent directions. The new Lean
noninterference reduction identifies exactly where the full intermediate
intervention premise is used to obtain native ancestry.

With `code` on `PYTHONPATH`:

```sh
python -m source_record_gluing.verify
python -m pytest -q code/source_record_gluing
```

To reproduce the compact capture and receipt:

```sh
python -m source_record_gluing.build
python -m source_record_gluing.verify --write-receipt
```

Source pins bind the scientific law, contract, code and Lean dependency
closure. Hashes provide custody; acceptance also requires live semantic
replay. Tests corrupt source preparations, clock grades, writers, payloads,
waits, dependencies, commits, intervention coverage and source custody.
The Lean audit covers every new theorem and accepts only standard axioms.
Full law interpretation and the native-order/count join are analytic proofs
in the scientific derivation, not claims made by passing JSON receipts.
