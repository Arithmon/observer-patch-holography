# Maintainer-style audit of the fermionic source realization

This is the author's audit, not an external maintainer approval. The review
uses the [contract](CONTRACT.md) and the [proof](../../extra/FERMIONIC_SOURCE_CLOCKS.md)
as the claim boundary. It checks the construction, the evidence independently
of its hashes, and the claims made about finite executions and resources.

## Findings corrected

1. **Weighted detector verification stopped at its effect.** The original
   producer assembled the three nonzero Kraus maps algebraically. A correct
   click probability does not prove a correct post-measurement state.
   The producer now executes the native QND helper circuit, its measured
   reset, the conditional M6 acceptance pulse and the second read. The
   independent checker contracts every branch and checks completeness and
   the actual output maps, including the zero empty-acceptance history.
   The detector composes these instruments in both parity blocks. Hostile
   controls reject an omitted reset, postselection, omitted maps, changed
   acceptance pulses and a changed output with the same click effect.
2. **Numerical closeness did not establish the direction of a bound.**
   In particular, subtracting the tiny positive trace-distance allowance
   from the swing in binary arithmetic could round upward. Receipt bounds
   now round conservatively. Independent 70-digit arithmetic checks both
   the primitive bounds and the composed exposure, trace-distance,
   accounting-error and swing inequalities. Mutation tests reject nearby
   false bounds inside the former relative replay tolerance. The displayed
   accounting-error bound is less than .0051, rather than an exact .005.
3. **Two resource descriptions were inaccurate.** The packet cutoff of 80
   is spatial, not momentum. Quantum hardware is O(q^3), but the conservative
   retained ledger is O(q^4) events and O(q^4 log q) central slots at fixed
   physical volume and observation horizon. The proof and README now agree
   with the actual ledger.
4. **The accounting and initialization boundaries needed precision.**
   The positive accounting operator uses the finite reflecting walk's own
   principal logarithm, not a compression of the infinite-lattice operator.
   Equation (12a) establishes its packet bound separately. As in the parent,
   classical program compilation/distribution is finite charged prehistory
   before quantum preparation; the explicit exposure bound does not promise
   a uniform runtime for arbitrary offline compilation. Dynamic messages,
   correction and reporting remain charged inside the displayed schedule.

## Contract review

| Deliverable | Evidence and disposition |
| --- | --- |
| Even CAR algebra at arbitrary occupation | Analytic signed-loop dimension/intertwining proof; independent occupation signs and exterior minors on both complete parity blocks. Reversed edge orientations and nondefault local orders have direct controls. Odd inter-parity coherence is explicitly outside the observable claim. |
| Deterministic preparation and native operations | Every finite branch is retained. The local cycle decoder is independently inverted over GF(2), proving its linear identity for every syndrome at each checked size. The analytic prefix construction supplies the all-size proof. Native gates, full QND and weighted instruments, and unknown-state code conversion are replayed. |
| Actual local clock and schedule | The independent reflecting-permutation formula is compared with every tested one-particle matrix entry and blank-bank column. Exterior functoriality proves the all-sector extension. Endpoint-star routes, source-cell coloring and bounded degree supply the all-size schedule; a separate collision test checks the implementation. |
| Reads, interaction and accounting | Complex resolved modes, the bounded full-Fock detector, and the non-Slater two-particle interaction output are independently reconstructed. The free-clock bound is not extended to interacting evolution. Finite accounting has its own spectral estimate. |
| Physical noise | The encoded-vacuum obstruction is checked against actual edge-error mixtures. The circuit bound counts data, visitors, helpers and processors, including intermediate off-code states. Central records are reliable within the declared quantum-noise model. Neither generic fault tolerance nor an optimal noise threshold is claimed. |
| Reproducibility and hostile rejection | Compact source-bound receipt, exact schemas, strict numeric types, independent producers-disabled replay, optimized-Python rejection and physics mutations with custody bypassed. The huge clock is explicitly analytic, not an executed huge lattice. |

The exit is a constructed parity-preserving fermionic process in the parent's
specified CCG source. The audit does not promote that model realization to
selection of statistics, parameters, coherence capability or laboratory
calibration by A1--A3. These were explicit contract exclusions, not omitted
steps in the encoding, preparation or detector construction.

## Reproduction

With `PYTHONPATH=code`, run the package build and verifier, then:

```text
python -m pytest -q code/m1_fermionic_source
python -m pytest -q code/m1_source_realization code/m1_operational_clocks code/m1_necessity
python -m pytest -q tools/test_public_surface_claims.py tools/test_premise_register.py tools/test_check_axiom_consistency.py
python tools/check_claim_registry.py
python tools/check_axiom_consistency.py
python tools/check_reader_style.py
```

Use a writable `--basetemp` when required by the host. The workflow also
runs the package on Ubuntu and Windows with the pinned dependency set.
Source custody includes this audit, the contract, proof, tests and workflow.
