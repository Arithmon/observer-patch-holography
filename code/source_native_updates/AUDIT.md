# Native writes and chain-scaling audit

Audit date: 2026-09-20. Baseline: #917 at `3a181999`. The user's larger request
to close the M1 gaps is **not fully satisfied** by this packet. The finite
write and general cleanup results below are genuine advances; no complete
q=13/q=21 native compiler or full-axiom selection/no-go theorem is claimed.

## Objective, deliverables and exit

The [contract](CONTRACT.md) defines the bounded construction and names every
remaining M1 obligation. The original #777 Objective/Deliverables/Exit and
#740's common-family boundaries were re-read, as were the canonical A1--A3
reference and the candidate scalar-law definitions. A smaller construction
has not been relabeled as the full issue's exit.

| Deliverable | Evidence | Audit conclusion |
| --- | --- | --- |
| Native arithmetic | `SourceNativeRecords.lean`: general halving, seven-mean sum, old-record scaling, scalar compilation | The new sum is computed from consumed source values; no payload-dependent preparation of C |
| Arbitrary chain cleanup | `SourceBusScaling.lean`: actual forward word, weighted contraction, repeated bounds, block halving, operation counts | Works for every n>=1; n=0 has its own exact clear result; archives outside the word are untouched |
| Captured write/read history | All 16 input pairs, 2,279 means per history, 64 reads | Each support edge and consumed writer is independently reconstructed; all decodings have positive margins |
| Error and resource accounting | Exact linear coefficients plus preparation, per-step, rounding and readout allowances | Minimum certified margin exceeds 0.0032996; noise allowances and excluded physical costs remain inputs |
| Serious negative controls | 136 new tests, real CLI failures, resealed mutations, alternative valid program, insufficient precision | Invalid or semantically changed evidence is rejected, including a null trace |
| Selection boundary | Exact two-word schedule control and raw-range theorem | Only the stated local-law/raw-value conclusions follow; neither is promoted to full A1--A3 non-derivability |

## Issues found and corrected during audit

- An explicit `null` representative trace initially bypassed the per-event
  comparison, because `None` also meant an intentionally compact history.
  The top-level verifier now requires exactly 2,279 representative events.
  Direct and real-CLI regressions reject `null`, empty and truncated traces.
- Merely importing the new audit through `Geometry.lean` would not ensure
  its execution on every dependency-only CI change: the per-change Lean
  workflow uses explicit roots and a changed-module selector. The new
  transitive audit is an unconditional root. Tests reject deleting,
  commenting out or replacing that root with a proof-only module.
- Windows exposed an implicit text-decoding assumption in the audit-coverage
  test. The source read now explicitly uses UTF-8. Temporary-directory access
  was restricted locally; Windows checks used a fresh workspace-local
  `--basetemp`, without altering the tests' scientific assertions.
- The resource review added record identifiers, publication flag and codebook
  entry widths alongside the scalar registers, arithmetic, program, port map
  and counters. The README names the representation and excluded physical
  costs, instead of equating the active registers with the full support.
- One new Lean linter warning was removed and the exact final proof and
  dependent artifact pins were rechecked. No old proof, tolerance or receipt
  was changed to obtain a pass.

## Independent checks

The producer expands a structured sequence of pair operations; the verifier
uses its own index/phase formula. They use different nearest-even arithmetic
and different receiver classification rules. Only serialization, source
identity and strict equality helpers are shared. A test disables every
producer execution entry point while verifying the retained artifact.

The verifier computes the entire ideal two-input linear map and checks the
write formula, clean initial bus, each receiver coefficient, uniform error
bound and subsequent cleanup. It reconstructs all sixteen complete fixed-point
histories and checks each committed state, every receiver, final state and
trace commitment. One full history remains directly inspectable; other
histories are deterministically expanded and semantically checked on every
run. The recorded hash is not a substitute for that replay.

Seventy event mutations change each of seven fields at ten positions spanning
the write, transfer and cleanup phases. Their tape commitment is recomputed.
The verifier rejects them. A second valid native program swaps the two
independent initial exports and still computes the correct sums; its complete
coherently regenerated evidence is rejected as a different controller.
Other controls cover input/record/version/scale substitutions, omitted and
duplicated cases, altered final states, source pins, floats, booleans,
noncanonical rationals, duplicate JSON keys and nonfinite numbers.

Real CLI tests reject optimistic margins/costs, missing histories, false
ancestry, forged commits, altered bytes and missing traces. Three full runs
add signed preparation, per-step and readout disturbances and remain within
the proved bounds. Q=2^8 both fails the certified margin and exhibits actual
readout errors. Matrix-basis controls include n=0, the smallest positive
chains, and larger chains; nonnegative coefficients check the weighted
absolute-value estimate. The universal claim comes from Lean, not those
finite samples.

All 39 new theorem declarations (27 scaling, 12 records/boundaries) undergo
the existing transitive `Lean.collectAxioms` gate. Only `propext`,
`Classical.choice` and `Quot.sound` are permitted. Its guarded `sorryAx` and
compiler-trust controls compile. A separate probe importing the new audit,
declaring `axiom injected : False`, proving an indirect theorem from it and
auditing that theorem failed with the expected rejection of `injected`.
The new source files contain no `sorry`, `axiom` or `native_decide` proofs.

## Executed validation

- Windows: **351 passed** across the new packet, reusable bus, encoded
  memory and mandatory/Lean-CI guards. The 136 new controls were rerun after
  the final proof-pin refresh: **136 passed**.
- Linux: the same focused set: **351 passed**, followed by successful
  standalone native-update verification.
- Linux frozen mandatory runner: standard shard index 0 of 6, all 24 stages
  passed, including external custody and FZ-11/FZ-12 history verification.
  This is not a claim that every mandatory shard was locally rerun.
- Lean: `lake build Geometry.SourceNativeUpdatesAxiomAudit` passed, and the
  full `lake build Geometry` passed (8,430 jobs). The final proof-only linter
  cleanup was followed by another successful audit build (8,254 jobs).
- Eleven toolchain, manifest, imported scalar/memory and new proof/audit
  inputs were byte-compared with the WSL build inputs after the final linter
  cleanup. All matched. The new theorem audits report only the standard axioms.
- `git diff --check` passed. The frozen runner remains SHA-256
  `04bab737fa2d1b7b51c5545242e40c20375eea98d1a41eb0645ba32238360291`.

The two retained JSON artifacts total 156,241 bytes and 3,636 lines. The
complete representative trace uses one scalar-array event per line. This
avoids repeating the earlier 200,000-line formatting overhead while retaining
the evidence needed for deterministic full replay.

## Remaining scientific obligations

The full native service still needs general write programs, bias supply,
many live versions, allocation, version publication/control and robust
precision/capacity scaling. General-chain existence on a supplied paired
embedding does not prove that every captured route has such an embedding.
The q=13/q=21 feedback histories of #910 are not relabeled as native means.

The schedule witness establishes ambiguity under supported scalar means; it
does not instantiate A1's complete carrier/federation/response clauses,
A2's accepted-domain/endogeneity clauses or A3's complete generated constraint
grammar. A complete qualifying countermodel or positive selection theorem
remains open. Declaring the desired metric menu as an A3 constraint would
assume the missing conclusion. M1 and #777/#779/#740 remain open.
