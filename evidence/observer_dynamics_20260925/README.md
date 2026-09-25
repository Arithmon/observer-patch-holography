# Observer dynamics evidence

`manifest.json` pins the immutable producer, specification, receipt and compact
input bytes. `sim-analysis/` preserves their original relative paths;
`oph-physics-sim/` contains the frozen import dependency closure needed by the
native repair producers. `lean/` contains the canonical compilation interface
and receipt for the corresponding RER theorem modules. The original
`NecessityCore.lean` and its historical receipt remain under
`sim-analysis/codex/necessity/formal/`.

Use `python3 code/observer_dynamics/verify.py` from the repository root.
Add `--replay-refinement --write-receipt` to reconstruct all 61 fresh-noise RNG
streams and record a fresh canonical verification. The command reads the
archive and never overwrites its historical receipts. Reproduction details and
interpretation boundaries are in `../../code/observer_dynamics/README.md`.

| Evidence family | Archive path below `sim-analysis/` |
| --- | --- |
| Native occupation histories | `codex/inflation/native_history/` |
| Weighted and nonlinear readbacks | `codex/necessity/universality/` |
| Incident attempted-event records | `codex/necessity/event_records/` |
| Finite-memory filter controls | `codex/necessity/temporal_filters/` |
| Absolute finite-window source covariance | `codex/source_derivation/native/` |
| Port transport, geometry and stationary controls | `codex/dynamics/` |
| All 48 schedule-conditioned scale readouts | `data/codex_audit_20260925/scale_ensemble.json` |
| All 61 refinement chains and nested covariances | `data/refine_ensemble/`, `codex/source_derivation/innovations/` |
| Layer and FLRW clocks; golden Fourier products | `data/codex_audit_20260925/`, `scripts/` |
| 120-case boosted-clock experiment | `codex/causal/` |
| Boundary/bulk counterexample certificates | `codex/necessity/boundary_bulk/` |
| Conditional CMB observer distribution and measured comparison | `codex/observer_cmb/`, `codex/boltzmann/`, `codex/primordial/` |
| Carrier orbit census and resolved large graph gaps | `data/census_L*.json`, `data/spectral3_L*.json` |

These are observer-like self-reading systems with local state, ports, repair
moves and records. The native covariance results require declared stationary
laws, clocks and readback functionals. The source-net population and read law
are supplied. No receipt here identifies primordial curvature or derives the
observed CMB amplitude, tilt or cosmological background. Large raw arrays are
not duplicated in this compact package.

The earlier two-step refinement pilot is retained separately at
`sim-analysis/data/codex_audit_20260925/inputs/data/refine/refine_settle_L6_to_L8.json`.
Its `fresh_records=uniform5`, `fresh_law_note` and preserved producer
`sim-analysis/scripts/refine_settle_cocycle.py` specify uniform fresh loads on
`{0,...,4}`. The historical receipt's generic `definition` string still says
`{0,...,5}`; that description does not override the recorded configuration or
producer. `sim-analysis/data/refine/refine_uniform5.log` records the successful
pilot. The separate `refine_uniform6_stalled.log` ends at sweep 4,752 with
positive excess 80,064 for the `{0,...,5}` attempt; it records an unfinished
attempt, not a nontermination proof. All original bytes remain unchanged.
The refinement ensemble contains 61 completed-chain receipts. Earlier catalogue
wording mentioning 64 chains is not a completed-run count and does not establish
that three other executions were launched, failed or lost; no 64-entry launch
roster accompanies these retained receipts.
