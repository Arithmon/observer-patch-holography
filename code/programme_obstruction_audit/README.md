# Fixed-formulation scientific audit for issue #1078

This is supplemental evidence for the [programme obstruction audit](../../docs/OPH_PROGRAMME_OBSTRUCTION_AUDIT_1078.md), using existing claim IDs. The scientific target is `0510578e317de5da0d485d09792dc8bcad2e2f3c`. It is not a security scan, new claim registry, full Lean build, physical experiment, or exhaustive model search. The programme-wide exclusion verdict is **inconclusive**; narrower attachment/no-entailment results and the complete analytic abstract scalar-seam witness are retained with their interpretations and escape routes.

## Run the seven offline controls

From this directory, with an ordinary Python interpreter and assertions enabled:

```sh
python3 run_all.py
```

Only the Python standard library and packaged files are required. No network or files outside this package are read. The entry point uses `sys.executable` for all seven subprocesses, validates package SHA-256 hashes, checks the immutable response fixture against its hardcoded SHA-256 and Git blob SHA-1, validates the fixed target/source provenance, and stops at the first failed assertion, process failure or 120-second per-script timeout. It prints a concise verdict and saves all individual JSON receipts plus a full `summary.json` in `run_receipts/`. Use `--output-dir relative-directory` to choose another receipt destination. Do not use `python -O`; these scientific controls rely on assertions. The publication run used **Python 3.14.7**, exact rational arithmetic and explicit Q(sqrt5) pairs; other supported Python versions can differ only in runtime metadata and their reported script/environment fields, not the exact mathematical results.

| Script | Prescribed inputs and controls | Meaning and stopping boundary |
| --- | --- | --- |
| `channel_author.py` | Rational 8x8 binary seam channel, four Kraus operators, private-controlled and nonnormal-exchange rotations; all 64 matrix units; completion n=2,3; two private inputs | Exact no-feedback control, CPTP/completion and feedback difference 144/625. Toy channel class only; no OPH source membership. Stop after these cases. |
| `channel_independent.py` | Independent private-first tensor order, weighted Kraus reconstruction, all 64 matrix units, completion projector identities and 64x4 private commutator equations | Independent exact feedback/restricted-class validation; scalar-only private commutant excludes every single-product exponential without logarithm-branch assumptions. |
| `joint_author.py` | Fixed twelve-vector Q(sqrt5) response fixture, at most60 stored paths,66 binary configurations, four oscillator cutoffs and three clock settings; s=2 units and idle controls | Corrected24 finite gates for carrier/triplet/CCR/clock attachments. Does not certify full source membership. |
| `joint_independent.py` | Independent D formula/rank, all 20 Cayley factors and60 path words,3960 event branches,431 A3 omission-coordinate controls, posterior and three current-unit controls | 907 prescribed gates; overlaps author/event checks and is not an extra907 scientific findings. Separate reconstruction imports no author implementation. |
| `kinetic_control.py` | A=I+path3 stiffness; e=0,+1/16,-1/16,-1/4; gamma=1,3/2 | Eight rational cases. Singular -1/4 is rejected. Same q-Gibbs marginal is an analytic integration result; full reference/vacuum changes and same-capacity inertia premise fails. |
| `attempt_independent.py` | Incidence reconstructed separately;66 states x30 seams x2 coins; s=2, idle half and unequal edge-weight controls | Exact named-process components B,nu,kappa and mean law. The kappa-weighted ratio needs unit atomic current; 3B/nu is current-unit invariant. No full-core independence assertion. |
| `compact_classifier.py` | Exact partitions of11,12 using imported simple dimensions3,8,10 | Unmodified upstream tiny classifier; finite arithmetic only, not reproof of compact-simple classification or physical gauge realization. |

The matrix-unit controls cover arbitrary complex operators by complex linearity despite real rational basis entries. Conjugated completion at n=2,3 is a regression check; the all-n identity follows analytically by unitary conjugation of the native decomposition. The unit/scheduler controls deliberately change their declared convention/law and do not falsify the frozen relation. No stochastic search, fitting or absence-of-witness inference occurs.

Individual commands are also portable:

```sh
python3 channel_author.py > author-result.json
python3 channel_independent.py > independent-result.json
python3 joint_author.py data/response.json > joint-result.json
python3 joint_independent.py > independent-joint-result.json
python3 kinetic_control.py > kinetic-result.json
python3 attempt_independent.py > attempt-result.json
python3 compact_classifier.py > classifier-result.json
```

`receipts/` contains the final publication execution and two retained upstream replays. `provenance/package_sha256.json` covers every package file except itself. The publication file manifest additionally covers that manifest, the report and the exact maintained-document delta. `provenance/packaging_adaptations.json` preserves original script hashes and describes the only adaptations: author joint receipt uses a stable fixture label plus actual-input SHA256; independent joint fixture path becomes package relative; kinetic control emits stdout without a scratch sidecar and reports its completed independent-review status. Arithmetic, cases and scientific scope are unchanged.

## Reproduce the separately scoped upstream components

The retained `upstream_components.json` and `upstream_galois.json` were executed on fixed source bytes with **Python 3.14.7, SymPy 1.14.0, mpmath 1.3.0**. They are not counts to add to the offline reconstructions. The seven components equal the upstream retained receipt: response, algebra, constraint grammar, central channels, four-level support, records and refinement. Both structural Galois laws pass response/observable rank12, restricted rank11, even/odd rank6,66 brackets,20 factors,60 actions,720 covariance pairs and3600 compositions. Whole-law conjugation includes the quadrupole and fixes i; it does not map arbitrary PSD preparations/effects.

To rerun those components, provide a **separate exact fixed-source fixture or checkout**, and an environment with those pinned dependencies. This optional step has explicit external-source/dependency requirements; the offline command above does not. The helper validates the source paths/hashes listed in `fixed_sources.json` and `additional_source_pins.json`, then uses the same upstream component calls and canonical JSON comparison as the original wrapper. It does not fetch anything or alter the source fixture. From a repository root, where `FIXED_SOURCE` names that separately prepared fixture:

```sh
python3 code/programme_obstruction_audit/replay_upstream_components.py --source-root "$FIXED_SOURCE" > upstream-component-replay.json
python3 "$FIXED_SOURCE/code/a5_closure/galois_source_response_control.py" > upstream-galois-replay.json
```

The exact Galois script has Git blob SHA-1 `ecb5f9af4b9ab9d3ccde131e7228e174d4385c11`, SHA-256 `71cfa71dce4080febc61e09b3675f8fc8822294acdbacfbd1f87012e9e0e618d`. It validates its four source pins internally. The packaged response snapshot has Git blob SHA-1 `e56b96a9429ac618aea9de7cd5ffa53a13354ff6`, SHA-256 `226ab001894f357d3732d8a898fbae80e28b9e16ca41104a6f293190b441e053`. Replaying a newer checkout by assumption is invalid; use the frozen hashes. The packaged optional component helper was also executed successfully on that fixed fixture; its receipt is `receipts/upstream_components_portable.json`. No full upstream package or installed dependency tree is vendored here.

These replays **did not run** top-level recursive source custody, producer regeneration, the full hostile-input suite, N4's separate 187-check numerical packet, scalar closure rebuilding or Lean compilation. The optional helper verifies the finite fixture paths it actually consumes; that is not a claim that the omitted top-level upstream verifier ran. All-level tower/operational model membership rests on the separately challenged analytic proof, not four finite support levels or a passing receipt. Historical source receipt bytes remain untouched. Runner integration controls additionally reject a modified fixture even when its manifest is edited consistently, and reject optimized Python before executing the scientific controls; `receipts/runner_refusal_controls.json` records those justified refusals.

## Independent review and evidence boundary

`provenance/review_status.json` records final methods, reviewers, note hashes and accepted corrections. Author replays establish reproducibility. Independent proofs use a separate implementation, tensor order or algebraic reconstruction, and source clause review challenges admissibility/quantifiers. The initial audit's current-unit overclaim was corrected; the repository already states the scoped rate relation correctly. The scalar-seam source is an already claimed full **abstract** model under its actual typed interpretation, rather than an absent common abstract witness. Neither correction is a repository theorem defect.

The exact current publication base is separately identified in `publication_base.json`; whole current-main scientific interpretation is not reaudited. A subsequent local Git-object custody check matched all 38 pinned source blobs against both Git SHA-1 and SHA-256; `provenance/local_object_verification.json` records that later verification without asserting a clean full filesystem status or rerunning scientific results. Existing generated premise-map rows, registry entries, theorem statements, live upstream producers and frozen evidence are unchanged. The maintained-document edit is an audit-evidence link in `docs/MODEL_SEARCH.md`. Because the new report changes the scanner's active surfaces, `claims/active_surface_inventory.json` is regenerated with unmodified `tools/check_axiom_consistency.py --inventory`; this generated integration inventory does not change the adopted basis, claim statuses or premise discharge. The portable archive omits this repository-wide inventory, just as it omits the maintained model-search document; the publication manifest covers both. All new generated audit receipts are produced by these packaged scripts. There are no measured observations, laboratory preparations, resource lower bounds or programme-wide decisions encoded in this package.
