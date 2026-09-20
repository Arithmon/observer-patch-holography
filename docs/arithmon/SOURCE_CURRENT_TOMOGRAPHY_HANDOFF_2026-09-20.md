# SOURCE-CURRENT-TOMOGRAPHY-0 handoff — 2026-09-20

## Repository state

- Base branch SHA: `56c77645a80b06f96c2d99088864c9c945f803f3`
- Synchronized fork-main SHA: `2b03a95caf5030272f7b426b964e816f650153f8`
- Final branch SHA: the commit containing this report; recorded in the external
  handoff because a Git commit cannot contain its own SHA.
- Scientific commits: 2 (Stage 1 syntax repair; Stage 2 inventory and claim
  integration), followed by the fork-main synchronization and custody fix.

## Primary scientific verdict

`INSUFFICIENT_ORDER_SENSITIVE_SOURCE_DATA`

The Stage 2 bounded inventory found 21 candidate packets and zero qualifying
source objects. No existing single packet jointly supplies twelve reversible
source-native perturbation families on the same W12 carrier, both mixed
composition orders as canonical raw histories, target freedom, and refinement
provenance. The exact machine verdict is
`SOURCE_CURRENT_ORDER_SENSITIVE_OBJECT_NOT_PRESENT`.

The inventory also commits a canonical 432-path snapshot covering every
non-cache file below the 13 audited directories. Both verification lanes
recompute the list and fail if a file is added, removed, or renamed, preventing
a new candidate surface from remaining silently outside the 21-row inventory.

The strongest registered source packet is reversible but its response-word
algebra has exact dimension four, is commutative, and contains only the
identity among 60 proper rechartings; it has no raw runtime histories. The
richest logged W12/refinement packet is append-only and strict-descent, and its
commutative addition law erases order. Combining their complementary fields
would create a new producer and was therefore rejected.

## Stage table

| Stage | Status | Result |
|---|---|---|
| Stage 0 | PASS / unchanged | Fail-closed baseline remains `INSUFFICIENT_SOURCE_DATA`. |
| Stage 1 | PASS / revised | Three malformed literal newline escapes repaired; contract and seven tests pass. |
| Stage 2 | NEGATIVE | 21 candidates, 0 qualifiers; pinned certificate, independent verifier, and mutation gates pass. |
| Stage 3 | NOT RUN | Not authorized by the Stage 2 stop condition. |
| Stage 4 | NOT RUN | Not authorized. |
| Stage 5 | NOT RUN | Not authorized. |
| Stage 6 | NOT RUN | Not authorized. |
| Stage 7 | NOT RUN | No reconstructed packet exists to test for naturality. |
| Stage 8 | NOT APPLICABLE | The data are insufficient; no common source packet supports an inequivalent-solution analysis. |
| Stage 9 | NOT RUN | Post-reconstruction fixture comparison is forbidden without a reconstruction. |
| Lean | NOT REQUIRED | No stable new positive reconstruction theorem exists to formalize. |

## Files added

- `code/a5_closure/source_current_order_sensitive_inventory.py`
- `code/a5_closure/verify_source_current_order_sensitive_inventory.py`
- `code/a5_closure/manifests/source_current_order_sensitive_inventory.json`
- `code/a5_closure/tests/test_source_current_order_sensitive_inventory.py`
- `code/a5_closure/source_current_tomography_stage2.md`
- `docs/arithmon/SOURCE_CURRENT_TOMOGRAPHY_HANDOFF_2026-09-20.md`

## Files changed

- `code/a5_closure/source_current_tomography_stage1_contract.py`
- `code/a5_closure/tests/test_source_current_tomography_stage1_contract.py`
- `code/a5_closure/README.md`
- `claims/claim_registry.yaml`
- `claims/selection_ledger.json`
- `claims/physical_identification_registry.json`
- `docs/SELECTION_LEDGER.md` (generated)

## Validation performed

Passed:

- `python3 code/a5_closure/source_current_tomography_stage0.py verify`
- `python3 code/a5_closure/tests/test_source_current_tomography_stage0.py`
  — 6 tests.
- `python3 code/a5_closure/source_current_tomography_stage1_contract.py verify`
- `python3 code/a5_closure/tests/test_source_current_tomography_stage1_contract.py`
  — 7 tests.
- `python3 code/a5_closure/source_current_capability_certificate.py verify
  --projection code/a5_closure/manifests/source_current_capability_projection.json
  --receipt code/a5_closure/receipts/source_current_capability.receipt.json`
- `python3 code/a5_closure/verify_source_current_capability_independent.py
  --projection code/a5_closure/manifests/source_current_capability_projection.json
  --receipt code/a5_closure/receipts/source_current_capability.receipt.json`
- `python3 code/a5_closure/source_current_order_sensitive_inventory.py all`,
  followed by byte-for-byte comparison with the pre-replay manifest.
- `python3 code/a5_closure/source_current_order_sensitive_inventory.py verify
  --inventory code/a5_closure/manifests/source_current_order_sensitive_inventory.json`
- `python3 code/a5_closure/verify_source_current_order_sensitive_inventory.py
  --inventory code/a5_closure/manifests/source_current_order_sensitive_inventory.json`
- `python3 code/a5_closure/tests/test_source_current_order_sensitive_inventory.py`
  — 11 tests, including an audited-directory staleness mutation.
- `python3 code/a5_closure/issue_566_bracket_space_stage1/test_stage1.py`
  — 6 tests.
- `python3 code/a5_closure/issue_566_bracket_space_stage2/test_stage2.py`
  — 4 tests.
- `python3 tools/check_claim_registry.py`
- `python3 tools/build_selection_ledger.py` and `--check`
- JSON parsing for all three modified registries, Python byte compilation, and
  `git diff --check`.

Unavailable in the current environment:

- Focused pytest: `python3 -m pytest ...` fails before collection with
  `No module named pytest`.
- `python3 code/a5_closure/test_audit.py`: 6 of 9 tests pass; 3 error on import
  because `sympy` is absent.
- Mandatory standard: attempted and fails on its first step because `pytest`
  is absent.
- Certificate smoke: attempted and fails for the same missing module.
- Mandatory full and certificate suites: not run because the same prerequisite
  failure blocks their command lists.
- `lake build OPHScreen`, full `lake build`, and CI: not run; Lean changes were
  not required and CI is external.

No dependency was installed or duplicated during this run.

The initially proposed Stage 0/1/2 entries were removed from
`tools/run_mandatory_suite.py`: that file is a pinned `control_artifact` in the
invariant-mining pre-generation freeze. The freeze was neither regenerated nor
repinned. Mandatory-runner integration is intentionally deferred to a separate
upstream-owned custody change.

## Claim boundary

- Physical current source bridge attained: **no**.
- Conditional `PORT-CURRENT-INNER` fixture used as reconstruction oracle:
  **no**. It appears only as an explicitly rejected downstream candidate and
  as an unchanged comparison boundary.
- New downstream physical claim introduced: **no**.
- Non-identifiability: not established or resolved. The current failure occurs
  earlier, because the required common order-sensitive source packet is absent.
- Abstract forced Lie-type theorem and all existing conditional fixture results
  remain unchanged.

## Remaining debt and next justified step

The only positive continuation is a new upstream, independently reviewed
source producer containing twelve reversible families `P_0,...,P_11`, all 132
ordered mixed histories, canonical raw-history hashes, common W12 and
refinement ancestry, two-sided reversibility witnesses, and a strict downstream
target firewall. Only after that object exists may Stages 3–7 be reopened.

Environment debt is separate: restore the repository's supported `pytest` and
`sympy` environment, then rerun focused pytest, standard mandatory, full
mandatory, and certificate suites. No scientific conclusion above depends on
treating those unavailable runs as passing.
