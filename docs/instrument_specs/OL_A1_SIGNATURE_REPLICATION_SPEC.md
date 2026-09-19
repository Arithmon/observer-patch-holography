# OL-A1 signature replication instrument: scouted specification

> **HISTORICAL SCOUT.** Scouted 2026-08-12 for V3 issue #737 against #728's
> OL-A1 row. INS-01 was subsequently preregistered and executed; this draft no
> longer describes current run state and authorizes no execution. The open,
> unexecuted factorial follow-up design is
> [OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md](OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md),
> with live work state owned by issue #737.

## Load-bearing findings from the scout

- F1: the archived "held-out inertia" is fitted on the parity-0 half and
  its eigenvalues give the inertia; only the cone margin uses the
  held-out half. The instrument adds a split-half inertia concordance
  readout (fit on the complementary half, compare inertia).
- F2 (critical): the (1,3) verdict rides on one eigenvalue at relative
  magnitude about 3e-5 (threshold 1e-12). At a relative threshold
  1e-3 * max|lambda| the robust inertia of the retained rungs is (1,2)
  with one degenerate direction, matching the independently measured
  (1,2) receipt in code/geometry/realized_event_receipts.py. A fresh
  seed can flip (1,3) to (2,2) with no change in physical content. The
  observation-ledger row and any public sentence must carry this.
- F3: the held-out split is stride-structured per rung (strides 8/35/141),
  so train/held-out parity structure varies between replicates; record
  and stratify.

## Instrument (Tier A, mandatory, laptop-scale)

Arms: A1 = 16384/128/96, A2 = 65536/256/96, C1 = 16384/128/6 (matched
support density control). Five replicates each via one fresh seed and
five declared replicate_ids. Observables: O1 threshold inertia; O2
robust inertia (tau = 1e-3 relative), reference (1,2)/(1,2)/(2,1); O3
degeneracy ratio (reference band 2e-5..4e-5); O4 cone margin and rung
ratio (retained 0.5725, declared band 0.35..0.80); O5 structural
diagnostics; O6 split-half concordance. Controls: C-SUPPORT (support-6
density match; fired = mean signature degradation plus halved cross
edges) and C-ANCESTRY (ancestry-permutation null must destroy (1,3) on
at least 4 of 5 replicates; zero extra runtime). Decision rule
REPLICATED/FAILED/INCONCLUSIVE as scouted, frozen before execution;
FAILED demotes the OL-A1 emergent rung with equal prominence.

## Blockers before freeze

1. New driver script in oph-physics-sim (the ladder script hard-codes
   seed 20260751 and requires exactly four rungs); pinned oph_fpe/ tree
   stays read-only.
2. Reconcile the evidence/einstein_convergence schema v1 (RER) vs v2
   (sim) manifest discrepancy.
3. The verdict needs a machine-readable OL-A1 row to write back to
   (issue #726 artifact) and an instrument register separate from the
   frozen-prediction ladder (the ladder is reserved for physical
   predictions against external data).

## Effort

S: Tier A about one day of work plus about 16 minutes of compute and
under 1 GB memory. M: Tier B replicates the original 262k support
96/384 control (about 2 hours compute, 8 GB). Runtime environment must
be single and recorded; N=5 resolves seed fragility at the 1-in-5
level only, and the report must say so.
