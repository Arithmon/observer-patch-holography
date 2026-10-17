# OL-A1 signature replication instrument

**INS-01 has status `FAILED` and remains the controlling completed instrument
for OL-A1.** The [instrument register](../registers/OBSERVATION_LEDGER_V3.md#ins-01-ol-a1-13-signature-replication-instrument)
owns its frozen decision rule, seed policy, external artifacts, verdict
receipts and reproducibility limits. This page summarizes that instrument;
it does not authorize another execution or alter its frozen protocol.

The [INS-02 factorial design](OL_A1_FACTORIAL_FOLLOWUP_DESIGN.md) is
`SPECIFIED` and unfrozen. It does not supersede INS-01 or restore OL-A1.

## Measurement and controls

The instrument compares three carrier/observer/support arms:
`16384/128/96`, `65536/256/96`, and the matched support-density control
`16384/128/6`, with five replicates each. Its observables are threshold
inertia, robust inertia at relative threshold `1e-3`, the degeneracy ratio,
cone margin and rung ratio, structural diagnostics, and split-half inertia
concordance. The registered controls test support density and ancestry
permutation. Full thresholds and the `REPLICATED`/`FAILED`/`INCONCLUSIVE`
decision rule remain in the frozen record.

## Scientific limits

- The `(1,3)` threshold verdict depends on an eigenvalue of relative magnitude
  about `3e-5` against a `1e-12` threshold. At relative threshold `1e-3`, the
  retained rungs have robust inertia `(1,2)` and one degenerate direction;
  a fresh seed can change `(1,3)` to `(2,2)` without changing the physical content.
- Inertia is fitted on the parity-zero half; only the cone margin uses the
  reserved half. Split-half inertia concordance is a separate observable.
  Strides `8/35/141` make the parity structure vary between rungs, so the
  split requires explicit stratification.
- Unequal coordinate scales make the signature fragile. Coordinate
  standardization resolves that numerical defect but does not make the
  ancestry control discriminate: the null reproduces data inertia in four
  of five permutations, and the reserved cone margin remains negative.
- The retained manifest, summary and fifteen per-arm receipts lack the raw
  feature matrices and fit captures required for independent observable
  reconstruction. Validate-only receipt replay is not that reconstruction.
  Five replicates resolve seed fragility only at the one-in-five level.

No successful replication or observed physical spacetime follows from this
instrument. A replacement needs a separately frozen, discriminating protocol
and an explicit change to the ledger's controlling instrument.
