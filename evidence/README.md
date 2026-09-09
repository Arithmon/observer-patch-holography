# Computational evidence archives

This directory contains scientific data and receipt packages cited by the
papers. Each package states its own replay and claim boundary.

- [`exact_federation_L6_canonical_20260909/`](exact_federation_L6_canonical_20260909/)
  is the level-six exact-federation run: 81,920 twelve-port carriers on the
  geodesic icosahedral tower under the canonical seam-mean law and the
  declared port gluing. It carries the primitive seam arrays, the initial and
  terminal loads, the per-schedule ledgers of `V = sum N_i^2`, the quotient
  hashes of sixteen schedules (one hash, equal to the component multiset),
  the isolated controls, the budgeted float schedule, the per-carrier
  response kernels, and a standalone verifier that imports no simulator code.
- [`closure_loop/`](closure_loop/) is the closure loop at carrier scale: the
  event logs of the canonical, integer, overwrite, tetrahedral and octahedral
  sources, the specifications recovered from those logs alone, the invariant
  vectors of the inhabited structure and the constructed realization, and a
  producer-free checker.
- [`source_net_causal_poset/`](source_net_causal_poset/) is the expected OPH
  causal poset: the source-record family receipts at `q <= 55`, the carrier
  realization receipts and event logs, a property table with the status of
  every property (proved in Lean, computed, or declared), a standard-library
  generator that rebuilds the poset from its definitions, and a checker.
- [`local_domain/`](local_domain/) is the complete current 11-file local-domain
  family used by the spacetime, particle, and cosmology papers. It contains
  the compressed stage-one arrays, all four staged receipts, the five
  downstream finite-interface receipts, the binding manifest, and a
  standard-library archive checker.
- [`source_causal_history_family/`](source_causal_history_family/) mirrors the
  simulator's 24-to-384-event source-history custody receipt with its
  theorem-level publication projection and checker. Its causal readout is
  superseded by the causal-poset package; the package stays as the custody
  input pinned by the particle-side source projections.
- [`particle_simulation_receipts/`](particle_simulation_receipts/) publishes
  the hash-bound report receipts for the calibration-null and direct
  permutation-transport assays cited by the particle paper. Their upstream
  raw arrays are not included.

Evidence directories are not interchangeable. A carrier count shared by two
packages does not make one run evidence for the observable tested by another.
