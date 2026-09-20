# Source-family count/volume acceptance audit (#782)

The acceptance contract is the supplied golden population and complete-neighbour
read law, with coordinate volume `dt d³x` and layer duration `Delta = a/c`.
The results below prove family convergence, beyond the existing finite poset
facts and polynomial reference constants. No weak-measure convergence,
count limit, geometric pair-volume identity or integral limit is assumed.

## Objective and deliverables

| Requirement | Checked result | Inputs and scope |
| --- | --- | --- |
| Weighted Alexandrov volume with the explicit error | `SourceNetVolumeError.weighted_alexandrov_error`, `weighted_alexandrov_tendsto`; `GoldenSourceVolumeLimit.golden_weighted_error`, `alignedVolume_tendsto`, `alignedCount_tendsto` | A convex h-covered population, all-neighbour read radius a>2h, positive supplied c, finite measurable cells with almost-everywhere disjointness/coverage and assignment error H, and the interior buffer. General unequal-mass families converge under a,H,h/a -> 0 and converging aligned durations. The golden specialization derives its spatial facts from the actual fractional golden orbit. |
| Golden quadrature estimate | Existing `GoldenSourceAssignment.golden_quadrature`, `golden_quadrature_tendsto` | Actual golden equal-cell masses L³/q³ and assignment bound 2 sqrt(3)L/q. No replacement regular-grid population. |
| Null-boundary and moving interval counts | `SourceCausalBoundary.diamond_frontier_null`; `GoldenSourceCausalLimit.generated_interval_count_tendsto`; `FlatDiamondNormalization.generated_interval_count_tendsto_volume` | Actual sampled tips converge in coordinates. Read-path agreement follows from the covering construction. The first limit retains the clipped window; the full formula additionally requires a timelike limiting diamond contained in the window. |
| Actual strict pairs and ordering-fraction limit | `GoldenSourcePairLimit.generated_strict_pair_count_tendsto`; `FlatDiamondNormalization.ordering_fraction_tendsto_one_tenth` | Both interval predicates and pair precedence use the generated order. For a positive-volume timelike limiting diamond in the window, 2C/[N(N−1)] tends to 1/10. |
| Proof index and paper citations | `PROOF_INDEX.md`; `SOURCE_NET_CAUSAL_LIMIT.tex`; `SOURCE_POPULATION_QUADRATURE.tex`; `CAUSAL_MANIFOLD_OBSERVATIONS.tex` | The cited source propositions and formalization status match the theorem statements. |

The finite vertical error is exactly

`|V − pi c³T⁴/24| <= 4 pi (T+Delta)(cT/2+H)²(H+cTh/a) + pi c³T³Delta/2`,

with `T = K Delta`, including all K+1 layers. `layerMass_sandwich` derives
both bounding balls from generated paths and actual assigned-cell unions;
the per-layer error is not supplied as a hypothesis. The aligned golden
limit allows arbitrary layer counts K_n whose durations K_n Delta_n tend to
T, and an eventually valid interior buffer. It derives the inner-radius
condition from the vanishing golden assignment/read-radius ratio.
The general `weighted_alexandrov_tendsto` theorem also covers unequal cell
masses, varying populations and varying convex windows. It derives the
eventual inner-radius condition and the vanishing total error from
`a_n,H_n,h_n/a_n -> 0`; a volume/count limit is never supplied.

## Geometric integral and normalization

`FlatDiamondVolume` connects actual Euclidean sections to Lebesgue volume.
`FlatLorentzVolume` constructs the boost, checks its determinant and causal
preservation, proves measure preservation, and obtains translated causal
subdiamond volumes, including null tips. `FlatDiamondPairIntegral` uses
Fubini and radial integration to identify the actual pair measure with the
polynomial integral in `OrderingFractionFourDimensional`. Thus V²/20 is a
conclusion about the geometric set, not a renamed input constant.

`FlatDiamondNormalization` proves the Jacobian factors for arbitrary c>0
and their cancellation in the ratio. The strict and closed pair regions
have equal measure because the pairwise null relation is proved null.
The vanishing event weight controls the distinction between N² and
N(N−1); growth of N is not assumed.

## Audit corrections and boundary cases

- Fixed continuity-set counting is not substituted for generated intervals:
  eventual membership is proved from actual paths for moving sampled tips.
- Representatives need not belong to their cells. The existing two-cell-width
  assignment and its counterexample to a one-cell-width claim are retained.
- The final clipped time cell is retained, including aligned endpoints and
  T=0. Its single-count discrepancy is at most L³Delta. For arbitrary pair
  masks the discrepancy is at most L⁶Delta(2T+Delta).
- Coincident and null-separated pairs are handled through proved nullity.
  The geometric volume formula also covers null subdiamond tips.
- A clipped diamond is never assigned the full flat volume automatically.
  Full-volume and 1/10 statements explicitly require containment in
  `[0,T] × [0,L)^3`; the difference from `[0,T)` is proved measure-null.
- The finite spatial error separates covering radius h from assignment H.
  Zero-hop endpoints are proved before using the inner cone construction.

The separate general moving-tip *finite error inequality* in
`eq:source-net-general-volume-error` remains analytic. The requested explicit
vertical Alexandrov inequality, the moving-tip count convergence, and the
moving-tip ordering-fraction limit are proved. This distinction is stated in
the paper; no obstruction or unproved limit is used to satisfy #782.

## Exit and validation

`Geometry.SourceCountLimitAxiomAudit` imports the complete development and
checks every one of its 142 new public theorems, plus the existing golden
quadrature and site-injectivity interfaces. The `audit_source_axioms`
command fails unless the entire transitive axiom set is a subset of
`{propext, Classical.choice, Quot.sound}`, and prints that set. All 145
receipts are covered. Two expected-error controls verify rejection of
`sorryAx` and compiler-trusted reduction (`Lean.ofReduceBool` and
`Lean.trustCompiler`); the controls introduce no axioms. The pinned environment is Lean
4.29.1 / Mathlib revision `5e932f97dd25535344f80f9dd8da3aab83df0fe6`.

Replay from `Lean/`:

```sh
lake build Geometry.SourceCountLimitAxiomAudit
```

## Second audit corrections

Re-reading the issue and the paper exposed three omissions in the initial PR:

1. The general weighted proposition had a finite error theorem but its named
   convergence theorem only specialized to golden equal masses. The new
   `weighted_alexandrov_tendsto` closes that formal API gap.
2. The paper changes had not been reflected in the mathematical claim registry,
   which caused the required registry-review CI job to fail. The golden claim,
   its parent mathematical component, assumption dictionary, and novelty and
   falsification rows now match the proved limits. Physical premise
   classifications, reviewed dependency projection and physical statuses retain
   their existing values; no validation rule was relaxed.
3. The previously unfinished local paper build exposed four new underfull
   warnings in the citation paragraph. Reflow and ordinary break opportunities
   in long module names resolve them without changing the warning allowlist.

The local spacetime paper build and existing warning gate pass: four existing
allowed underfull boxes, no unexplained boxes, no overfull boxes and no
reference/citation/glyph/font problems. Claim-registry validation and all 46
registry gate tests pass. The theorem-count floor and release-manifest checks
also pass. These local results do not assert completion of hosted CI.

The population, read law, time scale and coordinate measure remain supplied.
Additional instrument/routing operations require their own measure comparison.
The result neither selects a physical clock/measure nor proves the M1
selection obligations in other issues. No claim-registry, frozen-prediction
or physical-observation verdict is promoted by these mathematical proofs.
