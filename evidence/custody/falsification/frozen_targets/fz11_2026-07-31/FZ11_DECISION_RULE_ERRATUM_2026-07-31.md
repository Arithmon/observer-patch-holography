# FZ-11 append-only decision-rule erratum

Created: 2026-07-31T19:02:42Z  
Original target: `frozen_target_spin_six_primitive_port_prediction_2026-07-31.md`  
Original target SHA-256: `1e6dfe17cb917a66f729ae21061704137da65e5269cf925ac3ab2f6ae754d737`  
Original prediction receipt SHA-256: `8ac97d7c46199717ed031610efdda65c40f6a251e78715d6bc05888d598e66d8`

## Reason for the erratum

The frozen equations predict more than the first written failure rule scored.
They fix `C4 < 0`, set every intrinsic anisotropic coefficient at angular ranks
one through five to zero, and link `B0`, `B6`, and the rotated `I6` vector to
`C4`. The original rule triggered only on a negative `C4` and then scored the
sixth-order terms. It therefore left a resolved positive `C4` and a resolved
lower-rank intrinsic anisotropy without a failure verdict, even though either
contradicts the frozen branch.

This erratum corrects the decision rule before any eligible comparison. It
does not alter the operator, coefficient values, ratios, angular template,
branch premises, exposure exclusions, or prediction receipt. No qualifying
comparison data were read while it was prepared.

## Operator-scope clarification

The frozen formula defines a real, reciprocal, finite-range cosine kinetic
branch. Through the displayed order, the complete primitive twelve-port orbit
is its only hop support and the operator contains no independent isotropic
fourth- or sixth-order term. Proper-carrier covariance fixes equal weights on
that orbit, and continuum normalization fixes their common coefficient. This
is an explicit branch definition. A derivation of this kinetic choice from the
repair architecture remains the open physical theorem owned by issue #655.

## Corrected prospective decision rule

Every verdict first requires a dataset-specific contract that fixes one
post-freeze release, the joint likelihood or full covariance, the physical
sector and carrier frame, the `SO(3)/A5` orientation profile, the boost law,
source, medium, gravitational and instrumental nuisance models, trials
accounting, the sensitivity floor, and calibrated joint coverage. The carrier
contribution must be isolated. A photon test also requires equal action on
both transverse polarizations.

### Fail

The primitive twelve-port physical propagation branch fails if any one of the
following is established at five or more standard deviations under the
eligible contract:

1. the isolated intrinsic `C4` coefficient is positive;
2. an isolated intrinsic anisotropic coefficient at angular rank one through
   five is nonzero;
3. for a resolved negative `C4` with the required sixth-order sensitivity,
   the linked `B0`, `B6`, or rigid rotated `I6` vector is excluded after the
   fixed orientation profile;
4. the calibrated joint likelihood excludes the complete branch manifold
   `C4 < 0`, `B0/C4^2 = 10/21`, `B6/C4^2 = 32/315`, and the rotated `I6`
   coefficient vector.

### Support

Support requires all of the following:

1. the zero-coefficient minimal locally Lorentz-invariant Standard Model plus
   General Relativity baseline is excluded at five or more standard
   deviations;
2. the complete linked branch manifold agrees within two standard deviations;
3. the named systematic alternatives are rejected;
4. an independent eligible release replicates the result.

### Inconclusive

A null result is inconclusive because the branch supplies no positive lower
bound on the carrier scale. Insufficient sixth-order sensitivity after a
negative `C4`, incomplete covariance, an unresolved frame, polarization
splitting, or failure to isolate the carrier contribution is also
inconclusive.

A fail rejects the primitive twelve-port physical propagation branch. It
rejects OPH as a whole only if issue #655 proves that the branch is forced and
exclusive.

## Exposure boundary

The 17 July 2026 WMAP internal-linear-combination template campaign, its
cosmic-microwave-background likelihood class, and every data product examined
there remain excluded. Its template-only null has no role in this rule. The
linked coefficient manifold has no qualifying physical comparison.
