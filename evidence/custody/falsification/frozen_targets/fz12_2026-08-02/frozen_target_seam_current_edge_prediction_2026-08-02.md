# FZ-12 frozen target: source-seam edge propagation branch

Freeze time: 2026-08-02T11:52:27Z

Source repository: `FloatingPragma/observer-patch-holography`

Source commit: `bc5595f8dbb2d2886e2a64ddf447f69fbb00eb3f`

Canonical receipt SHA-256: `0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915`

## Scope

This target freezes a prospective conditional physical prediction of the
source-seam edge branch. The exact source theorem fixes the finite carrier
ray. It does not prove that the carrier ray is a physical kinetic operator.
Issue #666 owns that attachment.

The source result is exact. The canonical thirty-seam incidence current maps
onto the even-sum lattice

```text
D6 = {z in Z^6 : sum_i z_i is even}.
```

Its carrier differences have one common norm. Their complete projective
multiset is the thirty-direction edge orbit, including multiplicity. Their
normalized moments are

```text
M2 = 10 r^2,
M4 = 6 r^4,
M6 = (30/7) r^6 - (2/7) I6.
```

The conditional physical branch fixes the following premises before any
comparison:

1. The seam boundary acts as a displacement in the response-selected carrier.
2. The action extends homogeneously at every carrier point.
3. The complete thirty-direction edge orbit is the sole direct kinetic support
   through the displayed order. No second range, generator power, onsite
   correction, or independent isotropic term is added.
4. Source incidence multiplicity supplies equal weights on that orbit.
5. The dense `D6` action extends to the selected continuous field
   representation.
6. The tested sector carries the same scalar operator. A photon reading also
   requires the same operator on both transverse polarizations.
7. Local actions glue across scale and observer charts without changing the
   coefficient ray.
8. The carrier scale `a` is finite and strictly positive. A null verdict also
   requires a source-derived positive lower bound that makes the test
   sensitive to the branch.
9. One carrier frame, one orientation in `SO(3)/A5`, and one boost law are
   transported into the comparison frame. Source, medium, gravity,
   instrument, polarization, and orientation effects are isolated or
   profiled under a frozen nuisance model.
10. A time or clock bridge is required only if the spatial eigenvalue
    `Lambda` is interpreted as `omega^2`.

## Frozen prediction

For the thirty unit seam directions `w_j`, the branch fixes the spatial
kinetic eigenvalue

```text
Lambda_a(k,n) = (1/(5 a^2)) sum_{j=1}^{30}
                [1 - cos(a k w_j.n)].
```

Its long-wavelength expansion is

```text
Lambda_a = k^2 - (a^2/20) k^4 + (a^4/840) k^6
           - (a^4/12600) k^6 I6(n) + O(a^6 k^8).
```

Writing the displayed corrections as
`C4 k^4 + B0 k^6 + B6 k^6 I6`, the frozen relations are

```text
C4 = -a^2/20 < 0,
B0 =  a^4/840 > 0,
B6 = -a^4/12600 < 0,
B0/C4^2 = 10/21,
B6/C4^2 = -2/63,
B6/B0 = -1/15.
```

Intrinsic anisotropic coefficients at angular ranks one through five vanish.
The rank-six coefficient has one rotated `I6` shape. Once `C4` fixes the
scale, the remaining carrier freedom is one spatial orientation modulo the
proper icosahedral group.

## Prospective decision rule

The comparison remains unarmed. An eligible later release must provide a
joint likelihood or full covariance for same-sector `C4`, isotropic `B0`, and
the complete rank-six coefficient vector. Every physical premise above and a
dataset-specific nuisance, coverage, trials, sensitivity, and exposure
contract must be fixed before the data are opened.

- `FZ12-R01 FAIL`: isolated intrinsic `C4` is positive at five or more
  standard deviations.
- `FZ12-R02 FAIL`: an isolated intrinsic anisotropic coefficient at angular
  rank one through five is nonzero at five or more standard deviations.
- `FZ12-R03 FAIL`: after a negative `C4` trigger and adequate sixth-order
  sensitivity, the linked `B0`, negative `B6`, or rotated `I6` vector is
  excluded at five or more standard deviations with calibrated joint
  coverage.
- `FZ12-R04 FAIL`: the calibrated joint likelihood excludes the complete
  linked branch manifold at five or more standard deviations.
- `FZ12-R05 SUPPORT`: the zero-coefficient minimal locally Lorentz-invariant
  Standard Model plus General Relativity baseline is excluded at five or more
  standard deviations, the linked branch agrees within two standard
  deviations, named systematic alternatives are rejected, and an independent
  eligible release replicates the result.
- `INCONCLUSIVE`: every null, underpowered, incomplete-covariance,
  unresolved-frame, polarization-split, or non-isolated outcome. A null cannot
  decide the branch without the source-derived lower bound on `a`.

A fail rejects the seam-current edge propagation branch. It rejects OPH as a
whole only if a separate theorem proves that this branch is forced and
exclusive.

The minimal locally Lorentz-invariant Standard Model plus General Relativity
has no intrinsic local-vacuum coefficient of this form. Nonminimal effective
operators or another icosahedral medium can imitate the relation, so support
would distinguish this branch from the minimal baseline without identifying
OPH uniquely.

## Exposure boundary

FZ-11 remains the immutable primitive-vertex branch. It has
`B6/C4^2 = 32/315`, while this edge branch has `B6/C4^2 = -2/63`. The two
targets may not be pooled or selected after exposure.

The 2026-07-17 WMAP campaign, its CMB likelihood class, every data product
inspected for FZ-11, and all comparison data seen before this freeze are
ineligible for an FZ-12 verdict. No new target or comparison data were read
while this target was prepared.
