# FZ-11 frozen target: primitive twelve-port propagation prediction

Freeze time: 2026-07-31T18:13:34Z  
Source repository: `FloatingPragma/observer-patch-holography`  
Source commit: `66176656dc1143f9ec50ba1a6e409c403545857f`  
Canonical receipt SHA-256: `8ac97d7c46199717ed031610efdda65c40f6a251e78715d6bc05888d598e66d8`

## Scope

This is a prospective conditional physical prediction of the OPH primitive
twelve-port propagation branch. It is not a claim that A1 through A3 have
proved the physical photon bridge. Issue #655 owns the stronger derivation or
rejection of that bridge.

The branch premises are fixed:

1. The complete primitive twelve-port orbit is the sole intrinsic directional
   support through the displayed derivative order.
2. Proper-carrier covariance acts transitively on the twelve ports, which
   forces one equal coefficient on the orbit.
3. The continuum quadratic term is normalized to `k^2`.
4. The carrier scale `a` is finite and strictly positive.
5. The tested physical sector realizes the scalar symbol. A photon test also
   requires equal action on both transverse polarizations.
6. One carrier rest frame is transported coherently over the experiment. Its
   constant orientation in `SO(3)/A5` is profiled.
7. Source, medium, gravitational, instrumental and other intrinsic terms are
   modeled well enough to isolate the carrier coefficient vector.

“Spin six” means spherical-harmonic rank `j = 6`. It does not mean particle
spin.

## Frozen prediction

For the twelve unit carrier directions `u_i`, the branch fixes

```text
omega^2(k,n) = (1/(2 a^2)) sum_{i=1}^{12} [1 - cos(a k u_i.n)].
```

Its long-wavelength expansion is

```text
omega^2 = k^2 - (a^2/20) k^4 + (a^4/840) k^6
          + (2 a^4/7875) k^6 I6(n) + O(a^6 k^8),
```

where

```text
I6(n) = (25/132) sum_i P6(u_i.n)
```

is normalized to one on the twelve vertex directions. Write the three
displayed corrections as `C4 k^4 + B0 k^6 + B6 k^6 I6`. The frozen relations
are

```text
C4 = -a^2/20 < 0,
B0 =  a^4/840 > 0,
B6 =  2 a^4/7875 > 0,
B6/C4^2 = 32/315,
B0/C4^2 = 10/21,
B6/B0 = 16/75.
```

All intrinsic anisotropic coefficients at angular ranks `j = 1,...,5` vanish.
The `j = 6` vector has the unique rotated `I6` shape. Once `C4` is measured,
the scale and both sixth-order amplitudes are fixed. The only fitted carrier
freedom is one spatial orientation modulo the proper icosahedral group.
Binary refinement gives `B6(a/2) = B6(a)/16` at fixed physical momentum.

## Prospective decision rule

The physical comparison remains unarmed until a dataset-specific contract is
timestamped. An eligible release must provide a joint likelihood or complete
covariance for same-sector `C4`, isotropic `B0`, and the `j = 6` coefficient
vector, together with declared source, medium, gravity, instrument, boost and
orientation treatment.

- Trigger: `C4` is negative at five or more standard deviations and the
  experiment has enough sensitivity to resolve the two linked sixth-order
  terms.
- Fail: after the declared `SO(3)/A5` orientation profile, the linked `B0` and
  `B6` prediction is excluded at five or more standard deviations with
  calibrated joint coverage.
- Support: the zero-coefficient baseline is excluded at five or more standard
  deviations, the linked prediction agrees within two standard deviations,
  named systematic alternatives are rejected, and an independent release
  replicates the result.
- Inconclusive: `C4` is not detected, the linked terms are below sensitivity,
  the covariance is incomplete, or the carrier contribution cannot be
  isolated.

A fail rejects the primitive twelve-port physical propagation branch. It
rejects OPH as a whole only if the open #655 derivation proves that branch is
forced and exclusive.

The comparison baseline is minimal locally Lorentz-invariant Standard Model
physics plus General Relativity in local vacuum, for which the three intrinsic
coefficients vanish. A nonminimal effective theory or another icosahedral
medium can imitate the pattern, so a match distinguishes the branch from that
baseline without identifying OPH uniquely.

## Prior related exposure

The same `A5` angular template was used in a preregistered WMAP ILC sky search
on 2026-07-17. That campaign returned a null family-wide `p = 0.64`. Its
preregistration SHA-256 was
`2b83a001f75d4aa9f5a631b50d0fe8ad51950ae63146118ae369e8cfb80e84b2`.
The WMAP search, its CMB likelihood class, and every data product inspected in
that campaign are excluded from a new FZ-11 verdict. The linked `C4`, `B0` and
`B6` coefficient relations frozen here have not been compared with a
qualifying physical dataset.

No new comparison data were read while this target was prepared.
