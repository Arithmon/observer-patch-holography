# One source for dynamics and lensing

The [derivation](../../extra/DARK_SOURCE_DYNAMICS_LENSING.md) answers the
bounded dynamics/lensing target of #751. Even with the Einstein equation and
dominant energy, the released scalar anomalous-charge interface does not
fix a joint prediction. The certificate includes sharp flat-annulus stress
and lensing bounds, an exact joint inverse, and a regular finite-mass pair
with identical rotation and exterior but strictly different bending.

The global separation is at least `81/351232000` radians for the declared
dimensionless strength `epsilon=1/100`. This is a mathematical source pair,
not an observed galaxy or a microscopic OPH source realization. The finite
certificate also encloses the full bending difference on both sides, using
an exact rational logarithm bound at two source strengths. The finite
annulus also admits a continuum of nonnegative-pressure sources with the
same circular speed. Supplying the complete physical density and central
mass would remove the freedom; a fitted Newtonian mass is not that input.

These are classical continuum stress/metric models. They supply no new
bounded observer-patch implementation, quantum instrument or simulator
history. The public receipt is computational evidence for the mathematical
comparison, not an observer's physical readout.

From the repository root, with `PYTHONPATH=code` and pinned requirements:

```text
python -m dark_source_lensing.build
python -m dark_source_lensing.verify
python -m pytest -q code/dark_source_lensing
```

The producer uses closed-form pressures and ray integrals. The verifier
rebuilds the Einstein tensor identities and directly integrates a different
ray variable and the Abel inverse. Global admissibility and ray separation
are proved with exact rational bounds; decimal ray and transform controls
are numerical comparisons, not interval proofs. See the
[objective and exit](CONTRACT.md) for the physical boundary. No natural data,
new acceleration scale or fitted lensing amplitude enters this package.
