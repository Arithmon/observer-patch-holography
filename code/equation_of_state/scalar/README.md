# Scalar pressure readouts

Run from the scientific repository root:

```bash
python3 code/equation_of_state/scalar/scalar_eos.py
python3 code/equation_of_state/scalar/verify_scalar_eos.py --out code/equation_of_state/scalar/scalar_eos_verification.json
python3 code/equation_of_state/scalar/continuous_pressure.py
python3 code/equation_of_state/scalar/verify_continuous_pressure.py
python3 -m pytest -q code/equation_of_state/scalar
```

`scalar_eos_receipt.json` authenticates the existing 64-site golden scalar
execution and retains all 88 field-state readouts (four histories, 22 times).
For original field coordinates `q` and canonical momenta `p=Mv`, the supplied
three-dimensional dilation is

```text
H(s) = K/s^3 + s G + s^3 U,    V(s) = s^3,
p V = K - G/3 - U,            rho = (K+G+U)/V.
```

Here `s` dilates the full reference Dirichlet box relative to fixed length,
time and energy units. The action mass stays fixed in those units. The sum
of interior dual masses is not the box volume. Boundary contributions remain
in the gradient stiffness. The native scalar-mean simulator supplies none of
these identifications.

The canonical velocity decoded at step `n` is
`(q_n-q_(n-1))/tau - tau*A*q_n/2`. The split integrator conserves the modified
energy `H - tau^2 ||Aq||_M^2/8`; it does not exactly conserve the original
action energy `H=K+G+U`. The JSON retains both energies and signed pressure.
Transient `p/rho` can exceed `1/3` or be negative. It is not an equilibrium EoS.

`continuous_pressure_receipt.json` independently encloses continuous evolution
of the same action at all 22 recorded model times, using exact 81-term
operator Taylor polynomials and rational spectral/norm tail bounds. It retains
every polynomial state and the signed difference from the split pressure.
This is a mathematical reference calculation, not a new observed history.
The verifier uses a separate quadratic-field representation and Horner
evaluation. Outward rational pressure intervals include both polynomial
rounding and truncation error; no floating eigensolver certifies them.

The separate `thermal_cases` array supplies canonical Bose equilibrium,
zero chemical potential and `hbar=kB=c=1`. It retains all 64 modes at each of
36 declared `(s,m,T)` combinations. The original action has `m=1`; `m=0,4`
are explicitly supplied diagnostic variants. For spatial eigenvalues
`lambda_j`, it evaluates

```text
omega_j^2 = m^2 + lambda_j/s^2,
E_j = omega_j / (exp(omega_j/T)-1),
F = T sum_j log(1-exp(-omega_j/T)),
w = sum_j [E_j lambda_j/(s^2 omega_j^2)] / (3 sum_j E_j).
```

Thermal equipartition between kinetic and potential oscillator energies gives
`K_j = G_j + U_j`. Therefore the same mechanical formula simplifies to
`pV = 2 sum_j G_j / 3`. The evaluator uses this positive sum; subtracting
`K_j - G_j/3 - U_j` after rounding can erase a small massive-mode pressure.
This simplification applies to the supplied thermal ensemble. The signed
instantaneous pressure of the recorded classical states retains its full
formula.

The receipt includes oscillator energies, occupations and stress; the
verification report retains raw nearby-volume free energies. Independent verification reconstructs the
one-dimensional golden operator at 70 digits, forms its tensor spectrum and
differentiates Helmholtz free energy. Those thermal checks are numerical;
they are distinct from the exact transient and rigorous continuous enclosures.

The thermal prescription subtracts zero-point energy at every volume. It
predicts no cosmological vacuum or Casimir stress. Fixing 64 oscillators does
not bound their occupation or Hilbert-space dimension. No finite observer
capacity, thermalization, preparation law, physical clock calibration,
thermodynamic limit or continuum cutoff removal is supplied by this package.
At a fixed finite cutoff, neither the massive low-temperature limit `w=0`
nor the massive high-temperature limit `w=1/3` holds in general.

## Numerical contract and retained audit controls

`thermal_point` accepts a nonempty one-dimensional real spatial spectrum,
nonnegative eigenvalues and mass, and positive scale and temperature. Every
frequency must be positive: a massive zero-gradient mode is allowed, while a
massless zero-frequency oscillator has no normalizable Bose Gibbs state.
Inputs must be finite and exactly representable in the binary64 interface;
masked entries, Booleans, complex data and lossy conversions are refused.
Mixed sequences are checked element by element before NumPy chooses a common
type: `[2**53+1, 1.0]` must be refused, while `[2**53, 1.0]` remains valid.
Checking only the coerced array would miss the first sequence's lost integer
precision. Missing elements are refused before conversion can replace them
with `NaN`.

All intermediate mode quantities and sums use a private 90-digit context.
The free-energy logarithm uses complementary `expm1`/`log1p` expressions so
neither small `omega/T` nor exponentially small free energy is rounded away.
Every reported nonzero component must survive conversion to binary64 with
relative error at most `1e-12`, including occupations and individual `K,G,U`.
Resolved subnormal values remain available. An unrepresentable or insufficiently
resolved component raises `ValueError` naming the field; a positive value is
never replaced by zero, infinity or an occupation cutoff. These checks do not
constitute rigorous bounds on the high-precision transcendental evaluation.

The audit under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
reproduced these defects on main `de60560b`:

| Supplied one-mode input `(lambda,s,m,T)` | Old output | Independent result |
| --- | --- | --- |
| `(1,1,1e8,1e8)` | Pressure `0` | Pressure approximately `1.939922356231088e-9` |
| `(1,1,0,0.02)` | Free energy `0` | Free energy approximately `-3.857499695927840e-24` |

The old thermal verifier also accepted replacing all 377 positive entries
below `1e-11` in the retained occupation/energy/`K,G,U` arrays by `-1e-12`.
Its absolute tolerance concealed those corruptions. Thermal verification now
checks each nonzero value with relative tolerance `3e-10` against the independent
70-digit golden spectrum calculation, and checks structural zeros exactly.
The tolerance includes the supplied binary64 eigenspectrum's error; it is
not a permission to drop a small positive mode. The original receipt passes
these stronger numerical checks, so the retained 36-case sweep's scientific
conclusions do not change.

[`test_scalar_thermal.py`](test_scalar_thermal.py) independently differentiates
the oscillator partition function from the original inputs for energy,
pressure and mass-potential energy. It checks each component, changes energy
units, exercises normal/subnormal precision boundaries, and rejects zero,
negative and rescaled replacements of retained positive evidence. The live
scalar receipt and its verification report are regenerated by the commands
above. Parent execution, all 88 exact classical readouts, continuous-pressure
enclosures, frozen registrations, paper statements and claim payloads keep
their existing values and scope.
