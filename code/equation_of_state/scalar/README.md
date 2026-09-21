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
