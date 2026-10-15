# Controlling occupation removal by the full path defect

The nearest result is the contractive a posteriori estimate of
[Etienney, Robin and Rouchon, Lemma 1 and Section 6](https://arxiv.org/html/2501.09607v6).
We reuse it for a piecewise polynomial approximation `Y`, including its
endpoint jumps. For a conservative reference propagator on Hermitian
trace-class operators,

```
||rho(T)-Y(T)||_1/2 <= ||rho(0)-Y(0)||_1/2
  + integral ||Y'(t)-L Y(t)||_1/2 dt
  + sum ||Y(next_start)-Y(previous_end)||_1/2.
```

This packet adds a specialized, outward-enclosed evaluation of that defect
for the actual supplied spin-oscillator pulse family. Positivity of the
polynomial path is unnecessary. Hermiticity is checked exactly. This avoids
assuming that sparse endpoint samples bound an intervening trajectory.

## Target, domains and the meaning of the integration parameter

The oscillator Hilbert space is `ell^2(N)` with finite-dimensional retained
spin and, during detection, pointer factors. During a sideband the target is

```
H_s = delta*N - s*Omega*X_spin*(a+a†)/2,
L rho = -i[H_s,rho] - gamma*[N,[N,rho]]/2.
```

Here the parameter is supplied drive elapsed duration. Carrier gaps have
their finite duration, commuting number drift and dephasing. Instantaneous
spin rotations, preparation displacement and detector interaction use
normalized operation parameters on `[0,1]`; these are not additional claimed
physical dwell times. Published nominal duration remains `4*k*tau`.

There is a concrete conservative unrestricted propagator. The Hermitian
number noise has phase-diffusion unraveling with phase
`delta*t+sqrt(gamma)*W_t`. In each spin-X eigensector, going to this rotating
frame leaves a linear oscillator force with a continuous, bounded complex
coefficient. Its propagator is a Weyl displacement times a scalar phase;
the displacement magnitude is bounded by the integral of the force
magnitude. Each noise realization therefore gives a unitary, and its
expectation is trace preserving and completely positive. On finite Fock
matrix units the Ito calculation gives the displayed generator. Finite
spin rotations preserve the oscillator domain; concatenating finitely many
pulses gives the time-dependent propagator used here. Preparation and pointer
coupling have the same linear-force form. The initial thermal state has all
polynomial number moments finite. These facts supply the conservative
propagator and generator domain needed for the finite-support defect
calculation; finite-dimensional
complete positivity alone would not supply them.

All numeric apparatus parameters denote the exact binary rational values
recorded in `results.json`, including the stored values computed using
`math.pi`. They are not exact symbolic multiples of mathematical pi.
Carrier gap duration is the declared binary value `abs(phi)/omega_carrier`.
Products with that duration are then exact rational target coefficients.
Finite ladder square roots are enclosed rather than rounded into the target.
No apparatus uncertainties are included.

## Occupation boundary and polynomial integration

Let `P_d` retain occupations `0,...,d-1`. Number drift, number dephasing and
spin-only operations preserve this subspace. A linear force has only one
crossing edge, with `B=(I-P_d)HP_d`. For a Hermitian finite matrix `Y`,

```
||(L-L_d)Y||_1/2 = |force|*sqrt(d)*||R_d Y||_1,
R_d = I_sectors tensor <d-1|.
```

This identity holds for every Hermitian `Y`, including nonpositive polynomial
coefficients. The spin factor multiplying the force is unitary (`I`, `X`,
or the pointer `Z`) and does not change the norm. We use the always valid
upper estimate `||A||_1 <= sqrt(rank_bound)*||A||_F`. The boundary row has
at most two rows before detection and four during detection.

On a normalized cell write `Y(s)=sum_{k=0}^m c_k*s^k`, `0<=s<=1`, with
dimensionless generator `G=h*L_d`. Every stored coefficient is interpreted
as its exact binary matrix. The finite defect coefficients are

```
b_k=(k+1)*c_(k+1)-G*c_k     (k<m),
b_m=-G*c_m.
```

Triangle inequalities give integrated bounds

```
finite defect <= (1/2) sum_k ||b_k||_1/(k+1),
occupation defect <= h*|force|*sqrt(d) sum_k ||R_d*c_k||_1/(k+1).
```

Both cover the entire continuous polynomial. Taylor coefficients are merely
a way to propose a small-defect path; no claim that a floating Taylor step
is exact is made. The next stored endpoint may differ from the exact sum of
its coefficients, so that jump is separately enclosed and charged.

The interval kernel widens every elementary addition, subtraction,
multiplication, division used in a bound and square root. Positive reductions
use outward pairwise sums. Sparse Hamiltonian actions use no BLAS kernel in
the certificate. Numerical coefficient generation may use ordinary floating
arithmetic because the defect of the resulting binary coefficients is then
enclosed independently. The arithmetic contract is IEEE-754 binary64
round-to-nearest, gradual underflow and correctly rounded elementary
operations/square root; overflow and nonfinite norm bounds fail closed.
This is an executable numerical enclosure, not a Lean-checked proof of the
Python runtime or a certificate for platforms violating that contract.

## Preparation, detector and terminal selection

For exact thermal ratio `r=nbar/(1+nbar)`, occupations below `d-1` have
weights `(1-r)*r^n`; the last retained occupation gets `r^(d-1)`. The exact
half-trace distance to the unrestricted thermal state is `r^d`. Exact
rational comparisons additionally charge rounding of the stored weights.
Preparation applies `D(alpha)` through its own polynomial path and boundary
budget. Thus the displaced thermal tail is not assumed to equal `r^d`.

Detection starts a pointer in `|+>`, applies opposite half-displacements
`exp(+-i*h*(a+a†)/2)`, and reads `-Y_pointer`. Its off-diagonal trace gives
the imaginary characteristic function at literal probe `beta=i*h`. The
detector's occupation, time-integration and arithmetic defects are included;
the rejected spin branch is retained. This is the supplied instrument's
semantic dilation, not execution of the native compiler's source words.

The final approximate matrix is Hermitian but is not assumed positive or
exactly normalized. If its exact accepted mass is `q`, and the complete
half-trace error is `epsilon`, then the reference accepted mass is positive
when `q>2*epsilon`. The conservative conditional probability error is
`2*epsilon/q`; mean error is `4*epsilon/q`. For `q<=2*epsilon`, the packet
returns no conditional certificate. Approximate detector probabilities must
be nonnegative before a probability error is capped at one. Exact rational
readout sums, final ratio rounding and the outward-enclosed `sqrt(2)*h`
denominator account for scalar evaluation error as well.

For later composition, the same polynomial path also bounds the exactly
truncated finite apparatus: omit the occupation boundary and thermal tail.
Adding that finite-reference bound to the unrestricted-reference bound gives
`unrestricted_vs_exact_finite_half_trace_error_upper`. This is a controlled
occupation-removal comparison that can be combined with a separately valid
native splitting/compiler budget. It is not itself a certificate of the
frozen native compiler's floating execution or physical hardware.
