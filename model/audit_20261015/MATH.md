# Mathematical audit of the executed anharmonic model

Audit scope: the frozen `anharmonic_transfer` packet and the native primitives
it imports, at commit `93e59bfb`. This is an independent audit of the specified
finite computation and of useful theorem attachments. It does not replace any
retained run or change its prediction protocol.

The finite environment, cutoff-aware splitting estimate and terminal
conditioning estimate survive the algebraic checks below. The most useful
missing attachment is an existing a posteriori occupation-truncation theorem.
Two native words can also be shortened exactly. An attractive forced-oscillator
factorization needs particular care: it is exact before occupation truncation,
but false as an identity between the finite matrices now executed.

## 1. Checks of the existing guarantees

**Gaussian phase diffusion.** For the specified Hermitian noise operator
`sqrt(gamma)*N`, the exact channel is the expectation of conjugation by
`exp(-i*sqrt(gamma*t)*Z*N)`, where `Z` is standard normal. The positive
Gauss–Hermite approximation therefore is CPTP, including with retained spin
and reference systems. Hermite interpolation gives the claimed scalar
remainder `M!*(gamma*t*(n-m)^2)^M/(2*M)!`. Positivity gives
`||rho_nm||_1 <= sqrt(p_n*p_m)` and the sum is at most `d`, so
`environment.py:44–75` supplies a valid ideal half-diamond bound. No extra
factor from the spin dimension is required. Diagonal entries have zero error,
so the factor `d` can be sharpened to `d-1`; this would have no practical
effect on the current tiny remainder.

**Strang splitting.** `bounds.py:44–101` correctly uses positive-time CPTP
contractions. The double-commutator coefficients `1/24` and `1/12` match the
outer-noise, inner-Hamiltonian order. This is an application of the bounded
contraction case of [Iserles and Kropielnicka, Theorem 2.3 and §3.3](https://arxiv.org/html/2401.06635v1).
The finite relation `[a,a†]=I-d|d-1><d-1|` is retained. In particular the term
`Omega^2*d/2` in the bound for the commutator with `L=-i[H,-i[H,N]]` is real;
replacing it by the infinite-oscillator scalar would invalidate a uniform
finite-dimension estimate. The existing counterexample tests the dangerous
high/low occupation coherence, not just a low-energy Gaussian input.

**Conditioning.** `bounds.py:125–154` correctly bounds conditional TV by
`epsilon/q` when the rejected branch is retained. For accepted substates
`P,Q`, their masses `p,q`, and complete trace distance at most `epsilon`,

```
||P-Q||_1 + |p-q| <= 2 epsilon,
||P/p-Q/q||_1 / 2 <= epsilon/q.
```

This can be sharp and does not require an extra factor two. The code also
correctly identifies `q<=epsilon` as insufficient to guarantee nonzero
reference acceptance. These are mathematical statements about ideal
channels; numerical acceptance values are not certified lower endpoints.

**Present limitations are genuine.** The certificate excludes compiler
roundoff, node/weight rounding, preparation and detector approximation,
occupation removal and apparatus uncertainty. The packet discloses these
exclusions. A tiny ideal quadrature remainder cannot replace any of those
missing bounds. No empirical residual is explained away by this audit.

## 2. Two exact native-word simplifications

At `source.py:53–56`, free number drift uses the complete Walsh expansion.
For `d=2^b`, however,

```
N = sum_r 2^r (I-Z_r)/2,
exp(-i*u*N) = exp(-i*u*(d-1)*I/2) product_r exp(i*u*2^r*Z_r/2).
```

The environment already uses this identity. It also applies to carrier-gap
drift. At dimension 64 the current generic word has **838 native pulses**;
the exact specialized word has **7**, with the same full global phase.
The supplied angle, clock and apparatus law are unchanged.

The imported `anharmonic_observations/source.py:75–79` similarly compiles
controlled oscillator parity as a generic diagonal. Its action is simply
controlled `Z` between the spin bit and the oscillator's least significant
bit. The existing native `J(-pi)` pulse realizes minus controlled `Z`; adding
`I(pi)` supplies the scalar phase. At dimension 64 this takes **2 pulses**
instead of **2,054**. Both reductions were checked on the complete source
matrix in [test_math_probes.py](test_math_probes.py), including high occupations.

These are prospective compilation improvements. Changing the words changes
their hashes and resource receipts, so adopt them in a new packet rather than
rewriting the frozen evidence. They reduce declared service applications as
well as host work, but leave the much larger dense basis words in place.

## 3. A directly applicable occupation-truncation theorem

[Etienney, Robin and Rouchon, Lemma 1, Remark 2 and §3](https://arxiv.org/html/2501.09607v6)
give the a posteriori trace-norm estimate

```
||rho(T)-rho_d(T)||_1
 <= ||rho(0)-rho_d(0)||_1
    + integral_0^T ||(L-L_d)rho_d(t)||_1 dt.
```

It applies to a conservative Lindblad propagator and its Galerkin truncation;
the time-dependent version permits the piecewise pulse schedule. It needs the
whole truncated trajectory and a controlled integral, not agreement between
two cutoff runs. The paper provides polynomial-boson residual constructions.

For this particular model there is an especially simple application. Let
`P_d` select occupations `0,...,d-1`, with spin/reference identities implicit,
and embed the finite state `rho_d` into the unrestricted space. Number
dephasing, harmonic drift and spin-only rotations preserve `P_d`, so their
truncation residual vanishes. During a sideband only the edge
`|d-1> <-> |d>` crosses the boundary:

```
B = (I-P_d) H P_d = -s*Omega*sqrt(d)/2 * X ⊗ |d><d-1|,
(L-L_d)rho_d = -i*(B rho_d-rho_d B†).
```

The two terms are opposite off-diagonal blocks. Define
`R_d=I_spin,reference ⊗ <d-1|`. Their singular values give the exact formula

```
||(L-L_d)rho_d||_1 / 2
 = Omega*sqrt(d)/2 * ||R_d rho_d||_1
 <= Omega*sqrt(d)/2 * sqrt(tr(R_d rho_d R_d†)).
```

The bound holds for mixed states and arbitrary retained reference systems.
Only one extra occupation level is needed to evaluate the exact residual;
the probe independently compares the full generator difference with this
boundary-row expression. The parameter integrated here is supplied drive
elapsed duration, including its piecewise schedule; it is not host runtime.

This supplies a practical route to the missing occupation certificate, with
four remaining obligations:

1. Establish a conservative, well-defined unrestricted propagator for the
   chosen linear force plus number dephasing, including domains. The explicit
   phase-noise unraveling and forced-oscillator solution offer such a route;
   finite-dimensional complete positivity alone is not that proof.
2. Bound the residual integral along the continuous finite dynamics. Sparse
   endpoint samples or an unvalidated quadrature are only diagnostics. The
   existing splitting bound can help enclose numerical trajectories, or a
   validated time integrator can directly enclose the residual integral.
3. Include initial-state error. The current state is a normalized finite
   thermal mixture acted on by a finite displacement. It is not automatically
   the compression of an unrestricted displaced thermal state. The
   undisplaced thermal tail is `(nbar/(1+nbar))^d`, but displacement needs its
   own comparison.
4. Include readout error. The present finite displacement/characteristic
   instrument must also be compared with the unrestricted instrument before
   applying the existing terminal-conditioning estimate.

The population-only bound is convenient but can be looser than the two-row
trace norm. A global first-moment estimate is less promising: the elementary
bound `sqrt(<N>(t)) <= sqrt(<N>(0))+Omega*t/2` grows to order `10^2` over the
longest drive duration, so it cannot justify a 32- or 64-level cutoff here.
Actual boundary occupancy and its integral are the relevant quantities.

## 4. Exact finite dephasing is possible, but not automatically cheaper

The matrix `C_nm=exp[-gamma*t*(n-m)^2/2]` is positive semidefinite Toeplitz,
with diagonal one. The classical Vandermonde decomposition therefore gives
positive weights summing to one and at most `d` phases with

```
C_nm = sum_j w_j exp[-i*phi_j*(n-m)].
```

This is [Yang and Xie, Theorem 1](https://arxiv.org/html/1605.02431), applied to
the Gaussian moment matrix. Each phase is already executable by the same
`b+1` native pulses. For positive variance the matrix is strictly positive
definite, so an exact representation needs `d` rank-one phase atoms; at zero
variance one atom suffices. A general Schur-channel Kraus factorization also
exists, but its Kraus factors need not be unitary and would require a different
execution construction.

The existing 16-node rule is approximate and cheaper to enumerate at `d=32`
or `64`. Since its ideal error is already far below floating uncertainty,
an exact finite decomposition is a conceptual alternative rather than a
priority optimization. Computing its near-singular moment decomposition
would introduce a numerical conditioning problem to certify.

Even 16 nodes are unnecessary for the current splitting budget. At dimension
64, 256 subdivisions and 50 steps, a declared 4-node rule for sideband half
channels and 6-node rule for carrier gaps has a complete ideal environment
radius below `6.973e-6`, versus splitting radius `3.069e-4`. This is an
outcome-independent bound from the existing formula, not a rerun or improved
empirical result. Fewer nodes reduce all-branch host enumeration; the native
pulse count of each realized phase branch remains seven.

Do not combine adjacent finite half channels by silently replacing
`Q_M(t/2)^2` with `Q_M(t)`: Gaussian diffusion has a semigroup law, but a
finite Gaussian quadrature rule generally does not. Exact branch convolution
or a separately budgeted replacement would be required.

## 5. The forced-oscillator simplification and its cutoff trap

In a spin-`X` eigensector, write `beta=s*x*Omega/(2*delta)`, `x=+1` or `-1`.
For the unrestricted oscillator, completing the square and using the Weyl
relations gives

```
exp[-i*t*(delta*N-s*x*Omega*(a+a†)/2)]
 = exp[i*beta^2*(delta*t-sin(delta*t))]
   D(beta*(1-exp(-i*delta*t))) exp(-i*delta*t*N).
```

This is useful as an independent analytic oracle for individual sidebands,
and potentially as an alternate source compilation with a new occupation
comparison. The spin can be retained through its two conditional sectors.
Spin rotations between sidebands generally create superpositions of those
sectors, so the entire multi-step experiment is not one Gaussian orbit.

Substituting truncated `a_d` into the right-hand side does **not** produce
the exponential of the current finite Hamiltonian. The probe gives a concrete
dimension-eight counterexample: low-occupation action can agree to under
`1e-8` while full-operator disagreement exceeds `1e-2`. The finite CCR
boundary accounts for the difference. Moreover the current native
displacement itself uses dense basis words, so a large resource reduction
does not follow just from rewriting the formula. Keep this as an independent
comparison route until its source compilation and truncation budget are
actually demonstrated.

## Reproduction

From the repository root:

```bash
python -m pytest -q model/audit_20261015/test_math_probes.py
python -m pytest -q model/anharmonic_transfer/test_environment.py model/anharmonic_transfer/test_bounds.py
```

The audit probes use complete native matrices for the exact word reductions,
independent generator construction for the boundary residual, and an explicit
counterexample for the tempting cutoff identity. They neither load empirical
outcomes nor modify the frozen runtime or its receipts.
