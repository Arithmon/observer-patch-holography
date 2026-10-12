# Occupation refinement of the interacting retained history

The finite action in `scalar_coupling.py` is exactly the rectangular
occupation compression of an infinite two-oscillator Hamiltonian. Write
`P_d=1_{N0<d,N1<d}` and `N=N0+N1` in the original unit-frequency Fock basis.
On the four-dimensional matter/Gauss code the action is

```text
H = H_matter + (1+kappa/2)(N+1)
    + (kappa/4) sum_i (a_i^2+a_i*^2)
    - (kappa/2)(a_0+a_0*)(a_1+a_1*)
    + g(a_0+a_0*)B,
kappa=9/32, |g|<=1/4, B=n_L0-n_R0, ||B||=1.
```

The numerical finite onsite square correction is essential: it makes the
finite action `P_d H P_d`, rather than the square of a truncated canonical
coordinate. The scalar quadratic form is positive; the bounded matter
term and the relatively form-bounded linear force define its selfadjoint
interacting extension. The initial scalar state is the existing coupled
free vacuum, not a claimed interacting ground state. Matter can start in
any state entangled with an untouched reference.

## A number moment that does not assume Gaussian interaction dynamics

Passively rotate to common and difference annihilators. This leaves `N`
unchanged. The only number-changing quadratic term becomes
`(kappa/2)(a_-^2+a_-*^2)`; the remaining quadratic and matter terms commute
with `N`. Put `M(s,t)=<exp(sN)>` and `phi=exp(sN/2)psi`.

The weighted commutator and Cauchy inequality give

```text
|quadratic contribution to dM/dt|
    <= 2*kappa*sinh(s) <phi|N|phi> + kappa*sinh(s) ||phi||^2,
|linear contribution to dM/dt|
    <= 4*|g|*sinh(s/2) sqrt(<phi|N|phi>) ||phi||.
```

For the quadratic term, use
`|<phi|a_-^2|phi>| <= sqrt(<N_->(<N_->+||phi||^2))`
and bound that square root by `<N_->+||phi||^2/2`. For the force use
`||a_0 B phi|| <= sqrt(<phi|N|phi>)`. These inequalities hold for
matter/reference-entangled vectors; no scalar-state factorization is used.
Young's inequality with `epsilon=1/2` and
`sinh(s/2)<=sinh(s)/2` therefore yields

```text
dM/dt <= a*sinh(s)*dM/ds + b*sinh(s)*M,
a=2*kappa+(1/4)*epsilon=11/16,
b=kappa+(1/4)/epsilon=25/32.
```

The same constants include the disconnected control. Along the decreasing
characteristic

```text
tanh(s(t)/2) = exp(-a*t) tanh(s0/2),
M(s(t),t) <= M(s0,0) exp[(b/a)(s0-s(t))],
```

because `ds/dt=-a*sinh(s)`. The initial free vacuum has difference-mode
squeezing parameter `1/9`, so its exact initial moment is
`M(log(q0),0)=sqrt(80/(81-q0^2))` for `1<q0<9`. This formula is used only
at initialization. Taking `q0=8` gives `M(s0,0)=sqrt(80/17)`.

The calculation can first be made on total-number compressions
`R_L=1_{N<=L}`, with initial vectors `R_L psi` left unnormalized. Their
initial moment is no greater than the infinite initial moment, and the
compressed weighted commutators satisfy the same inequalities. The
finite-particle vectors form an essential selfadjoint core: acting `k`
times with this degree-two Hamiltonian is bounded on an initial
`N<=n` vector by a constant to the power `k` times
`product_{j=1}^k (n+2j+1)`. Thus these vectors are analytic for a
nonzero time radius. The selfadjoint operators `R_L H R_L`, extended by
zero on the orthogonal complement, eventually agree with `H` on each
core vector, hence converge in strong resolvent sense and their unitary
groups converge strongly. Bounded matter isometries and controlled
unitaries preserve this convergence through the measurement and
feedback. Lower semicontinuity of the positive number moment then
preserves its bound. Differentiability of an initially unknown infinite
exponential moment is therefore not assumed.

`occupation.number_moment_bound` uses a rational upper bound on
`exp(a*T)` to choose a rational `q<=exp(s(T))`. For `T=3/2` it obtains
`q=113/64`. Uniformly at every earlier prefix it bounds the moment by

```text
M <= sqrt(80/17) * R * min(1+(3/22)(R-1),
                         (121/100)[1+(3/88)(R-4)]), R=8/q.
```

The two factors bound `R^(3/22)` by concavity tangents at 1 and 4;
the second uses the exact rational inequality `(121/100)^22>=4^3`.
All tangent factors are positive for the relevant `R>=1`. Roots and
exponentials are rounded upwards rationally. At `T=3/2`, exact rational
comparisons certify `q>=7/4` and `M<13` (the implemented sharper moment
bound is approximately `12.108447`). No decimal approximation is used
to establish either inequality.

The complete first matter measurement can be retained coherently as
`V=sum_a |a> P_a`. It obeys `V*N=N*V`; the controlled matter feedback has
the same property. Thus both preserve the exponential moment, including
its sum over all unnormalized recorded branches. The argument continues
across the first interval, matter observation, feedback and second
interval with total duration `1+1/2`. Later classical dephasing of the
record is a common channel. No branch probability is divided out.

## Rectangular projection and complete retained-state error

Outside `P_d`, total occupation is at least `d`. Consequently
`||(1-P_d)psi(t)||^2 <= tau_d=M/q^d` uniformly through the history.
Only explicit lowering terms can cross from outside the rectangle to
inside. Their norms give

```text
||P_d H (1-P_d)|| <= C_d,
C_d = (kappa/2)*sqrt(d(d+1)) + (3*kappa/2)*d + (1/4)*sqrt(d).
```

The two onsite lowering-pair terms each contribute
`(kappa/4)*sqrt(d(d+1))`. Of the four bond monomials, double creation
cannot enter the rectangle; double annihilation and the two mixed terms
each have cross-norm at most `d` and coefficient `kappa/2`. The force
contributes `(1/4)*sqrt(d)`. Matter-only operations commute with the
projection and contribute zero.

The actual finite initial scalar preparation retains total occupations
`0,2,...,d-2` of the coupled vacuum and normalizes that projection.
Its omitted probability `v_d` is enclosed by the explicit rational
squeezed-vacuum recurrence. Its vector error is at most `sqrt(2*v_d)`.
Projected Duhamel comparison, continued through the commuting matter
measurement and feedback, gives

```text
projected_error <= sqrt(2*v_d) + T*C_d*sqrt(tau_d),
full_vector_error <= sqrt(projected_error^2 + tau_d).
```

The last sum is orthogonal because the first difference is inside `P_d`
and the final tail is outside. Trace distance between the normalized
retained histories is bounded by the same full vector error. This also
bounds all events of any **common** subsequent channel. Purification
extends the result uniformly to any initial matter/reference state with
the specified fixed scalar vacuum. It is not an unrestricted diamond
bound over arbitrary scalar inputs.

For fixed time the certified rate is `O(d*q^(-d/2))`; the initial
preparation tail decays faster. The rational rounding scales grow with
`d` and introduce no positive error floor. `certificate_table` evaluates
actual cutoffs 32, 64 and 128. Those are mathematical finite-dimensional
comparisons, not source execution at these dimensions.

The returned exact rational certificates are strictly below these
outward decimal bounds for total model time `3/2`:

| Onsite dimension | Retained-history trace-distance bound |
| --- | --- |
| 32 | 0.011408 |
| 64 | 0.000002498 |
| 128 | 0.00000000000006175 |

These comparisons include arbitrary matter/reference inputs and all
first measurement records. Independent tests verify the weighted
commutator inequality as a positive-matrix inequality on a finite
compression, check the boundary norm against a larger Fock matrix,
and identify the compressed action and preparation with those already
implemented by `scalar_coupling.py` and `field.py`. The tests are checks
of the implementation; the inequalities above provide the proof.

The bound stops immediately before the final scalar-coordinate sign
measurement. The finite sign projector and the infinite sign projector
are different maps, so contraction alone does not compare their terminal
instruments. Their separate readout comparison must be added explicitly.
Finite product splitting, floating coefficient/execution errors, physical
services and any spatial limit are also outside this occupation bound.
