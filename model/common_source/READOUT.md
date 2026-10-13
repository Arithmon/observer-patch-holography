# Infinite coordinate-sign readout

The finite source pointer measures `B_d = 1_(Q_d > 0)`. An occupation-limit
comparison must also account for its difference from the infinite oscillator
projector `P = 1_(Q > 0)`: functional calculus does not commute with compression.
Already at `d=2`, `B_2=(I+X)/2`, whereas
`V_2* P V_2=I/2+X/sqrt(2*pi)`. Their operator-norm difference is
`1/2-1/sqrt(2*pi)`, approximately `0.10105772`.

`halfline.py` evaluates the original infinite projector directly on finite
polynomial approximants. For normalized Hermite functions `h_n`, integration
of the oscillator differential equations gives, for unequal indices,

```text
integral_0^infinity h_m(x) h_n(x) dx
  = [h_m(0) h_n'(0)-h_m'(0) h_n(0)]/[2(n-m)].
```

The diagonal is `1/2`; other equal-parity entries vanish. For `m=2j,n=2k+1`
the entry is

```text
(-1)^(j+k) sqrt((2j)! (2k+1)!)/
[sqrt(2*pi) 2^(j+k) j! k! (2k+1-2j)].
```

The standard Hermite recurrence, derivative and polynomial formulas are
tabulated in [NIST DLMF §18.9](https://dlmf.nist.gov/18.9) and
[§18.5](https://dlmf.nist.gov/18.5). The displayed half-line formula follows
by the boundary integration above.

In the unnormalized Bargmann basis `|n>_B=(a†)^n|0>`, multiplying by
`sqrt(m! n!)` removes every factorial radical. Thus a finite vector with
Gaussian-rational coefficients has positive weight `S/2+T/sqrt(2*pi)`:
`S` is its factorial-weighted norm and `T` is an exactly computed rational.
The untouched second oscillator contributes its factorial norm; different
spectator and gauge labels do not interfere. A separately supplied rational
initial normalization divides both terms. Current branch mass is never used
as a normalization denominator.

The implementation encloses pi by Machin's identity
`pi=16*atan(1/5)-4*atan(1/239)` and the signed remainder of each alternating
Taylor series. Integer square roots then enclose `1/sqrt(2*pi)`. No floating
quadrature, target detector probability, or observed convergence supplies
these bounds. Both signs are retained, and their exact weights sum to `S`.
Unnormalized inputs of mass greater than one retain that mass.

The complete postread vector is `P_b u`, represented by multiplying the
polynomial Gaussian wavefunction by the appropriate half-line indicator.
Its density kernel is `1_b(x0) 1_b(y0) u(x) u(y)*`. This retains all
coherences within each sign subspace and all earlier classical records.
The resulting vector generally has infinite occupation support. The routine
does not project it back into the polynomial support. Consequently a proved
trace-distance bound for the complete pre-read history passes through this
same infinite projective measurement by contractivity.

For a separate comparison with the executed finite pointer, let `A=V_d*P V_d`
and `B=B_d`. For any finite unnormalized pure input `u` (including arbitrary
untouched reference factors), put `p_b=<u|A_b|u>`, `q_b=<u|B_b|u>` and
`z_b=<u|A_b B_b|u>`. The exact half trace norm between the two complete
postread outputs is

```text
(1/2) sum_b sqrt((p_b+q_b)^2-4*abs(z_b)^2).
```

Here `A_1=A`, `A_0=I-A`, and likewise for `B`. Sum this expression over
earlier committed records. This identity includes the infinite output tail;
checking only the positive probability does not certify that quantum output.
It follows from the trace norm of a difference of two rank-one operators.

## A lower bound for every output confined to a finite code

The infinite postread tail also obstructs approximation of the complete
quantum output by any state confined to a specified finite code. Let
`u_a` be the unnormalized pre-read vector for retained matter record `a`,
and put `v_ab=P_b u_a`. For the rectangular occupation projection `Pi_d`,
define

```text
p_ab = ||v_ab||^2,
a_ab = ||Pi_d v_ab||^2,
t_ab = ||(I-Pi_d) v_ab||^2 = p_ab-a_ab.
```

The inside and outside components are orthogonal. When both are nonzero,
the operator exchanging their normalized directions has operator norm one,
has expectation `2 sqrt(a_ab t_ab)` on `|v_ab><v_ab|`, and has zero
expectation on every operator supported inside the code. Taking the direct
sum of these witnesses over all old records and both signs proves

```text
(1/2)||rho_sign-sigma_d||_1 >= sum_(a,b) sqrt(a_ab t_ab),
```

for every code-supported output `sigma_d`. Zero components contribute zero.
The witness does not assume that the comparison state has the same branch
probabilities, is pure, or was produced by a particular finite detector.
All retained records contribute; none is conditioned on success.

This bound is evaluated on the actual interacting polynomial histories in
[sharp_output.py](sharp_output.py), with the approximation error retained.
If their complete recorded output `rho_poly` obeys
`(1/2)||rho_true-rho_poly||_1 <= epsilon_poly`, the certified lower bound
for the true output is

```text
max(0, sum_(a,b) sqrt(a_ab t_ab) - epsilon_poly).
```

Here the polynomial need not have unit mass. Its one initial normalization
is the same `S0` used in [infinite_history.py](infinite_history.py); no current
branch norm is used to renormalize it. The inherited error
`epsilon_poly=tau+epsilon+epsilon^2/2` includes initial preparation,
Hamiltonian Taylor remainders, discarded support and coefficient rounding.

The tail calculation is finite despite the infinite output. If `u_n` are
Bargmann coefficients with spectator occupation and gauge index held fixed,
the coefficient of `|m>_B` in `P_b u` is

```text
c_m = (1/m!) sum_n <m_B|P_b|n_B> u_n.
```

Sum `m! n1! |c_m|^2/S0` over `m,n1<d` and gauge to obtain `a_ab`.
Every matrix element is rational plus a rational multiple of
`alpha=1/sqrt(2*pi)`. Thus `a_ab` is quadratic in alpha and `p_ab` is
linear. Form the polynomial `t_ab=p_ab-a_ab` before interval evaluation;
the known inequalities `0<=a_ab,t_ab<=p_ab` permit safe intersections.
This accounts for the full omitted tail without truncating the postread
wavefunction or summing an assumed asymptotic series.

The vacuum is an independent control, not the interacting-history claim.
For `d=2m`, its omitted mass in each sign branch is exactly

```text
T_d = 1/4 - (1/(2*pi)) sum_(k=0)^(m-1)
                         binomial(2k,k)/(4^k(2k+1)).
```

Its two-outcome witness is `2 sqrt(T_d(1/2-T_d))`. Matching the sign
probabilities, both one half in this control, does not remove this quantum
output discrepancy.

There is an elementary all-cutoff resource bound for this vacuum control.
Write `c_k=binomial(2k,k)/4^k`. Its recurrence gives
`c_k-c_(k+1)=c_k/(2k+2)`, so telescoping the positive tail yields

```text
c_m/(2*pi) <= T_(2m)
           <= ((2m+2)/(2m+1)) c_m/(2*pi).
```

Induction gives `1/(2 sqrt(m)) <= c_m <= 1/sqrt(3m+1)` for `m>=1`.
The lower induction step reduces to
`(2m+1)^2-4m(m+1)=1`; the upper step reduces to
`4(m+1)^2(3m+1)-(2m+1)^2(3m+4)=m`. In particular,

```text
1/(4*pi*sqrt(m)) <= T_(2m) <= 2/(3*pi*sqrt(3m)).
```

The tail decreases and `T_2=1/4-1/(2*pi)<1/8`. The distance from the
vacuum's complete sign output to every normalized code-supported state is
therefore at least

```text
sqrt(3/(8*pi)) m^(-1/4).
```

A constructive upper bound uses a different, explicitly specified
comparison state: project each true sign branch into the code, normalize
its direction, and restore its original mass one half. The rank-one
formula gives distance exactly `sqrt(2 T_(2m))`, at most
`sqrt(4/(3*pi*sqrt(3))) m^(-1/4)`. This surrogate is not the output of
the finite `B_d` pointer. Together the two bounds show that the best
code-supported approximation has order `d^(-1/4)`. Achieving distance at
most epsilon consequently requires
`d >= 9/(32*pi^2*epsilon^4)` under this vacuum comparison contract.
These bounds do not assert the same rate for the interacting histories;
their finite-cutoff bounds are certified separately below.

## Complete comparisons with the actual finite pointer

The rank-one formula also compares different pre-read histories. For an
infinite-history polynomial `u_a/sqrt(S0)` and any polynomial finite-output
branch `w_ab`, use

```text
p_ab = <u_a|P_b|u_a>/S0,
q_ab = ||w_ab||^2,
|z_ab|^2 = |<u_a|P_b|w_ab>|^2/S0,
D_poly = (1/2) sum_(a,b) sqrt((p_ab+q_ab)^2-4|z_ab|^2).
```

The overlap is a finite Hermite sum, even though `P_b u_a` has infinite
support. It is linear in alpha before taking its squared modulus. The
formula therefore retains the infinite tail, relative phases, both signs
and every old record. It does not replace the quantum output by its outcome
probabilities. Its use of squared overlap avoids introducing an unnecessary
radical for the initial normalization.

For `d=2`, [finite_signal.py](finite_signal.py) supplies independently
certified Gaussian-rational vectors for the actual finite history, including
its finite sign pointer. Its two Taylor operator errors `epsilon1,epsilon2`
give a complete flagged-vector error
`eta=epsilon1+epsilon2(1+epsilon1)`: the intervening complete matter record,
feedback, and final complete pointer are isometries before record dephasing.
The resulting complete-output half trace-norm error is at most
`eta+eta^2/2`. Adding this error and `epsilon_poly` to the two sides of the
polynomial comparison gives an enclosure for the complete finite-versus-
infinite interacting outputs. This comparison includes their different
declared finite and infinite initial scalar preparations.

Larger even cutoffs use the following independent certificate for the
one-site sign projector. It verifies an approximate eigenbasis by exact
residuals; the eigensolver's convergence or printed eigenvalues are not
assumed correct. Put `A=a+a^dagger=sqrt(2) Q_d`. In the Bargmann basis A
has integer entries and is self-adjoint for the metric
`G=diag(0!,1!,..., (d-1)!)`. Let U be a supplied rational approximate
eigenbasis and Lambda a real rational diagonal matrix with no zero entries.
Use exact rational squared norms and directed roots to certify

```text
||U* G U-I||_F <= epsilon < 1,
||G^(1/2)(A U-U Lambda)||_F <= r,
M >= ||A||,    L >= max_j |Lambda_jj|,
gamma <= min_j |Lambda_jj|,    gamma > 0.
```

The A norm uses the G inner product; `M=2d` is a sufficient elementary
bound. The Gram estimate gives a G-unitary polar factor V with
`||G^(1/2)(U-V)||_F<=epsilon`. Hence, in orthonormal coordinates,

```text
||A-V Lambda V*||_F <= delta = r+(M+L)epsilon.
```

For Hermitian A and C=`V Lambda V*`, write their positive projectors in
their respective eigenbases. In the mixed eigenbases the entries satisfy

```text
(P_A-P_C)_ij =
  (1_(a_i>0)-1_(lambda_j>0)) (A-C)_ij/(a_i-lambda_j)
```

when the numerator is nonzero; all other projector-difference entries are
zero. Opposite signs make the denominator's magnitude at least gamma.
Thus `||P_A-P_C||_F<=delta/gamma`, without assuming an unverified spectral
gap for A. Even d ensures the finite Jacobi matrix has no zero eigenvalue,
as established in [the occupation proof](OCCUPATION_REFINEMENT.md).

For the computable rational seed `S=U 1_(Lambda>0) U* G`, the polar-factor
error adds at most `epsilon(2+epsilon)`. Therefore

```text
||B_d-S|| <= eta_d = delta/gamma + epsilon(2+epsilon).
```

[sharp_pointer.py](sharp_pointer.py) applies this certificate using
factorial-weighted rational arithmetic. The seed S is a computational
approximation, not an asserted exact projector or physical instrument.
On a polynomial input of total mass `M_poly`, replace its two actual
finite branches by `S u` and `(I-S)u`. Their complete density-operator
half trace-norm error is at most

```text
E_seed = M_poly (2 eta_d + eta_d^2).
```

Indeed each branch-vector error is at most `eta_d ||u||`; the complete
two-branch vector error is at most `sqrt(2) eta_d ||u||`. The displayed
bound is a conservative consequence of the rank-one density norm
inequality. Sum over the orthogonal old records using their total mass;
no factor equal to the number of histories is needed.

Finally suppose d contains the polynomial's entire pre-read support. Let
`E_pre(d)` be the independently certified distance between the true finite
and infinite histories immediately before the sign readout, from
[occupation.py](occupation.py). Evaluate the rank-one formula above with
the same polynomial under the infinite sign map and the rational finite
seed, obtaining `[D_low,D_high]`. The distance between the actual complete
finite and infinite histories is then enclosed by

```text
R = E_pre(d) + 2 epsilon_poly + E_seed,
max(0,D_low-R) <= D_actual(d) <= min(1,D_high+R).
```

The two polynomial errors compare each of the two readout maps with its
polynomial input. The finite map may be extended to a CPTP map on the
ambient Hilbert space when applying contractivity; its action on the code
is unchanged. The pre-read comparison and this explicit sign correction
are both required. These are finite-cutoff complete-output bounds for the
declared zero-current interacting histories, not an unrestricted detector
diamond bound or an arbitrary-matter uniform extension of the polynomial
certificate. Mathematical comparisons at larger d do not claim larger
native source words were executed. Finite splitting, numerical source
execution and physical service errors remain separate budgets.
