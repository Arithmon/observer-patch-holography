# The interacting feedback signal survives infinite occupation

For the declared common-source action and zero native preparation current,
the unconditional positive-coordinate probability changes by

```text
[-0.024045670143, -0.024043450955]
```

when the recorded matter bit enables the existing feedback pulse. This is
an exact enclosure for **two infinite oscillator Fock spaces**, their original
coupled free vacuum, and the actual coordinate half-line readout. The
interaction, intervals 1 and 1/2, matter measurement and feedback are those
already specified in [the finite construction](DERIVATION.md). No measured
probability, finite-cutoff response or floating simulation supplies an input
to this certificate. Removing the interaction gives exactly zero contrast.

The separate [occupation theorem](OCCUPATION.md) establishes convergence of
the rectangular finite-Hamiltonian family through its retained matter read
and feedback, uniformly over initial matter states and untouched references.
The present state-specific calculation also includes the final readout and
resolves the signal. Neither result is a spatial limit or a physical
calibration of the declared action and model-coordinate time.

The qualitative sharp finite-detector limit is also proved in
[OCCUPATION_REFINEMENT.md](OCCUPATION_REFINEMENT.md); this calculation
additionally certifies the infinite-history signal. The subsequent
[complete-output comparison](READOUT.md#complete-comparisons-with-the-actual-finite-pointer)
uses these same retained polynomial states to enclose actual finite sharp
instrument errors and prove unavoidable occupation-code output errors.
Those bounds are separate from this signal enclosure. The separate Ramsey
instrument in [SMOOTH_REFINEMENT.md](SMOOTH_REFINEMENT.md) has effects
`(I +/- sin(Q0))/2`, a complete-output cutoff bound and finite-energy output.
A sharp half-line projection can instead produce infinite kinetic energy;
the trace-distance bounds here do not bound the postread energy.

## A fixed observable across occupation cutoffs

The earlier finite detector is the spectral projector of the finite matrix
`Q_d`. The infinite detector is multiplication by `1_(x0>0)`. Compression
and taking a spectral sign do not commute. Already at d=2,

```text
1_(Q_2>0) = (I+X)/2,
V_2* 1_(Q>0) V_2 = I/2 + X/sqrt(2*pi).
```

Their operator difference has norm about 0.101058. Therefore the earlier
finite detector probability cannot be promoted merely by appending a
dynamical occupation-error bound. Its finite result remains valid.
[The readout calculation](READOUT.md) supplies exact Hermite matrix elements
for the infinite half-line projector, preserving the complete postread
wavefunction, including its infinite occupation tail.

The d=2 finite continuous experiment had contrast about -0.03166148. With
the infinite half-line effect but still d=2 dynamics, the exact comparison
in `test_halfline.py` gives about -0.02526221. Removing the dynamical cutoff
as well gives the interval above. Thus the signal survives, while its size
is materially affected by both approximations.

## The same unbounded Hamiltonian in exact coordinates

Write unnormalized oscillator monomials as
`|n0,n1,g>_B=(a0*)^n0(a1*)^n1|0,0,g>`, with gauge index `g=0,1,2,3`.
Their squared norms are `n0! n1!`. Creation shifts an index; annihilation
also multiplies its coefficient by that index. In these coordinates the
complete interacting Hamiltonian has rational coefficients:

```text
H = X0+X1+(-b0+b1)^2 + (73/64)(N0+N1+1)
    + (9/128) sum_i(a_i^2+a_i*^2)
    - (9/64)(a0+a0*)(a1+a1*)
    + (1/4)(a0+a0*)Z0.
```

This is the original projected-square convention before projection. In
particular no boundary term from squaring a truncated quadrature is lost.
Tests independently convert these integer actions to normalized Fock
matrices and compare with the existing d=2 and d=4 Hamiltonians.

Let `K` be the positive coupled scalar oscillator Hamiltonian. Its quadratic
form dominates `N0+N1+1` and `Q0^2/2`. Thus
`||Q0 v||^2 <= 2 <v,Kv> <= 2 ||v|| ||Kv||`, making the linear force
arbitrarily small relative to K, with a compensating norm term. The matter
Hamiltonian is bounded. The resulting matrix-valued oscillator H is
self-adjoint on `D(K)`. Finite Hermite polynomials are in every `D(H^k)`;
applying H changes total occupation by at most two. These are sufficient
domain facts for each finite Taylor remainder below. No convergence of
the infinite Taylor series at time 1 is assumed.

## Exact local residuals, with every discarded amplitude charged

The infinite free coupled vacuum is proportional to

```text
exp[-(a0* - a1*)^2/36] |0,0>.
```

Its squared normalization factor is `4 sqrt(5)/9`. Truncating this initial
polynomial at total occupation 12 gives a rational squared norm `S0`.
An integer-square-root lower bound for sqrt(5) supplies the exact omitted
probability bound `t0 >= 1-(4 sqrt(5)/9) S0`. The normalized projected
preparation differs from the original vacuum by trace distance at most
`sqrt(t0)`. Initial gauge state is the same blank code used by the finite
zero-current experiment. Subsequent polynomial coefficients are never
renormalized; all weights are divided by this **one initial S0**.

For each time step delta=1/16 and supplied finite polynomial v, form

```text
T_p(v) = sum_(k=0)^p (-i delta H)^k v/k!,       p=16.
```

Every power uses the full H on expanding support. There is no intermediate
occupation projection. The spectral theorem and the real-variable integral
Taylor remainder give

```text
||exp(-i delta H)v - T_p(v)||
    <= delta^(p+1) ||H^(p+1)v||/(p+1)!.
```

The next term is computed exactly, including all factorial norm weights.
After this polynomial is formed, components above total occupation 16 are
discarded, and retained real and imaginary coefficients are rounded to a
grid of spacing `2^-112`. The weighted norms of **both actual differences**
are computed exactly and added to the step's remainder. Directed integer
square roots bound all three norms. This is a certified numerical
approximation of the unbounded H, not a time-stepped finite Hamiltonian
with an assumed cutoff error. The largest expanded occupation is 50.

If v_j is the current approximate state, the true unitary's contraction
gives `e_(j+1) <= e_j + r_j`. No norm-growth factor for the Taylor polynomial
is needed: it is applied to the approximate input v_j, and its local
residual is evaluated on that actual input. Sixteen steps reach the first
matter observation. Its complete projectors split the state into the two
retained record blocks. On record one, feedback is exactly `(3I-4iX0)/5`;
on record zero it is identity. This record-preserving isometry retains the
incoming vector-error bound. Eight further steps evolve each unnormalized
block. Their additional error bounds combine by a sum of squares because
the record flags are orthogonal. Dephasing or duplicating a flag into an
environment makes the actual committed classical record and cannot increase
the resulting density-operator trace norm.

For each pre-readout polynomial, the [half-line formula](READOUT.md) gives
an exact probability expression `r+s/sqrt(2*pi)`. Machin's identity,
alternating arctangent tails and directed integer square roots enclose this
without any binary floating arithmetic. Both final outcomes are retained.
Their quantum outputs are the original polynomial Gaussian wavefunctions
restricted to the corresponding half line. They are not projected back
into finite occupation or replaced by classical scalar weights.

## What each error number controls

Let epsilon bound the final flagged-vector error relative to the normalized
projected preparation, and let tau bound the initial vacuum trace distance.
Both final sign projectors act on the same infinite Hilbert space. Their
complete recorded measurement contracts trace distance. Hence the complete
classical/quantum output has half trace-norm error at most

```text
tau + epsilon + epsilon^2/2.
```

The approximate polynomial need not have norm one. Its error for any one
event, including the **unconditional** positive-coordinate event, is bounded
by `tau+2 epsilon+epsilon^2`. Using only the half trace norm for this
unnormalized probability comparison would omit a trace defect. This bound
already includes both first records; it is not multiplied by their count.
The enabled-minus-disabled contrast uses the sum of the two probability
budgets. Exact rounding outwards to twelve decimal places gives the stated
interval. The wider total-occupation-12 control overlaps it, while purposely
insufficient calculations return an interval containing zero.

With coupling zero, the gauge dynamics, matter measurement and feedback act
on a factor disconnected from the scalar. The scalar starts in its even
coupled vacuum and remains there up to phase. Each scalar half line has
probability 1/2 and the contrast is exactly zero. The computed disconnected
control must enclose those values; its small approximation error is retained.

## The finite sign instruments have this same limit

The detector correction above is essential at finite d, but it does not
require replacing the mathematical finite family to obtain a limit. Take
`d=2^k` and embed its rectangular projection `P_d` into the original Fock
space. Let `A_d=P_d Q0 P_d`, extended by zero outside this rectangle.
Each A_d is selfadjoint. On every finite Hermite polynomial, A_d eventually
agrees exactly with Q0. Those polynomials form a core for Q0. In particular,
on the dense set of vectors `(Q0-i)u` with u in that core,

```text
(A_d-i)^(-1)(Q0-i)u - u = (A_d-i)^(-1)(Q0-A_d)u -> 0.
```

The resolvents have norm at most one, so convergence extends to the entire
Hilbert space. This is strong resolvent convergence. The spectral measure
of Q0, multiplication by x0, has no atom at zero. Approximating the half-line
indicator by continuous spectral functions, with a shrinking neighborhood
of zero and increasing outer cutoffs, therefore gives strong convergence
`B_d=1_(A_d>0) -> P_+=1_(Q0>0)`.

On the d-dimensional coordinate code, the tridiagonal determinant at zero
obeys

```text
det Q_(2m) = (-1)^m product_(j=1)^m ((2j-1)/2) != 0.
```

Thus the embedded finite negative projector is exactly `P_d-B_d`, and
converges strongly to P_-. Write the complete sign isometries as
`W_d=|+>B_d+|->(P_d-B_d)` and `W=|+>P_++|->P_-`. W_d is a contraction
on the ambient space and an isometry on `ran P_d`. The occupation theorem
gives convergence of the retained pre-readout vectors `psi_d -> psi`, so

```text
||W_d psi_d-W psi||
    <= ||psi_d-psi|| + ||(W_d-W)psi|| -> 0.
```

Keeping the first matter flag, then dephasing the final flag, proves
trace-distance convergence of the **complete finite sign instruments**,
including their terminal quantum states. With the fixed pure scalar
preparation, the initial matter input space has dimension four. Strong
convergence of the history isometries on this finite input space is uniform
in operator norm and remains so after tensoring an untouched reference.
It is not an operator-norm sign-projector bound on arbitrary scalar inputs.

Define also the mathematical exact palindromic factors by the same onsite,
bond, interaction, hopping and electric formulas in the finite construction,
at each d. The existing temporal proof gives
`S_d=3+(25/16)(2d-1)+(1/4)sqrt(4d-2)=O(d)` and combined error
`3 S_d^3/(8 n^2)` for intervals 1 and 1/2 with n steps each. Choosing
`n=d^2` makes this error `O(1/d)`. Both exact-factor feedback histories
therefore converge to the infinite histories certified here. Their contrast
tends to a value inside the displayed strictly negative interval and is
strictly negative for all sufficiently large d along this sequence.

This proves convergence, not a numerical rate for the discontinuous readout
or an explicit sufficient cutoff. The mathematical family uses exact
preparation, feedback and factor coefficients. The current native compiler
has deliberate size limits, and the common-source execution remains d=2,4;
this corollary does not claim larger native words were emitted or run.
Floating compilation and physical services still need their own budgets.

## Reproduction and retained evidence

```sh
python3 -m model.common_source.infinite_history --check model/common_source/infinite_results.json
python3 -m pytest -q -W error model/common_source
```

The receipt binds exact errors, all four readout branches per configuration,
occupation comparisons and hashes of each pre-readout polynomial. It stores
compact certificates instead of repeated full state arrays. Regenerate all
fine polynomial coefficients with

```sh
python3 -m model.common_source.infinite_history --polynomials --write PATH.json
```

Writes refuse overwrite. Replay compares every exact field and rejects
duplicate JSON keys. Polynomial rows are `[n0,n1,g,real,imaginary]` in sorted
order, with rational coefficients and the separately recorded S0. Applying
the two half-line indicators reconstructs the full approximate terminal
states and their coherences.

The emitted d=2/d=4 source words and their old receipts retain their finite
scope. This computation certifies the same declared action at infinite
occupation; it does not emit an infinite apparatus, establish its physical
services or source a spatial continuum, clock calibration, or empirical
coupling. Those existing physical identification questions are distinct
from the occupation and readout errors resolved here.
