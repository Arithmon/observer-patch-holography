# Tangent first laws and finite entropy changes

Maintainer follow-up (#1048): Gibbs tangent evaluation now conservatively
refuses subnormal thermal populations before converting them to exact
fractions. Exact arithmetic on rounded populations cannot undo their loss
of relative precision. For `T=diag(0,1e150)` and `lambda=7.45e-148`, the old
variance was `4.940656458412465e-24` instead of the independent partition
value `2.822350730472012e-24`; the entropy derivative had the same 75% error,
so the derivative ratio did not detect it. The new 500-digit scalar controls
cover both signs, ordinary and rescaled units, and gaps 700 through 746.
This is an explicit numerical-domain restriction, not failure of the entropy
identity. Ordinary receipt inputs, paper statements and frozen evidence are
unchanged by this correction.

This repair belongs to the standing evidence audit
[#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
It repairs the historical Einstein-closure receipts and their claim-registry
wording, building on the shared state, stable relative-entropy and Gibbs
repairs (#1029, #1041, #1046). It neither constructs a physical model nor
supplies the source, energy or temperature identifications in #1025/#1026.

## Reproduced failures

The first commit contains eight failing tests before the implementation
change. They reproduce on main `a11b9ab4` (and its predecessor `07067bf1`).

| Input to the old receipt | Old result | Required interpretation |
| --- | --- | --- |
| Wrong central weights `[0,0]`, `eps=-1e-6` | Edge defect `-0.4054651079332672` | An absolute error cannot be negative; the derivative mismatch is `log(3/2)` |
| Same wrong weights, `eps=1e-20` | All measured errors zero | The finite variation rounded away; this cannot certify normalization |
| Wrong weights at the default step | `split_identity_defect=NaN` | Measure the mismatch on the supplied split |
| Gibbs multiplier `1e-20` | Slope `1.539627470061857e-11` | Entropy subtraction noise replaced a nonzero derivative |
| Gibbs secant step zero or `1e-20` at multiplier `1.3` | `ZeroDivisionError` | Reject an undefined or unresolved secant explicitly |
| Boolean dimension or NaN central weight | Accepted | Reject malformed algebra/observable data |

There is also a sampling limitation independent of step size: a path with
fixed sector probabilities cannot see a wrong central normalization. The
new complete center diagnostic detects `log(3/2)` even on that path.

## One identity, two uses

Use ordinary matrix trace, natural logarithms and a faithful reference
`sigma`. Put `K_sigma=-log(sigma)`. For any normalized positive state `rho`,
expanding the definition of Umegaki relative entropy gives

```text
S(rho)-S(sigma) = Tr[K_sigma (rho-sigma)] - D(rho || sigma).
```

This is an exact **finite** identity, including noncommuting states. It is
not the equality of the two first-order terms. For a Hermitian traceless
tangent `D`, differentiation at the faithful reference instead gives

```text
dS_sigma[D] = -Tr[D (log(sigma)+I)] = Tr[K_sigma D].
```

The trace term vanishes on normalized tangents. The relative-entropy
remainder is second order, and is strictly positive for a distinct
normalized state. For example, from `sigma=I/2` to
`rho=I/2+epsilon X`, the modular change is zero while the entropy change
is `-2 epsilon^2 + O(epsilon^4)`. Omitting that remainder cannot be an exact
finite first law. This is the standard distinction used in
[Blanco, Casini, Hung and Myers](https://arxiv.org/abs/1305.3182), not a new
OPH law.

The existing `thermodynamics/conditional_repair_certificate.py` already
implements the classical cap identity with this remainder at 60 digits,
and its Lean heat-bath/Clausius bindings retain that scope. The new shared
`quantum_information/entropy_response.py` brings the older quantum receipts
into agreement. `finite_entropy_balance` **evaluates** the identity;
separate scalar spectra and high-precision entropy derivatives test it.
Subtracting two order-one entropies is no longer the production method.

## Joint counting stationarity and its geometric boundary

The existing [screen entropy identity](derive_screen_shock_spectrum.py)
has a joint state–capacity consequence. For positive integer sector counts
`d_i`, normalized probabilities `p_i`, and `M=sum_i d_i`, put `q_i=d_i/M`.
Then

```text
S_count(p,d) = H(p)+sum_i p_i log d_i = log M-D(p||q) <= log M.
```

Equality holds exactly at `p=q`. This compares every admissible finite
integer partition at the same total count `M`; it requires no derivative of
integer dimensions. The corresponding equality at dimension weights is
[proved in Lean](../../Lean/ObserverPatchHolography/DeSitterCapacityShock.lean).

In the **positive-real analytic relaxation**, take differentiable joint
paths `p(t),d(t)` with `sum p=1`, fixed `M`, and `p(0)=d(0)/M`. Writing
`u=p'(0)` and `v=d'(0)`, direct differentiation gives

```text
S_count'(0) = (log M-1) sum_i u_i + sum_i v_i/M = 0,
S_count''(0) = -sum_i (u_i-v_i/M)^2/p_i <= 0.
```

The second identity assumes twice differentiable constrained paths. Their
acceleration terms cancel by the same constraints. Equality directions
follow the moving maximizer; this is a nonstrict joint maximum, while the
area term alone can have positive fixed-capacity curvature. These derivatives
do not construct an operation that changes integer sector dimensions.

The count convention matters. For actual blocks
`rho=direct_sum_i p_i(rho_i tensor I_(d_i)/d_i)` with bulk dimension `b_i`,
the full entropy contains `sum p_i S(rho_i)`. Its exact deficit is

```text
S(rho) = log M-D(p||w/M)-sum_i p_i D(rho_i||I_(b_i)/b_i),
w_i=b_i d_i, M=sum_i w_i.
```

Thus full saturation also requires maximally mixed bulk blocks and uses
`w_i`, not edge counts alone. The simpler formula has one-dimensional bulk
blocks, or counts the entire sector as `d_i`.

Uniform microscopic repair statistics produce `q_i=d_i/M` when `d_i`
counts the accessible microstates in each observed sector. The existing
[native fixed-total calculation](../../model/physical_space_time/occupation_field.py)
checks this for the uniform stationary law on all 66 binary total-two
states. The finite preparation below supplies controlled approximate
stationarity from arbitrary binary inputs. Consensus or relative-entropy minimization against an
arbitrary reference does not select uniform counting: on two equally sized
sectors a reset kernel with stationary law `(1/4,3/4)` minimizes divergence
to that law, but the counting entropy derivative in direction `(1,-1)` is
`log 3`, not zero.

### Stationarity after finite native preparation

The [native warm-up theorem](../../model/record_laws/README.md#preparation-by-finite-native-warm-up)
removes exact stationary preparation on its binary total-two domain. After
`n` repair attempts, discard the complete warm-up tape and controller. For
every admitted record law, including the theorem's restricted law feedback,
the source marginal satisfies `dist(rho_n,I/66)<=epsilon_n`. The seam schedule
is exogenous, placement is fair, and `dist` is half trace norm. Under uniform
seam scheduling one may use `epsilon_n=min(1,33*(19079/20000)^n)`.
The parameter `n` counts repairs, including zero currents; it is not elapsed
physical time. An untouched quantum reference remains within the inherited
bound, but no decoupling of the discarded tape is asserted.

Choose **any complete partition** into nonempty microscopic sectors of
counts `d_i`, with `M=66`, `q_i=d_i/M`, and `p_i=Tr(P_i rho_n)`. Data
processing gives `TV(p,q)<=epsilon_n`. Let `eta=min q_i` and
`theta=epsilon_n/eta<1`. For every joint differentiable positive-real path
at this point, with `sum p=1`, fixed `M`, `u=p'`, and `w=d'/M`, direct
differentiation in the variation parameter (not repair count or a clock) gives

```text
S_count' = -sum_i u_i log(p_i/q_i) + sum_i w_i (p_i/q_i-1),
|S_count'| <= (||u||_1/2) log((1+theta)/(1-theta)) + ||w||_1 theta
           <= ||u||_1 theta/(1-theta) + ||w||_1 theta.       (E0)
```

Each coordinate differs from its reference by at most TV, so
`1-theta<=p_i/q_i<=1+theta`. For a zero-sum tangent, subtracting the midpoint
of a coefficient's range bounds its pairing by half the range times the
L1 norm. This proves (E0) uniformly over **all** such tangents, not only a
chosen path. The logarithmic bound is sharp: on singleton sectors, shift
two probabilities by `+epsilon,-epsilon`, take `u=(-1/2,1/2)`, and
`w=(1/2,-1/2)`. Both terms attain their bounds together. Capacity derivatives
are analytic variations, not native operations changing integer counts.

For the full quantum source density on this fixed 66-dimensional shell,
`||rho_n-I/66||_op<=epsilon_n`: a traceless Hermitian difference has positive
and negative parts each of trace `epsilon_n`. Thus `66*epsilon_n<1` gives
a faithful state. The same range argument applied to `log(rho_n)` gives
the complete trace-one tangent bound

```text
|dS_rho_n[X]| <= (||X||_1/2) log((1+66 epsilon_n)/(1-66 epsilon_n)),
X=X*, Tr X=0.
```

This includes noncommuting tangents. Also `log r<=r-1` and the same
positive/negative-part bound give
`log 66-S(rho_n)<=66 Tr[(rho_n-I/66)^2]<=132 epsilon_n^2`.
For the counting entropy the corresponding deficit is at most
`2 epsilon_n^2/eta`. The existing receipt computes finite attempt budgets
for derivative tolerances `1/100` and `1/1000`, uniformly over partitions
(`eta>=1/66`) and unit state/count tangent norms. Native transition columns,
complete tangent dual norms, coherent sharp controls and malformed-floor
tests check the attachment independently.

Record access and the conserved sector are essential. Even from the exact
uniform binary source, conditioning on current `+1` at seam `(0,1)` has
probability `5/66` and yields a uniform law on only ten states. Its entropy
deficit is `log(33/5)`, so the unread-state stationarity bound cannot be
assigned to that record-conditioned state. No postselected branch is used
in the preparation certificate. The full 78-state repair channel is
nonunital: one step sends its maximally mixed state to a nonuniform state
and strictly decreases entropy. The fixed-shell derivative bound therefore
does not extend to all 78-state tangents via the full-carrier warm-up bound.

This discharges joint counting stationarity on its stated finite/relaxed
domain. Identifying fixed `M` with fixed physical small-ball volume,
`sum p_i log d_i` with the area term, and the remaining block entropy with
the same continuum quantum fields still requires a common source dictionary
and controlled variations. [Jacobson's construction](https://arxiv.org/abs/1505.04753)
uses simultaneous geometry and field-state variations at fixed geodesic-ball
volume. The counting identity does not identify that geometric variation
family or its area coefficient.

### Native reference entropy and region boundaries

The [full native-return reference](../../model/observer_time_fields/REFINEMENT.md#the-full-reference-state-and-its-quantum-fluctuation-boundary)
has an absolute entropy bound, stronger than trace-distance convergence in
its growing dimension. Let `N=240*4^l`, `k=N/3`, `D=binom(N,k)` and
`tau=I/D`. Its transfer satisfies
`kappa=Tr(K)/D=1-chi`, `chi>=2/9`, and
`v=Tr((K-kappa I)^2)/D<=1/(16N)`.
For the ideal function `f(x)=exp(-acos(x))` or any finite return polynomial
`g_h(x)`, normalized to `sigma_h=g_h(K)/Tr(g_h(K))`,

```text
Delta_h := log D-S(sigma_h) = D(sigma_h||tau)
 <= Tr[(D*sigma_h-I)^2]/D <= 25*v/chi^2 <= 2025/(64N).      (E1)
```

For the ideal reference use `g_infinity=f`. All coefficients are positive,
so `|g_h(x)-g_h(kappa)|<=|f(x)-f(kappa)|<=|x-kappa|/chi` on `[0,1]`,
and the constant coefficient is `exp(-pi/2)>1/5`. A common scalar factor,
including the implemented finite-sampler normalization, cancels on forming
`sigma_h`. Thus `Var_tau(g_h)/tau(g_h)^2<=25v/chi^2`. The first inequality
follows from `log r<=r-1`, averaged with weights `r/D`, where
`r=D*sigma_h`. This is uniform over every cap, including zero, and proves
absolute entropy saturation without a dimension-dependent continuity bound.
Data processing also bounds the relative entropy of every admitted regional
partial trace by (E1), against the corresponding marginal of `tau`. It does
not turn that marginal into the identity on an unconstrained tensor factor.

The regional split gives a more precise physical boundary. For any set `R`
of `m` native ports, write

```text
a_j=binom(m,j), b_j=binom(N-m,k-j),
q_j=a_j*b_j/D, j=max(0,k-N+m),...,min(m,k).
```

Direct partial trace of `tau` gives
`S(tau_R)=H(q)+sum q_j log a_j` and the complementary formula with `b_j`.
Consequently its quantum mutual information is exactly
`I_tau(R:Rc)=H(q)`. These are full-state entropies in nats; they are not
correctable record capacities. For any, possibly coherent, fixed-`k` state
`rho`, put `Delta=log D-S(rho)` and let `p_j` be its regional-number law.
The marginals are block diagonal in `j` and `k-j`; their conditional
entropies are bounded by `log a_j` and `log b_j`. Hence

```text
I_rho(R:Rc) <= H(p)-D(p||q)+Delta <= log J+Delta,
J=min(m,k,N-m,N-k)+1.                                    (E2)
```

No diagonality of `rho` was assumed: coherences between distinct regional
charges disappear in either partial trace. Combining (E1) and (E2) bounds
the full coherent return reference by `log J+2025/(64N)` for every cap.

For the stationary native state, a separate exact bound is
`m*h(k/N)-log(N+1)<=S(tau_R)<=m*h(k/N)`, where `h` is binary entropy.
The upper bound is subadditivity. For the lower bound, compare with the
product Bernoulli state of density `k/N`: its relative entropy from `tau`
is `N*h(k/N)-log D<=log(N+1)`, since `k` is a binomial mode and has
probability at least `1/(N+1)`. Partial-trace contraction gives the bound.

This permits source-based controls on a proposed geometric dictionary.
On the twelve-port total-two carrier, the connected regions `(0,1,5,7)`
and `(7,8,10,11)` have **10 and 14 crossing seams**, respectively, but
identical sector counts `(28,32,6)` and identical reduced entropies and
mutual information. Across source refinement, a complete base-face patch
contains `m=12*4^l=N/20` ports and has `3*2^l` crossing seams: each of its
three sides is bisected at every refinement. Its stationary reduced entropy
grows proportionally to `4^l`, while the full return reference's mutual
information is at most `log(12*4^l+1)+2025/(64N)`. The latter divided by
crossing-seam count tends to zero. Exact native partial traces, coherent
controls, return spectra and the retained patch seam census test these
attachments independently.

These controls rule out assigning either regional reduced entropy or mutual
information a nonzero universal coefficient times **boundary-seam count
on this two-dimensional source**. They do not rule out a gravitational
area law: a whole `S2` screen's cell count could instead represent the area
of a three-dimensional region's boundary. That physical region map, its
bulk/edge split and the fixed-volume variation family still need construction.
The full return reference is a reconstructed state, not a demonstrated
native preparation of a spacetime vacuum.

## A complete generator test and its witness

Supply an independent Hermitian candidate `K`. Let

```text
R = K + log(sigma) - Tr[K + log(sigma)] I / n.
```

For every Hermitian traceless `D`, the response error is `Tr(R D)`. By
Cauchy--Schwarz its maximum absolute value on `||D||_HS=1` is `||R||_HS`;
if `R != 0`, the tangent `R/||R||_HS` attains it. Every such tangent is
admissible for sufficiently small positive and negative steps at a faithful
reference. Thus the first law holds on **all** normalized tangents if and
only if `K=-log(sigma)+c I`. The zero-residual case has no nonzero witness.
This proves the criterion, its completeness and exactly what it cannot
identify: the scalar energy origin. The code reports the norm and witness,
including witness trace/norm roundoff, without a Boolean theorem verdict.

This test accepts an independently supplied generator; its negative
controls use wrong coherent generators, including a variation of size
`1e-310` beside a scalar origin of `1e100`. Centering each operator before
adding them prevents the origin from erasing the response. Complete
independent Hermitian tangent bases and scalar characteristic-polynomial
derivatives check the norm and the attained response. A fixture that first
sets `K=-log(sigma)` only checks its own arithmetic consistency; the receipt
and registry now state that distinction.

## The center has the same scalar ambiguity

For the supplied direct-sum family
`rho = direct_sum_a p_a (rho_a tensor I_edge/d_a)`, spectral calculus gives

```text
S_bulk = H(p) + sum_a p_a S(rho_a),
S_edge = sum_a p_a log(d_a),
S = S_bulk + S_edge.
```

If `Z=direct_sum_a z_a I_a`, its tangent mismatch against edge entropy is
`sum_a (z_a-log(d_a)) dot(p_a)`. Write `e_a=z_a-log(d_a)`. The mismatch
vanishes for every normalized probability tangent if and only if all
`e_a` are equal: choosing `dot(p)=unit_a-unit_b` proves necessity, and
`sum dot(p)=0` proves sufficiency. The largest mismatch for transfers with
positive and negative masses each at most one is `max(e)-min(e)`, attained
by moving mass from a minimizing sector to a maximizing one. The code
returns that explicit witness. A one-sector center is correctly vacuous.

Consequently, a common shift of `Z` preserves the split's tangent laws;
choosing the particular representative `z_a=log(d_a)` is a convention.
Holding `p` fixed proves nothing about even the relative coefficients.
This connection to the full generator criterion removes a redundant
absolute-normalization claim without weakening the detectable mismatch.

At a finite comparison state `rho1`, the remainder also decomposes:

```text
D(rho1 || rho0) = D(p1 || p0) + sum_a p1_a D(rho1_a || rho0_a).
```

The equal edge factors cancel. This follows by taking the logarithm on
each direct-sum block. The receipts compare a full-matrix calculation with
that center/conditional-state decomposition, so moving probabilities are
covered rather than accidentally omitted.

## Gibbs slopes and numerical boundaries

For the supplied one-constraint family `rho(lambda)=exp(-lambda T)/Z`,
differentiating its energy eigenvalue probabilities gives

```text
dt/dlambda = -Var_rho(T),
dS/dlambda = -lambda Var_rho(T),
dS/dt = lambda, when Var_rho(T) > 0.
```

This is an analytic consequence of the family, not a finite-difference
measurement or a selection of physical energy units. A scalar `T` has no
energy coordinate and no such ratio. Setting `lambda=2pi` needs the stated
branch normalization. The complete generator theorem likewise identifies
`K=beta H+cI` only after a candidate `H` and `beta` are supplied.

The implementation reuses the validated Gibbs observable decomposition
and support checks. It computes variance from the positive pairwise sum
`sum_(i<j) p_i p_j (E_i-E_j)^2`, using exact rational accumulation of the
binary64 spectral data (with their total mass accounted for). It avoids
subtracting large moments or rotating a nearly pure state back to recover
tiny masses. Exact accumulation does not make the eigensolver exact.
Independent 160-digit partition-function derivatives check both tangent
components, including zero/tiny multipliers and a diagonal low-temperature
case. The separately labelled finite secant retains its truncation error.

Inputs must be finite numeric data, with no missing/masked values, Boolean
observables, numeric strings or lossy conversion. Tangents are Hermitian;
trace roundoff up to `64 eps ||D||_HS` is projected out explicitly, while
larger relative trace defects raise. State traces are never normalized
away: the finite balance reports their actual difference. For PSD inputs
of slightly unequal accepted trace, `D(rho||sigma)-Tr(rho-sigma)` is the
nonnegative Bregman remainder. Unresolved support, a rounded-away finite
comparison, an unrepresentable nonzero response or an unresolved remainder
raises instead of issuing a zero-error report. These remain numerical
diagnostics, not interval certificates of exact rank or exact equality.

## Impact and replay

`first_law_receipt` and `maxent_multiplier_receipt` retain their signatures
and legacy defect keys, but mark **schema 2**: the first-law defects and
primary Gibbs slope now concern analytic tangents. Finite changes and
secants have separate keys. No committed JSON receipt consumes these two
functions; their callers are the tests and mandatory gravity suite. The
existing gravity-suite assertions are tightened, and the dedicated
Windows/Linux quantum CI lane collects all new regression files. The
mandatory runner is byte-frozen by the invariant-mining campaign and stays
byte-identical to main; no freeze is regenerated to add test-list entries.

The entropy-bridge row `OPH-GR-D4D-ENTROPY-BRIDGE` and the shared derivation
paper fragment now distinguish the tangent theorem, finite remainder and
physical normalization. The active spacetime introduction already states
the cap remainder correctly. The two papers including the shared fragment
are rebuilt; the Standard Model branch skips the affected Einstein section.
Its PDF remains byte-identical. The changed Einstein PDF and the active
publication manifest are refreshed in preview mode; no release is bumped.
The recent book calibration correction from #1047 is integrated from main.
The live M1 necessity receipt binds the complete entropy-bridge registry
row: its canonical producer refreshes that one row hash. Physical outputs,
dependency roles, source pins and frozen evidence remain unchanged. No
book equation, Lean theorem, gravity-ladder status, D10 calibration result,
frozen registration, pinned historical evidence or source-selection claim
is promoted by this repair.

Replay with the repository's pinned requirements and `PYTHONPATH=code`:

```text
python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
python tools/check_claim_registry.py
python tools/build_gravity_ladder.py --check
```

The targeted independent controls live in `test_entropy_response.py` and
`test_first_law_audit.py`; maintainer-audit extensions are in
`test_entropy_response_audit.py`. They include malformed-input rejection, wrong
generator and normalization witnesses, noncommuting finite changes,
330-digit near-equilibrium scalar spectra, pure-source/faithful-reference
boundaries and a finite secant that demonstrably differs from its tangent.

The initial affected suite passed 872 tests on Linux and 866 on Windows, with
six pre-existing extended-precision skips on Windows. All 52 publication
manifest tests passed. Both paper builds passed the warning gate; the
changed Einstein pages 75 and 85 and continuation page 76 were rendered
and visually inspected.

Eight isolated implementation mutations were tested against the unchanged
74 shared-response controls. Each temporary copy contained its own package
and tests, avoiding the repository's import-path injection of the original
producer. Every mutation caused assertion failures:

| Deliberate error | Failing controls |
| --- | ---: |
| Omit the finite relative-entropy remainder | 16 |
| Reverse the remainder's sign | 16 |
| Report zero for every generator defect | 15 |
| Reverse the attaining witness | 15 |
| Suppress the central normalization spread | 1 |
| Use first rather than second power of energy units in the variance | 10 |
| Reverse the Gibbs entropy derivative | 9 |
| Accept traceful tangents | 1 |

## Maintainer audit of the response implementation

Ten further failing regressions were committed before the follow-up fix.
They exposed two numerical weaknesses in the first implementation:

- For `sigma=I/2+epsilon X`, a binary64 eigensolve erased coherences such as
  `epsilon=1e-20` and reported both zero tangent entropy and zero mismatch
  against a constant generator. The correct tangent on `X` is
  `-log((1+2 epsilon)/(1-2 epsilon))`. The finite change from epsilon to
  twice epsilon was consequently reported as `-2 epsilon^2` rather than
  `-6 epsilon^2+O(epsilon^4)`.
- At observable scale `1e-160`, separately rounded subnormal variance and
  entropy derivative produced slope `1.299901185770751` at multiplier `1.3`.
  A nonzero stored component alone does not establish adequate precision.

The response layer now diagonalizes the scalar-subtracted reference in a
private multiprecision context and uses `log1p` before restoring the scalar
origin. It preserves exact coordinate blocks using nonzero adjacency, with
no tolerance that could erase a small coupling. The output precision has
a 400-decimal-digit floor and grows with the component range. Merely
resolving first- and second-order input scales was insufficient: on a
five-site path with coupling `1e-80`, an intermediate implementation
reported a fourth-order response near `2.48e-300` instead of `3.125e-318`.
The retained path-series control catches this, as well as the quadratic
diagonal response on a three-site path.

An underflowing log entry is not automatically replaced by zero. If it
vanishes on every centered power `I,A,...,A^(n-1)`, exact rational arithmetic
certifies that its real or imaginary component vanishes for the centered
logarithm too: spectral interpolation represents `log(A)` by a real
polynomial of degree below `n`. Otherwise an unrepresentable nonzero entry
raises. This handles exact zeros in connected matrices, not just separate
blocks, and distinguishes them from small added coherences. The controls
use an exact `Q^2=2I` example and a nonzero response below binary64 range.
Support remains subject to the original faithful-state validator; this
does not infer exact rank or make the numerical log an interval proof.

The Gibbs tangent components must retain their exact-accumulation values
within eight binary64 epsilons relatively, and their ratio must agree with
the analytic slope within sixteen epsilons. Insufficient subnormal
precision raises explicitly. Independent unequal-weight complex-coherence
controls extend through amplitude `1e-300`, and the private numerical
context leaves the caller's global mpmath precision unchanged.

The CI audit also detected two binding defects. The live M1 row hash is
regenerated as described above; the independently replayed graph results
do not change. The invariant-mining failure came from adding entries to
its byte-frozen mandatory runner. Restoring that runner preserves the
original freeze and its independent verifier, while the existing quantum
CI workflow already executes the additional tests on both platforms.

After these repairs, the full affected suite passed **899 tests on Linux**
and **893 on Windows**, with the same six extended-precision skips. The
M1 necessity and invariant-mining controls passed all **116 tests**, and
the independent freeze verifier reported `GENERATION_LOCK_VALID`.
The follow-up adds 27 audit cases. All 13 isolated mutants were caught
against the 101 shared-response and audit controls: the original eight,
plus reinstating the binary64 reference log, reducing log precision,
accepting poorly resolved Gibbs components, flooring log underflow, and
discarding small reference coherences. The latter five caused 24, 2, 1, 1
and 20 assertion failures respectively.

A final center-boundary audit reproduced four more failures: distinct edge
dimensions `2^60` and `2^60+1` (also at powers 500 and 1000) had identical
rounded logarithms and falsely passed at `z=(0,0)`. Their mismatch is
`log1p(2^-power)`. At power 1100 the true mismatch is unrepresentable and
must raise, not pass as zero. The diagnostic now subtracts the common
central-energy origin first and evaluates relative dimension logarithms
in a private context whose precision follows the integer bit lengths.
A fifth control uses a ratio of large integers approaching `exp(1)` and
checks its small residual against an independent exponential construction.

The final totals are **904 Linux tests passed**, **898 Windows tests passed
with six existing skips**, and **32 maintainer-audit cases**. All **15**
isolated mutants fail against the 106 shared-response/audit controls.
Reinstating separately rounded dimension logs causes five failures;
reducing the relative-dimension precision to 80 digits causes one.
