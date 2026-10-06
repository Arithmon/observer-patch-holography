# Recovery evidence must describe the supplied state and counts

This repair under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
was developed on main `3d848b7f` and rebased onto `555852f8` after #1050/#1051
merged. It extends the finite-state validation
work and the squared-fidelity correction in #1038 to the remaining Stage 1
recovery and tomography pipeline. It supplies no physical model, recovery
instrument, or evidence distinguishing OPH from ordinary quantum mechanics.

## Reproduced defects

The first commit, `4d129953`, records **27 failing cases and four passing
existing controls**, before any implementation change.

| Input | Previous result | Correct treatment |
| --- | --- | --- |
| Zero matrix compared with `I/8` | Fidelity 1, obtained by manufacturing `I/8` | Reject: a zero matrix is not a normalized state |
| `diag(1-epsilon,0,0,0,0,0,0,epsilon)`, `epsilon=2^-40` | The recovered rare sector was zero | Preserve both resolved sectors; this classical Markov state recovers exactly |
| `diag(1,1e-100)` within the declared trace tolerance | Zero entropy | Retain the positive eigenvalue contribution, without normalizing the input |
| Empty tomography data | `I/8`, with apparently perfect recovery | Reject missing measurements |
| Empty, negative, fractional, Boolean or malformed counts | Zero or otherwise plausible expectations | Require matching binary outcomes and nonnegative integer counts, with positive shots per setting |
| Three estimates of XI with signed shot totals `(1,-100,1)` and shot counts `(1,100,1)` | `1/3` | The pooled estimate is `-98/102` |

The old inverse-root cutoff `1e-10` and subsequent output normalization were
separate errors. For the Markov example with `0<epsilon<=1e-10`, deleting the
rare B direction leaves a trace-deficient recovered operator with trace
`1-epsilon`. Renormalizing it yields `|000><000|`. Its trace distance from the
true state is epsilon, its squared fidelity is `1-epsilon`, and the Z read
changes by `2 epsilon`. A state normalization check alone cannot detect this
loss. The conditional sector is erased completely.

## One finite channel, with its support domain made explicit

Let tau be a normalized positive operator on BC, `t=Tr_C(tau)`, and Pi the
support projector of t. Let T be the Moore--Penrose inverse square root of t,
with every positive eigenvalue inverted. Define

```text
K_c = sqrt(tau) (I_B tensor |c>) T,
R_0(X) = sum_c K_c X K_c^*.
```

Direct multiplication gives `sum_c K_c^* K_c = T t T = Pi`. Thus R_0 is
completely positive and trace preserving on the supported input algebra;
the pseudoinverse alone does not make it trace preserving on all B inputs.

The new `petz_recovery_kraus` adds the explicit completion

```text
R(X) = R_0(X) + Tr[(I-Pi) X] tau.
```

For an orthonormal kernel basis `{v_r}`, its additional Kraus operators are
`L_jr = sqrt(tau)|j><v_r|`. Their adjoint products sum to `I-Pi`, because
`Tr(tau)=1`. Consequently the returned fixed map is completely positive and
trace preserving on the full B algebra. There is no division by the trace
of its output and no selection depending on the state being processed.
The replacement on the unused kernel is a disclosed completion choice,
not a derived physical law.

Positivity and `Tr[t(I-Pi)]=0` imply that tau is supported in
`Pi tensor I_C`: a positive operator with zero trace on a subspace
annihilates that subspace. Hence `R(t)=tau`. If rho_AB has marginal t,
the same support argument shows that its kernel-completion term vanishes.
Thus using tau=rho_BC recovers a normalized ABC state with the original
A marginal and BC marginal, and agrees with the usual supported Petz formula.
This proves the exact mathematical statements for finite states, including
singular references. The executable implementation evaluates them numerically.

The tests check every Choi entry through matrix-unit action, trace preservation,
reference recovery, arbitrary positive inputs on the completed kernel, and
65-digit independent matrix-function evaluations in dimensions 2x2, 2x3 and
3x2. A symbolic Bell-reference calculation proves that its channel replaces
any B input X with `Tr(X) P_Bell`; a valid recovery channel need not leave its
input marginal unchanged.

The general channel and fidelity definitions are the standard ones in
[Watrous, chapters 2 and 3](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).
They are used here as software controls, not new OPH assumptions.

## Complete classical conclusion, and the quantum distinction

For any classical distribution p_abc, the recovered distribution is

```text
q_abc = p_ab p_bc / p_b     when p_b > 0,
q_abc = 0                 when p_b = 0.
```

Summing the expression proves normalization and preservation of both AB
and BC. It is Markov, and a second reconstruction gives q again. For any
candidate Markov distribution `r_abc=s_b u_a|b v_c|b`, expansion of the
logarithm gives

```text
D(p || r) = I_p(A:C|B) + D(p_B || s_B)
            + sum_b p_b D(p_A|b || u_A|b)
            + sum_b p_b D(p_C|b || v_C|b).
```

The identity uses the usual support convention: a positive numerator over a
zero denominator gives infinity, and zero-p_b terms contribute zero.
Nonnegativity of the three added divergences proves that q minimizes
`D(p||r)` over the entire classical Markov family, with minimum
`D(p||q)=I_p(A:C|B)`. Equality fixes the candidate distribution q, including
its zero-weight sectors. This does not turn the generally nonconvex Markov
family into an affine constraint set. Nor is optimization over the second
argument of D the A3 optimization over its first argument.

An exact rational distribution is constructed on **each of the 255 nonempty
supports** of three classical bits. Independent rational marginal arithmetic
checks the recovered q; high-precision scalar logarithms check the entropy
identity. Separate Markov controls retain rare masses down to `1e-310`.

The corresponding quantum map need not preserve AB, be idempotent, or return
a Markov state. Noncommuting controls retain all three failures while checking
the A and BC preservation identities. The Fawzi--Renner theorem guarantees
the existence of a suitable recovery channel, not a fidelity promise for this
particular unrotated map. Its squared-fidelity bound in bits is `2^-I`.
[Fawzi and Renner](https://arxiv.org/abs/1410.0664)

The regression also replays the full-rank example obtained from Table I and
equation (20) of [Zhang's published counterexample](https://yuxuanzhang1995.github.io/assets/pdf/agentic/petz-cmi.pdf),
with maximally mixed admixture `10^-4`. A separate 65-digit implementation
reconstructs the exact decimal data, contracts indices and evaluates the
matrix functions. It gives `I + log2(F_Petz)` approximately `-0.00180240517`
bits. The ordinary Petz fidelity is below the optimal-recovery bound even
for this strictly positive state. The benchmark preserves that result rather
than adjusting its fidelity to satisfy an inapplicable inequality. This is
an attributed external control, not a new OPH result; the regression is a
high-precision numerical replay, not a replacement for that work's certificate.

## Missing measurements can hide one full bit

The tomography contract is a complete set of local Pauli settings, with
strict binary outcomes in Qiskit's count order. Repeated estimates of one
observable pool signed counts and shots before division. With equal shots,
this is the old mean; unequal shots retain their actual weight.

Completeness is essential for this unconstrained linear inversion. The states

```text
rho_0 = I/8,
rho_1 = (I + Z tensor Z tensor Z)/8
```

give identical uniform outcome probabilities in every local-Pauli setting
except ZZZ. Nevertheless `I(A:C|B)` is zero for rho_0 and one bit for rho_1.
For rho_1 every proper marginal is maximally mixed, while its full entropy
is two bits. An independent product-projector calculation checks all 26
remaining settings. Filling the missing ZZZ expectation with zero therefore
chooses a Markov answer that the data cannot distinguish from a non-Markov
state. Rejecting incomplete data closes that ambiguity; it does not imply
that tomography with a separately declared restricted model is impossible.

Shot noise can make complete linear inversion nonpositive. The existing
normalized-positive-part estimator remains an explicit tomography operation.
It is neither a channel nor a claim to be the nearest trace-one PSD matrix.
New reports expose its raw minimum eigenvalue, removed negative spectral
mass, Frobenius correction, per-setting shots and estimator name. They also
retain the complete count dictionaries for later replay. They do not report
those corrections as confidence intervals or certify statistical error.
The analysis and state-distance APIs never call this estimator.

## Numerical and downstream boundaries

Strict numeric conversion rejects masked data, mixed Booleans and lossy
binary64 conversions before the shared density validator. State normalization,
Hermiticity and negative spectral roundoff retain the shared `1e-12`
convention. Inputs are not rescaled and positive eigenvalues are not floored.
Inverse roots reject unresolved dense support instead of guessing rank;
coordinate-diagonal tiny positive entries remain explicitly supported.
Square roots and fidelities remain floating calculations and can have
spectral roundoff near singular states. Their values are not outward-rounded
certificates and are not clipped to make them satisfy a theorem.

The new code is shared with the finite quantum-information package. The
Stage 1 wrappers retain bit units and the legacy squared-fidelity result
fields. The existing optional-hardware imports remain lazy: numerical tests
need no Qiskit or provider access. No quantum jobs are submitted by this repair.

The frozen Stage 1 summaries contain neither complete counts nor reconstructed
states, so their measured recovery metrics cannot be recomputed from this
archive. Their bytes remain unchanged, and their reading notice records the
additional limitations. The structured circuit reference family can be rebuilt
analytically: its CMI and trace-distance values are retained to roundoff;
fidelity changes by about `2.1e-8` and `2.4e-8` at theta=.6 and 1, respectively,
consistent with the old inner-square-root evaluation's rank sensitivity.
This does not revalidate the unreplayable measured metrics.

The separately frozen issue-509 benchmark, other receipts, registry payloads,
papers, book, simulator, physical model gates #1025/#1026 and byte-pinned
mandatory runner are unchanged. Existing papers already distinguish supported
Petz action from a full channel and recovery existence from a chosen map.
The full-domain completion here makes that existing distinction executable.

## Reproduction and hostile controls

```text
PYTHONPATH=code python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent -W error
python -m pytest -q --confcutdir=code/ibm_quantum_cloud/tests code/ibm_quantum_cloud/tests/test_stage1_markov_fingerprint.py -W error
```

The existing quantum-information workflow runs the first command on both
platforms; mandatory shard zero already runs the Stage 1 file on both.
The affected core suite passes 927 tests on Linux and 921 on Windows, with
six existing extended-precision skips on Windows. The Stage 1 file passes
49 tests on each platform, for **95 added cases** across the two suites.
Warnings are treated as errors.
The tests include independent rational and symbolic controls, high-precision
matrix functions, complete Born-count reconstruction, invalid input rejection
under optimized Python, and a guard forbidding state repair during analysis.

Twelve isolated deliberate faults were detected: restored positive support
cutoff, omitted kernel completion, wrong inverse power, lost complex adjoint,
root fidelity reported as squared, constant perfect fidelity, zero-state
replacement, equal setting weights, old incomplete tomography, negative
counts, fabricated zero-shot data, and fidelity forced to meet the
optimal-recovery bound. Each ran in a separate copied package tree; collection
or import errors were not counted as detection. These controls exercise
incorrect scientific answers, not merely an implementation's success flag.
