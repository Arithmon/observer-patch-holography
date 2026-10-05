# Tangent first laws and finite entropy changes

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
new tests are part of that suite and the Windows/Linux quantum CI lane.

The entropy-bridge row `OPH-GR-D4D-ENTROPY-BRIDGE` and the shared derivation
paper fragment now distinguish the tangent theorem, finite remainder and
physical normalization. The active spacetime introduction already states
the cap remainder correctly. The two papers including the shared fragment
are rebuilt; the Standard Model branch skips the affected Einstein section.
Its PDF remains byte-identical. The changed Einstein PDF and the active
publication manifest are refreshed in preview mode; no release is bumped.
The recent book calibration correction is separate upstream work. No
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
`test_first_law_audit.py`; they include malformed-input rejection, wrong
generator and normalization witnesses, noncommuting finite changes,
330-digit near-equilibrium scalar spectra, pure-source/faithful-reference
boundaries and a finite secant that demonstrably differs from its tangent.

The full affected suite passed 872 tests on Linux and 866 on Windows, with
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
