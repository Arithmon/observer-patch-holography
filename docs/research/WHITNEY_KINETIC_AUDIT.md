# Whitney kinetic reduction: positive action and numerical precision

This repairs existing finite evidence under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
It does not change the supplied action, select a physical source or clock, or
advance the physical-model decision in #1025/#1026. The exact reduced metric
in the paper remains positive; its former binary64 evaluation did not always
preserve that fact.

## Reproduced defect

Reviewed baseline: `7521d5c3a829ceebb0f86df1c34c70de7a52cb25`. Set all edge
coefficients to zero, all 13 complex matter coefficients to the real value
`T=1e9`, charge to `1/4`, and the scalar imaginary velocity to
`x=(12,-1,...,-1)`. The remaining velocity components vanish. The original
source action reduces independently to the one-variable minimization

\[
 \min_z \{Az^2+B(1+eTz)^2\}
 =\frac{AB}{A+B e^2T^2},\qquad
 A=169(30-10\sqrt5),\quad B=\frac{114}{5}(10+10\sqrt5/3).
\]

The exact positive answer is approximately `2.0656721888405686e-14`.
The previous dense reduction returned `-1.6542323066914832e-13`; the public
packet path also accepted a negative momentum/velocity contraction. A
positive determinant sign did not establish a positive metric. Even at
`T=1e6`, where density evaluation is practical, the old Gaussian log amplitude
was wrong by about `4.6e-6`. Both independent receipt verifiers repeated the
same binary64 Schur subtraction and inherited its cancellation.

The failure controls were committed before the implementation repair in
`13baf783`. [The retained original-input controls](../../code/electromagnetism/test_whitney_kinetic_precision.py)
derive the scalar modes and absolute determinant from the original cone,
simplex moments and icosahedral adjacency spectrum. They do not use producer
quadrature or a rounded reduced matrix as the expected answer.

## Reduction and domain

Write the positive quadrature action as `||H v + V eta||²`, where `H` contains
the 56 slice columns and `V` the 12 mean-zero gauge columns. Normalize each
column of `[V,H]`, apply Householder QR, and undo the column scales. The
trailing triangular block is a factor `F` of the minimized action; the leading
block gives its minimizing gauge velocity. This replaces subtraction of
nearly equal Gram matrices. Gaussian density uses `2 sum(log(abs(diag(F))))`
directly. A density-only calculation omits the potential, so an unused
quartic overflow does not reject a resolved kinetic density.

The scalar edge derivative now pairs its endpoint values before multiplying
by barycentric weights. Uniform matter at zero edge field therefore has
exactly zero edge derivative, even at large amplitude. The vertical scalar
factor uses the Ward identity `i e Psi lambda` directly.

The same constant-mode cancellation affected the spatial gradient. At
`a=0`, uniform neutral matter `T=1e80` and zero potential couplings, the old
gradient calculation invented energy `9.640784334323916e128` instead of zero.
Subtracting the transported constant mode before differentiation fixes it.
A nearby field `psi=(1e15+1,1e15,...,1e15)` independently checks the nonzero
energy `30-10 sqrt(5)`; it was previously reported as `9.088123195321897`.
Zero mass/quartic terms are omitted before evaluating powers, and genuinely
required energies outside binary64 range are explicitly refused.

The implementation checks both the scale of the horizontal factor before
projection and the conditioning of the normalized vertical columns. Its
binary64 resolution indicator is

\[
 68\epsilon\,\frac{\|H\|_2}{s_{\min}(F)}
 \kappa_2(V_{\rm normalized}).
\]

Returning a dense metric additionally requires
`56 epsilon ||abs(F).T abs(F)||₂ / s_min(F)² <= 1e-7` and successful Cholesky
factorization. Both indicators must be at most `1e-7`. These are numerical
resolution policies, **not certified forward-error or quadrature bounds**.
There is no eigenvalue floor or projection onto positive matrices.

The distinction matters in actual controls: the uniform `T=1e6` factor and
Gaussian density are evaluated, while the unresolved dense `T=1e9` metric and
cotangent are refused. With matter `(1e12,1,...,1)`, testing only `cond(F)`
would accept a projection already wrong by about `2e-5`; the combined test
refuses it. The nearby `(1e6,1,...,1)` configuration is retained and checked
against original degree-four moments.

Binary64 entry points validate each original scalar before array coercion.
Booleans, strings, missing values and nonrepresentable exact inputs are
refused. Exactly representable integers, fractions, decimals and subnormal
configuration values remain valid inputs. This is a numerical evaluator's
contract, not a restriction on the analytic configuration space.

## Independent replay and audit controls

The two receipt verifiers now share one original-input monomial assembly.
They lift supplied scalars and vertices into a private mpmath context before
products, the gauge rechart or constrained reduction. Their Cholesky-checked
moment calculation is independent of the producer's quadrature/QR algorithm.
Only final reported values are converted to binary64; unrepresentable reports
are refused. The context includes 80 guard digits beyond a scale-based working
budget. This remains a finite numerical replay, not an interval certificate.

The controls exercise each action component, all 13 scalar modes, inertia,
minimizer, absolute density and public cotangent. They include real/imaginary
uniform fields, nonuniform complex fields, signed and zero charge, changes of
coordinate frame, concentrated fields, and large neutral or scaled-charge
fields whose unused potential overflows. The independent replay also rejects
indefinite reduced metrics with positive determinant and checks nonuniform
dressing monomials separately. Representation tests challenge mixed numeric
containers before conversion. The Windows/Linux workflow runs these controls
with warnings treated as errors.

The final peer audit also challenged a tiny physical tangent hidden behind a
large pure-gauge tangent. For uniform matter one and gauge vector
`xi=(0,1,-1,0,...,0)`, supply edge velocity `D xi`, imaginary scalar velocity
`xi/4`, and an additional center velocity `1e-100`. The exact recharted
velocity is just that last entry. Binary64 gauge solving instead left
`O(1e-16)` noise; a working-precision budget based only on the largest input
also lost it. Public binary64 packet evaluation now refuses that unresolved
cancellation, while direct tiny velocities through `1e-200` remain valid.
This additional normwise gate includes the gauge solve's condition number;
it is not a componentwise error certificate. Independent replay budgets from
the full original exponent span and numerator/denominator bit lengths, before
any products, and retains the `1e-100` signal. It also preserves near-one
Fraction/Decimal differences. Exact-zero comparisons use an absolute
roundoff bound rather than declaring numerical solve residuals exactly zero.

## Evidence impact

The state, trial-history, packet and real-continuum receipts are live artifacts
and must be regenerated with their canonical producers. The state and packet
floating diagnostics can change at rounding scale for their ordinary retained
samples. Canonical regeneration of the independently computed real-sector
receipt also changes floating values at roundoff scale. Exact state moments,
history bounds, packet projection formulas and the real-sector equations and
initial data are unchanged. History pins
the refreshed state receipt, so generation order matters. The source-current
inventory verifies unchanged: it inventories the Lean Whitney modules, but
does not include these numerical electromagnetic scripts.

The postdiction, observation and premise ledgers are checked against the
refreshed evidence. Classical parent trajectories, frozen registrations,
claim payloads, Lean proofs and paper/book sources are unchanged. The papers
already distinguish exact positivity and analytic domain results from finite
quadrature diagnostics; the implementation repair does not strengthen those
scientific claims.

## Validation record

The final nine Whitney test modules pass with `-W error`: 569 tests on Linux
(Python 3.12), and 568 on Windows (Python 3.13) with one expected skip because
Windows `longdouble` has no precision beyond binary64. They comprise the
interacting-quantum, state, history, packet, neutral-packet and real-continuum
suites plus the three kinetic audit modules. Each canonical receipt also
passes its standalone independent verifier. Postdiction ledger parity,
observation/premise register checks and both source-current inventory readers
pass without changing their payloads.

Isolated mutations restore the faulty Schur subtraction, remove precision
guards, condition only the final factor, omit gauge reduction, corrupt its
minimizer, use a dense determinant, restore unpaired derivatives, evaluate an
unused potential or reject all valid inputs. The corresponding controls fail.
Seven separate verifier mutations, including the largest-input-only precision
budget, also fail their intended assertions. A numerically equivalent rounded
`J @ R` vertical-factor variant survives the retained resolved cases; no
inaccurate accepted example was established, and it is not counted as a
detected defect. These mutation checks are evidence against specified
failure mechanisms, not a completeness claim.

Canonical regeneration changes retained state/packet diagnostics by at most
`4.3e-14`/`5.0e-14`, respectively, and real-continuum floating values by at
most `3.6e-15`. The trial-history payload changes only source pins. Source
hashes have been checked against Git's normalized committed bytes as well as
the local working files.
