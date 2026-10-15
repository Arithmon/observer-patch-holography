# Local action recovery while the entire field participates

The [previous interacting inverse](../interacting_reconstruction/DERIVATION.md)
identified a supplied local quartic scalar action by isolating each onsite
and pair generator. This construction removes that isolation capability.
Every calibration trial now runs the complete eight-site onsite/bond source
schedule; only the initial number preparation and terminal local detector
change. The inverse still receives detector laws, opaque labels and separate
phase-clock laws, without backend matrices or coordinate labels.

The mathematical addition is a derivative certificate for a fixed family of
source products. It identifies the sum generator from local marginal
observations even though the short finite source schedule is not its exact
continuous exponential. The supplied canonical action family, oscillator
basis, preparation and measurement services, time scaling, and topology gap
remain premises.

## Local probes of the complete action

The same finite scalar family is used:

```text
H = sum_i [p_i^2/2 + A_i q_i^2/2 + g_i q_i^4]
    - sum_(i<j) kappa_ij q_i q_j,
A_i = m_i^2 + sum_j kappa_ij.
```

All eight sites retain four occupation levels during calibration. The
polynomials `p^2`, `q^2` and `q^4` are formed in the infinite oscillator
before compression. The supplied basis frequency is `Omega=2`. Nine real
product preparations suffice: the product blank, and one occupation-one
excitation at each labelled site. The source contracts them as independent
columns to reduce host work. The resource ledger still counts a fresh
preparation and full field schedule for each binary instrument law; the
columns are not one coherently superposed preparation.

For an onsite probe at `i`, the local bit-one effect projects onto

```text
|e_±> = (|j> ± i |k>)/sqrt(2),
(j,k) = (0,2) or (1,3).
```

The local instrument acts as identity on every spectator. Both outcomes
are retained, and both phase signs use fresh preparations. For a pair
probe the bit-one effect projects onto `(|00> ± i|11>)/sqrt(2)` within
the full sixteen-dimensional two-site register; its complement retains
all other occupations. Spectators are traced out exactly in computing
these laws. This is local marginal observation, not postselection of the
other sites into their starting states.

Write `s(t)=[p_-(t)-p_+(t)]/2` for the difference of the bit-one laws.
At zero duration, the first derivative of the onsite probe is

```text
s_i,jk'(0) = <k|p_i^2/2 + A_i q_i^2/2 + g_i q_i^4|j>.
```

Spectator onsite terms have no first-order local coherence. A bond incident
on the probed site also gives no first-order contribution here because
each spectator number state has zero expectation of `q`. The initial
occupation-one state is real and has the same vanishing coordinate mean.
For the pair probe from the full product blank,

```text
s_ij'(0) = -kappa_ij/(2 Omega).
```

Only the direct pair term connects the initial pair `|00>` to `|11>` at
first order. Indirect paths can generate nonzero finite-duration pair
coherences even when the direct edge is absent. Consequently topology is
inferred from controlled derivatives, not from whether a finite-time
correlation is exactly zero.

## The source family is a product, not an assumed semigroup

At every probe duration the source executes the same symmetric schedule:
all onsite half steps, every bond step, and all onsite half steps again.
All generators participate, including those outside the measured site or
pair. The retained calibration uses one full schedule at each duration.
In general, write its ideal matrix as

```text
U(t) = product_l exp(-i t alpha_l H_l).
```

The chronological factor sequence, matrices `H_l`, and scaling coefficients
`alpha_l` are fixed independently of `t`. Each `H_l` is real symmetric in
the supplied occupation basis. The common parameter `t` scales the model
duration of each factor; it is not host execution time, a pulse counter, or
the sum of all service durations. A separate supplied phase clock reports
this parameter on its declared unwrapped phase window.

The first derivative is

```text
U'(0) = -i sum_l alpha_l H_l = -i H.
```

At finite `t`, generally `U(t) != exp(-itH)` and
`U(t+u) != U(t)U(u)`. Those identities are not needed for calibration.
Termwise complex conjugation preserves the factor order, giving
`U(-t)=conjugate(U(t))`. It follows that the imaginary local coherence
`s(t)` is odd for each real initial preparation. This argument does not
replace `U(-t)` by `U(t)^dagger`; that replacement would be invalid for a
general ordered product.

For each factor remove an arbitrary scalar `c_l I`, which changes only a
global phase and leaves detector laws invariant. Suppose

```text
S >= sum_l |alpha_l| ||H_l-c_l I||.
```

Applying the product rule and summing its multinomial coefficients yields
`||U^(r)(t)|| <= S^r` after this phase removal. The coherence observable
`O=(P_- - P_+)/2`, tensored with spectator identity, has norm `1/2`.
Applying the product rule again to `<psi|U(t)^dagger O U(t)|psi>` gives

```text
|s^(5)(t)| <= (1/2) sum_(r=0)^5 binom(5,r) S^r S^(5-r)
            = (2S)^5/2.
```

This is a bound on the scaled factor family. A bound merely on
`||sum_l alpha_l H_l||` is insufficient: three noncommuting real factors
with generators `X`, `Z`, and `-(X+Z)` have zero sum, but their product has
nonzero higher derivatives and nonzero finite-duration coherence.
The tests retain this cancellation counterexample.

## A certified derivative from two complete field executions

Oddness and the fifth-derivative bound imply

```text
s(t) = h t + b t^3 + r(t),
h = s'(0),
|r(t)| <= (2S)^5 |t|^5/240.
```

For two decoded clock readings `0<u<v`, the generalized Richardson estimate

```text
h_hat = [v^2 s(u)/u - u^2 s(v)/v]/(v^2-u^2)
```

cancels the cubic coefficient of the actual product family. Its error is
bounded by

```text
|h_hat-h| <= (2S)^5 u^2 v^2 (u^2+v^2) / [240(v^2-u^2)].
```

For `u=v/2`, this reduces to `S^5 v^4/18`. The retained protocol uses
`u=1/16384`, `v=1/8192`, and the declared norm-sum cap `S=96`.
The inverse processes the two decoded clock readings individually; it does
not replace them by an exact ratio two. The previous exact source bound
gives a sum below `87` for the default action. That fact is checked only
in the later backend audit and is not used to tune the declared inverse cap.

No separate probability error for approximating `exp(-itH)` is required in
this calibration argument. The target being estimated is the derivative
of `U(t)` at zero, which is exactly the desired sum generator. Its cubic
product error is canceled together with the other cubic terms. Adding a
global splitting error divided by `t` would be unnecessary for this target
and much less effective numerically.

Freezing the schedule is essential. If the smaller duration uses two split
steps while the larger uses one, their cubic coefficients generally differ
and the cancellation fails. The tests give a two-generator example in which
this mixed schedule violates the proposed fourth-order radius despite equal
generator sums and equal norm-sum bounds. Both individual fixed schedules
satisfy the bound.

An absolute error at most `delta` in each probability law contributes the
additional radius

```text
delta [v^2/u + u^2/v]/(v^2-u^2),
```

which is `3 delta/v` in the dyadic case. The ideal Taylor radius decreases
with shorter durations while probability sensitivity increases. Clock
uncertainty, source floating errors and physical service errors are separate;
they are not automatically included by this deterministic probability term.

## Recovering the action and predicting continued records

The resulting derivatives enter the existing finite inverse. With
`u0=H_20/sqrt(2)` and `u1=H_31/sqrt(6)`,

```text
g = Omega^2 (u1-u0),
A = Omega^2 + Omega(10u0-6u1),
kappa_ij = -2 Omega s_ij'(0).
```

Their interval propagation, four-level identifiability, and the two-level
onsite degeneracy are unchanged from the
[previous derivation](../interacting_reconstruction/DERIVATION.md).
The family premise that a pair weight is either zero or exceeds a positive
threshold remains necessary for exact nonedge classification. After all
28 pairs are probed, the inverse tests homogeneity, a common positive bare
mass squared, and the cube graph. Its unit-edge Gram matrix describes the
tested finite spatial relations up to isometry; it does not establish an
emergent metric for arbitrary fields.

The reconstructed action then predicts new coupled two-read histories with
disturbance and record-dependent feedback. These predictions use the inferred
continuous generator and are formed before the corresponding withheld source
executions. For those finite-duration predictions, the source splitting
budget is still needed. The calibration derivative argument does not remove
it. The inherited comparisons remain

```text
inferred finite continuous ↔ executed source:
    TV <= parameter uncertainty + withheld splitting error;
executed source ↔ unrestricted oscillators on the same finite cube:
    TV <= withheld splitting error + occupation error;
inferred finite continuous ↔ unrestricted oscillators on the same finite cube:
    TV <= parameter uncertainty + occupation error.
```

The occupation certificate applies to the particular withheld preparation,
instruments and feedback, not to these new excited calibration experiments.
The low generator matrix elements nevertheless agree across the relevant
compressions within the supplied polynomial family. This is progress from
isolated-service identification to observations of a fully participating
finite field. Endogenous action and clock selection, controlled interacting
spatial and energy limits, relativistic covariance, and independent empirical
validation remain outside this result.
