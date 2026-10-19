# Observer restriction across the source federation

The elementary object stays a twelve-port carrier. The retained source joins
20, 80 and 320 copies on successive cell subdivisions, giving 240, 960 and
3,840 ports. This audit constructs a compatible **static site-observable
family**, executes its occupation readout from retained currents, and measures
the failure of a proposed dynamical refinement map. It supplies exact finite
premises and a finite-horizon error bound. It does not establish a continuum
limit or exclude a limit on another, suitably controlled observable domain.

Run the subsection with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.federated_observer.refinement
python -m pytest -q -W error model/federated_observer/test_refinement.py
```

## The actual parent map and its domain

The pinned builder's `cell_ancestor_ids` determines which coarse cell contains
each fine cell. The executable checks its four-child ordering. Extending that
ancestry to ports by retaining each local port label is a **declared readout
choice**:

\[
 p(v)=12\,\mathrm{ancestor}(\lfloor v/12\rfloor)+(v\bmod12).
\]

It is not a derived identification of local tangent frames. The source uses
edge-neighbour cells, not the different shared-vertex graph in
`code/source_routing_refinement/hierarchy.py`. Every fine edge either stays in
one parent cell or maps to an actual coarse cell edge. The stronger statement
for individual ports fails: fine seam `(1,38)` maps to coarse pair `(1,2)`,
which is not a coarse seam. The executable checks all source edges, without
replacing the source graph by a more convenient hierarchy.

| Refinement | Fine edges | Edges mapping to a non-edge of the coarse port graph | Edges mapping onto a single coarse port | Coarse edges lacking a direct fine-edge lift |
|---|---:|---:|---:|---:|
| 0 to 1 | 2,520 | 74 | 2 | 0 |
| 1 to 2 | 10,080 | 280 | 4 | 16 |
| 0 to 2 | 10,080 | 400 | 4 | 8 |

These port-edge failures leave the cell subdivision intact. They identify the
additional routing compatibility a port-level dynamical tower would need.

## A compatible static observable family

Let \(A_l=\mathbb C^{V_l}\) be the algebra of functions on the **port labels**.
For a refinement with \(b=4^{l_f-l_c}\) descendants per parent port, define

\[
 (Jf)(v)=f(p(v)),\qquad
 (Rg)(u)=\frac1b\sum_{p(v)=u}g(v),\qquad E=JR.
\]

Disjoint, equally sized fibres prove that `J` is a unital injective
star-homomorphism, `R J = I`, and `E` is a positive unital idempotent
conditional expectation. Its bimodule identity is
\(R[(Jf)g]=f\,Rg\). The uniform site state
\(\omega_l(f)=|V_l|^{-1}\sum_v f(v)\) satisfies
\(\omega_f(Jf)=\omega_c(f)\); the maps compose through the intermediate
level. These statements are certified for every fibre, with exact rational
weights. They use the uniform site measure, not the builder's spherical-area
measure.

This commutative site algebra has trivial modular flow; its uniform reference
is tracial. The compatible maps do not attach the nontrivial reference-state
modular flows on history-transfer matrix algebras, or establish a quantum
local-net limit.

This is not the complete algebra or probability law of the fixed-two-particle
experiment. Fine block occupations take values `0,...,b`, whereas a single
binary coarse port takes values `0,1`. Under the uniform fine two-particle
state, the probability that both particles lie in one parent block is
\((b-1)/(|V_f|-1)\), which is positive; the corresponding single-port event is
impossible in a coarse binary configuration. In particular, static site-state
compatibility must not be substituted for a full field-state refinement law.

## What the dynamics actually do

Write \(L_l\) for the integer graph Laplacian and \(M_l\) for the number of
source seams. The native fair-coin update gives the exact first-moment kernel

\[
 P_l=I-\frac{L_l}{2M_l}.
\]

This identity holds on every fixed-number binary sector; it does not assume
independent occupations. An auxiliary convention in which each seam attempts
at rate one has density semigroup \(T_l(t)=e^{-tL_l/2}\), with expected global
attempt count \(M_l t\). Here `t` is **per-edge attempt intensity**. It is not
a calibrated observer clock, modular parameter or physical duration.

For a declared coarse intensity multiplier \(c\ge0\), the audit computes
every sparse entry of

\[
 D_c=\frac12(L_fJ-cJL_c).
\]

All three level comparisons give
\(\|D_1\|_{\infty\to\infty}=1\); the fourfold coarse-rate control gives
\(\|D_4\|_{\infty\to\infty}=19\). Thus copying port observables does not
intertwine these source generators.

A more direct witness does not depend on any chosen coarse generator. Fine
states `(0,1)` and `(0,13)` have the same two-particle block counts under
the adjacent-level map. Port `1` has an edge to parent block `2`; port `13`
does not. The next-attempt mean occupation of that block differs by
`1/5040` at level 1 and `1/20160` at level 2. The same difference is `1/2`
in its initial per-edge intensity derivative. Multiplying every coarse rate
by a scalar cannot distinguish two identical coarse initial states. The
projection is therefore not an autonomous Markov description on this domain.
The smaller raw-attempt discrepancy at the finer level reflects its larger
number of seams.

There is nevertheless a controlled finite-horizon comparison. Duhamel's
identity gives

\[
 T_f(t)J-JT_c(ct)
 =-\int_0^t T_f(t-s)D_cT_c(cs)\,ds.
\]

The semigroups contract the supremum norm. Residual rows have zero sum, so
an indicator, or any observable between zero and one, costs half the row
absolute sum. Consequently the total variation distance between the
projected fine one-particle transition law and its coarse comparison is at
most

\[
 \min\{1,\tfrac12t\|D_c\|_\infty\}.
\]

For `c=1` this is `min(1,t/2)`. It is a response-kernel bound; no full
many-particle history-distance bound is asserted. The receipt also computes
exact indicator errors at 1, 2, 4 and 8 steps of the common uniformization
kernel \(Q_l=I-L_l/12\), whose step intensity is `1/6`. The first error is
`1/12`, attaining its telescope bound. These comparison steps are not single
native source attempts.

The nonvanishing worst-case defect blocks an *exact* tower with this map and
this domain. It does not rule out convergence for smoother observables,
different coarse variables, effective memory, or a suitably rescaled weak
limit. Such a result needs a declared test-observable class and a vanishing,
controlled error on that class; three finite levels alone do not supply it.

## Occupation readout survives the failure of Markov closure

For a parent block \(B\), its actual occupation after a retained history is

\[
 N_B(n)=N_B(0)+\sum_{r=1}^{n}j_r
       \bigl[1_B(b_r)-1_B(a_r)\bigr].
\]

Here `(a,b)` is the declared oriented seam and `j` the recorded signed
native current. The reader consumes only the disclosed initial block
counts and the complete authentic crossing records. It does not query the
hidden instantaneous configuration. Internal currents cancel; the formula
also restricts consistently to any union of the blocks.

For each comparison the executable checks all eight endpoint/coin cases on
every actual seam. It then executes all 64 fair-coin words of a six-attempt
fixed-two-particle history, including 238 zero-current events, verifies all
384 prefix readouts, and resumes each history from its block-count checkpoint
and event cursor. An omitted boundary-current control produces a wrong
readout. Retaining currents restores exact reconstruction of the ongoing
observation; it does not make the coarse counts alone a Markov state.

## Which existing theorems this connects

| Existing source | Premise supplied or checked here | Boundary |
|---|---|---|
| [`ConsensusTower`](../../Lean/Tower/ConsensusTower.lean) | Composing site embeddings and compatible uniform site states | The full matrix-algebra/observer interfaces are not instantiated; its `generator_natural` premise fails on the tested site subalgebra. |
| [`CPRestrictionNet`](../../Lean/QFT/CPRestrictionNet.lean) | A concrete finite conditional-expectation pattern with positivity, state preservation and tower composition | Our commutative site algebra is not the full regional matrix-factor net. The Lean interface's complete positivity is certified by Kraus witnesses, not just its positivity field. |
| [`SourceRegionalNet`](../../Lean/QFT/SourceRegionalNet.lean) | Separates actual source labels/records from declared algebra and region adapters | Its split-fibre noncommutative adapter is a different finite witness; it does not certify this federation's tensor factors or refinement. |
| [`RegionalContinuity`](../../Lean/Screen/RegionalContinuity.lean) and [`DiscreteGauss`](../../Lean/Screen/DiscreteGauss.lean) | Actual signed endpoint incidence makes the regional balance law an executable observer readout | Initial preparation and complete current access are premises; a continuum field/current identification is separate. |
| [`ObservableNormalForms/Refinement`](../../Lean/ObservableNormalForms/ObservableNormalForms/Refinement.lean) | An explicit finite defect and telescoping error bound illustrate the required quantitative receipts | This dynamics comparison is not an instantiated normalizer theorem, and no summable refinement-tail or inverse-observation hypothesis is supplied. |
| [`LocallyCovariantLimit`](../../Lean/QFT/LocallyCovariantLimit.lean) and [`LimitRestriction`](../../Lean/QFT/LimitRestriction.lean) | The required compatibility squares are tested against the actual source | Their algebraic colimit, regional and time-slice conclusions retain their tower, region, restriction and eventual-inverse premises. They do not infer Lorentzian geometry or a continuum dynamics from a connected graph. |

The concrete next step is a compatible effective observable/history family
for these joined carriers, with a controlled scaling error and the relevant
state, locality and continuation bounds. No enlargement of the elementary
twelve-port object is required by that task.

## A reusable exact modular-compatible federation construction

The paper [*Screen microphysics and observer synchronization*](../../paper/screen_microphysics_and_observer_synchronization.tex)
provides the **Reference multiresolution carrier** definition and
**Fixed-cutoff regulator certificate** theorem. Added detail factors carry
finite algebras \(D_j\) and faithful states \(\tau_j\); a finite-depth local
presentation circuit \(W_r\) expresses their product algebra and reference
state in patch coordinates. In the product coordinates, refinement is
\(A\mapsto A\otimes I\), and restriction is
\(\mathrm{id}\otimes\tau_{s\setminus r}\). Conjugating both maps by the
declared circuits gives injective embeddings, state-preserving completely
positive conditional expectations and exact composition. It also gives

\[
 \alpha_u^{\widehat\omega_s}(\iota_{rs}(A))
 =\iota_{rs}(\alpha_u^{\widehat\omega_r}(A)),
\]

where `u` is the reference state's modular parameter. All \(D_j\) may be
copies of the same finite algebra: larger regulator algebras can consist
entirely of federated fixed-size factors. Faithfulness and product-state
factorization supply this modular compatibility automatically. The theorem
does not require enlarging the primitive carrier, nor does it assign an
independent physical factor merely because a screen cell exists.

Applying this construction to the executed native observer requires an
explicit factor/interface identification, the actual reference-state choice
and presentation circuits, and a match of source histories, allowed
instruments and retained records under refinement. The operator and
conditional-expectation identities alone do not provide that match. In
particular, selecting \(W_r\), selecting a physical reference state, proving
generator/instrument compatibility, and controlling continuum observation
tails are additional obligations; our measured direct-parent defect remains
valid. The paper's **Reference state, ground state, and vacuum** and
**Regulator-to-continuum identification rule** remarks explicitly require
positive-transfer/reflected-Gram evidence and the corresponding limit
identification. The exact reference construction is a concrete route to
reuse, not a substitute for those source attachments.

Unitary presentation circuits preserve the reference spectrum; they cannot
convert an arbitrary native history state into a product reference or replace
its transfer-logarithm dispersion by a square-root dispersion.
