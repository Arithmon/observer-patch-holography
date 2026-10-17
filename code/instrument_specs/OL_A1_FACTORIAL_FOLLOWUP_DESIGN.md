# OL-A1 follow-up: carrier count by absolute support size

This is an unfrozen campaign
design, not a reverse-engineering-reality result document.

Status: design only. This document does not freeze a campaign, draw a seed,
authorize a run, or change the `FAILED` verdict of INS-01. Any execution needs
a new immutable preregistration and a separate freeze event under the owning
simulation-instrument lane #737.

## What the first campaign did and did not show

The conformant INS-01 campaign returned `FAILED` overall, but its component
results should not be collapsed into one undifferentiated statement:

- The 65,536-carrier A2 arm reproduced the threshold `(1,3)` form in five of
  five replicates and matched its preregistered robust reference in five of
  five replicates.
- The 16,384-carrier A1 arm reproduced the threshold form in two of five
  replicates and matched its robust reference in one of five. This is direct
  evidence of small-rung or seed sensitivity, not evidence that the small rung
  is known to be the cause.
- The A2/A1 cone-margin ratio missed the declared band in all five paired
  replicates.
- The ancestry-permutation null did not discriminate reliably. The current
  estimator therefore did not attribute the form to ancestry structure.
- The low-support control greatly reduced cross edges, but the frozen firing
  rule did not fire. That control did not isolate a clean support-size
  mechanism.

It is therefore plausible that the smallest rung was too small, but INS-01
cannot establish that explanation. Carrier count changed together with
observer count, and only one low-support control was included. A factorial
follow-up is needed to separate finite-size, support, and ancestry effects.

## Candidate factorial ladder

The next freeze should vary carrier count and absolute observer support size
independently, while holding observer count and every other producer and
analysis choice fixed within each paired seed block. "Support size" here means
the integer number of carriers read by one observer. The derived support
fraction is support size divided by carrier count; because that fraction
changes along the carrier-count axis, this document does not call the second
axis "support density." A reasonable design envelope is:

| Axis | Candidate levels | Purpose |
| --- | --- | --- |
| carrier count | 16,384; 65,536; 131,072; 262,144 | retains the seed-sensitive rung and adds three larger rungs, so the small-rung hypothesis is directly testable |
| absolute support size | 48; 96; 192 | low, reference, and high absolute support at every carrier count; support fraction is reported as a derived quantity |
| observer count | one value fixed at the future preregistration, provisionally 256 | prevents observer-count scaling from being confused with carrier scaling |
| fresh replicate blocks | determined prospectively before freeze | chosen from declared minimum effects, confidence-interval widths, multiplicity, and attrition rather than assuming that eight or any other convenient count is sufficient |

These numbers are design candidates, not frozen values. Before a freeze, a
runtime-only pilot may measure cost but must not expose signature, margin, or
control outcomes. A target-blind prospective calculation must set the number
of replicate blocks using predeclared minimum scientifically relevant main
effects, interactions, ancestry contrasts, equivalence widths, familywise
error or false-discovery control, desired power, and an attrition allowance.
It must publish the calculation and its inputs; it may not use INS-02 outcome
estimates. The final preregistration must either keep the complete factorial
grid or declare a target-blind resource reduction before any outcome is
inspected. This ladder carries no authorization to execute.

Each fresh master seed forms one block spanning every retained
carrier-count-by-support-size cell and every matched control. The freeze must
pin a deterministic, collision-free derivation from the master seed to each
cell seed, randomize execution order within each block, and analyze the paired
contrasts at block level. Seed redraw, cell substitution, optional stopping,
and analysis that treats paired cells as independent are forbidden.

## Stronger ancestry controls

A mere label permutation may preserve the structures used by the estimator.
The next campaign should preregister all of the following:

1. An ancestry-destroying, degree-preserving rewiring that preserves carrier
   count, degree sequence, support counts, and one-point record marginals while
   breaking parent-child and shared-lineage correlations.
2. A depth-stratified lineage shuffle that preserves the population at every
   declared depth while breaking cross-depth family identity.
3. A sham relabeling that preserves the complete ancestry graph. The analysis
   must remain invariant under this negative control.
4. A sensitivity control with a synthetic ancestry-dependent perturbation of
   known sign and predeclared magnitude. It checks whether the fixed estimator
   can recover an ancestry effect when one is present; it has no promotion
   authority.

Every control must run through the same observable code path and the same
fixed estimator as its matched main cell. A separate structural verifier may
check the invariants of a transform, but it may not replace the estimator in
an inferential gate. Before freeze, the preregistration must define one
continuous ancestry-sensitive estimator output and quantitative thresholds for
all of the following gates:

- **destruction:** each destructive transform reduces the declared lineage
  association metric by at least a fixed amount while its preserved graph and
  marginal quantities remain within fixed tolerances;
- **estimator contrast:** the paired actual-minus-destroyed contrast has the
  declared sign and its multiplicity-adjusted interval clears a minimum effect;
- **sham equivalence:** the paired actual-minus-sham interval lies wholly
  inside a predeclared equivalence band;
- **synthetic sensitivity:** the same estimator recovers the injected effect's
  sign and clears a predeclared minimum fraction of its known magnitude.

The thresholds, interval construction, multiplicity family, estimator code,
and missing-data rule must all be frozen. Failure of destruction, sham, or
synthetic-sensitivity gates prevents an ancestry-specific positive conclusion,
even if a threshold signature appears in the main cells.

## Candidate discriminating estimator

INS-01's ancestry control was a binary indicator: the permutation null either
destroyed the (1,3) threshold form or it did not, and the fourth eigenvalue
behind that form sits at relative magnitude near 3e-5 against a 1e-12 cut. A
statistic of that shape cannot separate the architecture arm from its ancestry
null when both arms land on the same side of a fragile threshold, and that is
the recorded INS-01 failure mode: the null did not discriminate. The follow-up
therefore requires one continuous, within-run-calibrated estimator whose
separation property is stated before freeze. The candidate defined here is a
design object. Every constant in it is a named placeholder until a
target-blind freeze binds it, and nothing in this section draws a seed,
authorizes a run, or changes the INS-01 verdict.

### Estimator definition

Fix the predeclared retained dimension `M_RUNGS` (candidate value 4, the rung
count behind the (1,3) form; the freeze binds the value). For one run, let
lambda_1, ..., lambda_n be the eigenvalues of the fitted bilinear form, the
same matrix whose eigenvalue signs give the threshold inertia, sorted by
absolute value in descending order. The eigengap-ratio statistic is

    G = log( |lambda_{M_RUNGS}| / max( |lambda_{M_RUNGS + 1}|, EPS_EIG * |lambda_1| ) ),

with `EPS_EIG` a predeclared positive relative floor that keeps G finite on a
rank-deficient fit whose tail eigenvalue is exactly zero. The floor is a
numerical guard applied identically to the observed run and to every null
draw; it is not an inferential cut.

G is invariant under overall rescaling of the spectrum, consults no absolute
eigenvalue threshold, and measures exactly the separation whose fragility
drives the INS-01 eigenvalue caution: retained rungs against discarded tail.
A verdict built on G cannot flip because a single eigenvalue crosses a fixed
cut; it moves continuously with the spectrum.

Producer-side schema requirement, named as a freeze precondition: the
committed INS-01 receipt schema emits exactly four eigenvalues per fitted
form, which equals `M_RUNGS` at the candidate value, so lambda_{M_RUNGS + 1}
does not exist in the present producer output and G cannot be carried by the
present schema. Before the freeze, the INS-02 producer must either emit at
least `M_RUNGS` + 1 eigenvalues per fitted form in the receipts or retain the
raw feature matrices the custody section requires, so that G is computable
and independently recomputable. This requirement is explicit here rather than
left as an implicit consequence of the custody section.

### Within-run permutation-null calibration

For the same run, draw `N_PERM` independent realizations of the strongest
declared destructive transform (the ancestry-destroying, degree-preserving
rewiring, control class 1 in this document), recompute the fitted form through
the identical observable code path, and obtain G_null_1 through
G_null_{N_PERM}. Two calibrated outputs per run:

- the standardized score
  T = (G_obs - median_j G_null_j) / max(MAD_j G_null_j, `EPS_SCALE`),
  with the median absolute deviation scaled for consistency at the normal
  distribution and `EPS_SCALE` a predeclared positive floor;
- the within-run exceedance
  p = (1 + #{ j : G_null_j >= G_obs }) / (1 + N_PERM).

The null distribution is computed inside the run it calibrates, so the
comparison never crosses seeds. The seed-to-seed spectrum drift that defeated
the cross-run margin-ratio band has no path into this null.

Null-dispersion validity gate: if MAD_j G_null_j falls at or below the
`EPS_SCALE` floor, the null is degenerate for that run. In that regime T
divides numerical noise by the floor and can grow without bound while the
exceedance p simultaneously saturates under ties, and the two candidate
inference rules below can disagree on identical data. A run whose null MAD
attains the floor is therefore voided for inference, and a cell is
ESTIMATOR-INVALID when the fraction of its blocks voided this way exceeds
the predeclared `MAX_VOID_FRAC`. This gate is part of the design and its
placeholder freezes with the others; it exists precisely because a null that
preserves the fitted form's spectrum nearly exactly is the fragile-spectrum
failure class this estimator is meant to eliminate, and such a null must
surface as an invalid instrument rather than as a false positive or a false
negative.

The sham relabeling (control class 3) preserves the complete ancestry graph.
Under a label-invariant implementation of the observable path, G_sham = G_obs
is the expected outcome; the analysis asserts that equality to the predeclared
tolerance `EQ_BAND` as a pipeline negative control. The gate is not vacuous:
it does real work exactly when the implementation is label-sensitive through
iteration order, tie-breaking, or hashing, and a sham failure is therefore an
implementation-defect finding that voids attribution, not a scientific
negative result.

### Paired-contrast form across factorial cells

Block b and cell c yield T(b, c). Cell-level inference uses the paired block
structure of the factorial ladder: a one-sided location test of the block
values T(., c) against zero (candidate: Wilcoxon signed-rank with a
Hodges-Lehmann interval), with familywise control at `ALPHA_FAMILY` across
retained cells (candidate: Holm). Combination of the within-run exceedances
(candidate: Fisher) is the declared alternative form. Exactly one of the two
rules survives to the freeze, chosen target-blind.

### Decision rule per cell

- DISCRIMINATES: the multiplicity-adjusted lower confidence bound for the
  block location of T(., c) clears `DELTA_MIN` and the sham-equality,
  synthetic-sensitivity, and null-dispersion gates all pass.
- NO-DISCRIMINATION: the multiplicity-adjusted upper bound lies below
  `DELTA_MIN` while the sham-equality, synthetic-sensitivity, and
  null-dispersion gates all pass.
- ESTIMATOR-INVALID: the sham or synthetic gate fails, or the voided-block
  fraction under the null-dispersion gate exceeds `MAX_VOID_FRAC`; the cell
  issues no attribution verdict in either direction.
- INCONCLUSIVE: none of the conditions is met.

### The discriminating property, stated exactly

Two probability spaces carry the statement, and both must be named because
the observable path is deterministic once a run capture and a transform draw
are fixed. Within one fixed run capture M, the randomness is the transform
draw pi over the declared class (the pi-law): G(pi(M)) is a random variable
under the pi-law while G(M) is a constant. Across the campaign, the
randomness is the run law induced by the frozen seed derivation, under which
M itself is random.

The property this design tests is per-run dominance, not distribution-level
stochastic dominance between pooled arms: for run-law almost every capture M
at a retained cell, G(M) exceeds the median of the pi-law of G(pi(M)) by a
positive margin delta(M), and the block location of delta is bounded below by
a positive constant. This is a stronger statement than comparing the marginal
law of G under the architecture arm with a pooled null law; the
distribution-level comparison does not by itself deliver the paired
within-run consequences without a coupling assumption, and it is not what
the pipeline computes. Under per-run dominance the block location of T is
strictly positive, the within-run exceedances p are stochastically smaller
than uniform under the pi-law, and the paired block test gains power with
block count; `DELTA_MIN` is then bound as the minimum scientifically relevant
separation, not as a detectability limit.

Failure routes are split by detector coverage, and the coverage claim is
correspondingly narrow. The total-insensitivity route, where the fitted
form's spectrum depends only on quantities the rewiring preserves, collapses
the pi-law onto the point G(M); the null-dispersion gate voids such runs, and
the synthetic-sensitivity control is the predeclared detector for this route:
an injected ancestry-dependent perturbation of known sign and predeclared
magnitude must move T by a predeclared minimum fraction of that magnitude
through the same code path. The partial-failure route, where the pi-law
overlaps G(M) with a positive location below `DELTA_MIN`, is not detected by
the sham or synthetic gates; it surfaces at the cell as NO-DISCRIMINATION or
INCONCLUSIVE, never as a gate failure. A DISCRIMINATES verdict is therefore
statistical evidence, at the frozen confidence level, for the per-run
dominance property at the tested cells at the frozen margin `DELTA_MIN`, and
nothing stronger; it says nothing about untested cells or about margins below
the frozen one.

## Prospective precision and power sizing

The committed campaign retains the manifest, summary, and fifteen per-arm
receipt files and no raw feature matrices, so G and T cannot be recomputed
from committed artifacts. The pilot quantities are therefore restricted to
derived receipt fields, except for cost, which no committed field records and
which only the runtime-only cost pilot may measure:

- between-replicate dispersion of the continuous receipt endpoints: per arm,
  the cone margin and degeneracy ratio carried by each run receipt; across
  arms, the paired A2-over-A1 margin ratios stored at summary level
  (`margin_ratios_a2_over_a1` in the campaign summary), which is a cross-arm
  paired quantity rather than a per-arm receipt field. Together these inform
  the variance-inflation placeholder `KAPPA`;
- replicate-level completion patterns per arm, which inform the attrition
  placeholder `ATTRITION` only insofar as runs failed for runtime reasons;
- per-cell cost and hence `N_PERM` feasibility, which no committed field can
  inform: the committed manifest and receipts record no wall-clock, memory,
  or other cost fields, so these two quantities come only from the
  runtime-only cost pilot permitted by the ladder section.

A runtime-only pilot may refine cost and `N_PERM` feasibility; per the ladder
section, it must not expose signature, margin, or control outcomes.

Effect-size parameterization: the effect is the block location delta of T at a
cell. The within-run calibration standardizes the null scale of T near unity
by construction, so the sizing is parameterized by delta alone rather than by
an estimated null variance. `DELTA_MIN` is a scientific declaration of the
minimum relevant separation, fixed target-blind at freeze; it is not estimated
from any run, committed or future.

Sizing formula per cell: with alpha' the Holm worst-case per-cell level
`ALPHA_FAMILY` / n_cells and `POWER_TARGET` = 1 - beta, the required block
count before attrition is

    B >= KAPPA * (z_{1-alpha'} + z_{1-beta})^2 / DELTA_MIN^2,

and the drawn block count is ceil(B * (1 + ATTRITION)). Exceedance resolution
requires N_PERM >= ceil(1 / alpha') - 1, so the minimum attainable within-run
p lies at or below alpha', with equality exactly when 1 / alpha' is an
integer. The numbers that await the pilot analysis rather than receiving
invented values here: `KAPPA`, `ATTRITION`, per-cell cost, and the feasible
`N_PERM`. The numbers that are target-blind scientific declarations, not
estimates: `DELTA_MIN`, `ALPHA_FAMILY`, `POWER_TARGET`, `M_RUNGS`,
`EPS_SCALE`, `EPS_EIG`, `MAX_VOID_FRAC`, and the sham equivalence tolerance
`EQ_BAND`. This document binds none of them to numeric values.

## Kill band template

The kill band for the discriminating-estimator endpoint, stated as a template
with named placeholders:

> At the largest retained carrier rung, with the sham-equality,
> synthetic-sensitivity, and null-dispersion gates all passing, the
> multiplicity-adjusted upper confidence bound for the block location of T
> lies below `DELTA_MIN` at familywise level `ALPHA_FAMILY`. Endpoint
> verdict: FAILED, no ancestry-specific stabilization, regardless of
> threshold-form reproduction in the same cells.

The complementary bands use the same placeholders: the endpoint returns
DISCRIMINATES when the adjusted lower bound clears `DELTA_MIN` at every
retained carrier rung at or above the declared reference rung,
ESTIMATOR-INVALID when a sham, synthetic, or null-dispersion gate fails
anywhere, and INCONCLUSIVE otherwise. Binding rule: every placeholder
(`M_RUNGS`, `N_PERM`, `EPS_SCALE`, `EPS_EIG`, `MAX_VOID_FRAC`, `DELTA_MIN`,
`EQ_BAND`, `ALPHA_FAMILY`, `POWER_TARGET`, `KAPPA`,
`ATTRITION`, and the drawn block count) freezes to a numeric value before any
seed draw. A seed drawn while any placeholder is unbound voids the campaign
under the owning lane's freeze discipline, and a frozen band admits no
post-draw revision. A threshold-form reproduction in the main cells cannot
rescue a kill-band FAILED on this endpoint, and a kill-band FAILED on this
endpoint cannot erase a separately reported component pass, consistent with
the outcome separation stated in this document.

## Design-only status of the estimator sections

The estimator, sizing, and kill-band sections are design only, on the same
footing as the rest of this document: no seed is drawn, no run is executed or
authorized, no freeze exists, and no placeholder is bound. INS-01's FAILED
verdict stands as the controlling completed verdict, OL-A1 remains `owed`, and
nothing in these sections modifies either state. Execution requires an
immutable preregistration and a separate freeze event under the owning
simulation-instrument lane #737.

## Separate questions and outcomes

The follow-up must report the following separately before applying any overall
verdict:

- threshold-signature stability at each predeclared threshold;
- the continuous normalized eigenvalue spectrum and degeneracy margins;
- carrier-count main effect, absolute-support-size main effect, and their
  interaction, with the derived support fraction reported but not relabeled as
  an independently varied axis;
- ancestry-destroying effects versus sham effects;
- the discriminating-estimator endpoint (block location of T per cell,
  sham equality, synthetic recovery, null-dispersion validity) as its own
  endpoint, distinct from threshold-signature stability;
- cone-margin scaling and the old ratio-band result as a distinct endpoint;
- cell-level pass, fail, or inconclusive status with multiplicity and missing-
  data rules fixed in advance.

The future frozen decision rule should distinguish at least these scientific
outcomes: size-supported stabilization, support-sensitive stabilization,
ancestry-specific stabilization, no discrimination, and unresolved. A large
rung passing cannot overwrite a failed ancestry control or failed margin
endpoint. Conversely, a failed small rung cannot erase a reproducible large-
rung component result. Only a separately frozen overall rule may return
`REPLICATED`, `FAILED`, or `INCONCLUSIVE` for the emergent OL-A1 claim. Even
`REPLICATED` is not an automatic ledger promotion: it can support OL-A1 only
for the exact frozen run configuration, after an independent verifier
recomputes the declared outcomes from retained captures, and after the
scientific register explicitly selects the successor as controlling. A
controlling `FAILED` verdict blocks or demotes OL-A1 to `owed`.

## Custody and stopping rule

The final campaign must retain raw feature matrices or sufficient source
captures for independent observable recomputation; derived receipt fields
alone are insufficient. It must pin the exact simulator commit, manifest inventory, code,
configuration, seeds, environment, captures, receipts, and independent
verifier before comparison. No adaptive seed redraw, optional stopping,
post-hoc threshold, or selective cell omission is allowed.

Until that freeze and execution happen, INS-02 remains `SPECIFIED`, OL-A1
remains `owed`, and INS-01 remains the controlling completed verdict. INS-02
does not supersede INS-01 merely by existing as a design.
