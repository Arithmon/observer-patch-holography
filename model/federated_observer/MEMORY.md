# Exact coarse fields with memory

The declared parent restriction does not give a Markov process on the coarse
fields: two preparations with the same initial block counts can predict
different next counts. Retaining memory yields an exact finite prediction
from the same native source. This is a first-moment result. It neither
recovers all joint occupation probabilities nor replaces the quantum
instrument by a classical channel.

## The process and the observable

The fine federation has 80 elementary twelve-port carriers, 960 ports and
2,520 seams. One step selects a seam uniformly and executes the existing
native binary endpoint operation with a fair coin. For a mean occupation
profile `v`, its exact one-step transfer is

```text
P = I − L_f / 5040.
```

The parameter `n` counts these global native attempts. It is not an observer
clock, a duration in seconds or the auxiliary uniformized step used in
[the refinement comparison](REFINEMENT.md). Symmetry makes the same matrix
act on site observables and on mean occupation profiles, but their
normalizations still matter.

The actual retained cell ancestry gives four fine descendants for each of
240 parent port addresses. `J` copies a coarse value to those descendants;
`R` averages their values. Thus `RJ=I`, `E=JR` is a projection, and
`x=Rv` is the **mean occupation per descendant**. The expected particle
count in a block is `4x`. Uniform site inner products, with weights `1/960`
and `1/240`, give `R=J*`. These site weights are not an assertion that
the full fixed-two configuration distributions commute with refinement.

Occupation fields become observer records when an initial preparation is
disclosed and all signed incident currents are retained, as in the source
tape and [block-current readback](REFINEMENT.md). The memory calculation
predicts expectations of those fields; it does not measure an evolving
hidden configuration to supply its recurrence.

## Elimination of unresolved fields

Write `F=I−E`, `x_n=R v_n` and `y_n=F v_n`. All four blocks come from the
fine native transfer:

```text
A = RPJ,       B = RP restricted to ker R,
C = FPJ,       D = FP restricted to ker R.

x_(n+1) = A x_n + B y_n,
y_(n+1) = C x_n + D y_n.
```

Substitution of the second equation, by finite induction, gives

```text
x_(n+1) = A x_n + B D^n y_0
          + Σ_(j=0)^(n−1) B D^(n−1−j) C x_j.
```

This identity holds for every finite initial profile and every integer
`n≥0`. The executable checks horizons through eight and rejects requests
beyond sixteen to bound its resource use. It evaluates the convolution
using sparse rational vectors. The kernel factors are calculated from the
fixed graph, the initial unresolved profile and earlier coarse fields;
they never query the future fine trajectory. A separate fine propagation
compares every one of the 240 output coordinates exactly.

`A` is a stochastic coarse transfer obtained from **this fine graph**. It
is not identified with the independent level-zero native transfer. The
memory terms may have signed entries and are not additional transition
probabilities. This predictor still stores fine vectors to apply its
kernel factors; exact closure does not imply a computational cost
independent of the fine graph or a finite sufficient memory window.

The graph Laplacian satisfies `0≤L_f≤2 d_max I`, with `d_max=6`. Hence
`419/420≤P≤I`. Orthogonal compression gives `0≤D≤I`, and self-adjointness
gives `B=C*` in the weighted inner products above. Each memory kernel
`B D^k C = C* D^k C` is positive semidefinite as a quadratic form. This
does not mean that its individual entries are nonnegative or that its
decay is uniform under refinement.

If `y_0` is not disclosed, it cannot simply be inferred from `x_0`.
The exact equation retains the forcing term `B D^n y_0`. Contraction
gives two distinct bounds:

```text
||B D^n y_0||_coarse ≤ ||y_0||_fine,
||R P^n(v_0−J R v_0)||_coarse ≤ ||y_0||_fine.
```

The second bounds the complete forecast error from replacing the initial
profile by its conditional block average while retaining the memory
kernel. It is not a bound for additionally dropping the memory kernel.
The norms are normalized root mean squares. A measured preparation rule
or a separately justified uncertainty bound is still necessary to use
these statements as experimental predictions.

## Executed preparations and results

The main preparation chooses uniformly among the six fixed-two
configurations within one four-port descendant block. Its initial mean
is `1/2` on those four ports and zero elsewhere, so `y_0=0` without
assuming a four-particle preparation. All arithmetic is rational.

| Calculation | Exact result |
| --- | --- |
| Full memory prediction versus fine evolution, attempts 0 through 8 | Every coarse coordinate agrees exactly |
| Memoryless `A^n x_0`, one attempt | Agrees |
| Memoryless prediction, two attempts | Maximum mean-occupation error `1/135475200` |
| Memoryless prediction, eight attempts | Maximum error `350025677730509483674478911/1705313536899783721589145600000000` |
| Total expected particles | Exactly two at every checked horizon |

The errors are small on this attempt scale because only one of 2,520
seams is selected at each attempt. Their nonzero exact values establish
the failure of this memoryless reduction; they do not establish a large
physical effect or an empirical discrepancy.

A second control uses the fixed-two states `{0,1}` and `{0,13}`. Their
coarse initial profiles coincide. After one native attempt, the expected
count in parent block 2 differs by `1/5040`; its mean occupation differs
by `1/20160`. Both memory predictions agree exactly with the fine
propagation through four attempts when their initial unresolved terms
are retained. Dropping that term loses the distinction immediately.

The tests independently enumerate all native seam/coin outcomes of the
six-state preparation, check the convolution indices against explicit
powers, and use an independently assembled edge-difference propagation.
They also check an equal-level control with identically zero memory,
weighted adjoints, hidden-domain guards and exact input types.

## Attachment and remaining limit problem

The nearest existing block-elimination material is the Schur construction
in [the local-tick packet](../local_tick/README.md) and
[`SourceLaw.markov_reader`](../../Lean/SourceLaw/IdentificationLadder.lean).
The latter identifies a reversible **Gaussian** chain precision under
its stated invertibility hypotheses; it does not identify this
non-Gaussian exclusion process. The recurrence above is a direct finite
block identity with its short induction supplied here, rather than an
unproved attachment of the Gaussian theorem or a newly formalized Lean
result.

This construction gives an exact finite coarse predictive description
despite nonlumpability. A continuum description still needs compatible
multi-level observables and preparations, control of the memory kernels
and unresolved terms in a declared scaling regime, and convergence of
the relevant multi-time laws. Relativistic field locality and a physical
spacetime dimension do not follow from this first-moment identity.
