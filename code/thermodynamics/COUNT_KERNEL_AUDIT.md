# Exact equilibrium and structure of the retained count kernels

The collar realization probe derives finite Markov kernels from the supplied
weighted count table. It checks the unique stationary law, communicating
classes, periods, reversibility and coordinate lumpability exactly. This fixes
four reproductions on main `b8714c98`, retained before the implementation repair:

| Input | Previous result | Correct result |
| --- | --- | --- |
| Two-state transition probabilities `epsilon` and `2 epsilon`, `epsilon=2^-48` | After 2,048 steps, an almost uniform vector passes the `1e-12` stationary residual test | The stationary law is `(2/3,1/3)`; a residual must be interpreted relative to the mixing scale |
| The same chain, initially uniform | Relative entropy to the purported equilibrium is about numerical zero | It is `log(9/8)/2`, approximately `0.0588915` |
| Loopless chain `0 -> 1`, `1 -> {0,2}` equally, `2 -> 0` | Not aperiodic because it has no self-loop | Return walks of lengths two and three give period one |
| Integer masses `[[2^54-1,1],[2,2^54-2]]` | Aggregation first rounds `2^54-1` to `2^54` | Preserve the original masses and normalize their exact sums |

`exact_count_chain.py` owns the finite algebra. The probe delegates its graph
operations and equilibrium solve to it; fixed-count power iteration and its
absolute residual gate are removed. Fifteen syntactic coordinate maps reuse
the calculations for their four distinct partitions.

## Numerical and source contract

- Rational and integer inputs retain their exact values. Binary floating-point
  weights retain their represented values, including subnormals. Validate each
  original scalar before mixed containers can round integers, coerce Booleans
  or discard masked data. Weights must be finite and nonnegative; matrices are
  nonempty and square. Empty rows require an explicit absorbing convention.
- The retained `counts` are **accumulated binary64 weighted masses**, with weight
  field `transition_history_mean_modal_mass`. They are not exact integer event
  counts. Rational conversion cannot recover rounding during acquisition or
  accumulation, nor does it establish sampling uncertainty or a physical law.
- The two pinned exported matrices must exactly reproduce their original
  binary64 constructions: row-normalized `C` and `(C+C.T)/2`, including the
  declared absorbing completion. The subsequent rational kernels are explicitly
  normalized models of the stored count masses. They are not assertions that
  rounded exported rows sum exactly to one. Export support must agree with the
  count kernels; a positive edge erased by export underflow is refused.
- Source hashes, binary64 representations, state labels and their ordering are
  checked. Exact stationary laws solve the original balance equations and
  normalization. A unique closed class permits a unique law, including zero
  mass on transient states; multiple closed classes require separately supplied
  class masses and are refused by the unique-law solver.
- Period comes from the GCD of directed support-cycle lengths, evaluated using
  breadth-first levels. A self-loop is sufficient for period one but is not
  necessary. The [NetworkX reference](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.dag.is_aperiodic.html)
  describes this standard finite-graph criterion. Reducible chains have no
  single reported global period; closed-class periods are reported separately.
- Reversibility and strong lumpability compare exact rational fluxes and block
  sums. Rounded display fields and the retained legacy tolerance diagnostic do
  not decide these properties. Exact rational fields accompany their displays.
- Entropy readouts use 90-digit numerical logarithms on exact finite iterates
  in a private arithmetic context. A positive-term KL formula retains small
  differences; direct probability ratios preserve near-boundary populations.
  Contraction loss is evaluated as the relative entropy between the exact joint
  laws `A_ij=p_i K_ij` and `B_ij=(pK)_j pi_i K_ij/pi_j`. For a faithful stationary
  `pi`, the chain rule gives `D(A||B)=D(p||pi)-D(pK||pi)` without subtracting two
  nearly equal numerical totals. This requires no detailed balance and retains
  genuine zero loss in periodic chains. A nonstationary reference is refused.
  These samples are not interval proofs. The archived spectral-gap estimate is
  identified as inherited and is not used to establish equilibrium. Oversized
  rational receipt strings remain subject to the interpreter's decimal-digit
  limit and raise a serialization error; no process-wide limit is disabled.

## Retained result and scientific impact

The live probe receipt advances to schema v4. Its exact normalized count models
preserve the previous qualitative result: the eight-state repair-load
aggregation is irreducible, period one and nonreversible, and it fails the
fine-chain strong-lumpability test. The fine chain has one closed class, the
absorbing state at index 12. The record-family charge is constant. This does
not identify a physical common reference, energy or clock, or close
THERMO-REALIZATION. Other quotient classes remain outside this bounded audit.

Numerical summaries now use the exact count model; tiny rounding-level changes
in the retained stationary and entropy diagnostics are expected. No pinned
source bytes, frozen prediction rules, paper claims, Lean premises or registry
payloads change. The separate conditional-repair thermodynamics certificate
and protected-memory witness remain unchanged.

## Reproduction and independent controls

```sh
python -W error -m pytest -q code/thermodynamics
python code/thermodynamics/collar_matrix_realization_probe.py --write
python -m pytest -q tools/test_mandatory_suite_workflow.py tools/test_post_r2029_audit_surfaces.py
```

The independent controls enumerate directed spanning trees instead of solving
the producer's balance equations. They cover all 144 strongly connected
three-state support graphs, with periods obtained by explicit simple-cycle
enumeration; periodic and slowly mixing chains; unique and nonunique stationary
laws; arbitrarily small flux/lumpability defects; and supplied-weight precision.
For the retained table, symmetric weighted degrees independently determine the
reversibilized law, and integer Laplacian cofactors evaluated by the Leibniz
determinant give the eight-state stationary law. Source-export mutations include
one-ulp changes and added minimum-subnormal edges. Valid sources and nonzero
small effects are positive controls alongside the refusal cases.

The new controls execute on both operating systems in standard and manual-full
CI, and on Linux in the nightly workflow. The audit additionally retains slow
positive entropy losses at `epsilon=2^-400` and `2^-1000`, checked against
independent 500-digit closed-form trajectories. Those losses were reported as
zero on the first PR head despite its correct equilibrium. Paired public
classification controls also reject tolerance substitutions that survived the
initial helper-only tests. This is one finite-evidence repair under the standing audit
[issue #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
