# Finite-state entropy and the role of sector labels

This shared implementation replaces duplicated density-matrix, entropy,
partial-trace and modular-splitting operations in the MaxEnt, collar-alignment,
null-net and Einstein-closure evidence. It uses natural logarithms and the
ordinary matrix trace. It implements finite numerical diagnostics, not a
physical state selection or a proof of exact rank from approximate matrices.

The [entropy-response audit](../geometry/ENTROPY_FIRST_LAW_AUDIT.md) separates
tangent first laws from finite entropy changes with relative-entropy
remainders. `entropy_response` tests independently supplied generators on
all normalized tangents, supplies central probability-transfer witnesses,
and computes Gibbs tangent slopes separately from finite secants.

The follow-up [algebra and repair audit](ALGEBRAS_AND_REPAIR.md) fixes regional
separation and basis normalization, replaces dense standardness matrices with
reduced spectra, and classifies exactly when finite relative-entropy
completion is an affine quantum repair. It includes independent recognition
of supplied conditional expectations and a noncommuting repair convergence
proof.

`expectations.repair_generator` accepts rates in any common clock unit.
Its identity and spectral-resolution checks use the generator divided by
`max(rates)`; the returned generator, gap and defect magnitudes retain the
input rate units. The report's `rate_scale` is this divisor, and `tolerance`
applies to defect magnitudes divided by it. Unresolved relative gaps, lost
positive rates, and unrepresentable rescaled generators still raise; a small
absolute rate alone does not imply poor numerical resolution.

## Audit findings and corrections

The audit used main commit `0f7aa44259f53e4f5a1f4cd4904758039975755e`.

| Finding in the previous implementation | Correction and retained regression |
| --- | --- |
| `markov_chain_state(2,2,2,2,[.5,.5],seed=3)` returned an ordinary mixture with `I(A:C\|B)=0.04714088603973687` nats, despite claiming a multi-sector Markov state. The purported central label was absent. | Construct the direct sum and retain its label in the collar. The corrected state has zero CMI to numerical precision. Tracing out that label reproduces the old non-Markov control. |
| Empty or zero-weight block families could report perfect alignment; NaNs and negative states were not validated. | Require a nonempty normalized block family, finite normalized positive states, and consistent dimensions before any criterion is evaluated. |
| `takesaki_defect(...,n_probes=0)` returned zero for a nonaligned state whose modular defect exceeded 2.12. Random probes could not justify the claimed all-observable criterion. | Remove the sampling arguments. Test every matrix unit of the left algebra; linearity supplies completeness. |
| Relative entropy of orthogonal pure states was reported as `690.7755278982137` rather than infinity. A negative density also returned a finite result. | Validate the states and preserve reference zeros. Detected support escape returns infinity. Full-algebra modular logarithms reject singular states. |
| Gibbs construction silently ignored unmatched constraints, multipliers or central-sector energies; the projection routine could return unconverged parameters or nonunique multipliers. | Validate complete pairings, independent constraints modulo identity, faithful targets, and convergence. Use one common Gibbs energy shift across sectors and reject spectral underflow instead of flooring it. |

Positive eigenvalues are no longer discarded below `1e-14` or replaced by a
`1e-300` floor. The shared `1e-12` tolerance only admits floating-point error
in Hermiticity, normalization and negative eigenvalues. Accepted Hermiticity
roundoff is symmetrized and negative spectral roundoff is set to zero for
spectral calculations. Inputs are never implicitly normalized. No positive
support leakage is rounded away to make relative entropy finite. Relative
entropy rejects unresolved support: a dense reference with positive
eigenvalues on the eigensolver roundoff scale, a computed kernel that does
not annihilate its reference, or leakage indistinguishable from matrix-product
cancellation. It neither assigns false infinity nor drops that leakage to
force a finite answer. A diagonal reference supplies its spectrum directly;
explicit positive entries as small as `1e-310` remain supported. These are
numerical guards, not a rank certificate: exact-support claims still require
exact or separately certified support data.

The follow-up audit of this refactor retains four additional controls:

- A real Hadamard change of basis applied to spectra `(.5,.25,.25,0)` and
  `(.25,.5,.25,0)` preserves common support and gives `D=.25*log(2)`.
  Kernel roundoff previously returned infinity. The routine now returns the
  correct value when support is resolved or raises an explicit support error.
  A companion state with mass in the missing direction must give infinity
  or the same explicit uncertainty error, never a finite divergence.
- A smallest positive binary64 sector weight multiplying `I_2/2` vanished.
  Direct sums now reject component underflow; explicit zero sectors and
  representable tiny sectors remain valid.
- Complex central energies or inverse temperatures produced complex weights
  that were silently cast to real sector probabilities. The scalar validators
  now reject complex, Boolean, nonfinite and nonscalar inputs, including
  convergence and alignment tolerances.
- For `H=30 X`, both Gibbs eigenweights are positive, but reconstruction of
  the dense matrix erases the smaller direction. Gibbs states, sector blocks
  and Duhamel covariance now validate the reconstructed faithful state too.
  Diagonal controls preserve the same small probability when representable.

## One normal form for the finite A3 objective

The existing [A3 specification](../../docs/AXIOM_REFERENCE.md) supplies fixed
positive weights `w_P`, local states `rho_P`, references `tau_P`, a compatible
feasible family K, and declared local traces. Set `W=sum_P w_P`, `p_P=w_P/W`,
and form the labelled direct sums

```text
R(rho) = direct_sum_P p_P rho_P,
T      = direct_sum_P p_P tau_P.
```

Use each declared local trace on its summand. On a positive block,
`log(p_P rho_P)=log(p_P) I+log(rho_P)` on its support. Thus

```text
sum_P w_P D(rho_P || tau_P) = W D(R(rho) || T).
```

Both sides are infinite in exactly the same support-escape cases. This is
an identity for the existing finite A3 data; it adds no new physical premise.
For a nonempty cover W is positive. If an empty cover is permitted, its
injectivity already forces K to be a singleton and the objective is zero;
there is then no remaining optimization to encode.
R is affine, so optimizing on K is equivalent to optimizing on its image
R(K). The A3 cover's injectivity makes the correspondence one-to-one.
The weights, references and compatible feasible image must still be retained;
optimization over the entire direct-sum state space is a different problem.

The label is the cover index, not a constructed joint state on the physical
union of patches. Three pair distributions that make AB, BC and CA each
perfectly anticorrelated have compatible uniform singleton marginals but no
joint classical ABC extension. Their labelled direct sum is nevertheless a
valid encoding of the A3 objective. The test suite checks this obstruction
by complete enumeration of the eight assignments.

This normal form also explains the bulk/edge accounting. For normalized
sector probabilities p,

```text
S(direct_sum_j p_j rho_j) = H(p) + sum_j p_j S(rho_j).
D(direct_sum_j p_j rho_j || direct_sum_j q_j tau_j)
  = D(p || q) + sum_j p_j D(rho_j || tau_j).
```

The second formula has the usual support convention: a positive p block
with zero q has infinite divergence. When p=q, its classical term vanishes.
When p varies, omitting that term changes the selection problem. For
`rho_j=rho_bulk,j tensor I_edge,j/d_j`, the first formula adds
`sum_j p_j log d_j`, with H(p) counted once in the bulk entropy. This is the
identity already used by the Einstein-closure code. An auxiliary cover label
does not acquire physical edge entropy merely by being introduced in software.

These identities follow by evaluating block logarithms and traces; they use
the finite entropy definitions in [Watrous, section 5.2](https://cs.uwaterloo.ca/~watrous/TQI/TQI.pdf).
They reorganize existing OPH mathematics without selecting its physical inputs.

## Why a central label cannot be replaced by a mixture

The [quantum Markov structure theorem](https://arxiv.org/abs/quant-ph/0304007)
uses a decomposition of the conditioning system B:

```text
H_B = direct_sum_j (H_bL,j tensor H_bR,j),
rho_ABC = direct_sum_j p_j (rho_A,bL,j tensor rho_bR,C,j).
```

The sector j is retained in B. Every entropy in
`I(A:C|B)=S(AB)+S(BC)-S(B)-S(ABC)` contains the same H(p); these terms cancel.
The within-sector product terms cancel too, giving zero. Replacing the direct
sum by an unlabelled mixture does not preserve this argument.

The smallest retained control is an equal classical mixture of `A=J=C=0`
and `A=J=C=1`. With J in the collar, `I(A:C|J)=0`. After J is erased,
`I(A:C)=log(2)`. No approximation, simulation fit or optimizer decides this
difference. In multiple sectors, modular splitting must likewise be tested
inside each central block; zero CMI does not imply an unconditioned product
across one globally chosen left/right tensor cut. The numerical cross-lane
test verifies both assertions on the same corrected state.

The API consequence is explicit: `markov_chain_state(da,dl,dr,dc,weights)`
now returns tensor dimensions `(da,len(weights)*dl,dr,dc)`, ordered
`A,(sector,bL),bR,C`. The single-sector interface is unchanged. An application
without access to that sector label must use the erased state and test its
actual CMI; it cannot assume the label exists as a physical repair.

## Downstream scope and reproduction

| Existing evidence | Audit outcome |
| --- | --- |
| Null-net standardness, #524 | Fixes the general multi-sector constructor. The published single-sector modular-locality and nonlocal-Gibbs controls still pass. |
| Collar alignment, #543 | Replaces sampled commutators with a complete basis and rejects vacuous/invalid families. The Bell counterexample remains Markov but not aligned. |
| MaxEnt closure, #539 | Preserves the multiplier counts, nonclosure example and exact product-subfamily control while enforcing their input and convergence premises. The live receipt is regenerated; its numerical changes are at roundoff scale. |
| Einstein closure, #526–#528/#503/#578 | Reuses the same central-label entropy identity. The first-law and normalization-countermodel tests still pass. |

The shared helpers are one implementation, not four independent calculations.
The follow-up [MaxEnt coordinate audit](../maxent/PROJECTION_COORDINATES.md)
fixes dependence on observable units, basis and scalar energy origins, and
connects its invariant stopping error to a bound on the true optimum. It
retains the nonzero coarse-graining closure defect.

The [sector Gibbs audit](../collar_alignment/GIBBS_SECTOR_AUDIT.md) extends
that energy-origin repair to collars. It shares observable decomposition
and thermal diagonalization with MaxEnt, retains relative sector partition
functions, and prevents a large scalar origin from manufacturing alignment.
The audit proves the finite interaction/log-state identity and separates
center-probability error from conditional-state error in relative entropy.
Independent controls use explicit index contraction, scalar closed forms,
SciPy matrix logarithms, complete classical enumeration and known product/Bell
states. Tests include malformed matrices, missing sectors, invalid tensor
partitions, singular references, tiny positive eigenvalues, actual label
erasure and execution under `python -O`.

From the repository root with the pinned requirements and `PYTHONPATH=code`:

```text
python -m pytest -q code/quantum_information code/collar_alignment code/geometry code/maxent
```

The corrections concern these finite evidence implementations. They do not
amend A1–A3, prove a physical Markov collar, promote the gravity branch, or
complete the separate model-selection issues #1025/#1026. Frozen evidence and
archived simulator copies retain their historical bytes; this shared module
governs the live consumers named above.
