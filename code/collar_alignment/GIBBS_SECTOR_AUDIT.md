# Sector Gibbs states must preserve energy differences

This is a focused repair under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033),
building on the sector-label audit (#1029) and the MaxEnt coordinate audits
(#1036/#1037). The remaining collar constructor still diagonalized scalar
energy origins together with interactions. A changed origin could therefore
turn a nonaligned Gibbs state into an apparently aligned one. Its absolute
Hermiticity tolerance also depended on energy units.

## Reproductions before the repair

The first commit adds eight failing tests against main `7f5805a2`; the branch
then incorporates the null-net integration at `07067bf1`. Write
`X=[[0,1],[1,0]]`, `J=X tensor X` on two qubits, and use ordinary matrix trace.

| Input | Previous result | Independent answer |
| --- | --- | --- |
| `H=J`, central energy `1e20`, `beta=.4` | Approximately `I/4`; alignment accepted | `(I-tanh(.4) J)/4`; off-diagonal magnitude `0.09498724056380622` |
| `H=1e20 I+J`, central energy zero, `beta=.4` | `I/4`; alignment accepted | The same interacting state |
| Sectors `diag(0,1)` and `(0)`, both central energies `1e20` | Sector probabilities `(.66666667,.33333333)` | `(.57768120,.42231880)` |
| Sectors `diag(0,2)` and `1e20 I+X`, central energies `(0,-1e20)` | Probability ratio `.56766764` | `exp(-1)=.36787944` |
| `H=1e-200 X`, central energy `1e200`, `beta=1e200` | Energy multiplication overflows | `(I-tanh(1) X)/2`, a well-resolved normalized state |
| `H=[[0,1e-200],[0,0]]`, `beta=1e200` | Accepted as `I/2` | Reject: non-Hermitian |
| Mixed Boolean or numeric-string matrix entries | Silently coerced to numbers | Reject malformed observable input |

For the first two rows the true mutual information is
`beta tanh(beta)-log cosh(beta)=.0740260995142574` nats. The old result was
zero to roundoff. The failure was in the supplied state, so four mutually
consistent diagnostics of that wrong state could not validate the physics.
Follow-up controls also reject integer or extended-precision conversion
that erases an energy gap, such as the literal integers `2**53,2**53+1`.

## One Gibbs problem, with the center retained

Let the declared finite algebra be the direct sum of matrix algebras on
spaces of dimensions `d_a`. Supply Hermitian `H_a`, real central energies
`e_a` and a finite real inverse temperature `beta`. Define

```text
Z_a   = Tr exp(-beta H_a)
tau_a = exp(-beta H_a)/Z_a
Z     = sum_a exp(-beta e_a) Z_a
pi_a  = exp(-beta e_a) Z_a / Z
tau   = direct_sum_a pi_a tau_a.
```

This follows by exponentiating each summand of
`H=direct_sum_a(H_a+e_a I_a)`. Discarding the local partition functions,
or shifting each block without restoring its normalization factor,
changes the center probabilities. At `beta=0`, `pi_a=d_a/sum_b d_b`.
Finite negative beta is also mathematically allowed; no unbounded-spectrum
or zero-temperature claim is made here.

For any labelled state `rho=direct_sum_a p_a rho_a`, spectral functional
calculus gives `log(p_a rho_a)=log(p_a) I_a+log(rho_a)` on its support.
Consequently, with the convention that a zero-weight term contributes zero,

```text
D(rho || tau) = sum_a p_a log(p_a/pi_a) + sum_a p_a D(rho_a || tau_a)
             = [beta Tr(rho H)-S(rho)] - [beta Tr(tau H)-S(tau)].
```

Both terms in the first line are nonnegative, and equality requires the
same center distribution and the same conditional states in every sector.
Thus the constructed state is the unique minimizer of the stated finite
dimensionless free-energy functional. This is the usual Gibbs variational
principle, not a new OPH physical law; see, for example,
[Chowdhury, Low and Wiebe](https://arxiv.org/abs/2002.00055).
The displayed proof also distinguishes two independent errors: wrong
sector probabilities and wrong conditional states. Agreement between
several alignment diagnostics says nothing about the first error.

If the center probabilities are prescribed instead, minimization only
selects the conditional `tau_a`; it does not replace those probabilities
by `pi_a`. In particular, A3's declared cover weights are not additional
variables silently optimized by this calculation.

## An exact interaction test for alignment

Within a declared collar sector write the cut as
`L=(A,bL)` and `R=(bR,D)`. The Hilbert--Schmidt projection onto one-sided
operators is

```text
Pi(T) = (Tr_R T / d_R) tensor I_R
      + I_L tensor (Tr_L T / d_L)
      - (Tr T / (d_L d_R)) I.
```

Because `Pi(I)=I`, the Gibbs logarithm obeys the exact identity

```text
(1-Pi) log(tau_a) = -beta (1-Pi) H_a.
```

Including the central factor `log(pi_a) I` does not change that identity.
For `beta != 0`, the following are therefore equivalent within this finite
Gibbs family: the conditional state is a product across the declared cut;
its logarithm is a sum of one-sided operators; and the Hamiltonian is a sum
of one-sided operators. To prove the first equivalence, take the logarithm
of a faithful product state in one direction, and exponentiate the two
commuting tensor terms in the other. The second equivalence is the displayed
identity. It also gives an exact norm relation:

```text
||(1-Pi) log(tau_a)||_op = |beta| ||(1-Pi) H_a||_op.
```

A finite central-energy penalty changes the sector's probability but never
eliminates its positive weight or its conditional interaction. Exact
blockwise alignment cannot be obtained by making a bad sector rare. At
`beta=0`, every conditional state is maximally mixed and the converse about
the Hamiltonian is false. The tests retain that exception explicitly.
The identity does not assert that zero collar CMI implies alignment; the
existing Bell counterexample still disproves that implication.

## Shared implementation and numerical scope

`quantum_information/gibbs.py` now owns the observable decomposition and
thermal eigensolve previously private to MaxEnt. The collar constructor
validates its four subsystem dimensions and delegates to `gibbs_sectors`.
MaxEnt keeps its public functions and optimizer, using the same numerical
primitives. These two consumers are not independent verification oracles.

1. Subtract a diagonal scalar anchor before scaling and diagonalizing the
   bounded traceless variation. Check Hermiticity relative to the variation,
   so neither large identity terms nor small energy units hide a bad input.
2. Retain the anchor, remaining mean and central energy as separate scalars.
   Form each sector log weight from its conditional log partition function
   and those offsets. Accumulate the supplied binary64 scalar offsets and
   their beta products as exact rational numbers; subtract one common log
   weight before conversion back to floating point. This handles common
   offsets whose products overflow binary64 without erasing relative mass.
3. Normalize all sectors together. Check conditional support, center mass,
   and every product of center mass with conditional spectral mass. Local
   checks alone would miss a positive direction lost at global normalization.

Input conversion must preserve every numeric component exactly in binary64;
mixed Boolean/string inputs and lossy integer or extended-precision inputs
raise. This strengthens the shared MaxEnt validator as well. No precision
rule can recover a diagonal splitting the caller already rounded away while
forming `1e20 I+diag(1,-1)`; the routine evaluates the actual supplied array.
Underflow and unresolved dense-state support raise instead of clipping
eigenvalues, deleting sectors or inserting a probability floor. Operations
whose centered differences or coefficients exceed the working range may
also be rejected. Exponentials, logarithms, eigensolves and reconstructed
states still have floating-point error. The exact scalar accumulator is
not an interval certificate for the whole computation.

## Controls and downstream impact

The controls use analytic Pauli and diagonal spectra, SciPy's separate
scaling-and-squaring exponential on the complete direct-sum matrix, and
separate matrix logarithms for the relative-entropy/free-energy identity.
They cover unequal dimensions, signed temperature, energy-unit changes,
central-energy changes, compensated offsets, sector regrouping and
permutation, local basis changes, complete alignment criteria, and the
beta-zero exception. Adversarial cases include non-Hermitian or malformed
inputs, lossy conversions, positive mass lost only at joint normalization,
and validation with Python assertions disabled. No external simulation,
new empirical fit or generated bulk dataset is needed.

| Consumer | Impact |
| --- | --- |
| Collar Gibbs/alignment evidence (#543) | Corrects false alignment and incorrect center probabilities for the reproduced inputs; one-sided and cross-cut controls retain their intended distinction. |
| MaxEnt closure (#539, #1036/#1037) | Removes duplicated thermal machinery, strengthens conversion checks, retains the optimizer, covariance and closure-defect contracts. The live acceptance computation is replayed without rewriting pinned evidence. |
| Gravity premise ladder, rung 3 | The Bell countermodel and independence-limited state-alignment status remain valid; the correction supplies no new axiom-side alignment derivation. |
| Geometry, null-net and Einstein-closure consumers of the shared package | Their affected regression suites are replayed. This PR introduces no change to their input models or physical normalization. |
| Papers, flagship, book, public claims | No new physical conclusion or target changes. The active spacetime paper's section on the energy-normalization hinge still requires a declared beta and physical energy identification. Historical `prop:msachar`/`cor:msareduction` labels from #543 no longer exist in the active paper; code comments now identify their historical origin. |
| Receipts, registries, custody and frozen registrations | No bytes changed. This branch preserves the merged null-net receipt and all frozen targets. |

Run the focused controls with `PYTHONPATH=code`:

```text
python -m pytest -q code/quantum_information/test_gibbs.py code/collar_alignment/test_gibbs_energy_origins.py -W error
```

The existing two-platform finite-quantum-information CI workflow discovers
both files and also runs `code/quantum_information`, `code/collar_alignment`,
`code/geometry` and `code/maxent` together.

Local validation after incorporating main `07067bf1`: all 773 affected
tests pass on Linux; Windows passes 767 with six skips for unavailable
extended precision. Both runs treat warnings as errors. Claim registry,
axiom consistency, reader style, the generated gravity ladder, clean-export
public surfaces and the Linux Lean proof inventory also pass.

The live MaxEnt acceptance passes and exactly matches a separate replay of
`origin/main` on the same Windows environment. Both differ from the stored
receipt in six floating-point fields: the largest entropy difference is
`8.59e-16` nats and the largest Pinsker-bound difference is `1.47e-8` near
zero (the bound takes a square root). Those pre-existing platform differences
are not attributed to this refactor; all nonnumeric fields agree. The
stored receipt remains unchanged.
