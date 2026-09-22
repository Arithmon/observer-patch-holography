# Conditional occupation selection

This package checks finite witnesses for the analytic theorem in
`paper/tex_fragments/PAULI_STABILITY_SELECTION.tex`. Its starting relation is
`a a† - q a† a = I`, with one unknown real scalar `q`. It does not assume CAR.
The normalized annihilation vacuum, positive Hilbert metric, invariant operator
domain, and exact phase generator `N=a†a`, `[N,a†]=a†` are supplied premises.
The short-word identity is `q²-1=0`. Thus the declared scalar family has two
candidates. A fixed signed generator with an accessible negative mode and a
global lower energy bound excludes `q=1`; `q=-1` then implies square-zero
creation. Linearity and the common mode relation give CAR by polarization.
Alternatively, finite dimension of the entire nonzero matter Hilbert space
excludes `q=1` by the trace of a commutator. These are separate sufficient
conditions, both within the same supplied scalar algebra and number law.

From the research repository root:

```sh
python3 code/pauli_stability/build.py --write
python3 code/pauli_stability/verify.py
python3 -m pytest -q code/pauli_stability/test_pauli_stability.py
```

The producer normal-orders finite words in the unknown `q`; the verifier uses
an independent ladder recurrence. The edge matrix comes from the current
package's `spin_action`; independent Pauli-matrix assembly checks four exact
orthonormal eigenmodes with energies `-1,+1,-1,+1`. It is a supplied single
hopping factor, not the complete source Hamiltonian. The five parent multiplets
have 15 internal channels and total integer hypercharge zero. Placing `n`
bosons in one negative mode of every channel gives charge zero and energy
`-15n`. The arbitrary-`n` proof is analytic; retained samples are finite
witnesses. No operator Gauss constraint or gauge energy is included.

The receipt also checks a two-mode CAR sea and a neutral fully filled
negative-mode edge sea. Its positive three-dimensional truncated-boson control
has nonzero double creation and a bounded negative energy, and obeys the exact
number-generator law. It violates every constant-`q` relation at the top of
the ladder. Stability or finite dimension alone therefore does not select CAR.
Fixed particle number and a uniform number-energy shift are separate escapes;
the latter preserves fixed-number continuous one-body readouts, but does not
automatically preserve the original same-step Cayley update.

This selector is spin blind. It is not a physical spin--statistics theorem or
a derivation from the observer axioms. Finite response capacity does not
identify the full matter Hilbert space. An OPH attachment must derive
the algebra and number law on bounded self-reading patches with local state,
ports, readback, records and feedback or repair, and publish evidence binding
that realization to the energy or finite-capacity premise. The present public
evidence bundle establishes the conditional mathematical step only.
