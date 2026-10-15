# Native quantum transfer and executed coherent measurements

The same native C/E rules now determine a positive transfer operator on
the **entire 4,356-dimensional binary operator space**, extending the
[66-state classical reconstruction](../primitive_observer/DERIVATION_TRANSFER.md).
Its unique logarithm predicts the executed coherent readout, including
complete quantum outputs and untouched reference correlations. No oscillator
Hamiltonian or fitted measurement curve enters this calculation.

At four native attempts the declared phase experiment predicts

| Candidate | Plus probability | Minus probability |
| --- | ---: | ---: |
| C | `1921/2592 = 0.7411265432` | `671/2592 = 0.2588734568` |
| E | `12401/20000 = 0.62005` | `7599/20000 = 0.37995` |

The difference is `12259/101250 = 0.1210765432`. The readout attains the
full output-state trace distance for this preparation, and four is the
best integer attempt count for this experiment. This makes the prior
[analytic phase-memory witness](../observer/QUANTUM_MIXING.md) executable.
These are conditional native-model predictions, with no laboratory dataset
or empirical agreement claimed.

## What is supplied, derived and executed

The supplied instance is the existing twelve-port carrier, its binary
total-two subspace, independent uniform choices among thirty seams, fair
placement coins, and one complete repair candidate, C or E. The coherent
preparation and final two-outcome apparatus are also declared. Selecting
these services from a broader interface remains an obligation in
[source selection](../SOURCE_SELECTION.md).

From those inputs the construction derives each quantum channel, its
operator-space logarithm, the exact phase decay and both measurement
probabilities. There is no independently supplied transfer Hamiltonian or
decay constant. C and E remain distinguishable complete candidate laws;
their common classical action does not choose between them.

For each candidate and each count `n=0,1,2,4,8`, the computation:

1. Constructs the two-particle binary basis from the native load domain.
2. Calls the actual native rule on every configuration, seam and coin;
   constructs Kraus maps by equality of the candidate's committed record.
3. Prepares the stated uniform mixture of ten coherent phase pairs.
4. Runs `n` native channels through the shared
   [primitive observer runtime](../primitive_observer/runtime.py), summing
   every forgotten seam, coin and repair label as a CP operation.
5. Executes the complete plus/minus measurement and retains both
   unnormalized 66-dimensional terminal density matrices.
6. Compares both probabilities and every terminal matrix entry with the
   independently constructed exact block transfer. It also compares
   `exp(-n H_L)` with the unread runtime density.

Both preparation phases are executed. A second experiment uses an
entangled reference qubit and compares all 132-dimensional terminal
matrices. A serialized checkpoint after two native attempts is independently
replayed, then continued through attempts three and four and the readout.
No normalized branch or probability is discarded or postselected.

This is native **channel semantics** execution. It does not emit compiled
source pulse words. Reported attempt and readout counts are logical
resources. Numerical residuals check this finite implementation; they are
not outward error bounds or physical service certificates.

## Full operator-space source attachment

Let `X` be the 66 binary configurations and let `Phi_L` be the channel with
candidate `L`, after the seam, coin and native records are forgotten.
Use the Hilbert–Schmidt inner product `Tr(A* B)` on `M66`.

For a matrix unit `|x><y|`, native record equality preserves the signed
difference `z=x-y`. The existing signed-difference theorem identifies each
block exactly as

```text
Phi_L[z] = (1-kappa_L(z)) I - L_exclusion[z]/60,
kappa_C(z) = sum_seams |z_a-z_b|/60,
kappa_E(z) = sum_seams 1[z_a != z_b]/30.
```

The common particles move on the zero-difference sites; positive and
negative difference sites remain fixed. The implementation checks equality
between the actual native record enumeration and this independent
exclusion construction for **every** binary matrix unit. The block sizes
are one block of size 66, 132 blocks of size ten, and 2,970 scalar blocks.
This avoids constructing or diagonalizing one dense `4356 x 4356` matrix.

All blocks are real symmetric, entrywise nonnegative, with row sums at
most one. Exact rational diagonal-dominance bounds give

```text
(2/3) I_HS <= Phi_C <= I_HS,
(11/30) I_HS <= Phi_E <= I_HS.
```

These are **operator-space** inequalities, distinct from positivity or
complete positivity of maps on density operators. A CP map need not be
positive as a Hilbert–Schmidt operator: conjugation by Pauli `Z` sends
Pauli `X` to `-X`, giving a negative Hilbert–Schmidt quadratic form.

Strict operator-space positivity supplies the unique self-adjoint log

```text
H_L = -log(Phi_L) >= 0,      Phi_L^n = exp(-n H_L).
```

The zero-difference block is the connected native exclusion kernel.
Every other block has strictly subunit row sum. Thus the ground direction
is uniquely `I66/sqrt(66)`. The existing exclusion-gap attachment gives
the common finite gap

```text
-log((55+sqrt(5))/60) = 0.04716030350797...
```

Its unit is inverse **Euclidean native attempt spacing**, with that
spacing set to one. The ground operator is a representation of the
stationary tracial state, not a newly prepared coherent native ket.

The complete matrix-unit identity extends by complex linearity to every
operator, not just density matrices or chosen observables. Tensoring with
an arbitrary untouched reference preserves the identity. In particular,
for all complex observables `A,B`,

```text
Tr(A* Phi_L^n(B))/66 = <A, exp(-n H_L) B>_HS/66.
```

For any finite family of operators `A_i` and nonnegative integer counts
`n_i`, the bond form is nonnegative:

```text
sum_ij Tr(A_i* Phi_L^(n_i+n_j+1)(A_j))/66
    = <V, Phi_L V>_HS/66 >= 0,
V = sum_i Phi_L^n_i(A_i).
```

This is the positive-transfer algebra used by
[`SourceLaw.ReflectionPositivity`](../../Lean/SourceLaw/ReflectionPositivity.lean),
now attached to the full native operator space. It is a two-point/bond
Gram statement. It does not establish all noncommutative history
reflection axioms or replace the actual CP instrument composition for
intervening measurements. The new attachment and its scope are analytic
and executable; no new Lean theorem is claimed.

## Exact coherent prediction and discrimination

Choose the first native seam `(u,v)`. With `|uw>` denoting the configuration
occupied at ports `u,w`, define

```text
A = sum_(w outside {u,v}) |uw><vw|,       S = A+A*,
rho_+/- = (1/20) sum_w (|uw> +/- |vw>)(<uw| +/- <vw|),
F_+/- = (I +/- S)/2.
```

The twenty involved configurations are distinct. `S` has eigenvalues
`+1,-1,0`, so the two effects are positive and sum to the identity.
The apparatus uses the positive square roots of these effects as its
postmeasurement maps. That choice is declared; specifying effects alone
does not select an instrument.

The exact native blocks give `Phi_L(A)=r_L A`, with
`r_C=5/6` and `r_E=7/10`. Therefore

```text
p_L(plus | rho_+, n) = (1+r_L^n)/2,
p_L(plus | rho_-, n) = (1-r_L^n)/2,
dist(Phi_L^n(rho_+), Phi_L^n(rho_-)) = r_L^n.
```

For the same plus preparation, the candidates' diagonal populations are
identical and their difference is `(r_C^n-r_E^n) S/20`. Consequently the
trace distance is `Delta_n=(r_C^n-r_E^n)/2`, exactly attained by the binary
readout. Its equal-prior, single-trial discrimination success is
`(1+Delta_n)/2`, about `0.56053827` at four attempts. This is not a claim
about the best possible preparation, ancilla, adaptive test or many-trial
channel-discrimination protocol.

The sign of `Delta_(n+1)-Delta_n` changes when `(25/21)^n` crosses `9/5`.
The exact inequality `(25/21)^3 < 9/5 < (25/21)^4` proves the unique
integer maximum at four. Waiting longer is not always a stronger test:
by eight attempts the separation has fallen to about `0.08746001`.

Dephasing the preparation removes this discriminator: the two channels
then have the same complete outputs and both pointer probabilities are
one half. The separate reference test deliberately has the same pointer
marginal under both laws, so checking its surviving reference correlations
is necessary; a matching marginal alone would be insufficient.

## An exact boundary on the logarithm

Neither `exp(-s H_C)` nor `exp(-s H_E)` is a positive quantum map for every
nonnegative real `s`. Both preserve the diagonal subalgebra with the same
classical kernel `K`. Their functional-calculus square roots therefore
restrict to `sqrt(K)` on that subalgebra.

The construction rechecks the
[exact diagonal half-step obstruction](../primitive_observer/DERIVATION_TRANSFER.md#an-exact-fractional-step-obstruction)
using rational arithmetic. Applying the putative half-step to the positive
configuration state occupied at ports `(2,4)` gives a population at `(3,11)`
bounded above by

```text
-4008679803811/42998169600000000 < 0.
```

A negative population already disproves positivity, without constructing
a large Choi matrix. Thus this logarithm does not generate a continuous
CP semigroup. Integer powers remain the original valid native channels.
Nor is `H_L` a Hamiltonian acting on the 66-dimensional carrier kets:
it acts on their 4,356-dimensional operator space. A physical unitary,
clock interpretation, local continuum field family and experimental
realization require further attachments.

## Reproduce

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_quantum_history --check model/native_quantum_history/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/native_quantum_history
```

The [receipt](results.json) binds the source files and records exact
fractions, portable numerical diagnostics, complete outcome laws and
checkpoint/reference checks. The tests also compare the actual Kraus
execution with block transfer on a dense complex operator tensored with
a reference, detect corrupted native records, and retain the dephased
preparation and fractional-step controls.
