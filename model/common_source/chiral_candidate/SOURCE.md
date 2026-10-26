# Finite Gaussian CAR execution through the native source

The quadratic candidate can be executed without constructing its exponentially
large Fock vector. This uses an exact representation theorem for finite
quadratic fermion dynamics. It does not substitute a fitted spectral matrix
for execution of the declared source controls.

The conventions are

\[
 H=\sum_{ij}h_{ij}c_i^\dagger c_j+
 \frac12\sum_{ij}(\Delta_{ij}c_i^\dagger c_j^\dagger+
 \overline{\Delta}_{ij}c_jc_i),\qquad
 h=h^\dagger,\quad\Delta=-\Delta^T,
\]

\[
 w_{2j}=c_j+c_j^\dagger,\quad
 w_{2j+1}=-i(c_j-c_j^\dagger),\quad
 \Gamma_{ab}=\frac i2\langle[w_a,w_b]\rangle.
\]

Thus an empty mode has the covariance block `[[0,-1],[1,0]]`, and
`n_j=(1+Gamma[2*j,2*j+1])/2`. The parameter in every exponential is the
programmed model/modular parameter. Neither host runtime nor pulse count is
identified with a calibrated physical clock.

## What actually executes

`local_gate(kind, coefficient, duration)` constructs a chronological native
word and executes **all four** local Fock input columns through the existing
six-level source compiler. It retains the resulting complete 4×4 matrix.
The kinds `hopping`, `pairing` and `onsite` mean respectively

\[
 zc_0^\dagger c_1+\bar zc_1^\dagger c_0,\qquad
 dc_0^\dagger c_1^\dagger+\bar d c_1c_0,\qquad
 e n_0.
\]

For real coefficients, the hopping word uses the commuting `XX+YY`
combination, and pairing uses `XX-YY`. For a complex coefficient, an actual
onsite phase conjugation supplies its argument. Scalar phases and the
zero-angle requests inside those words are retained. Onsite evolution is
also contracted on all four columns, with the second mode an identity
spectator in this local identity check.

For each executed matrix U, the adapter extracts
`R_ab=Re Tr(w_b U† w_a U)/4` and checks the entire operator identity
`U† w_a U=sum_b R_ab w_b`. It also checks unitarity, reality and orthogonality.
The completed source words are even Gaussian unitaries in exact arithmetic.
Their CAR action therefore determines their adjoint action on every CAR
polynomial. A Gaussian density operator is fully determined by its covariance
and Wick's rule, so this representation retains its full state, rather than
only its number expectations. Overall ket phase is absent from a covariance;
coherently controlled unitary interferometry must additionally retain the
local source words/phases. It cannot be inferred from covariance alone.

The implementation applies no polar, norm, covariance or probability repair.
The displayed arrays and residuals are floating diagnostics of these exact
identities, not outward arithmetic or physical-service error certificates.

## Why a nonadjacent pair is still fermionic

An ordinary qubit SWAP gives the wrong answer on two occupied fermion modes.
The source `fswap` word implements

```text
|00> -> |00>, |10> -> |01>, |01> -> |10>, |11> -> -|11>.
```

It is a real hopping rotation by π/2 followed by +π/2 number phases on both
modes. All four Fock columns, including the minus sign, are checked. This
Gaussian unitary swaps the two CAR generators and their Majorana pairs.

For a pair i<j, the adapter moves mode j to i+1 with adjacent FSWAPs, applies
the two-mode word, and reverses the route. This is the finite-CAR proof of
the required Jordan–Wigner signs. It is not a bosonic tensor embedding of a
two-qubit matrix onto nonadjacent sites. The accumulated transformation
composes the **actual source-derived** FSWAP maps, including their numerical
residuals; it does not replace them with perfect permutation matrices.

Every routing and unrouting FSWAP is charged in the source resource ledger.
Each block updates only two or four rows of the accumulated transformation.
The code never allocates the full many-mode Fock matrix. For 90 work modes
it stores a 180×180 Majorana transformation, then applies it once to the
covariance. Additional Gaussian reference modes and their cross-correlations
are retained, with identity evolution on the reference. Cached local source
contractions do not reduce the cost of the physical word being represented.

`evolve(h, delta, duration, steps=..., covariance=...)` uses a symmetric
forward/reverse product of the declared onsite, hopping and pairing terms.
The split count and term order are supplied architecture data. Comparison
with continuous quadratic dynamics needs its own splitting and numerical
error analysis; this adapter does not assert equality to that continuous
flow for noncommuting terms.

## Complete occupation records

`occupation_read(Gamma, mode, weight=...)` executes the actual native
occupation-pointer word on a field mode and a fresh **neutral record qubit**.
The pointer is not another fermion mode. Its two Kraus matrices are checked
on the complete two-dimensional field algebra: `P0=|0><0|` and
`P1=|1><1|`. Number is even, so this read needs no fermionic swap string.

Both outcomes are retained. A Gaussian branch carries a normalized covariance
descriptor and a separate probability weight; conditioning that descriptor
does not renormalize the retained classical branch law. The probabilities
are `p_n=(1+(2*n-1)*Gamma_ab)/2`. For unmeasured indices i,j the finite Wick
identity gives

\[
 \Gamma'_{ij}=\Gamma_{ij}+
 \frac{s}{1+s\Gamma_{ab}}
 (\Gamma_{ia}\Gamma_{bj}-\Gamma_{ib}\Gamma_{aj}),\qquad s=2n-1.
\]

The measured mode has `Gamma'_ab=s` and zero cross-correlations with the
remaining modes. All unmeasured and reference entries are updated by this
identity. An exactly zero-probability outcome remains present, with an
undefined conditional covariance represented by `None`. Probabilities are
not clipped. Ill-conditioned numerical conditionals are reported through
their CAR-bound residual, rather than silently repaired.

Each conditional branch remains Gaussian. A classical mixture of different
Gaussian branches need not be Gaussian: later evolution must retain the
branch tree, not replace it by its averaged covariance and assume Wick's
rule still holds. No work mode is re-prepared between a read and continuation.

## Scope

This realizes the declared finite quadratic CAR continuation through the
native local control identities. The coefficient map, Higgs background,
processor load/export and initialization services, and any clock apparatus
remain explicitly supplied. Pairing is parity even but does not conserve
the unbroken occupation number. Electromagnetic gauge invariance and the
appropriate charge assignment require the candidate action's separate
checks. No empirical mass, strong-interaction backend, interacting gauge
field, laboratory mass unit or Standard Model mass agreement is introduced
by this execution adapter.
