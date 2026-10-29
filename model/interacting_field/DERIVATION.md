# Interacting CAR field: source map, estimates and limitations

This is a CP-06 / [#1142](https://github.com/FloatingPragma/observer-patch-holography/issues/1142)
component, developed on draft [#1151](https://github.com/FloatingPragma/observer-patch-holography/pull/1151)
at `67aa7b6f4a904a1a90ea6d77bf2f66c741f2ba0c`. The action and apparatus below
are independently specified before executing their compiler. No native
repair residual is used as a definition of agreement with this field law.

The nearest existing constructions are the [Gaussian CAR source](../common_source/chiral_candidate/SOURCE.md)
and the [spatial Bose field](../spatial_interaction/DERIVATION.md). We reuse
complete local source controls, exact fermionic routing, symmetric splitting
and contractivity of retained instruments. The new contribution is their
mapping to a spatial **non-Gaussian fermion** action, an independently readable
connected four-point witness, complete continued records on growing cubes,
and a source-specific obstruction to identifying that action with the
record-derived product-reference modular generator. None of the general
CAR, product-formula or trace-distance facts is claimed as new.

## 1. Declared action and exact invariant space

Let `Lambda_R={-R,...,R}^3` with two supplied CAR modes per site. Mode order
is lexicographic site order followed by component 0,1. On each forward bond
`y=x+e_a`, put

```
h_(x alpha,y beta) = i c (sigma_a)_(alpha beta)/2,
h_(y beta,x alpha) = conjugate(h_(x alpha,y beta)),
H_R = sum_pq h_pq c_p^dagger c_q + U sum_x n_(x0)n_(x1),
c=1/8, U=1.
```

The local terms, preparation, number reads and feedback preserve total
number. Therefore `wedge^2 C^(2|Lambda_R|)` is an **exact invariant sector**
of every completed operation. No fermion mode, doubler or high momentum is
removed. The exact dimensions at R=1,2,3 are 1431, 31125, 234955. The code
never applies Gaussian/Wick evolution to the interacting state.

For a two-particle vector, use the antisymmetric matrix `A_pq=psi_(p,q)`
for p<q. Then `gamma=A A*`, `Tr gamma=2 ||psi||^2`. A local number-preserving
four-dimensional Fock matrix acts on local occupations 0,1,2. When a
spectator lies between two routed modes, the orientation of `c_p^dagger
c_r^dagger` supplies the CAR sign. `TwoFermions.apply_pair` applies that
sign to the single-occupation block and retains the vacuum and double
occupation scalars. Tests compare it to independently assembled full-Fock
matrices, also with a complex untouched reference.

The quadratic terms use the existing native hopping words. The contact
word is `diag(1,1,1,exp(-i U t))`, compiled through the existing source
diagonal compiler and contracted on all four Fock inputs. Nonadjacent
mode routing is the FSWAP conjugation identity, including its double
occupation minus sign. Its full source word is independently contracted;
all requested FSWAPs are counted. The large execution uses the exact CAR
identity instead of numerically contracting every repeated global route.
A wrong bosonic SWAP gives a different retained state in the control test.

The spatial lattice, spinor interpretation, CAR statistics, action, U,
initial occupations, code load/export, routing, blank pointers and readout
are supplied apparatus. The source compiler is reused; the endpoint-biased
repair law has not selected these services or this action.

## 2. A complete detector witness of non-Gaussian interaction

Use two sites a,b and prepare the normalized Slater state

```
psi_0 = ((c_a0^dagger+c_b0^dagger)/sqrt(2))
        ((c_a1^dagger+c_b1^dagger)/sqrt(2)) |0>.
```

Execute the on-site contact at both sites for t. The two double-occupation
coefficients acquire `exp(-iUt)` and the separated coefficients do not.
The four eigenvalues of the one-body density are `(1+cos(Ut))/2` twice and
`(1-cos(Ut))/2` twice. Hence

```
W = 2 - Tr(gamma^2) = sin^2(Ut),
<n_a0 n_b0> - (gamma_a0,a0 gamma_b0,b0 - |gamma_a0,b0|^2)
    = -sin^2(Ut)/4.
```

The actual same-component pair probability is zero. The nonzero Wick
prediction is a diagnostic of a Gaussian replacement, not a discarded
physical probability. At U=1,t=1/4,
`W >= (1/4-(1/4)^3/6)^2 = 9025/147456 > 0`, using the elementary sine
remainder bound. With U=0 the output remains Slater and the connected
four-point vanishes. The completed contact also violates the two-mode
Gaussian identity `det(U_one)=U_vacuum U_double` for this angle.

These are standard finite CAR/Wick criteria; see
[Surace and Tagliacozzo, Fermionic Gaussian states](https://arxiv.org/abs/2111.08343).
Here they are attached to actual compiled detector records. Four number
settings and twelve neutral phase/number settings reconstruct all entries
of gamma; a seventeenth prepared copy supplies the complete paired number
law. Every complement is retained. Incompatible tomography settings use
separate prepared copies. No assertion is made that one undisturbed specimen
provides all settings. The receipt counts all 17 preparations, contact words,
routing, 18 number pointers and 18 committed detector bits. Matrix-derived
gamma is used only as an independent reconstruction check.

This witness deliberately isolates the interaction by a supplied pulse
intervention. The following protocol keeps kinetic and contact terms together.

## 3. Continuing records and a positive interaction contrast

Let D occupy both components at the origin, and S occupy origin component
0 and neighbor `+e_x` component1. Prepare `(D+S)/sqrt(2)` with a source
hopping word. Evolve for 1/16, read the occupied component1 orbital
`(origin+i neighbor)/sqrt(2)`, retain both outcomes, apply phase pi/2 to
the neighbor component1 if the first bit is one, evolve again for 1/16,
and repeat the complete read. There is one field preparation and no reset
between reads. The state disturbed by the first read continues into the
second evolution. Four classical outcomes and their unnormalized quantum
output vectors remain in the execution. The JSON receipt retains their
complete probability law and cq-distance comparisons; the full large vectors
are regenerated by the replay rather than stored in the receipt. Number reads and rotations act on their full Fock
spaces, including configurations which have propagated outside the apparatus.

The U-deletion control changes U alone. The feedback-deletion control omits
the conditional phase alone. An independent sparse matrix exponential of
the entire H_R provides a reference, separate from the source split schedule.

A short additional protocol certifies that the interaction is observable
without relying on a floating contrast. For the same Slater preparation and
orbital read, compare H=K+UD to K at `t=1/(32B)`, where `B=6c+U=7/4`.
The interaction-minus-free probability derivative is U/2 at t=0.
For the adjoint generators, telescoping their nth powers gives

```
||ad_H^n(O)-ad_K^n(O)|| <= n (2B)^(n-1) 2U ||O||.
```

Consequently the higher-order remainder is bounded by
`2 U t (exp(2Bt)-1) <= 2Ut/15`, since `2Bt=1/16` and
`exp(x)<=1/(1-x)` for 0<=x<1. The exact Hamiltonian contrast is at least
`11Ut/30 = 11/1680`. Subtract the two ideal product-formula errors below.
At 16 steps the executed ideal contrast is at least
`141294691/21579694080 > 0.00654` in every box and on the infinite fixed
lattice. Both read outcomes remain; no outcome is selected or renormalized.

## 4. Volume-uniform complete-history bounds

Six spatial edge matchings partition the hopping terms. Each matching
also pairs spinor modes disjointly: each Pauli matrix has one nonzero entry
in each row and column. A matching's one-particle eigenvalues have absolute
value at most c/2; its two-particle norm is at most c. Since contact
double-occupation projectors are mutually orthogonal on N=2, the full
contact term has norm at most U. Thus the sum of the seven split-generator
norm bounds is `S=B=6c+U`, independently of volume. The same bound holds
on `wedge^2 l2(Z^3 x {0,1})`; that bounded sector Hamiltonian is self-adjoint.

The palindrome uses half steps of the six matching generators, the contact
generator twice at the center, and the reverse half steps. Its Taylor
coefficients through degree two agree with exp(-itH). The integral Taylor
remainders for each unitary-valued path are at most `S^3 |dt|^3/6`; adding
the exact-exponential and product remainders and telescoping n steps gives

```
||U_split(t)-exp(-itH)|| <= S^3 |t|^3/(3 n^2).
```

Purify an arbitrary input and dilate the fixed instruments, retaining all
classical records and the untouched reference. Intervening isometries are
contractions. Summing the evolution errors therefore bounds the half trace
distance of the **complete quantum-record output**, with no multiplier for
the number of outcomes. Adapted operations must be the same functions of
retained records in both compared protocols. For the two 1/16 segments at
n=16 the exact bound is `343/100663296`.

For the finite/infinite comparison, the initial occupied sites and every
intervention are confined to `{origin, origin+e_x}`. Instruments may change
occupations there but do not transport particles elsewhere. In the joint
Stinespring Taylor expansion, fewer than R hopping factors cannot produce
a difference with a larger cube: they have not reached a missing boundary
bond. Contact terms do not move particles. The expansion coefficients of
total Hamiltonian degree below R agree. Summing the two integral remainders
across fixed intervening isometries gives

```
||V_R - V_infinite||_(initial apparatus-supported subspace)
    <= 2 (B T)^R/R!,   T = sum_j |t_j|.
```

This is a bound on that initial subspace, **not all states on the finite
cube**. It covers its untouched references. It applies to ideal continuous
Hamiltonian evolutions; the temporal splitting bound is added separately.
At R=3,T=1/8 the spatial bound is `343/98304`; the combined ideal-source
bound is `351575/100663296`, approximately 0.003493. Increasing radius and
split count controls the observation law on the infinite lattice at fixed
spacing. This is not a spatial continuum or ultraviolet limit.

Native pulse matrices are floating contractions. Their full local target
residuals include number-sector compression residuals. Tiny forbidden
off-number and off-diagonal floating entries are compressed away in the
global sector replay; the literal floating local matrix is therefore not
asserted to preserve that sector exactly. Routing, preparation,
readout, feedback and every repeated evolution call are counted in numerical
residual sums; the diagnostic telescoping estimate uses `expm1(sum residuals)`
to allow norm gain. Those residuals are **not outward certified**. They are
not added to the rational theorem as if they were rigorous error bars.
Physical service error, angle synthesis, storage and initialization accuracy
remain separate. No state or output branch is rescaled to repair roundoff.

## 5. What the kinetic symbol establishes

Fourier transformation of the infinite one-body hopping gives
`h(k)=-c sum_a sigma_a sin(k_a)`. Its eight zeros have coordinates 0 or pi
and Jacobian chirality `-(-1)^(number of pi coordinates)`: four of each sign.
All are retained in the action and in the execution. This direct count is
consistent with the standard lattice chirality obstruction of
[Nielsen and Ninomiya](https://doi.org/10.1016/0370-2693(81)91026-1);
we do not import a general no-go without its hypotheses. The interaction
vanishes on the exact N=1 sector, so that sector's dispersion is unchanged.
This does not identify the excitations of an interacting finite-density vacuum.

For componentwise |k_a|<=kappa, the elementary sine remainder gives
`||h(k)+c sigma.k|| <= c kappa^3/2`. On the free N=2 momentum patch the
unitary error at model duration t is at most `c t kappa^3`. Local contact
scattering can leave that patch. The free bound therefore does not certify
an interacting low-energy truncation and does not justify throwing away
other nodes. The local interaction supplies a nontrivial four-point response;
control of a relativistic interacting EFT/continuum, a chiral regulator,
vacuum, current matter lift and matched renormalization remain required.

The empty Fock state is not asserted to be a ground state: the Weyl kinetic
spectrum has negative as well as positive energies. The execution concerns
a prepared finite-number state and a bounded invariant-sector generator.
It does not derive a positive physical vacuum, a relativistic rest mass,
confinement or a photon gap. Number conservation is the supplied CAR U(1),
not a derived electromagnetic current or a stress/Ward map to gravity.

## 6. Exact compatibility boundary with the recovered native clock

The stationary target reference of #1151 on its H66 occupation carrier is
`pi_x = product_i theta_i^(x_i)/Z`, sum x_i=2. Under the canonical
occupation-preserving CAR dictionary,

```
K_source = -log(pi) = -sum_i log(theta_i) n_i + log(Z) I.
```

Thus its modular unitary is quadratic and preserves every Slater input.
This holds for **every positive table**, including the weights actually
inferred from the retained record counts. The inferred table need not have
the exact target's three-level degeneracy. The receipt executes both sets
of one-mode phase controls and retains their zero-Slater-defect checks.
It does not rerun 31 million clock acquisition events or relabel the field
Hamiltonian coordinate as the recovered modular parameter.

The contact witness W>0 consequently rules out equality with this source
modular flow under the occupation/CAR-linear dictionary. Also the complete
H66 matrix algebra cannot be unitarily identified with the much larger
field N=2 algebra merely by renaming its basis. Neither statement excludes
an embedding with new degrees of freedom, a different joint reference,
nonlinear CAR identification, changed law, or a protocol with repair and
extra coherent interventions. Those require an actual compatibility map.
A clock reconstruction advance does not by itself select the interaction,
physical duration or energy normalization. This is a scoped CP-06/CP-12
obstruction with an explicit surviving path, not a whole-program no-go.
