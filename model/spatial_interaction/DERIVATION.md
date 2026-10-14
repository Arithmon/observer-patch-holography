# A spatial interacting observer family at fixed lattice spacing

This branch supplies a repulsive Bose--Hubbard action on growing cubic graphs.
It adds spatially extended interacting observations and a controlled
infinite-volume comparison. It does not establish a relativistic continuum
QFT, remove a spatial regulator, derive a Hamiltonian or geometry from native
repair, extend the existing charged gauge/scalar common action, or calibrate
physical duration. Every time below is the supplied Hamiltonian **model
coordinate**, with lattice spacing and hbar set to one by convention.

The action belongs to the established Bose--Hubbard family; see
[Jaksch et al., *Cold bosonic atoms in optical lattices*](https://arxiv.org/abs/cond-mat/9805329).
The additional result here is its realization as one spatial family of
recorded experiments, with an exact invariant execution sector, explicit
whole-history error budgets and an interaction witness uniform in volume.
The nearest existing project results are the
[two-site source splitting bound](../common_source/scalar_coupling.py), the
[continued instrument comparison](../common_source/COMPOSITION.md), and the
[free-field spatial observer comparison](../observer/scalar_response/DERIVATION.md).
We transfer their unitary-telescoping and retained-instrument reasoning to
this specified interacting family; no new general convergence principle or
new lattice action is claimed. Compatibility here is fixed independently of
the execution: the same preparations, local occupation/phase/Ramsey
instruments and record labels define each comparison.

## Action, sector and observed quantities

Let Lambda_R={-R,...,R}^3, R>=1, with induced nearest-neighbor edges and the
open-graph combinatorial Laplacian L_R. The full Fock Hamiltonian is

    H_R = m N + J dGamma(L_R) + (U/2) sum_x n_x(n_x-1),
    m>0, J>=0, U>=0.

The graph Laplacian includes its actual degree diagonal. This is an open
(discrete Neumann) graph boundary; it is not the constant diagonal of a
Dirichlet compression. Its edge form is

    dGamma(L_R) = sum_{edges xy} (a_x^dagger-a_y^dagger)(a_x-a_y).

All terms preserve total particle number. Execution in N=2 is exact for
the stated preparations and number-preserving instruments; it introduces
no occupation-cutoff error. The Hilbert dimension is V(V+1)/2 for V sites.
This compressed sector is not a tensor product of unrestricted local
matrix algebras. We use only the stated operations that act faithfully in
the invariant sector.

The apparatus occupies A={o,e_x}. The preparation may coherently combine
the normalized Fock vectors D=|2_o> and S=|1_o,1_e_x>. The history retains
a local Ramsey outcome, applies an outcome-controlled phase of pi/2 to
the split configuration when that committed bit is one, continues
Hamiltonian evolution, and reads the same local Ramsey instrument again.
The two Hamiltonian segments each have model-coordinate duration 1/16. Every
outcome remains in the joint law; there is no branch postselection. A
later fixed native readout or repair chosen from the committed records is
a common channel applied to that law and cannot enlarge its error.

Particle number, hopping strength J, contact coupling U, boundary
conditions and time convention are supplied. Observations can distinguish
interaction-dependent interference and transport in this model. At N=2,
mN contributes only a common phase. These recorded experiments therefore
do not estimate m or measure the vacuum gap.

## Uniform norm and source-splitting bound

For every finite cubic box and the infinite cubic lattice,
0<=L<=12I. Hence, on N=2 and after omitting 2mI,

    0 <= H'_R = H_R-2mI <= B I,    B=24J+U.

The contact term is U times the projector onto configurations with both
particles at the same site. Color the graph edges by coordinate axis and
the parity of the lower endpoint along that axis. These six color classes
are disjoint matchings. Each one-particle matching Laplacian has norm at
most 2, so its two-particle kinetic generator has norm at most 4J. The
onsite contact generator has norm at most U. Thus the sum of the seven
generator norm bounds is S=24J+U, independently of the number of sites.

A symmetric palindrome of the six matching half-steps around a full
onsite step agrees with exp(-ihH') through degree two in h. The third
derivative norm of each expression is at most S^3: for the palindrome,
the product-rule multinomial sum is the cube of the sum of the scaled
generator norms; all other factors are unitary. Integral Taylor
remainders give a one-step difference at most S^3 |h|^3/3. Telescoping n
unitary steps gives

    ||U_n(T)-exp(-iTH')|| <= S^3 |T|^3/(3 n^2).

The operator distance is also at most 2. This theorem concerns ideal
source substeps. Numerical contractions, polar corrections, compiler
roundoff and physical service errors need their separate accounting.

For a fixed history, purify the initial state and retain every instrument
outcome in a separate register. Every intervening instrument is then an
isometry; classical feedback is a block-controlled isometry. Telescoping
the complete history therefore bounds its Stinespring vector difference
by the **sum** of the segment operator bounds. For normalized pure
outputs, trace distance (half trace norm) is at most their vector
distance. Partial trace and record measurement contract trace distance.
Thus the same sum, capped at 1, bounds both the complete retained
instrument output and the full classical record-law total variation.
There is no branch-count multiplier and no division by a small branch
probability. The bound also holds with a spectator reference.

## Boundary comparison for the complete history

Compare Lambda_R with any larger centered box, or with Z^3, at the same
lattice spacing and with the same parameters. Embed finite occupation
states by retaining their site labels. The initial support lies in A.
An allowed instrument or intervention preserves every particle outside A
and can rearrange occupations inside A. It must preserve N=2. This
includes an origin-number read, local phase feedback, and repeated two-outcome
Ramsey instrument P_+, I-P_+, where

    P_+ = |phi_+><phi_+|,    phi_+=(D+iS)/sqrt(2).

Such operations do not enlarge the set of occupied sites beyond its
previous union with A. Each Hamiltonian factor can move a particle by at
most one edge; diagonal terms do not move particles. The distance of A
from the smaller box boundary is d=R-1. A product containing at most d
Hamiltonian factors, with any number of the fixed apparatus operations
inserted, agrees under embedding. In particular, reaching the boundary
on the last hopping factor is harmless: its changed degree or missing
outgoing edge can act only in a subsequent Hamiltonian factor. The first
potentially unequal Taylor degree is therefore q=d+1=R.

Let T=sum_j |t_j| over all Hamiltonian segments. Scale their durations
simultaneously by lambda, retain their record registers, and write the
resulting complete Stinespring maps as W_R(lambda), W_infinity(lambda).
For the localized initial state, their derivatives through degree q-1
agree under embedding. The qth derivative of either map has norm at
most (BT)^q by the product rule and isometry norms. Integral Taylor
remainders along real lambda give

    ||W_infinity(1) psi - E W_R(1) psi||
        <= 2 (BT)^q / q!,    q=R.

No exp(BT) factor is necessary because the derivative estimate uses
unitary evolution at real times. The same expression, capped at 1,
bounds full-history total variation and retained-output trace distance.
For adaptive protocols, a fixed sequence of stages with record-controlled
operations is admissible; use the maximum duration bound at each stage.
The elementary receipt API represents stages with fixed durations.

The proof applies on the infinite two-particle Hilbert space because H'
is bounded there by the same B. At fixed finite T the factorial bound
tends to zero as R grows, defining the infinite-volume observed history.
This is a growing-volume limit at **fixed spacing**, not a spatial
continuum theorem. Add the temporal splitting bound separately when
comparing an executed finite-box splitting with the infinite-lattice
continuous Hamiltonian experiment. Floating/service errors remain separate.

## A volume-uniform interaction witness

Use a separate short experiment with initial state psi=(D+S)/sqrt(2),
continuous Hamiltonian evolution and the final P_+ read just defined.
There is no intermediate occupation read: that would destroy the
coherence needed by this witness. Compare H_U=H_0+U C against H_0 with
everything except U fixed, where C projects onto all doublons. Put

    Delta(t) = p_U(P_+,t)-p_0(P_+,t).

Direct differentiation gives Delta(0)=0 and Delta'(0)=U/2. For the
Liouvillian L_H(rho)=-i[H,rho], ||L_H||_(1->1)<=2B and
||L_HU-L_H0||_(1->1)<=2U. Telescoping the kth powers bounds their
difference by 2U k (2B)^(k-1). Consequently

    |Delta(t)-Ut/2| <= 2Ut (exp(2Bt)-1),    t>=0.

Choose t=1/(32B), with U>0. Since exp(1/16)-1 <= 1/15,

    Delta(t) >= Ut/2 - 2Ut/15 = 11Ut/30 > 0.

This proof holds in every box containing the apparatus, as well as on
the infinite lattice. Its positive margin does not shrink with volume.
Subtract both interacting and free temporal error bounds before claiming
a positive ideal source-splitting contrast. Numerical execution alone
does not certify additional roundoff or physical service bounds.
Nonzero correlations without the matched U-deletion control would not
suffice: even the free theory has non-Gaussian Fock preparations.

## Vacuum and spectral gap

In every particle sector H>=mN. The vacuum has energy zero, while every
nonvacuum sector has energy at least m. On each finite open graph the
constant one-particle vector has kinetic energy zero and no contact
energy, so the vacuum gap is exactly m. On Z^3 the one-particle
dispersion is

    E(k) = m + 2J sum_{a=1}^3 (1-cos(k_a)).

Arbitrarily broad one-particle wave packets approach energy m. Hence
the infinite-volume excitation spectral infimum is m as well; it need
not be a normalizable eigenstate at that threshold. The direct sum of
the bounded, self-adjoint finite-particle Hamiltonians defines the
self-adjoint full Fock Hamiltonian on its usual direct-sum domain.

This is a stable spectral gap of the **supplied** repulsive lattice
action. It is not the Yang--Mills mass-gap problem, a relativistic
particle-mass derivation, a gap generated by repair, or an observed gap
in the fixed-two-particle history. Relativistic propagation, a physical
clock, identification with nature, further matter/gauge sectors and the
joint relativistic continuum construction remain separate obligations.
