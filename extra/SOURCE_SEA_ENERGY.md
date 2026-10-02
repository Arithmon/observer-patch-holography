# A filled source vacuum, positive excitation energy and readable clocks

The source now protects quantum outputs and actual public records at fixed
nonzero subthreshold error strength. Its positive observable
`dGamma(|i log U|/tau)` was deliberately kept separate from the generator
of its walk. This note closes that mathematical separation for a new,
explicitly prepared filled-vacuum experiment. It does not rename the old
empty-vacuum experiment or identify a Floquet branch with laboratory energy.

The starting points are the [massive walk](MASSIVE_OPERATIONAL_CLOCKS.md),
its [all-sector source compiler](FERMIONIC_SOURCE_CLOCKS.md), and
[fixed-strength protection](FIXED_STRENGTH_PUBLIC_RECORDS.md). The earlier
finite edge-sea example in `paper/tex_fragments/PAULI_STABILITY_SELECTION.tex`
already explains the elementary particle/hole construction. The work here
is its integration with the actual walk, reflecting boundaries, local
continuum vacuum, complete clock read and paid noisy preparation.

Particle/hole normal ordering and Slater determinants are standard. Relevant
precedents include [Gupta and Short](https://arxiv.org/abs/2412.03466) on the
Dirac sea and the modular-energy boundary in discrete time, and
[Kivlichan et al.](https://arxiv.org/abs/1711.04789) on preparation by adjacent
fermionic rotations. We give the bounds used here explicitly. In particular,
a free filled sea is not asserted stable under unspecified interactions.

## 1. Use the same walk, including its boundary

Write `tau=sqrt(3)a/c`, `mu=m tau`, `0<mu<pi/2`, and
`U=exp(-i mu beta) U_0`. In infinite volume,
`U_0=diag(F C,C F*)`, with the parent's involution C and opposite tetrahedral
flights F,F*. In a finite cube use exactly the parent's reflecting flight:
a flight that would leave the cube is sent to the opposite-chirality bank
at its starting cell. No periodic wrap is a short physical link.

The reflection permutation S obeys `beta S beta=S*`. The precoin and
postcoin are exchanged by beta, so in either volume

```
beta U_0 beta=U_0*,
Re U=cos(mu) Re U_0.                                       (1)
```

Consequently every eigenphase theta of U lies in
`[-pi+mu,-mu] union [mu,pi-mu]`. This proves both gaps, at +1 and -1,
for every finite reflecting cube, rather than inferring them from a torus.
Different masses have separate blocks and use their own mu. Let

```
h=(i/tau) Log U,    P=1_(h<0),    Q=I-P,
K=(U*-U)/(2i)=sin(tau h),    delta=sin(mu).
```

The principal logarithm is well-defined, `|h|>=m`, `|h|<pi/tau`,
`|K|>=delta`, and `P=(I-sign K)/2`. These are finite operator identities.
The gap persists at a boundary even when translation invariance does not.
The number of filled modes need not equal half the finite dimension. In the
one-cell reflecting cube all flights reflect, U_0=beta, and P=0. The finite
operator controls retain this case; the local limit below does not assume
periodic half filling at a reflecting boundary.

In infinite volume the parent gives the complete phases
`+/-omega,+/-(pi-omega)`, each twice, where
`cos omega=cos(mu) sqrt(1-sum sin(a k_i)^2/3)`. All eight valleys and both
bands remain present. We do not use the principal logarithm of U squared:
that would change the energy assigned to the high band.

## 2. A positive generator on the actual fermionic source

For a finite workspace, fill every negative eigenmode of h once. This gives
the Slater state Omega_P, with one-body covariance P. Set

```
E_vac=Tr(h P),       H_sea=dGamma(h)-E_vac I.                (2)
```

In an eigenbasis of h, a positive mode has particle annihilator c_j;
a negative mode has hole annihilator c_j*. CAR give

```
H_sea=sum_(e_j>0) e_j c_j* c_j
      +sum_(e_j<0) |e_j| c_j c_j* >=0.                    (3)
```

On the active logical Fock space, the filled state is its unique zero-energy
vector. Helper and transport banks have the parent's specified blank inputs
and return conditions; no unique ground state of an apparatus Hamiltonian is
asserted. The source does not preserve arbitrary odd inter-parity coherences
through a classical parity tag: the full even-observable instrument is the
comparison, as in the parent. The filled state is a ground state under all
finite unitary perturbations for this specified Hamiltonian. The full Fock
identity, including every parity sector and external spectator, is

```
exp(-i n tau H_sea)=exp(i n tau E_vac) Gamma(U)^n.           (4)
```

The scalar has no effect on the complete even-observable source instrument.
Particles and holes both have positive excitation costs |e_j|, but a hole
is an absence in a filled mode. Replacing h by |h| on the old occupation
space changes the dynamics of negative modes and is not (2).

The spectrum of H_sea is every subset sum of the positive numbers |e_j|.
Thus its finite excitation partition function equals the parent's
accounting partition function, while its vacuum and field interpretation
are different. In the parent's fixed-side even-torus thermodynamic
comparison the limit is still

```
log Z_a -> 32 sum_(n in Z^3) log(1+exp(-b sqrt(m^2+(c/3)^2|2 pi n/L|^2))).
                                                               (5)
```

This uses the already proved complete-spectrum dominated convergence, not
a discarded band or a truncated sum. The torus in (5) is a spectral
comparison; the implemented source below has reflecting boundaries.

## 3. The boundary cannot determine a distant vacuum read

We need more than finite-speed propagation: P is a spectral projection and
is not finite range. Here is a direct exponential estimate, using the
standard weighted-resolvent argument often called the Combes--Thomas method.

Measure separation in tetrahedral flight layers, or conservatively by
sup-norm cell distance. A flight changes the distance to a set by at most
one. For a bounded 1-Lipschitz distance function d and `W=exp(s d)`, onsite
coins commute with W and the weighted flight permutation differs by at most
`exp(s)-1` in operator norm. Hence

```
||W K W^-1-K|| <= exp(s)-1.
```

Choose `s=log(1+delta/2)`. The resolvent at it has norm at most
`1/sqrt(delta^2+t^2)`. Its weighted perturbation has norm at most delta/2;
the resolvent identity therefore bounds the weighted resolvent by
`2/sqrt(delta^2+t^2)`. Let B be the boundary layer where reflecting and
infinite flights differ, and let A be a tested region at distance r from B.
Decouple the exterior with the same reflection rule. Both pieces retain
(1), so the direct-sum comparison K_tilde has the same gap. The difference
K-K_tilde has norm at most 2 and is supported at B. The resolvent identity
and one off-diagonal weighted estimate give

```
||(P-P_tilde) 1_A|| <= (2/delta) exp(-s r)
                    <= (4/delta) exp(-s r).                (6)
```

Indeed the sign-difference integral is
`(1/pi) integral_R [(K-it)^-1-(K_tilde-it)^-1] dt`;
the extra factor 1/2 for P and the integral
`integral_R dt/(delta^2+t^2)=pi/delta` give the first bound.
Bounded distance weights can be increased to the actual separation after
the estimate, so no unbounded similarity is assumed.

For any target epsilon choose the *integer* margin

```
r=ceil(log(4/(delta epsilon))/log(1+delta/2)).                (7)
```

This controls the entire projected vector from A, not only a diagonal
density. At fixed m,c and a=L/q, an epsilon=q^-12 margin has physical width
`ar=O(log q)`. The cube also contains the whole native causal cone of the
named finite experiment. This costs O(q^3 log^3 q) active modes. It is a
proved local vacuum comparison; the two global sea vectors need not have
large fidelity and their Hilbert spaces need not be identified.

## 4. A local continuum vacuum, with the high band retained

Use the parent's `R=I_2 tensor C` and
`H_m(k)=m beta+(c/3) alpha.k`, with `H_m(k)^2=E(k)^2 I` and `[H_m,R]=0`.
Taylor expansion of the *one-step* walk gives

```
||(U*-U)/(2i tau)-R H_m(k)||
  <= (m tau+sqrt(3)a|k|)^2/(2 tau) =: z_a(k).               (8)
```

To see the cancellation, put `B=diag(k.s)` and
`A=diag(B,-CBC)`. The first derivative of U at a=0 is
`-i [m tau beta+a A] R`; its anti-Hermitian part averages A with RAR,
which is precisely tau R H_m. The second derivative of the product of
the mass and flight unitaries is bounded by `(m tau+sqrt(3)a|k|)^2`.
The integral Taylor remainder gives (8) without an exponential prefactor.

For self-adjoint X,Y gapped by g, the same sign resolvent integral gives
`||1_(X<0)-1_(Y<0)||<=||X-Y||/(2g)`. Here take
`g=min(m,sin(mu)/tau)>=2m/pi`. Therefore

```
||P_a(k)-P_0(k)|| <= min(1, pi z_a(k)/(4m)),
P_0(k)=(I-R H_m(k)/E(k))/2.                                (9)
```

This is a statement about all eight internal components. On R=+1 the limit
is the negative-energy Dirac projection. On R=-1 it is the positive-energy
projection: the high Floquet band has the opposite filling in the same
low-momentum chart. Calling the entire eight-component limit a single
ordinary Dirac vacuum would be false. At another valley the scalar sign
of U exchanges the R labels, with the same four low excitations per valley.

For normalized Fourier data f, the squared norm of
`(P_a-P_0)f` is bounded by
`[pi z_a(K)/(4m)]^2+4 ||1_(|k|>K) f||^2`, with the sampled alias error added
by contractivity. The parent's H2 sampling inequality supplies
`A_a(f)=O(a^2)|| |k|^2 fhat||`. Taking K=a^-1/4 gives norm error O(sqrt(a))
for each fixed H2 packet, and (6) transfers it to the reflecting source.
For normalized resolved kernels f,g this controls every covariance
`<f,Pg>`. Wick's determinant formula then controls each fixed finite
polynomial in smeared CAR fields. One elementary bound for an r by r
determinant is `r! r epsilon` when all entry errors are at most epsilon
and all entries have modulus at most one. The continuum state is the
quasi-free state with projection P_0; convergence is local in this stated
observable class, not global sea-vector convergence.

## 5. Local operations have a positive, finite energy cost

Let f be a normalized source mode, `p=<f,Pf>`, and apply the *existing*
resolved operation `exp(i theta N_f)`, where `N_f=c(f)*c(f)`. Its one-body
unitary is `V=I+(exp(i theta)-1)|f><f|`. The resulting covariance is VPV*.
Writing `e_+=<Qf,|h|Qf>` and `e_-=<Pf,|h|Pf>`, direct expansion gives

```
Delta E=Tr h(VPV*-P)
       =4 sin(theta/2)^2 [p e_+ +(1-p)e_-]
       <=4 <f,|h|f>.                                    (10)
```

Its expected particle-plus-hole number is
`8 sin(theta/2)^2 p(1-p)<=2`. The signs and constant vacuum subtraction
are essential. The source does not obtain free excitations by filling a sea.

For fixed compact H2 kernels in R=+1, the parent's bound on the *full*
|h| expectation, including high-band leakage, makes (10) uniform as a goes
to zero. The same statement holds in a reflecting workspace: for every
unitary, `|theta|<=(pi/2)|1-exp(i theta)|` yields
`<f,|h|f><=pi ||(I-U)f||/(2 tau)`. Away from the boundary, the one-step
source action on such a smooth R=+1 packet has `(I-U)f=O(a)`.
No nonlocal logarithm compression is substituted for the finite Hamiltonian.
For the additional normalized particle used below,
`chi=Q psi/sqrt(1-<psi,P psi>)`, its excitation energy also obeys
`<chi,h chi><=<psi,|h|psi>/(1-<psi,P psi>)`. The proved blocking bound
stays below one, so the clock's extra particle has uniformly finite energy.

The full CAR read/encoder commutator bounds of the parent are state-uniform.
They therefore hold in Omega_P and its excited states without alteration.
We retain compactly supported full CAR modes; projecting a field onto Q and
then pretending it is a compact local field would not preserve that proof.
General cell-scale operations and the R=-1 chart need not have bounded cost.

## 6. A clock that is read above the filled background

The old many-site `I-Gamma(I-F)` click detector cannot simply be reused:
occupied background modes can make its no-click probability tiny. A filled
mode within the accepted subspace already makes the click certain. This
failure is retained as a control. We use the already compiled *single
resolved-mode* QND instrument `{N_f,I-N_f}` instead; both outputs are real
source records, with no vacuum subtraction or postselection.

Take the parent's two masses, common velocity u, `v=c/3`,
`gamma=(1-|u|^2/v^2)^-1/2`, `k_j=gamma m_j u/v^2`, and common normalized
spinor xi in R=+1 with `H_mj(k_j)xi=gamma m_j xi`. Let g_j(t) be a normalized
Gaussian envelope of width sigma centred at ut, times `exp(i k_j.x) xi`,
in mass block j. Use the parent's finite smooth cutoff to radius P and
sample on the actual lattice, normalizing the resulting modes. The target
read mode is

```
f_theta(t)=(g_1(t)+exp(i theta) g_2(t))/sqrt(2).
```

The mode shapes contain translation and the fixed plane waves, but not the
unknown relative clock phase. Prepare the Slater state with covariance

```
C_0=P+|chi><chi|,    chi=Q f_0(0)/||Q f_0(0)||.             (11)
```

Section 7 constructs this deterministically; (11) is not a successful
creation branch selected out of a failed preparation. The initial blocked
fraction is `b_0=<f_0(0),P f_0(0)>`, and the one-particle trace distance
between chi and f_0(0) is exactly sqrt(b_0). At t=n tau the actual binary
read probability is

```
p_theta(t)=<f_theta(t),P f_theta(t)>
            +|<f_theta(t),U^n chi>|^2.                    (12)
```

It lies in [0,1] because C_t is a projection. Since the sea has no
cross-mass covariance, the first term is independent of theta. It is
nonetheless retained in every probability. At zero time and theta=0,
the probability is exactly one, even for nonzero b_0.

For uncut continuum reference envelopes, the parent packet estimate bounds
the difference from the rigid translated state by

```
e(T)=sqrt(3)v/(m_min sigma)+sqrt(15) T v^2/(8 m_min sigma^2).
```

The rigid relative phase is `Delta m t/gamma`. The reference probability
for this resolved-mode read is `(1+cos(Delta m t/gamma+theta))/2`.
It has unit swing; its detector is different from the parent's pointwise
weighted detector and we do not reuse that detector's visibility.

Let d be the parent's complete lattice dynamic error through T on |k|<=K,
A its sampling error, T_K its Fourier-tail norm, and b_cut a bound on each
normalized cutoff-mode error. Put
`rho=pi z_a(K)/(4m_min)+2 T_K+4 A+4 b_cut+epsilon_boundary+8 exp(-1000)`.
Here z_a uses m_max, and the positive alias term applies to the finite witness.
The Gaussian negative-energy weight is at most
`3 v^2/(4 m_min^2 sigma^2)`: the Dirac negative projector varies by at most
`v |k-k_j|/m_min` on the reference positive spinor. Equation (9) and the
triangle inequality then give the conservative uniform bound

```
b_0,b_t <= b=3v^2/(4m_min^2 sigma^2)+rho,
|p_theta(t)-(1+cos(Delta m t/gamma+theta))/2|
 <= b+sqrt(b)+4(e+d+4A+2T_K)+8 b_cut
       +epsilon_boundary+8 exp(-1000).                   (13)
```

For clarity, compare first with the unprojected packet: replacing its input
ray by chi changes any probability by at most sqrt(b_0), independently of
time. Its overlap with the normalized resolved probe has the parent's
continuum, sampled-dynamic and tail errors; the inequality
`||<f,x>|^2-|<g,y>|^2|<=2(||f-g||+||x-y||)` supplies the factors in (13).
The remaining positive background is at most b_t. Projection comparison
uses (9) on the Fourier ball and contractivity outside it. An H2 sampling
error A changes a normalized projection expectation by at most 4A.

Use the parent's product quintic cutoff, equal to one in the cube of
half-width P and zero beyond P+sigma, centered at the prescribed packet or
probe center. For every translation, summing Gaussian tails over lattice
cells bounds the discarded mass by `3 erfc((P-a)/(sqrt(2)sigma))` times
the normalization correction. Poisson summation bounds that correction by
positive exponential aliases uniformly in the translation. The witness has
`2 pi^2 sigma^2/a^2>1010`, so it is covered by the retained exp(-1000)
allowance. Normalizing a truncated unit vector changes its norm by at most
twice the discarded norm. Thus the convenient uniform bound
`b_cut=4 sqrt(3 erfc((P-a)/(sqrt(2)sigma)))` is conservative for both
preparation and moving probe; the polynomial cutoff has the same support
and no greater discarded norm. The constants in (13) overcount its use in
the two mass components and normalized source state.

The finite witness uses
`c=3, m=(100,100.1), |u|=.6, sigma=10, P=80, a=10^-9`,
`K=max |k_j|+1`, and a boundary margin from (7). The executable outward
interval check certifies probability error less than **0.038** and an
observable native-tick swing greater than **0.923**, through
`T=2 pi gamma/Delta m+tau`. The boundary comparison error is 10^-10 and
(7) gives the integer margin 1,423,037,275 cells. These are analytic finite-size bounds; the
astronomical corresponding three-dimensional apparatus is not simulated.

The input-time label and the availability-time label are different. The
existing returned mode-elimination/QND circuit has positive duration and
charged flights, and its result is exported to the live archive. Its spatial
extent is finite. Translating the predetermined read envelope does not
execute an instantaneous distant measurement or grant a free runtime clock.

## 7. Prepare the sea and pay for it

On any finite connected source encoding graph, choose orthonormal occupied
orbitals for P, and append chi when (11) is requested. Complete them to a
unitary Z. Preparation starts from the occupation product with precisely
these first r modes filled, in the appropriate parity representation.
The parent loop-preparation circuit retains and repairs every syndrome.
Products of edge A operators along tree paths prepare the desired product
occupation without an outcome filter. This costs at most O(M^2) edge
rotations for M modes, including the odd-root convention.

Adjacent-row QR decomposes Z into at most M(M-1)/2 SU(2) rotations and M
phases. All rotations and identity cases are retained. Each logical pair
can be moved together on a spanning tree by at most 2M fermionic swaps,
then restored; the actual fermionic signs and overwritten source buffers
are the parent's compiled operations. This deliberately conservative
construction has O(M^3) native logical rotations and depth O(M^3).
Blank routing banks can be included as unoccupied spectator modes in the
preparation unitary. They increase the graph size by a fixed factor, already
absorbed in the O(M^3) bound, and are restored after routing. The finite
routing table uses the total graph-mode count for its explicit integer bound.
It needs no new full-M6 reversible command. The executed small circuits
verify their final covariance, full occupation amplitudes and source gates;
an arbitrary change of occupied-orbital basis leaves the prepared density
unchanged and must pass verification.

This preparation is much slower than the parent's empty-vacuum preparation.
Its duration and its stored public history must be counted. With
`M=O(q^3 log^3 q)`, a conservative complete request/depth envelope is
`G=O(q^9 log^16(2q))`, including the later finite-horizon clock experiment,
mode reads and all routing. Here the seven extra logarithms beyond M^3
are a conservative allowance for the parent's less-than-four-power
Solovay--Kitaev synthesis, finite angle words, addresses and schedule
bookkeeping. The orbitals and rotation words are a known compiled program,
as in the parent's source capability. This is no unknown-input encoder and
no claim that the source autonomously discovers or selects its program.
Offline classical diagonalization is not identified with a native quantum
operation or assigned a physical work cost. Running the existing protection and
archive compiler with k levels and N-bit archive words gives the envelope

```
active inventory       O(G N^2 s^k),
scheduled depth        O(G N^2 s^k),
active locations       O(G^2 N^4 s^(2k)),
diagnostic inventory   O(G^2 N^4 s^(2k) log q),
storage-time volume    O(G^3 N^6 s^(3k) log q).              (14)
```

Fixed native factors include all edge-code helpers, leakage reduction,
complete resets, export entry corrections and the expander's large degree.
For fixed small eta with A eta<=1/4, take
`ell=ceil(log2(2q)), k=ceil(log2(64 ell)), N=m_arch^2`,
`m_arch=ceil(sqrt(4096 ell))`, and `kappa=ceil(log2 s)`.
Then protected bad strengths are at most A^-1(2q)^-128 and archive/export
bad-region norms at most (2q)^-128. The respective envelopes in (14) are

```
O(q^9 log^(18+kappa)(2q)), O(q^9 log^(18+kappa)(2q)),
O(q^18 log^(36+2kappa)(2q)), O(q^18 log^(37+2kappa)(2q)),
O(q^27 log^(55+3kappa)(2q)).                               (15)
```

Synthesize each requested variable rotation to q^-64; the added protection
and Boolean gates remain in the exact fixed library. The complete joint
instrument error is therefore `O(q^-55 log^(36+2kappa)(2q))`. Since
`0<=H_sea<=M pi/tau=O(q^4 log^3 q)`, the extra excitation-energy error
is `O(q^-51 log^(39+2kappa)(2q))`. This expectation is evaluated after the
parent's complete decoding map; an explicit failure/abort can be assigned the
same finite E_max. It is not the energy of a noisy raw processor or helper.
These estimates cover the long preparation,
its idles, every retained named outcome and subsequent clock, using the
parent's disjoint-region norm argument, not independent success conditioning.
Diagnostic storage noise is separately charged and unused in control.

Every finite program is thus protected at the same fixed positive admitted
noise strength. The ideal finite signal in (13) survives on a cofinal tail.
The enlarged prehistory is not a fixed-time preparation claim. After it,
the original clock's local schedule still has the parent's constant wall-time
dilation; its noisy implementation uses the same paid faster local drives.
No physical energy is assigned here to those drives or to entropy erasure.

## 8. What energy selection is still impossible to infer

For any integer-valued spectral function J commuting with h,
`h'=h+(2pi/tau)J` has exactly the same one-step U. Two explicit choices matter.
Taking J=P makes h' strictly positive and selects the empty ground state.
Taking J=sign(h) retains the *same filled sea* and the *same integer-time
channels* while increasing every particle and hole cost by 2pi/tau.
Both normal-ordered excitation Hamiltonians are positive. Thus positivity,
the chosen sea, and all integer-time channels together still do not fix
the energy scale of individual excitations. Their partition functions differ.

The new result is an exact positive **stroboscopic** generator of the source
experiment, a local vacuum limit and finite-energy, readable excitations.
It does not identify this interpolation with the apparatus's intra-step
drive Hamiltonian. That identification would have to include the supplied
pulses, controls and work stores. Modular-energy resonances under added
interactions are not ruled out by the free gap or ground-state theorem.
This boundary is demonstrated by the two branch controls, rather than
hidden in a generic appeal to a stable physical vacuum.

The [contract](../code/m1_sea_energy/CONTRACT.md) is met by these constructed
objects and quantitative estimates. A1--A3 source selection, an autonomous
clock and empirical energy calibration retain their existing status.
