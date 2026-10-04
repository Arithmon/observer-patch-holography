# Fermionic clocks and reads from the proper-code source

The coherent source construction in [COHERENT_SOURCE_CLOCKS.md](COHERENT_SOURCE_CLOCKS.md)
implemented the massive clock on vacuum and one-particle states. Ordinary
occupation swaps did not supply its fermionic many-particle extension. Here
an explicit local encoding, deterministic preparation, native circuits and
charged schedule close that construction gap. The resulting process implements
the full even CAR algebra, the massive walk on every occupation sector, the
resolved number instruments and the bounded Fock detector. It also implements
a quartic interaction with a non-Gaussian two-particle output.

This is a realization theorem for the **same declared coherent source with
CCG**, not selection of that source or of fermionic statistics by A1--A3.
CCG, the available drive norms and the specified clock parameters retain
their existing status. No additional transport axiom is introduced. In
particular, the encoding gates, loop constraints and read instruments below
are constructed; they are not extra services granted to the source.

The edge encoding is the Bravyi--Kitaev superfast encoding. Its generators,
loop constraints and bounded-degree locality are prior work:
[Bravyi and Kitaev, section 8](https://arxiv.org/abs/quant-ph/0003137),
[Setia and Whitfield](https://arxiv.org/abs/1712.00446), and
[Setia, Bravyi, Mezzacapo and Whitfield](https://arxiv.org/abs/1810.05274).
The contribution here is its explicit integration with the audited OPH
proper-code source: complete preparation and conversion instruments, an
all-sector clock/read compiler, causal preparation schedule, and the
physical-interface noise/accounting bounds. The standard encoding is not
an OPH discovery or a proof that physical particles must be fermions.

## 1. Source and observable algebra

Keep the parent's full primitive algebra `C12 tensor M6 tensor private
buffers`, complete twelve-port response, scalar static seams and within-fibre
graded process links. A qubit lives in a fixed two-dimensional proper code
of a private buffer. Each edge-qubit owner has its own M6 processor and visitor
buffer; each helper has its own processor and buffer. Two-qubit operations
use the existing four-dimensional pair code in M6. These are extra counted
interfaces, not simultaneous uses of one processor.

Every load, store, pack, separate, wait and flight is the parent's complete
instrument: code projection, complement failure flag, sender reset and
receiver overwrite. In particular, none is a reversible full-algebra
operation on an otherwise inert buffer. All basis changes below are actual
native M6 pulses. They do not add primitive response directions or alter the
port Lie algebra. The parent's support construction accommodates these finite
private factors and process links within one support fibre; its static seams,
degree-one spherical bridge and retained-factor refinement are unchanged.
This is an extension of its constructed support, not a claim about the
separate simulator's currently implemented instruction set.

Let G=(V,E) be a connected simple graph with one qubit per edge, a fixed
orientation and an ordering of edges incident to each vertex. For canonical
orientation i<j define

```
B_i = product_(e incident to i) Z_e,
A_ij = X_(ij) product_(e <_i ij) Z_e product_(e <_j ij) Z_e,
A_ji = -A_ij.
```

On fermionic Fock space write `b_i=1-2 c_i* c_i` and
`a_ij=-i gamma_i gamma_j`, `gamma_i=c_i+c_i*`. Both versions square to one;
all B commute; A_ij anticommutes with B_i and B_j and with precisely those
distinct A edges sharing one endpoint. Other pairs commute. For each closed
simple path C=(v_0,...,v_l=v_0), impose the signed relation

```
S_C = i^l A_(v_0 v_1) ... A_(v_(l-1) v_0) = I.                 (1)
```

The S_C are commuting Hermitian Pauli involutions and commute with every
A and B. Their X supports are the binary cycle vectors. A cycle basis has
r=|E|-|V|+1 independent constraints, so its common +1 eigenspace has
dimension `2^(|E|-r)=2^(|V|-1)`. Products of loops obey the cycle relations:
the signs follow by the displayed anticommutation rules, or by reducing a
walk one backtrack at a time with `A_ij A_ji=-I`. Thus a basis enforces (1)
for every loop, including its sign, not just its binary support.

Since `product B_i=I`, this code represents the even total-parity block.
For the odd block replace B at one root by -B. Then `product B_i=-I`;
all other relations remain valid. A fixed central parity tag chooses the
block, with its initialization and record counted. This is not a new physical
superselection postulate. For the specified parity-preserving operations,
deleting coherences between the two parity sectors leaves every transcript
probability unchanged. Both blocks are needed for the full even observable
algebra; odd fields are not asserted to be commuting local qubit observables.

For completeness an intertwiner can be built explicitly. The state fixed by
all S and the even B is unique. Associate it to the empty occupation vector,
or to the occupied root in the odd representation. Products of A along tree
paths from the root toggle any desired parity-compatible occupation pattern.
Fix their phases by the occupation rule
`c_i |s> = (-1)^(sum_(j<i) s_j) s_i |s-e_i>`.
The resulting vectors are orthonormal eigenvectors of all B. There are
`2^(|V|-1)` of them, and path independence follows from (1). They form an
isometry J_p with `A J_p=J_p a` and `B J_p=J_p b`. This also proves faithfulness
on each parity block: B distinguish every occupation and connected path
products connect every two occupations of the same parity.

Each A has weight at most 2D-1 and each B weight at most D when degree<=D.
Changing two adjacent edges in one vertex ordering conjugates the entire
representation by CZ on those two edge qubits. Reversing an edge orientation
conjugates by its Z. Adjacent exchanges generate all ordering changes.
States, loops and observables must all be conjugated; changing a state alone
is a fault, not a free rechart or a CCG basis choice.

## 2. Preparation, every outcome, and code extension

Reset all edge qubits to |0>. Measure an independent set of S_C by the QND
circuit in section 3, retaining each sign. Every nonempty product of basis
loops has nonzero X support, and hence zero expectation in |0>. Consequently
all `2^r` outcome strings have probability `2^-r`. Given a syndrome s, choose
Z corrections z with `C z=s` over GF(2), where C is the cycle/edge incidence
matrix. Such a correction changes precisely the required loop signs and
preserves every B. All corrected branches are the same code vacuum (up to
an immaterial branch phase). No outcome is rejected or resampled. In the
odd representation the same physical state has one root fermion.

For a fundamental-cycle basis, one Z on each negative chord suffices. Section
5 gives a local basis and a causal decoder, avoiding a serial traversal of
long fundamental cycles on a large workspace. Pure buffer resets are declared
preparation interventions as in the parent, not an assertion that unconstrained
A3 selects a pure ontic state.

There is also a state-preserving code conversion. Add a leaf vertex and its
edge by appending |0>; the new mode is empty and all old A,B intertwine. To
add an edge between existing vertices, place it last in both local orders,
append |0>, measure its new loop and apply Z_new to the negative outcome.
Writing V for append-|0>, the two corrected maps are

```
K_+ = (I+S_new)V/2,
K_- = Z_new(I-S_new)V/2 = K_+,
K_+* K_+ = K_-* K_- = I/2.                                  (2)
```

Here Z_new V=V and Z_new anticommutes with S_new. Old A are unchanged;
endpoint B become old B times Z_new, which agrees on V. Since the new loop
commutes with the new algebra, (2) intertwines the old algebra on the entire
old code, including any spectator entanglement. Repeated leaf and chord
extensions realize connected graph growth. This is a charged quantum
conversion, **not** raw partial-trace refinement of the encoding. Tracing a
new gauge qubit need not recover the old encoded density. The underlying
source's physical-factor refinement and this conversion are distinct maps.

## 3. Native fermionic gates and complete instruments

The parent's exact response sum is `diag(i I3,0)`. A pi pulse gives CZ on
the four-dimensional pair code, up to a scalar. Its native u(2) controls
give H, S and arbitrary two-dimensional rotations. Thus CNOT is H-CZ-H;
no black-box encoded two-qubit gate is assumed.

For a Hermitian Pauli P, rotate its X factors with H and its Y factors with
S* followed by H. Accumulate Z parity into a helper initially |0> with
CNOT(data,helper), apply `exp(i theta Z_helper)` with P's sign, undo the
CNOTs and return the data basis. The **whole helper unitary** is
`exp(i theta P tensor Z_helper)`, so the claimed `exp(i theta P)` follows
with a clean returned helper. The dense checker tests the whole unitary,
not just its action on |0...0>.

Instead reading the helper after parity accumulation, while returning the
data basis, yields exactly the Kraus maps `(I+P)/2` and `(I-P)/2`. A negative
Pauli sign uses a physical X on the helper. Both outcome maps and their
sum-of-squares identity are retained. Intermediate native pulses can leave
the gauge code; the complete logical circuit returns it. The noise proof
therefore cannot assume the system stays in a fixed occupation sector during
these pulses.

The algebra gives, without an occupation restriction,

```
n_i = (I-B_i)/2,
H_x = c_i* c_j+c_j* c_i = (i/2) A_ij(B_i-B_j),
H_y = -i c_i* c_j+i c_j* c_i = A_ij(I-B_i B_j)/2,
H_z = n_i-n_j = (B_j-B_i)/2.                                 (3)
```

For `g=[[a,b],[-b*,a*]]`, let theta=atan2(|b|,|a|),
xi=(arg a+arg b)/2, zeta=(arg a-arg b)/2. Then
`g=D(xi) exp(i theta sigma_y) D(zeta)`, where D(t)=diag(exp(it),exp(-it)).
The two Pauli terms of each factor in (3) commute. Its Fock lift Gamma(g)
therefore uses six exact Pauli rotations, not a Trotter approximation.
`exp(i alpha n_i)` uses one rotation and a scalar. The Fock lift of a mode
transposition uses Gamma(-i sigma_x) followed by phase pi/2 on both modes:
eight rotations. It includes the fermionic sign on double occupation and
the full occupation-order signs in the presence of other modes. An ordinary
qubit swap agrees on vacuum/one-particle vectors but is wrong in general.

A density interaction is equally explicit:

```
exp(-i g n_i n_j)
 = exp(-ig/4) exp(ig B_i/4) exp(ig B_j/4) exp(-ig B_i B_j/4).  (4)
```

It uses three rotations, with weight at most 2D on any chosen pair and at
most 2D-2 for an edge. Hence the clock's weight bound 2D-1 covers adjacent
interactions as well. Equations (3)--(4) preserve every loop and hold on
both parity blocks. Arbitrary pulse precision is the same declared exact
control idealization as in the parent; bounded errors are separately charged
by channel norm, not silently rounded to an exact identity.

An executed four-mode example begins with modes 0 and 2 occupied and splits
each over its dual rail. The state is
`(|0101>+|0110>+|1001>+|1010>)/2` in integer-bit occupation order. It can be
prepared from the even vacuum by a pi/2 A_02 rotation followed by two SU(2)
mixers; scalar phases are irrelevant. Applying (4) with g=pi to modes 1,3
changes the last amplitude's sign. Its one-body density is I4/2, so
Tr(gamma^2)=1 rather than the value 2 for a pure two-particle Slater state.
With the dual-rail even bilinears X_L,Z_L,X_R,Z_R,
`<sqrt(2)(Z_L X_R+X_L Z_R)>=2 sqrt(2)`.
This is a computed non-Gaussian entanglement witness, not a loophole-free
laboratory Bell experiment. The coupling is declared. The free clock's
dispersion estimates are not asserted for nonzero interaction.

## 4. Finite reflecting workspace and exact clock compiler

Take q^3 cells of spacing a in a finite cube. Each has 16 active modes A
(two masses, two chiralities, four tetrahedral directions) and 16 blank
modes B. The encoding graph contains active coin chains, chirality pairs,
mass-readout pairs, local A_j--B_j pairs, and flight pairs A_(x,j)--B_(y,j).
For an outgoing tetrahedral flight leaving the cube, replace its target by
B_(x,j xor 4), reversing chirality at the same cell. Every B has exactly one
flight predecessor: a missing interior predecessor is replaced by that
reflection. There is no long periodic wrap charged as a short flight.

Axis-neighbor cell hubs A_0 are additional **encoding-only** edges. They
make the graph connected and permit local preparation; no clock hopping is
inserted on them. The graph has

```
|V|=32q^3, |E|=63q^3-3q^2, r=31q^3-3q^2+1,
degree <= 12, interaction-edge length <= ell=sqrt(3)a.        (5)
```

The looser bound |E|<=192q^3 is used in the resource certificate. It avoids
relying on a favourable edge count for noise or memory sufficiency.

Decompose each four-mode coin into adjacent SU(2) Givens rotations and
phases. Apply the precoin on left chirality, swap each active mode to its
unique flight B, restore each B to its local A, apply the postcoin on right
chirality, and apply the local mass mixers. These are **fermionic** swaps.
The resulting one-particle matrix is exactly the finite reflecting version
of `U_m=exp(-im tau beta) diag(FC,CF*)`, with tau=ell/c, and the blank bank
returns to its vacuum. Functoriality `Gamma(UV)=Gamma(U)Gamma(V)` proves
the same statement for every occupied configuration and superposition within
either parity block, including arbitrary spectator entanglement. For general
Fock states it gives the same entire even-observable history; it does not
promise preservation of odd inter-parity coherences in a central classical
tag. The evidence separately
replays the QR/two-bank circuit against a direct reflecting-permutation
formula for all one-particle matrix entries at the tested sizes. The
all-size Fock conclusion uses (3) and functoriality, not a many-body matrix
simulation at astronomical q.

There are 80 logical layers per cell: 16 phases, 32 SU(2) mixers, 16 ferry
swaps and 16 restore swaps, or **464 Pauli rotations**. A Pauli support for
one operation lies in the endpoint stars. Place edge owners at their
geometric midpoints; co-located owners remain distinct interfaces. A helper
visits those owners and returns, using counted pair-code operations. A
closed tour through the stars, shortened past unused vertices, bounds the
required flight length. With w=23 a safe two-tour bound is 8w ell, at most
4w flights, 10w+1 native pulses and 28w+4 code events. Define

```
T_P(w) = 8w ell/c + (10w+1)pi/Omega + (28w+4)eta.             (6)
```

Here eta>0 is one code-event duration and Omega is the native drive-norm
bound. All chosen single-code pulses have logarithm norm<=pi. A weight
bound alone would not control travel for an arbitrary disconnected support;
(6) uses the explicit connected-star/loop routes in this construction.

Anchor each clock operation at its source cell (so each logical layer has
one operation per helper). Its edge-owner support lies in
`[anchor-1.5a,anchor+1.5a]^3`. Anchors with equal
coordinates modulo four have disjoint owner sets. Serializing the 64 colors
and the 464 rotations permits one helper per cell and one processor per
edge owner. Therefore

```
T_macro = 64*464*T_P(23).
```

For Omega=pi c/a and eta=a/c, `K=T_macro/tau=20534513.36899866` is independent
of q and a. This is a conservative existence bound, not a proposed efficient
device. Matter speed and clock beat frequency both divide by K; the native
classical record-flight speed c does not. The parent's time-dilation ratio
is unchanged for velocities and matter speed expressed in this common
wall time. No laboratory clock or energy calibration follows from K.

## 5. Local deterministic preparation in diameter time

Within each cell choose the 31-edge tree: B_j connects to A_j, A_j for j>=8
to A_(j-8), then j>=4 to j-4, and A_3--A_2--A_1--A_0. Every vertex is at
tree distance<=6 from its hub. Let H contain these trees and all axis hub
edges. Use the following cycle basis:

* all xy unit plaquettes in the z=0 plane;
* all xz and yz unit plaquettes at every permitted location;
* for each edge outside H, that edge plus the path in H going to the first
  hub, in coordinate order x,y,z to the other hub, and down its local tree.
  Cancel backtracks in paths within a single cell.

The grid part has `(2q+1)(q-1)^2` cycles, its exact cycle rank. Each extra
cycle contains a unique non-H edge. The following decoder proves independence
and completeness constructively. Set corrections on internal tree edges,
all z grid edges, y edges at z=0, and x edges at y=z=0 to zero. Denote the
remaining grid-edge Z-correction bits by x(x,y,z), y(x,y,z). Set

```
x(x,y+1,0)   = x(x,y,0)   XOR s_xy(x,y,0),
x(x,y,z+1)   = x(x,y,z)   XOR s_xz(x,y,z),
y(x,y,z+1)   = y(x,y,z)   XOR s_yz(x,y,z).                    (7)
```

Run the y prefix first, then the two z prefixes in parallel: 2(q-1)
neighbor-message rounds. For an extra edge e, set its bit to its own syndrome
XOR the corrections on its short H path. Only at most three grid-edge bits
enter, all within one unit cube. Disseminate those grid bits to the cube's
anchor in three coordinate sweeps, then send the conditional correction to
the edge owner in one additional round. Each sweep can carry the at most
24 oriented grid-edge slots of a cube as separate recorded bits; 1024
local send/receive/XOR events per cell per round exceeds even six directions
of both sends and receives, all 48 extra-cycle XORs and their bookkeeping.
No instantaneous global syndrome processor is used.

Equation (7), followed by the extra-edge assignment, gives C z=s for every
s. There are r cycle rows and exactly r non-tree correction variables, so
C restricted to those variables is invertible. This proves the decoder for
all `2^r` outcomes, rather than selecting successful measurement records.
The independent finite check reconstructs C and obtains its inverse by
Gaussian elimination, without calling (7).

Every extra loop has length<=6+3+6+1=16 and Pauli weight<=23*16=368.
Anchor each loop at its minimum cube corner. Its support lies in
`[anchor-.5a,anchor+1.5a]^3`, also separated by the 64-color schedule. Per anchor there
are at most 13 internal non-tree edges, 16 reflected flights, 16 diagonal
flights and three plaquettes: at most 48 loop instruments. Their local
support can be toured through at most D*16 edges, so the padded T_P(368)
also bounds each QND instrument and helper reset. The complete preparation
schedule is bounded by

```
T_vac <= 10 eta + pi/Omega + 64*48*T_P(368)
                  + (2q+2)(ell/c + 1024 eta).                (8)
```

This includes parallel blank-data and visitor initialization, all loop
outcomes, finite-speed classical decoding, parallel conditional Z pulses,
loads/stores and the parity tag. Ideal commuting measurements need not
be simultaneous; the coloring provides an executable schedule. Under the
drive/event scaling of section 4, T_vac=O(qa/c) at fixed physical size.
Thus deterministic preparation does not secretly require volume times a
long global loop length in elapsed time. Quantum hardware is O(q^3).
The conservative ledger retains O(q^4) event records, including padded
message rounds and idle slots, and O(q^4 log q) central record slots with
their identifiers at fixed physical volume and observation horizon.
This construction does not provide generic fault tolerance.

As in the parent, the fixed instruction table and local counter schedules
are compiled and distributed before blank quantum preparation. That
classical initialization has finite charged prehistory; no quantum coherence
is present then. The bounds here cover subsequent quantum exposure and
do not claim a uniform time bound for arbitrary offline table compilation.
Dynamic syndrome messages, corrections and the final report are included
in the displayed schedules. Static identifiers do not require a fresh
global broadcast at each quantum event.

## 6. Actual resolved reads and detector

For any normalized finite kernel f on the graph, eliminate its components
up a rooted spanning tree. If the current parent/child amplitudes are x,y,
use `[[x*,y*],[-y,x]]/sqrt(|x|^2+|y|^2)`; if both vanish use I. All zero
subtrees remain in the schedule. This sends f to the root vector. The
compiled Gamma(V_f), a QND read of root n, and its inverse therefore have
exact Kraus operators `N_f=c(f)*c(f)` and `I-N_f`, with no discarded branch.
A root phase and return similarly implement `exp(i theta N_f)`.

A spanning tree has height H<=3(q-1)+16. At each depth, serialize at most
D children, 32 vertex slots per cell and 64 colors. One forward Givens
pass costs at most `6*D*32*64*H*T_P(23)`; a returned read uses two passes
and its root QND instrument. This is O(qa) wall time with the stated
drive scaling, independent of particle number. The observable comparison
uses the input time of this explicitly scheduled instrument; it is not
instantaneous access to a spatially extended mode.

The [resolved-read theorem](MASSIVE_OPERATIONAL_CLOCKS.md)
now transfers by an exact isometry, not an assumed CAR lift. For its
compact H2 kernels f,g and separated supports it bounds the overlap
`|<g_a,U^(-n)f_a>| <= Q_a(f)/||S_a f||` and the induced N_f probability
change under a finite N_g phase. Norms and all ancilla extensions are
preserved by J_p. Hence its conclusion holds for **every encoded Fock
state**, including high occupations; the restriction is the stated
resolved operation class, not a low-energy state set. At macrostep endpoints
the continuum error and cone are the parent's with t_wall=K t. Intermediate
native pulses are not separately identified with that free field. Their
macro duration vanishes as a tends to zero. Finite boundaries are placed
outside the entire native stencil support cone of the experiment.

The clock readout is also realized, including its multi-particle meaning.
Within each cell, phase and mix each mass pair so that the requested mass
superposition is one output mode. For weight 0<=w<=1, first QND-read its
occupation, reset the measured helper to |0>, then apply its native rotation
with angle `asin(sqrt(w))` on the occupied branch and zero on the empty
branch. Read the helper again. Retain all four classical histories
`(occupation,acceptance)`: their data Kraus maps are respectively `I-n`,
`0`, `sqrt(1-w)n`, and `sqrt(w)n`. Thus the empty-acceptance history has zero
map without postselection. Both reads, the reset and the conditional pulse
have charged positive padded durations, including weights zero and one.
The verifier replays the full instrument, including its quantum outputs;
equality of the final click effect alone would not certify those outputs.
Compose the no-click adjoints, then take the complement. The effect
on the full Fock space is exactly

```
E_theta = I - Gamma(I-F_theta),                              (9)
```

where F_theta is the parent's one-particle effect. The product of the
commuting local no-click factors proves (9). Summing one-particle number
effects on Fock space would instead exceed one and is not this instrument.
The concrete four-mode check includes vacuum, double and higher occupation
in both parity blocks. A full clock read has at most 56 phase/mixer rotations
and 16 padded QND/acceptance slots per cell, hence `72*64*T_P(23)` elapsed
time. Local records are combined by a charged classical octree OR; no
global zero-latency report is assumed.

## 7. Physical noise: a vacuum obstruction and a valid bound

Independent computational dephasing on the **edge qubits** is not independent
fermion-occupation dephasing. After an idle time T with interface rate lambda,
an edge has phase-flip probability `p=(1-exp(-lambda T))/2`. Consider even the
encoded empty vacuum. In any fundamental-cycle basis, the r chord Z errors
have independent unit syndrome columns. Conditional on all other edge errors,
there is precisely one chord-bit string giving zero syndrome. Its probability
is at most `(1-p)^r`. B is unchanged, and B=+1 together with all S=+1 fixes
the vacuum uniquely, so

```
1 - <vac|rho(T)|vac> >= 1-(1-p)^r.                           (10)
```

Thus the parent's unencoded `2 N lambda T` sector estimate is false for
physical encoding noise even at N=0. With r from (5), fixed nonzero physical
volume and a nonzero wait, fidelity bounded away from zero requires
lambda=O(a^3). Fidelity tending to one requires lambda=o(a^3); this rate
also suffices by the all-edge union bound for that idle experiment.
An ideal terminal syndrome repair would restore this particular vacuum
error, but does not prove robustness of noisy preparation, general logical
states or noisy correction. No generic single-qubit error-correction claim
is made for this encoding.

For the actual entire circuit a weaker but valid uniform bound counts every
active noisy quantum interface. With E_up=192q^3 and H_aux=q^3+1, allocate

```
R = 3 E_up + 2 H_aux.                                      (11)
```

These are edge-data buffers, edge-visitor buffers, all their active M6
processors, helper buffers and helper processors. Treat fixed-basis
dephasing of an active M6 as one full-interface channel, including off-code
levels. Classical syndrome/control records are retained central data;
their finite operations are charged in (8). The noise model here acts on
the listed quantum interfaces, not as unmodeled corruption of central bits.

For each interface `||L||_diamond/2<=lambda`; Duhamel/telescoping and channel
contractivity give, through arbitrary gates, measurements, corrections and
spectator systems,

```
D_trace(actual,ideal) <= lambda R T_total.                  (12)
```

It also covers idle buffers and temporary excursions outside the gauge code.
For a fixed positive accounting observable with range E_max, its expectation
error is at most E_max times (12). Apply the parent's positive Floquet
prescription to the finite reflecting one-particle walk U_q: second-quantize
`|i log U_q|/tau` on the valid Fock code, and assign E_max to invalid-code/bank
outcomes. This is the finite walk's own observable, not a compression of the
infinite-lattice logarithm. This
is one fixed finite positive observable; bad syndromes are not dropped.
It remains a declared accounting reference, not a physical Hamiltonian or
proof of energy conservation during driven preparation.

The reflecting boundary does not invalidate finite ideal packet accounting,
although a principal logarithm is nonlocal and a causal support argument
alone would not prove that fact. For every unitary U and |theta|<=pi,
`|theta| <= (pi/2)|1-exp(i theta)|`. Spectral calculus and Cauchy--Schwarz
therefore give the boundary-independent estimate

```
<psi, |i log U| psi>/tau <= pi ||(I-U)psi||/(2 tau).          (12a)
```

For a compact prepared R=+1 packet whose first step is inside the workspace,
the reflecting and infinite-lattice `(I-U)psi` agree exactly. The coin acts
as identity on its R=+1 spinor at zero momentum; the flight estimate
`||F(k)-I||<=sqrt(3)a|k|` and the mass estimate
`||exp(-im tau beta)-I||<=m tau` imply
`||(I-U)psi||<=m tau+sqrt(3)a || |k| psi_hat ||`.
Thus its ideal accounting is at most
`(pi/2)(m+c || |k| psi_hat ||)`, uniformly bounded for the parent's compact
sampled C2 packets. Use the larger mass for a two-mass superposition. This
finite-workspace bound is less sharp than the translation-invariant energy
estimate, but needs no identification of their nonlocal logarithms. The
finite free walk conserves its own accounting observable exactly.

There is a cofinal consequence, not just an isolated finite witness. Keep
the physical cube length L=qa, the native observation horizon and read
resolution fixed. Equations (6)--(8), the tree-pass bound and the report
bound give a uniform `T_total=O(1)` under Omega=pi c/a, eta=a/c. Meanwhile
`R=O(q^3)` and the whole Fock accounting range is `E_max=O(q^4)`.
Consequently the sufficient bounds are

```
state error      = O(lambda q^3),
accounting error = O(lambda q^7).                            (13)
```

For example the explicitly chosen noisy-approximation family
`lambda_q=lambda_0 q^-8`, with fixed finite lambda_0>0, has state error
`O(q^-5)` and extra accounting error `O(q^-1)`. Each finite member has
positive event times, finite drive and strictly positive allowed noise;
its exact CCG reference is the zero-noise construction. Thus full occupation
and the encoded vacuum admit a controlled refinement family with that
drive/noise capability. The q^-7 sufficient accounting scale comes from
the full spectral range and is **not** a proved necessary threshold.
Neither it nor a decreasing hardware noise rate is selected by the axioms.
The necessary vacuum obstruction (10) remains q^-3.

## 8. Finite analytic clock witness and reproducibility

Use the parent's interval-certified packet: c=3, a=10^-9, masses 100 and
100.1, reference velocity .6, envelope sigma=10, spatial preparation cutoff 80, detector
scale 100 and extent 600. Its ideal clock swing exceeds .6753 and its beat
is .08 with dilation factor 1.25. Prepare it by the reversed tree pass in
section 6, starting with the odd code root. Take q=2^41 and
`n=ceil(25 pi/tau)=136034952318`. Place its compact initial support in the
cube interior: initial coordinate half-extent<=90 plus n*sqrt(3)a<236 is less than 1024, while the
half-width exceeds 1024. Therefore the finite reflecting boundary cannot
affect this experiment. It is not necessary to assume Gaussian tails vanish;
the parent already charges its finite packet cutoff.

The conservative charged ledger gives:

| Quantity | Bound or value in the declared units |
| --- | ---: |
| Gauge preparation, (8) | 1.504 x 10^6 |
| Packet preparation | 3.884 x 10^11 |
| Macrostep wall duration | .011855607 |
| Final instrument | .001839664 |
| Octree report latency | 1269.607 |
| Total exposure | 3.900 x 10^11 |
| R | 6.147 x 10^39 |
| E_max=32q^3 pi/tau | 1.852 x 10^48 |
| Sufficient positive lambda cap | 1.1265 x 10^-102 |
| Trace-distance allowance | 2.701 x 10^-51 |
| Extra accounting error | less than .0051 |
| Clock swing after this error | greater than .6752 |

The verifier recomputes the analytic quantities with independent 70-digit
arithmetic and replays the parent's interval bound. Binary receipt values
are rounded conservatively; one-sided checks enforce every reported upper
or lower bound and their composed exposure, trace-distance and accounting
inequalities. Relative numerical agreement alone does not certify a bound.
The report uses 41
octree stages, length at most sqrt(3)a(q-1) and 64 events per stage. Every
quantum operation, idle/color slot, setup message and report gets a finite
retained identifier and enough central record slots. No huge lattice state
was executed: this is an analytic finite-existence witness. These deliberately
loose hardware, time and accuracy requirements demonstrate consistency, not
practical performance or an empirical prediction. In particular, (10) and
the global accounting budget prevent a claim of population-independent
physical robustness.

The package [code/m1_fermionic_source](../../code/m1_fermionic_source/README.md)
contains the exact GF(2) decoder, source circuits, finite executions and
independent verifier. It checks fermions from occupation signs and exterior
minors; it does not trust producer Pauli arithmetic. Local-cycle preparation
is independently inverted over GF(2). The large witness is evaluated by
separate high-precision arithmetic. Small finite checks support the analytic
all-size proofs above; they are not substitutes for those proofs. The
receipt binds the parent source, proof, verifier, tests, workflow and the
three scoped claim rows, and hostile controls test actual mathematical
rejection as well as custody.

The closed gap is precise: the named coherent source now supplies **one
local, timed, all-sector fermionic process with actual reads and an
interacting extension**. Its statistics, mass parameters, interaction coupling,
coherence capability and laboratory energy/time identification are not
selected by this construction. Those limits do not leave the encoding,
preparation, many-particle signs or measurement implementation as future
premises of this theorem.
