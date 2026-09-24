# A sharp native-time budget for coherent quantum transport

## 1. The operation class and the clock

The physical comparison arena is the affine three-dimensional response space
E used in `code/source_selection_model/RECORD_GLUING.md`. Its operational
attachment is a proposed law, not a consequence of the internal gauge algebra.
Fix a native flight speed bound c>0. A flight of duration t>=0 has mutually
orthogonal internal projectors Pi_s, summing to the identity, and displacements
b_s satisfying |b_s|<=c t. Its Fourier matrix is

    F(k) = sum_s exp(-i k.b_s) Pi_s.                         (1)

It is unitary because these projectors are orthogonal. During the flight the
modes follow the straight segments b_s u/t, 0<=u<=t. Distinct modes move in
parallel. An onsite unitary changes internal states only when they coincide;
all its elapsed time, or a declared ideal zero-time limit, is counted. A word
consists of flights, onsite unitaries and waits. Its duration T is the sum of
its successive stage durations, not the length of the net endpoint displacement.
Independent branches are parallel, rather than serial software evaluations.
All intermediate positions, controller operations and records are retained.

This is coherent transport of modes, not copying an unknown quantum state.
On the fermionic Fock space a flight permutes modes and an onsite one-particle
unitary lifts by exterior powers. Classical occupation records can be carried;
classical copies require separately prepared registers and onsite controlled
operations. The scalar speed constraint bounds every actual mode trajectory.
Arbitrary concatenations therefore have no native influence outside |x-y|<=cT.
If arbitrary directions and finite local control are admitted, flights/waits
attain the round operational cone. A fixed finite direction menu instead has
its own polyhedral reachability cone. Neither its hull nor a continuum field
speed may be substituted for the charged elapsed time of a particular word.

In particular, independent record flights permitted by RG still have speed
c even if an emergent matter mode is slower. Restricting observation to that
matter mode does not remove those admitted faster record operations.

## 2. The exact first-order budget

Let U(k) be a word of duration T>0. At k=0 let an isometry V from a two-state
space select a degenerate eigenband U(0)V=exp(-i phi)V. The velocity matrices
of this band, per native elapsed time, are

    A_i = (i/T) V* U(0)* (partial_i U)(0) V.                 (2)

The derivative of the band projector contributes no first-order eigenvalue
term at a degenerate band. Differentiating the product in (2) gives a sum of
compressed, unitarily conjugated flight displacement operators. For flight r
of duration t_r>0 this has the form

    (t_r/T) sum_s v_rs,i E_rs,
    v_rs=b_rs/t_r,  |v_rs|<=c,
    E_rs=V* W_r* Pi_rs W_r V >=0,  sum_s E_rs=I_2.           (3)

Here W_r is the appropriate zero-momentum prefix. Onsite derivatives vanish.
Waits and onsite time can be included with v=0. Thus combining all stages
gives a positive operator-valued velocity resolution

    A_i=sum_a v_a,i E_a,   E_a>=0, sum_a E_a=I_2, |v_a|<=c. (4)

The weights t_r/T are essential. Quantum interference and additional internal
states are retained in E_a; (4) is not a classical mixture of output states.

Write E_a=w_a(I+r_a.sigma), with w_a>=0, |r_a|<=1, sum w_a=1 and
sum w_a r_a=0. Zero-weight effects contribute nothing. Decompose

    A_i=u_i I + sum_j M_ij sigma_j.

Then M=sum_a w_a v_a r_a^T. The nuclear norm (sum of singular values) obeys

    ||M||_* <= sum_a w_a |v_a| |r_a| <= c.                  (5)

This follows from the rank-one singular value |v_a||r_a| and the triangle
inequality; equivalently maximize Tr(O^T M) over orthogonal O. The same
trace witness for a d-component Clifford band gives d v<=c for an isotropic
velocity v: (sum_i n_i alpha_i)^2=|n|^2 I bounds every summand in (3).
The two-state three-dimensional case is

    M=v R, R orthogonal  =>  3 |v| <= c.                   (6)

An orientation reversal changes chirality but not the bound. A scalar tilt
u is allowed in necessity; (5) does not assert a full characterization of
the joint pair (u,M). For the untilted case u=0 it is exact.

### Sufficiency and saturation

Take a real M with singular decomposition M=sum_{j=1}^3 s_j u_j r_j^T,
s_j>=0, sum s_j<=c. During the fraction s_j/c of one period T, perform a
conditional flight with velocities +/-c u_j and projectors
(I+/-r_j.sigma)/2. Wait for the remaining fraction. Its first-order
generator is exactly k.M.sigma, and its native elapsed time is T. Every
stage is a declared lossless unitary flight. Hence (5) characterizes all
untilted first-order Weyl matrices in this class. For general real axes this
is a finite addressed experiment; periodic integer-lattice realizations use
commensurate directions or controlled rational approximations.

For M=(c/3)I choose three equal-duration orthogonal flights. Extra coin states,
more directions, any finite number of local scatterings, grouping whole words
into a longer period and changes of time units cannot make the ratio v/c
exceed 1/3. The theorem concerns the stated primitive grammar and a regular
degenerate-band linear limit. General quantum cellular automata need not
possess a lowering into that grammar with the proposed native time.

This is the obstruction relevant to OPH: classical flight completeness plus
this coherent extension cannot give an isotropic Weyl/Dirac sector whose
light speed equals the complete record-cone speed. It is a quantitative
decision of an operation class, not a claim that every quantum source fails.

### Saturation derives a geometric design and its minimal population

For M=(c/3)I, equality in the trace proof of (5) requires, for every
nonzero-weight effect, |v_a|=c, |r_a|=1 and v_a=c r_a. The POVM normalization
and isotropic velocity then give

    sum_a w_a n_a=0,  sum_a w_a n_a n_a^T=I/3,  |n_a|=1.   (6a)

Conversely these two moment identities with E_a=w_a(I+n_a.sigma) and
v_a=c n_a attain the bound. Thus the optimizers are precisely weighted
spherical two-designs (after the allowed orthogonal change of spin axes).
This is derived from saturation; it is not an assumed successful geometry.
In an actual pure flight the effects arise by compressing orthogonal
projectors, and every such finite POVM has a finite isometric dilation.
Distinct equal-velocity outcomes may be merged before counting directions.

The three coordinate vectors sqrt(w_a) n_a span a three-dimensional space
orthogonal to the nonzero vector sqrt(w_a), so at least four nonzero outcomes
are required. With exactly four, the square matrix whose rows are
sqrt(w_a)(1,sqrt(3)n_a) has orthonormal columns by (6a), hence orthonormal
rows. Its row norms force w_a=1/4, and its off-diagonal products force
n_a.n_b=-1/3. Therefore a regular tetrahedron with equal weights is the
unique minimal optimizer up to rotations, chirality and relabeling.

This minimum counts distinct velocity outcomes across the complete word.
For a single flight their nonzero orthogonal projectors require internal
dimension at least four. It is not a lower bound of four on the internal
dimension of a multi-stage word: the three-axis construction below reuses
two internal channels in three stages and has six velocity outcomes.

An explicit four-channel realization uses b_s=a s_s with
s_s=(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1). Define

    S = [[0,1,1,1], [1,0,i,-i], [1,-i,0,i], [1,i,-i,0]],
    C=S/sqrt(3),  P_+=(I+C)/2,  P_-=(I-C)/2.               (6b)

Here S*=S, S^2=3I, and each P has rank two. Apply the onsite unitary C,
then fly channel s by b_s at speed c. The period is T_a=sqrt(3)a/c and

    U_tet(theta)=diag(exp(-i theta.s_s)) C.                 (6c)

Let B_i=diag(s_s,i). Direct multiplication shows that
X_i^+=sqrt(3)P_+ B_i P_+ obey the Pauli algebra on P_+;
X_i^-=sqrt(3)P_- B_i P_- obey its opposite orientation on P_-.
The compressed effects are

    P_+ |s><s| P_+ = (P_+ + sum_i s_s,i X_i^+/sqrt(3))/4.

This realizes the uniquely minimal four-direction design as one charged
flight, with coherent onsite mixing and no hidden sequential axis motion.
Its low-energy speed is (a/sqrt(3))/T_a=c/3, as required by the theorem.
Four channels are quantum modes, not four of the source's twelve central ports.

The same process carries an exact native record at speed c, without invoking
an additional transport menu. In an empty-register background, encode zero
by leaving the source empty and one by the locally prepared one-particle
state C|s>. The onsite coin gives C^2|s>=|s>, then
the flight sends it to displacement a s in time sqrt(3)a/c. The receiver's
occupation effect has probabilities zero and one. Every direction s has
this lossless experiment. Thus the factor-three clock difference is an
attained observable response within this very walk, not only a hull bound.
These are declared finite input states; the filled quasienergy reference
for the thermal comparison is unchanged. The preparation uses the complete
finite channel family; no uniform bound
on its continuum quasienergy is claimed. Restricting admissible interventions
to a bounded-energy low band would change the complete-clock requirement.

The whole spectrum is accessible, rather than only the chosen band:

    det(lambda I-U_tet)=lambda^4
      -(2/3)(cos(2x)+cos(2y)+cos(2z))lambda^2+1.            (6d)

Writing cos(2epsilon)=(cos(2x)+cos(2y)+cos(2z))/3,
0<=epsilon<=pi/2, the four eigenvalues are +/-exp(+/-i epsilon).
The zero-quasienergy points are exactly {0,pi}^3; even parity selects P_+
and odd parity P_-, so four cones have each chirality. The same points
also contain the pi-quasienergy band. The points {pi/2,3pi/2}^3 are
degeneracies at +/-i, not additional zero-energy species.

Since sin^2 epsilon=(sum_i sin^2 theta_i)/3, distance delta to {0,pi}^3
obeys epsilon>=2 delta/(pi sqrt(3)). With the same declared logarithm and
zero-level convention as below, the excitation energies are epsilon/T_a
and (pi-epsilon)/T_a, each twice. All modes enter the partition function.
The first has the uniform offset bound E>=4c|n|/3; the second is at least
pi/(2T_a). Thus its entire thermal limit is also (11), with eight Weyl
cones at c/3. The pi band has no hidden low-energy contribution.

The whole four-channel dynamics also has an explicit estimate, without
preparing a momentum-dependent eigenband or discarding the high band.
For physical momentum k, set B=diag(s_s.k) and
X_i=X_i^++X_i^-=(sqrt(3)/2)(B_i+C B_i C). These X_i commute with C and
obey the three-component Clifford algebra on all four channels. Exactly,

    U_tet(a k)^2 = exp(-i a B) exp(-i a C B C).

The unitary Duhamel product bound is at most
(a^2/2)||[B,CBC]||<=3a^2|k|^2 per pair of periods. Its comparison generator
is B+CBC=(2/sqrt(3))k.X. An unpaired final period has error at most
a||(B-CBC)/2||<=sqrt(3)a|k|; final physical-time rounding costs at most
a|k|/sqrt(3). Thus for j=floor(t/T_a), |k|<=K and 0<=t<=T,

    ||U_tet(a k)^j - C^j exp(-it(c/3) k.X)||
      <= (sqrt(3)cT/2)aK^2 + (4/sqrt(3))aK.                (6d1)

At the other seven zero points replace C by eta C, eta=(-1)^parity.
On the actual zero band P_eta that carrier is the identity; the complementary
band retains its alternating phase. Dropping C^j would create order-one
error at odd ticks even at k=0. Grouping two ticks cannot justify treating
the pi sector as an additional zero-energy field or deleting its original
quasienergy. This estimate uses one common native time for every channel.

### A broader finite-stencil obstruction, without a flight decomposition

Consider any translation-invariant finite-range one-particle unitary

    U(k)=sum_{b in S} exp(-i k.b) K_b,                     (6e)

on a lattice, with finite internal dimension. It need not factor into our
flights. Its charged period is T, and exact native speed c requires
|b|<=cT for every nonzero coefficient. At any momentum, put
G(n)=i U(k)* partial_n U(k), l=min_b n.b and h=max_b n.b. Then

    l I <= G(n) <= h I.                                  (6f)

Here is an elementary proof. First take a lattice-rational direction and
rescale so every n.b is an integer. After shifting by l, the matrix-valued
function V(z)=sum_b z^(n.b-l) exp(-i k.b)K_b is a polynomial unitary on
the unit circle. The maximum principle, applied to every scalar matrix
element, gives ||V(z)||<=1 inside the disk. For any vector psi, differentiating
||V(r)psi||^2<=||psi||^2 from below at r=1 gives
<psi,(G(n)-lI)psi>>=0. The boundary unitary identity makes this matrix
Hermitian. Applying the same argument to the reversed polynomial gives
hI-G(n)>=0. Density of rational directions proves (6f) for every n.
This is the finite-polynomial causal-velocity bound; it does not assume
that the coefficients K_b are positive, mutually orthogonal or probabilities.

Compressing (6f) to a two-state degenerate band gives the exact necessary
convex-geometric condition

    u + M closed_unit_ball subset conv{b/T : b in S}.      (6g)

Indeed the expectation of A(n) in Bloch state r is n.(u+M r); all its
support-plane inequalities are (6f). In particular an untilted isotropic
speed v requires the ball of radius |v|T to lie inside conv S. A finite
polytope contained in the cT ball cannot contain that ball when cT>0.
Thus no such finite stencil, even outside the flight grammar, attains
an isotropic field speed c at a finite level.

There is also a uniform resource bound. If N is the number of nonzero
displacement terms and 0<|v|<=c, the caps
{n in S^2 : n.b>=|v|T} must cover the unit sphere. Each has area at most
2pi(1-|v|/c), since |b|<=cT. Subadditivity therefore gives

    N (1-|v|/c) >= 2.                                    (6h)

A fixed N excludes even a refinement limit with |v|/c->1. Reaching
|v|/c>=1-eta requires N>=2/eta. Zero displacement contributes no cap and
cannot be counted as an extra direction. Growing internal dimension alone
does not change the bound. Nor is exact finite-level isotropy needed for
the limiting exclusion: (6g) contains a centered ball of radius
s_min(M)-|u| whenever this quantity is strictly positive, where s_min is
the least singular value. Substitute that positive radius for |v| in (6h).
If it is nonpositive, this argument supplies no positive-radius cap bound;
the origin need not even belong to the translated ellipsoid. For example
U(k)=exp(-ik_x)I_2 with T=c=1 has N=1, u=(1,0,0) and M=0.
Its velocity set is the singleton {u}, so substituting a clamped zero radius
would incorrectly demand 1>=2. The stationary identity has N=0 and likewise
cannot be subjected to the positive-radius bound. Both are legitimate
unitaries, retained as exact scope controls. In an approximately luminal
sequence M->cR and u->0, the radius is eventually positive and tends to c,
so the conclusion N->infinity is unchanged.
For example an eight-corner cubic support has
inradius a and circumradius sqrt(3)a, so even an endpoint-only timing
assignment to (7) gives |v|/c<=1/sqrt(3); its actual intermediate flights
give the stronger 1/3 bound. The regular tetrahedral support has inradius
a/sqrt(3), so (6c) saturates both its geometric and its flight bound.

These are exclusions, not promises that a large stencil realizes a desired
quantum field. A growing direction set is necessary for a luminal isotropic
finite-range refinement, whereas every positively timed register-flight
lowering remains blocked by (6), even with infinitely growing directions.
Finite-stencil dynamics, finite autonomous Hamiltonians and the proposed
classical-flight lift therefore cannot be interchanged to hide the clock
problem. The supplied continuum coherent witness in Section 6 avoids both
finite-cutoff hypotheses; its source admission is a separate physical question.

## 3. An explicit unitary process attaining the bound

On the unit three-torus take q divisible by four, a=1/q and two fermionic
modes at every site. A period has three conditional flights of length a,
first along x with sigma_x projectors, then y with sigma_y projectors, then
z with sigma_z projectors. Each lasts a/c. Its total native duration is

    T_a=3a/c,    U(theta)=exp(-i theta_z sigma_z)
                         exp(-i theta_y sigma_y)
                         exp(-i theta_x sigma_x).           (7)

The physical momentum is k=theta/a. The factors are exact translations and
onsite basis changes. They are not Euler approximations to a Hamiltonian.
At every intermediate time modes travel at speed c; basis changes occur
onsite. For bounded k,

    U(a k)=I-i a k.sigma+O(a^2 |k|^2),
    U(a k)^floor(t/T_a) -> exp(-i t(c/3) k.sigma).           (8)

The product Taylor remainder is at most (a^2/2)(sum |k_i|)^2, and the
single exponential remainder is at most a^2 |k|^2/2. Thus the one-period
operator error is at most 2 a^2 K^2. Telescoping unitary products and charging
the uncompleted final period gives the explicit uniform estimate

    error <= (2 c T/3) a K^2 + a K,  |k|<=K, 0<=t<=T.      (8a)

Every axis step is charged. Assigning duration a/c to
the three-flight word would artificially multiply its field speed by three.
Net displacement, a simultaneous replay of dependent stages, or a reset
of the clock at an intermediate writer cannot justify that reassignment.

Positive onsite gate overhead only lowers the effective speed. If its total
per-period duration is eta_a T_a with eta_a->0, whole-horizon overhead is
O(eta_a T), so (8) is retained. A uniform bound on onsite generator strength
does not in general permit such a refinement for fixed-angle basis changes;
their local control resource is not free. The ideal word and its flight
budget are the object of (5)--(8), not a bounded-hardware implementation claim.

## 4. The whole spectrum, without selecting one good Weyl point

The complete 2 by 2 matrix has eigenvalues exp(+/-i epsilon), 0<=epsilon<=pi,
where direct multiplication gives

    cos epsilon = cos x cos y cos z + sin x sin y sin z.     (9)

All eight zero-quasienergy points on the momentum torus are retained:

* x,y,z in {0,pi}, with an even number of pi entries (four points);
* x,y,z in {pi/2,3pi/2}, with an even number of 3pi/2 entries (four).

There are also eight pi-quasienergy points with the opposite parity. To prove
completeness, Cauchy--Schwarz bounds the absolute right side of (9) by
sqrt(cos^2 x cos^2 y+sin^2 x sin^2 y)<=1. Equality requires both
cos x sin y=0 and sin x cos y=0, fixing these possibilities, and z must
align with the remaining unit two-vector. No numerical search establishes
this classification.

Near every zero, U is I-i sum_i (D k)_i sigma_i+O(|k|^2) with D orthogonal.
The first four points have positive determinant and the other four negative
determinant. Every cone has native speed c/3. Thus the full massless content
is eight Weyl cones, four of each chirality, rather than one chosen species.
The exact finite derivative matrices and chirality signs are independently
checked. The pi cones have energy of order 1/a in the declared observable.

Define the reference before making a comparison. On each finite torus let
h_a be the principal Hermitian logarithm of the actual U divided by T_a.
At eigenvalue -1 use a fixed branch convention. Fill all strictly negative
levels and leave zero levels empty. This defines a stationary Slater vacuum.
Normal order dGamma(h_a) relative to it. Its excitation energies are
E_a(theta)=epsilon(theta)/T_a, with two excitation choices per momentum,
including the declared zero-level degeneracy. The finite vacuum is preserved
exactly by the Fock lift of (7); no comparison substitutes a modified state.
This is a declared Floquet quasienergy observable. Its logarithm is generally
nonlocal and is not asserted to be an autonomous native Hamiltonian or a
derived OPH thermodynamic energy calibration.

For every beta>0 the entire finite partition function is

    log Z_a = 2 sum_theta log(1+exp(-beta E_a(theta))).      (10)

For q->infinity through multiples of four,

    log Z_a -> 16 sum_{n in Z^3} log(1+exp(-beta (2pi c/3)|n|)). (11)

The corresponding normal-ordered energy and Gibbs entropy converge too.
In fact the following global bound supplies an explicit summable envelope:

    epsilon(theta) >= dist(theta,Z)/(8 sqrt(3)).            (11a)

For epsilon>=pi/8 this follows from torus diameter sqrt(3) pi. Otherwise
Cauchy--Schwarz in (9) gives cos(2x) cos(2y)>=cos(2epsilon)>0. Both factors
have the same sign and absolute value at least cos(2epsilon), so x,y are
within epsilon of the same family, either multiples of pi or pi/2 modulo pi.
Write R exp(i phi)=cos x cos y+i sin x sin y. Then
cos epsilon=R cos(z-phi), R<=1, so |z-phi|<=epsilon modulo 2pi. The angle
phi is within arctan(tan^2 epsilon)<=epsilon of the corresponding z coordinate
of a genuine zero. Hence dist(theta,Z)<=sqrt(6) epsilon in this case.
Assign each lattice momentum to a nearest zero, breaking ties deterministically.
The zeros lie on the momentum grid when q is divisible by four, so
E_a>=pi c |n|/(12 sqrt(3)) for its integer offset. Each zero supplies at most one point
per offset. Exponential summability dominates the full catalog, while the
fixed-offset limit is (2pi c/3)|n| at each of the eight zeros. The functions
E/(exp(beta E)+1) and log(1+exp(-beta E))+beta E/(exp(beta E)+1) obey the
same domination. In particular there is no hidden extensive zero-energy
sector, and no mode is discarded to obtain (11).

Fixed-band dynamics about all eight points obey (8) with their actual D.
On at most K fermionic excitations, exterior-product telescoping multiplies
the one-particle operator error by at most K. Finite smeared CAR covariance
comparisons follow by spectral projector convergence away from each zero;
the explicitly specified zero-mode states are compared separately. This
does not claim global norm convergence between inequivalent infinite seas,
an interacting theory, or the scalar vacuum theorem of `code/m1_vacuum_fidelity`.

## 5. Exact delay and fixed finite Hamiltonian sites

Let H be Hermitian on a finite tensor product of full matrix site factors,
with distinct fixed site positions. Suppose every pair of local observables
A_x,B_y at distinct sites satisfies

    [exp(i t H) A_x exp(-i t H), B_y]=0
                       whenever 0<t<|x-y|/c.               (12)

Each matrix entry is analytic in t. Vanishing on a nonempty interval makes
it identically zero. Consequently the evolved algebra of x commutes with
every other full site factor for all t, and therefore stays in the factor
at x. Differentiation shows [H,A_x] belongs to that factor for every A_x.
Every derivation of a full matrix algebra is inner. Subtract its onsite
implementer at each site; the remainder commutes with the full tensor
algebra and is scalar. Thus H is a sum of onsite terms and a scalar.

Equivalently expand H in tensor products of the identity and traceless local
matrix bases. A coefficient with support on two or more sites produces a
nonzero double commutator with suitable local matrix units and is forbidden
by the derivative of (12). This also proves the conclusion directly.

For two qubits this obstruction has a quantitative, basis-complete witness.
If the interaction part is sum_{a,b=1}^3 h_ab sigma_a tensor sigma_b, then

    sum_{i,j=1}^3 ||[[H,sigma_i tensor I],I tensor sigma_j]||_HS^2
        = 256 sum_{a,b=1}^3 h_ab^2.                       (12a)

Onsite and scalar terms contribute zero. Pauli orthogonality eliminates cross
terms, and each interaction coefficient contributes four terms of size 64.
For H=sigma_x tensor sigma_x, prepare the first qubit along +y and the second
along +/-x. The first-qubit effect (I+sigma_z)/2 has probabilities
(1+/-sin(2t))/2. Thus an actual remote preparation, not just a matrix norm,
gives a nonzero response for arbitrarily small positive t.

The claim is about all local preparations/effects and exact zero influence,
not a thresholded numerical response. Finite oscillator quadratic systems
have the same conclusion at the canonical level: their finite-dimensional
linear propagator is an analytic matrix exponential. This is the standard
strict-causality versus Hamiltonian distinction, here applied to the complete
native-clock requirement; it is not attributed as a new general theorem.

Moving modes in (1) evade the fixed-site premise: their local algebra ownership
changes along their trajectories. A sampled quantum cellular automaton has
no continuous intermediate fixed-site dynamics satisfying (12). Continuum
local fields use infinite-dimensional local algebras and singular generators.
None of these cases is a counterexample to the finite-site statement.

## 6. A precise coherent escape, and its source boundary

The missing operation is more specific than "some quantum dynamics." The
positive classical-position velocity resolution (4) is incompatible with
three luminal Pauli components. A continuum operator can have those components
because the velocity operators do not commute.

Let E act on a fundamental two-component SU(2) module. A real-linear
Hermitian velocity map A:E->Herm(2) covariant under this SU(2) has the form
A(v)=gamma v.sigma. To see this, Herm(2)=R I plus the adjoint three-vector;
there is no invariant linear functional E->R, and an intertwiner between
the two irreducible real three-vector representations is a scalar. Thus,
after choosing handedness and time scale, the derivative part is fixed.
A covariant onsite term can add mu I; covariance does not determine that
energy offset. The witness chooses zero onsite offset and the Weyl equation

    i partial_t psi = -i c sigma.grad psi.                 (13)

Its squared generator is -c^2 Delta. Its density and current obey
partial_t rho+div j=0, rho=psi*psi, j=c psi*sigma psi,
and |j|<=c rho. Integrating this identity outside an expanding ball proves
exact propagation inside radius ct for compactly supported initial data.
Fourier modes have dispersion +/-c|k|. The field dynamics and its causal
metric therefore use the same c without a rescaled algorithmic tick.
The standard free CAR quantization supplies spacelike locality of the even
observable algebra and a positive excitation-energy Fock representation.
This continuum witness is standard free-field mathematics, not a new field
equation or a finite-source implementation.

For a four-component parity-paired module, an onsite mass matrix beta
anticommuting with all three alpha_i gives the massive Dirac equation and
the same current bound. By contrast no positive velocity resolution (4)
can represent c alpha_i: the trace witness gives 3c<=c. This states exactly
what a direct coherent transport primitive would have to supply beyond the
register-flight grammar.

An OPH source realization still must justify the fundamental module's
operational meaning, this coherent spatial generator, its physical time and
the refinement of source-local instruments. Internal SU(2) alone does not
identify a gauge response with spin or motion. Equation (13) is a concrete
candidate transport law, not an adopted axiom or an A1--A3 derivation. The
finite-site rigidity result forbids claiming exact cutoff-level (12) from
ordinary interacting finite Hamiltonians. A continuum-limit causal law and
a finite exact native-flight law are different physical proposals.

## 7. Relation to the preceding results and external mathematics

The sparse read-family comparison (`code/m1_necessity`) separates continuum
requirements from M1's finite menu. The interface repair (`code/m1_interfaces`)
removes microscopic partitions and unwanted scalar spectral modes. The vacuum
fidelity theorem (`code/m1_vacuum_fidelity`) controls the original scalar vacuum
but excludes its numerical clock. The present result studies
a different, exactly unitary finite-dimensional operation family, charges
its native intermediate motion, characterizes its entire Weyl speed budget
and retains all of its low-energy species. It does not replace their scalar
reference action or claim that a good low-energy mode proves native M1.

Standard quantum-walk Dirac limits are described by Mlodinow and Brun,
Physical Review A 97, 042131 (2018), arXiv:1802.03910, and by Bisio,
D'Ariano and Perinotti, arXiv:1608.02004. Their continuum constructions do
not by themselves charge a separately proposed native intermediate clock.
The tetrahedral qubit POVM is a standard SIC measurement; see Renes,
Blume-Kohout, Scott and Caves, Journal of Mathematical Physics 45, 2171
(2004), https://arxiv.org/abs/quant-ph/0310075. Tetrahedral Dirac walks have
also been studied by Nzongani and collaborators,
https://arxiv.org/abs/2404.09840. Neither the tetrahedron nor the existence
of a relativistic quantum-walk limit is claimed as a new general result.
The Hamiltonian/strict-causal distinction is discussed by Ranard, Walter
and Witteveen, Annales Henri Poincare (2022),
https://doi.org/10.1007/s00023-022-01193-x.
Positive operator-valued compressions, the nuclear norm, exterior powers,
Floquet logarithms and the free Weyl/Dirac current are standard inputs.
The contribution here is the complete native-time velocity budget, its
derived optimal geometry and explicit minimal realization, the broader
finite-stencil clock exclusion, full spectral and thermal accounting, and the precise
comparison of the resulting quantum transport with OPH's complete clock law.
