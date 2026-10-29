# Continuous equivariant degree: independent proof candidate

Date: 2026-10-29. Source base: `f94cde24b8b3486ad6e1e823498bb29587f026bb`.

Acceptance status: **FROZEN_PAIR_AND_ALGEBRAIC_CROSSCHECK_ONLY**. The original
Bernhard Mueller message, its §III, and the four scripts at their advertised
commit were not located. No execution or independent review of those scripts
is claimed. The argument below is an independently written proof candidate.
The separate numerical/finite verifier checks arithmetic and representation
typing; it does **not** independently audit the continuous regularization
argument. In particular, this branch does not emit
`CONTINUOUS_TWIN_DEGREE_FLOOR_PROVED`. The remaining acceptance obligation is
hostile mathematical review of Lemma 1, especially its relative approximation
and free-part transversality step, and the normalizer comparison in Lemma 2.
No failing arithmetic or counterexample to the proposed theorem was found.

## Objects and statement to review

Let K be the geometric realization of the twenty **ordered** faces parsed
from `CoreAxioms.orientedFaces`. Equip K with their integral fundamental cycle.
The convex plus radial realization identifies K with a sphere by an oriented
homeomorphism. This transports its sixty proper incidence permutations to a
faithful rotation action rho on an oriented S². The twin target has the
rotation action rho' reconstructed from the Galois-conjugate coordinates.
Target orientations use the *same* committed convention, ambient orientation
minus the standard xyz orientation; the frozen radial degrees are +1 and +7.

The claim is about **every continuous** f with
`f(rho(g)x) = rho'(g)f(x)` for every g in this *fixed marked* A5. No outer
automorphism is inserted in this equation. The untwisted comparison replaces
rho' by rho. No three-dimensional linear-map, quantum-channel, or observer
typing is implied.

Candidate conclusion: twin degrees lie in {7,13,17,23} modulo 30, so have
absolute value at least seven. Untwisted degrees necessarily lie in
{1,11,19,29} modulo 30, and degree one exists. These are **necessary residue
conditions**, not a claimed sufficient classification of all degrees.

## Lemma 1: cyclic local-index congruence without a derivative hypothesis

Let C_p, p prime, act faithfully by proper rotations on two oriented spheres.
Write x+,x- and y+,y- for their respective fixed poles. Suppose a continuous
equivariant f sends x+ to y+ and x- to y-. In oriented tangent charts at
x+ and y+, the generator acts by multiplication by exp(2πia/p) and
exp(2πib/p), with a,b nonzero modulo p. Then

    deg(f) = b a^(-1) (mod p).

Here is a relative regularization argument, including the critical-pole case.

1. Choose small disjoint invariant closed disks at the two domain poles.
   Continuity makes their images lie in small disjoint invariant geodesic
   disks at the corresponding target poles. Use oriented exponential charts
   there. The target action is exactly the linear tangent rotation in these
   charts. All convex interpolations below stay inside those target disks.

2. Choose the integer k in {1,...,p-1} with ak=b modulo p. On a smaller
   disk about x+, replace the chart expression of f by ε z^k. Interpolate
   to the old expression with a radial cutoff on the larger disk. Both
   expressions are equivariant, so the interpolation is an equivariant
   homotopy. Take ε positive and sufficiently small. The new map has an
   isolated zero at x+ of local index k; its derivative may still vanish.
   About x-, a similar cutoff makes the map constant at y-, which excludes
   y+ from that disk. Nothing is assumed about the original fibre at x+.

3. Keep a still smaller neighborhood of each pole fixed. Approximate the
   resulting continuous map smoothly and equivariantly relative to those
   neighborhoods. One construction views the map in R³, first uses smooth
   approximation relative to the closed smaller disks, and then averages
   over C_p:

       Fbar(x) = (1/p) sum_g rho_target(g)^(-1) F(rho_domain(g)x).

   The map being approximated is equivariant and has norm one. Uniform
   approximation within δ<1 makes Fbar nonzero everywhere. Normalize Fbar;
   normalization commutes with rotations. The normalized straight-line
   homotopy remains nonzero and preserves the protected disks. Relative
   smooth approximation here follows by ordinary Euclidean local smoothing
   and a partition of unity with cutoffs supported outside the protected
   disks; the map was already smooth on larger neighborhoods of those disks.

4. The complement of the protected disks is a compact surface with boundary
   on which C_p acts **freely**. Its boundary avoids y+. Make the smooth map
   transverse to y+ on this complement by an arbitrarily small equivariant
   perturbation supported away from that boundary. Explicitly, each point
   has a neighborhood with disjoint translates. Pick finitely many such
   neighborhoods covering the compact inverse image of a small closed
   target disk around y+. On each, vary the map by two independent tangent
   parameters with a bump positive on a smaller covering set, and extend
   the variation to its translates equivariantly. At every point mapping
   to y+, some such pair of parameters spans the target tangent plane.
   The evaluation family is therefore transverse to y+. Parametric
   transversality (equivalently Sard applied to this finite parameter family)
   gives a small parameter choice making y+ regular on the complement.
   Perturbations are small enough that no other points enter the fibre.

5. The new fibre consists of x+, of local index k, and finitely many free
   orbits of size p. Proper rotations on source and target give equal
   local indices at all members of each orbit. Each such orbit contributes
   p times a signed integer. Excision/additivity of the degree gives
   `deg(f)=k+p m`. Every alteration was a homotopy, so the original continuous
   map has that degree too. This proves the candidate lemma.

If x+ maps to the other target pole, choose that pole as y+; its tangent
weight changes sign. If the two domain poles map to the *same* target pole,
the lemma's hypothesis fails. A constant map to a fixed pole is an immediate
C_p-equivariant degree-zero control. Thus pole separation is essential.

## Lemma 2: A5 forces separation of the two poles

For a subgroup of order 3 or 5, its A5 normalizer is dihedral of order 6 or
10. An element outside the cyclic subgroup conjugates a generator to its
inverse. In a faithful proper rotation representation, such an element
interchanges its two fixed poles: fixing either would force its rotation
axis to coincide with that of the cyclic subgroup, hence commute with its
generator, a contradiction for p odd. Both rho and rho' are faithful, so
the same source normalizer element interchanges both pairs of poles.
Equivariance therefore gives f(x-)=the pole opposite f(x+).

For p=2 the normalizer is V4. Its three nonidentity elements are distinct
half-turns. Two commuting distinct half-turns in SO(3) have perpendicular
axes; consequently either other half-turn exchanges the chosen axis's
poles. This holds in both faithful target actions as well. Thus separation
also holds at order two. The p=2 conclusion does not follow just from the
formal equality g=g^(-1).

## Consequences, assuming these lemmas survive independent review

The chosen order-five source element is

    (0,2,4,1,6,8,3,5,10,7,9,11),

fixing ports 0 and 11. Its source trace is phi and its twin trace is 1-phi.
For source weight a=1 at a chosen pole, `tr=1+2 cos(2π/5)=phi`;
the twin has weight b=±2 since `1+2 cos(4π/5)=1-phi`. Thus Lemma 1 gives
degree ±2 modulo 5. The sign depends on which target pole receives the
chosen source pole. It is one congruence for the **global** degree, not an
independent sign per preimage or face.

At order three, source and target weights are units modulo 3, so degree is
±1 modulo 3. At order two both weights are 1, so degree is odd. Combining
these conditions gives exactly 7,13,17,23 modulo 30; their least absolute
integer representatives have absolute value seven. The existing exact
regular-value replay of radial_minus supplies an attaining +7 map.
The untwisted weights give ±1 modulo 5, hence residues 1,11,19,29; radial_plus
supplies an attaining +1 map.

A continuous injective map S²→S² is a homeomorphism, by invariance of domain
and compactness, and has degree ±1. Under the reviewed lemmas there can be
no injective twin-equivariant map in the specified category.

Coherent proper relabeling conjugates the action without changing degree.
Squaring the **same** marked generator on both sides multiplies a and b
by two and preserves b/a modulo 5. Squaring it on the target alone changes
the equivariance equation. An outer A5 automorphism exchanges the marked
order-five classes and changes the representation typing; it is not a
fixed-incidence recharting. Reversing the common target orientation changes
both radial degrees to -1 and -7 and preserves both absolute floors. A
branch-dependent orientation convention is not an allowed comparison.

## Independent algebraic cross-control and lift boundary

`equivariant_degree_crosscheck.py` constructs all nine irreducible
characters of 2I over Q(phi). Opposite lifts, including the central -1,
remain separate classes. Q is the defining SU(2) representation; Q' has
character sigma(chi_Q), with sigma(phi)=1-phi. Its symmetric square has the
twin SO(3) character. Character inner products give the McKay graph and the
two multiplicity sequences through degree 94. The independent verifier
starts with the affine E8 graph and uses `v_n=A v_(n-1)-v_(n-2)` instead.

The unique Q' in Sym^7(Q*) gives a nonzero covariant homogeneous polynomial
pair P of degree seven. Its common projective zero set, if nonempty, is
A5-invariant. A point stabilizer in the rotation action is cyclic of order
at most five, so an orbit has at least twelve points. Two degree-seven
polynomials cannot have twelve distinct common projective zeros. Hence
P has no base point and defines a degree-seven rational map. This existence
argument does not require matching HMS coordinates to the OPH labels.

Conversely a rational map with projective equivariance and a primitive
homogeneous pair of degree d satisfies `P(gv)=c(g) Q'(g)P(v)`. The multiplier
is constant: two primitive homogeneous pairs of the same degree defining
the same projective map differ by a nonzero constant, by coprimality in
C[x,y]. Composition makes c a character of 2I. For completeness its
abelianization vanishes: in the binary (2,3,5) presentation
`x²=y³=z⁵=xyz`, the abelian relation matrix

    [ 1 -1 -1 ]
    [-1  2 -1 ]
    [-1 -1  4 ]

has determinant -1. Thus c=1. The central -1 then forces d odd. This
justifies the homogeneous lift for this **rational** category; it says
nothing by itself about an arbitrary continuous map.

The HMS (3.18) pair has gcd one over Q and no common zero at infinity.
The producer uses rational Euclidean division; the independent verifier
uses coprimality modulo 101 with nonvanishing leading coefficients. This
checks base-point freeness of that explicit formula. Source-typed full
equivariance of that coordinate formula is not asserted by the gcd test.

## Checked literature and remaining custody

[Houghton–Manton–Sutcliffe, hep-th/9705151v2](https://arxiv.org/html/hep-th/9705151)
gives degree-seven formulas in (3.18), the family in (3.19), and the
two-dimensional target representation in (4.6). These are rational-map
statements, not a proof of Lemma 1 for continuous maps.

[Kostant, math/0411142v1, Example 1.7](https://arxiv.org/pdf/math/0411142)
lists the near/far two-dimensional E8 numerators with denominator degrees
12 and 20. The checked numerators are {1,11,19,29} and {7,13,17,23}.

[Doyle–McMullen, §5, Proposition 5.3](https://math.dartmouth.edu/~doyle/docs/icos/icos/node5.html)
classifies the four commuting rational maps of degree below 31, with degrees
1,11,19,29. It does not classify twin-equivariant continuous maps. The
category difference prevents importing it as the universal degree theorem.

Searches covered fetched OPH refs and paths, GitHub code/issues/commits in
FloatingPragma/observer-patch-holography and muellerberndt/oph-physics-sim
using Kostant, Galois, equivariant, degree and degree-seven. No original
§III or four-script commit was identified. Indexed simulator results pointed
to `90095f32e7b88dd120158130802cf9caccf7f3d5`. The complete recursive
simulator tree at that commit was then inspected (2,474 entries, not
truncated). Its `port_gram_equivariant_action.py` and independent verifier
concern finite twelve-port recharting of a conditional metric completion;
they are not the announced continuous-degree packet. No branch named for
degree was returned by branch search. This bounds the search and does not
exclude unpublished scripts or scripts with unrelated names.

Even if the continuous theorem is accepted, applying it to the A1 support
requires a source-selected relative map, the same oriented realized target
sphere, compatible refinement maps, and a proved equivariance equation.
The seed/degree receipts do not supply these premises. Any degree-one twin
exclusion is conditional on them; physical identification remains open.
