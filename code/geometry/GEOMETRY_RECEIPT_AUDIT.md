# Finite geometry: support, refinement and oriented caps

The finite geometry readout needs the full record-support complex, a
certificate for the actual refinement, and a side choice for each oriented
cap. None can be replaced by a lower-dimensional projection or a successful
fit. This audit corrects those substitutions in the executable evidence for
Theorem 4.3c in `paper/tex_fragments/PAPER.tex`.

The repair fixture admits one seeded nonzero bit and records of size at least
two. Its first effective update erases that bit: a larger record has a zero
majority, and a two-member record uses the disclosed zero-valued tie break.
Other updates change nothing, and a complete schedule necessarily visits a
record containing the defect. Every schedule therefore reaches the all-zero
state in one effective update. This proves the fixture's confluence on its
whole admitted family. Incomplete schedules and states with several defects
are rejected; a random schedule battery is unnecessary for this proof. The
scalar KMS comparison also validates its positive finite inputs and bounded
tolerance. Matching a supplied target does not select that target or its
physical clock normalization.

The findings are reproduced against main commit
`439a406e97de54d3ce12f99e9adafe0938dee051`. The theorem already distinguishes
simplicial refinement from topology-preserving subdivision. The corrections
bring its finite helpers into agreement with that distinction; the paper's
axioms and physical claims are unchanged.

## Reproduced failures

| Input to the audited implementation | Observed result | Correct interpretation |
| --- | --- | --- |
| One joint record `(0,1,2,3)` | `S2`, Euler characteristic 2 | The full support is a 3-simplex, a ball, with Euler characteristic 1. Its triangular boundary is a different complex. |
| One edge mapped identically to two vertices with no edge | Simplicial refinement accepted | The image of the edge is absent. Checking triangles alone is insufficient. Isolated vertices and higher simplices also require checks. |
| Four exact equatorial points | Cap normal approximately `(707114.6,-500005.5,-500005.5,0)` | A great circle has normals `(0,0,0,1)` and its negative. An inhomogeneous equation with right-hand side 1 cannot represent a plane through the origin. |
| Circle at angular radius 2.4, north pole intended inside | Normal approximately `(1.091686,0,0,-1.480466)` | The intended normal is its negative. The circle alone does not select its inside. |
| Function named `reconstruct_from_cross_ratios` | Consumed point coordinates and computed its own ratios | Fixture production and reconstruction from supplied scalar data are different operations. |

Malformed face data are rejected before a topological verdict. Empty
complexes and isolated vertices are valid finite complexes, but fail the
closed-surface predicate. Numerical cap inputs require finite unit vectors,
a resolved plane and spacelike normalization, an explicit side witness, and
a declared residual budget. Radius clipping and silent sign selection are
absent.

The maintainer-style audit of commit `2dff2c9c` also reproduced a conversion
failure: on platforms with extended-exponent floating types, the finite
coordinate `1e400` became infinity when cast to binary64 and satisfied the
infinity gauge check. A nonzero `1e-400` could similarly become a zero gauge
entry. Validation checks the original components and rejects conversion
overflow or erased nonzero components. Mixed Boolean coordinate lists are
rejected before array coercion. The wider-exponent regression executes on
Linux; Windows skips it when its `longdouble` has the binary64 range.

Record supports are validated whenever repair consumes them, including after
mutation. Support validation does not construct an exponentially larger face
closure merely to check labels. A caller-supplied initial-state argument is
rejected instead of silently overwritten by the fixture seed. These input
corrections preserve the theorem premises and the retained cyclic-tower
receipt, which reproduces exactly without artifact changes.
The defect count uses Python integers: a sum of 128 NumPy `int8` ones wraps
to -128 in fixed-width arithmetic and must not pass the one-defect gate.
Signed and unsigned overflow controls, and fresh subdivision labels beyond
the NumPy integer range, exercise the exact combinatorial arithmetic.

## The support-to-incidence theorem

For a finite family of nonempty record supports R, define

\[
K(R)=\{s\ne\varnothing:\ s\subseteq r\text{ for some }r\in R\}.
\]

This is the smallest abstract simplicial complex containing every record
support as a simplex. Closure under nonempty subsets is immediate; every
complex containing R must contain all those subsets. Consequently

\[
\dim K(R)=\max_{r\in R}|r|-1.
\]

For the empty family use dimension -1. Relabeling commutes with this
construction, and adding a record already contained in another changes
nothing. Every finite abstract complex has such a presentation: use its
maximal faces as records. Hence record-support incidence by itself imposes
no special dimension or topology. This statement concerns the finite
representation; it does not claim that every presentation is physically
admissible under A1-A3.

The implementation stores vertices, edges, triangles and all higher faces.
The Euler characteristic includes every dimension. For a support of size m,

\[
\chi(K)=\sum_{j=1}^m(-1)^{j-1}\binom mj=1.
\]

In particular, four separate triangular records describing the boundary of
a tetrahedron give a sphere. A single four-patch record includes the interior
3-simplex. Both have exactly the same vertices, edges and triangles, so no
test confined to the two-skeleton can distinguish them. Independent binary
boundary-matrix elimination verifies Betti numbers `(1,0,1)` for the boundary
and `(1,0,0,0)` for the solid simplex.

The surface predicate requires no higher simplex, two incident triangles per
edge, and a nonempty single-cycle link at every vertex. Connectivity and
Euler characteristic 2 then give the sphere by surface classification.
Coherent orientation is retained as an explicit output and checked as in the
existing six topology controls. Simplicial complexes and their homology use
the standard conventions of [Hatcher, Algebraic Topology, section 2.1](https://pi.math.cornell.edu/~hatcher/AT/AT.pdf).

## Simplicial maps and actual subdivisions

A vertex function is simplicial exactly when it is total and the image of
every nonempty simplex is a target simplex. Vertex identifications are
allowed. `refinement_is_simplicial` checks that definition in every stored
dimension. The identity vertex map from a solid tetrahedron to its boundary
fails even though all triangular faces have valid images.

A constant vertex function from a sphere to a point is simplicial. It also
defines a deterministic stochastic pushforward on classical vertex states:
its matrix has nonnegative entries and unit column sums. Such a channel has
an extension to a completely positive trace-preserving measure-and-prepare
map. Nevertheless it annihilates the sphere's second homology. Simplicial
validity and this classical channel property certify no topology
preservation. The simplicial helper therefore carries no Petz/CPTP label.

**Subdivision certificate.** For a finite complex of dimension at most two,
supply one distinct fresh midpoint label m_e per edge e. Replace each edge
`{a,b}` by `{a,m_ab}` and `{m_ab,b}`. Replace each triangle `{a,b,c}` by

```text
{a,m_ab,m_ac}, {b,m_ab,m_bc}, {c,m_ac,m_bc}, {m_ab,m_ac,m_bc}.
```

Keep all isolated vertices. The verifier requires exact equality of the
entire resulting complex with the supplied fine complex. Missing triangles,
extra isolated vertices, repeated midpoint labels and omitted edges fail.
It reads the supplied coarse/fine complexes and midpoint assignment, rather
than trusting a producer's Boolean verdict.

**Proof of topology preservation.** Realize an abstract coarse simplex in
barycentric coordinates. Send old vertices to themselves and m_ab to
`(a+b)/2`, extending affinely on each fine simplex. On a coarse edge the two
segments partition the edge. On a triangle the four triangles have disjoint
interiors and cover it; their oriented area determinants are all 1/4.
The maps agree on shared faces. Their union is a continuous bijection from
the compact fine realization to the Hausdorff coarse realization, hence a
homeomorphism. This proves the certificate for every finite complex in its
stated dimension, including boundary and disconnected components.

The endpoint-valued projection returned by `refinement_subdivide` is a
different map. It collapses some fine edges and is not a homeomorphism. Its
linear realization and the midpoint homeomorphism send each fine simplex
into the same coarse simplex; straight-line interpolation within that
simplex gives a compatible homotopy. Thus the endpoint projection is a
homotopy equivalence. An independent integral fundamental-chain calculation
checks unit degree on the tetrahedral sphere, in addition to the binary
homology check. The producer certifies its full subdivision before returning
the endpoint projection.

## Cross ratios are reconstruction data

Represent a complex coordinate z by `(z,1)` and infinity by `(1,0)`.
For homogeneous vectors u and v let `[u,v]` be their determinant. Then

\[
\operatorname{CR}(u_1,u_2;u_3,u_4)
=\frac{[u_1,u_3][u_2,u_4]}{[u_1,u_4][u_2,u_3]}.
\]

This formula handles exact infinity without subtracting infinities. An
invertible two-by-two matrix multiplies each bracket by its determinant,
which cancels from the ratio. Fix three distinct labelled gauge points
g1,g2,g3. The coordinates

\[
w_i=\operatorname{CR}(z_i,g2;g1,g3)
\]

send the gauge to 0,1,infinity. The corresponding fractional linear map is
invertible. Two distinct labelled configurations with identical w_i are
therefore related by one common Mobius transformation, and the converse holds.
This is a complete finite reconstruction up to that freedom, including
points at either stereographic pole.

`cross_ratio_receipts` is the coordinate-based fixture producer.
`reconstruct_from_cross_ratios` accepts only the scalar receipt vector and
its gauge labels. It validates the exact gauge entries, distinctness,
finite coordinates or canonical infinity, and numerical resolution. A test
disables both source-coordinate helpers while the consumer runs. These
scalar fixtures do not become measured modular data merely by receiving
the name "receipt".

The complete `cap_from_cross_ratios` consumer takes the same scalar vector,
boundary-point indices and a distinct inside-point index. In the chosen
gauge it maps a homogeneous coordinate `(u,v)` to
`(2 Re(u conjugate(v)), 2 Im(u conjugate(v)), |u|^2-|v|^2)` divided by
`|u|^2+|v|^2`, then applies the oriented-plane reconstruction below. It never
reads source point coordinates. The end-to-end tests disable both coordinate
receipt production and stereographic projection while this consumer runs.

## A circle, its two caps and the Lorentz normal

Let distinct unit vectors p_i lie on a nondegenerate circle. Introduce the
rows `(-1,p_i)` of a matrix A. A plane vector u=(t,v) satisfies

\[
Au=0,\qquad v\cdot p_i=t.
\]

**Theorem.** If A has rank three and its one-dimensional nullspace is
spacelike, `||v||^2-t^2>0`, the boundary determines exactly the two unit
Lorentz normals

\[
n=\pm\frac{(t,v)}{\sqrt{\|v\|^2-t^2}}.
\]

Rank three gives uniqueness up to real scale. Lorentz normalization fixes
its magnitude. Writing c=v/||v|| and
`cos(alpha)=t/||v||`, spacelikeness gives `0<alpha<pi` and

\[
n=(\cot\alpha,\csc\alpha\,c).
\]

The cap interior is `n_spatial dot p - n_time > 0`. A supplied unit witness
off the boundary selects exactly one sign. Flipping the sign exchanges
`(c,alpha)` with `(-c,pi-alpha)` and selects the complementary cap. Great
circles have t=0 and fit the same theorem. No choice based only on the
unoriented boundary can distinguish these two caps. Boundary framing or an
oriented order could supply the same information; this API uses a side
witness explicitly.

The implementation uses the least singular direction for noisy data and
checks every sample's normalized incidence residual. It returns that maximum,
the Euclidean plane residual, the third singular value, and the side margin.
Nonunit data are rejected rather than normalized implicitly. An unresolved
rank or Lorentz norm is an error. Overdetermined data require an explicit
residual budget; no outlier is removed. Tests cover circles on both sides of
a hemisphere, rotations, input reorderings and degenerate data.

**Error bound.** Let A0 be a true rank-three plane matrix with unit null
direction u0. Let the observed matrix satisfy `||A-A0||_2 <= epsilon`, and
let u be the computed Euclidean unit direction. Write `r=||Au||_2` and let
s3 be the observed third singular value. If `epsilon<s3`, singular-value
perturbation and the decomposition perpendicular to u0 give

\[
\sin\angle(u,u_0)
\le\min\left(1,\frac{r+\epsilon}{s_3-\epsilon}\right).
\]

Indeed `||A0u|| <= r+epsilon`, while the restriction of A0 to the orthogonal
complement of its kernel has smallest singular value at least
`s3-epsilon`. The error input is externally supplied; a small fit residual
does not estimate calibration uncertainty. This bound concerns the plane
direction before Lorentz normalization. Near a point circle the spacelike
norm approaches zero and normalization becomes ill-conditioned; no uniform
cap-normal error bound is asserted there. Floating diagnostics do not
certify exact rank or an interval enclosure of roundoff.

## Connection and limitation: conformal geometry versus angular scale

With Lorentz metric `diag(-1,1,1,1)`, a sphere point p is the future null ray
represented by k=(1,p). Its cap incidence is `<n,k>`. For any proper
orthochronous Lorentz transformation Lambda, normalize the transformed ray
by its positive time component: `k'=Lambda k/(Lambda k)_0`. Then

\[
\langle\Lambda n,k'\rangle
=\frac{\langle n,k\rangle}{(\Lambda k)_0}.
\]

The boundary and inside sign are preserved, and `Lambda n` is the transformed
unit cap normal. The associated action on the sphere is Mobius; see
[Oblak, From the Lorentz Group to the Celestial Sphere](https://arxiv.org/abs/1508.00920)
for the standard correspondence. The tests check a nontrivial boost using
independent four-vector arithmetic and compare all supplied cross ratios.
That boost changes the angular radius while preserving those ratios.

Thus conformal data determine the configuration up to Mobius transformations,
and oriented incidence transports correctly under those transformations.
An angular radius relative to a fixed round metric additionally requires a
chosen conformal frame. The arbitrary gauge triple in this finite fixture
does not identify a physical observer's metric or clock. The unit-radius
representation and the inside witness are disclosed inputs.

## Reproduction and integration

From the repository root with its pinned requirements:

```text
python -m pytest -q code/geometry/test_finite_incidence.py code/geometry/test_conformal_readout.py code/geometry/test_quotient_cap_readout.py
```

The geometry CI workflow runs the complete geometry suite on Windows and
Linux. Independent controls include exact rational area determinants and
nullspaces, binary boundary-matrix homology, an integer fundamental-chain
map, analytic cap normals and Lorentz boosts. Malformed faces, omitted
simplices, duplicate midpoint assignments, invalid gauge data, nonfinite
coordinates and unresolved circles are negative tests. Rejection guards
also execute under optimized Python.

The modules split combinatorial incidence (`finite_incidence.py`) from
projective/cap calculations (`conformal_readout.py`). The existing
`quotient_cap_readout.py` imports those operations and retains the repair
fixtures and surface classification. Existing triangular topology controls
and the selected cyclic tower remain supported. Historical frozen receipts
retain their bytes. This audit does not derive spherical support from the
axioms, certify an infinite refinement limit, or select a physical model for
issues #1025/#1026. It is independent of the pending entropy and algebra
audits #1029/#1031 and changes none of their implementation files.
