# Native exchange statistics

## Complete topology of the declared path domain

The existing collision-free native configuration complex has fundamental
group **Z/2**, strengthening its previously established first homology
`H_1(Z)=Z/2`. Therefore its finite-dimensional flat unitary exchange
transports, up to changes of local frame, are exactly direct sums of the
trivial and sign characters. This classifies the possibilities under the
declared local deformations. It does not choose a character for the native
quantum instrument.

The statement concerns the binary total-two carrier and exactly the local
cells in [defect_history.py](../clock_defect_energy/defect_history.py): a
single particle crossing a native triangular face while the other remains
outside that face, and commuting moves along vertex-disjoint seams. It
assumes neither a continuum limit nor additional native control words.

The triangular fillings are **declared deformation relations**, not a
consequence of allowing moves along their three edges. Omitting those 180
fillings while keeping the same movement graph and 315 disjoint squares
changes the domain. Its square-boundary rank over `F_2` is 215, so its first
cohomology has dimension `300-66+1-215=20`; no fundamental group for that
changed domain is claimed here. There is also a direct continuous family:
put phase `exp(i*theta)` on one directed one-particle seam and its inverse
on the reverse. Every disjoint square has zero exponent, while an affected
spectator triangle has exponent one. This is flat on the domain with only
squares for arbitrary theta, and violates the original triangle relation
unless theta is a multiple of `2*pi`. The exhaustive two-character result
therefore uses the declared face fillings essentially. This negative
control leaves the predecessor's deformation contract unchanged.

### Exact finite identification

Let K be the native triangular surface. Its incidence is checked directly:
12 vertices, 30 edges, 20 distinct triangles, connected one-skeleton, exactly
two faces on every edge, connected circular vertex links, and a coherent
orientation of all triangles. Euler characteristic is two. These are the
hypotheses of the classification of connected closed surfaces, which gives
`|K|` homeomorphic to `S^2`. Euler characteristic alone would not suffice;
the validator separately rejects an open surface and a closed torus.

The ordered simplicial deleted product is

```
D(K) = union of |sigma| x |tau| over vertex-disjoint simplices sigma,tau.
```

Particle swap acts freely. A product of oriented simplices has the tensor
boundary `d(sigma x tau) = d(sigma) x tau + (-1)^dim(sigma) sigma x d(tau)`.
Swapping the two factors has orientation sign
`(-1)^(dim(sigma)*dim(tau))`. The implementation retains this sign when it
projects onto an unordered representative. In particular, reversing two
edge factors reverses the square orientation.

| Dimension | Unordered product cells | Count |
| --- | --- | ---: |
| 0 | vertex x vertex | 66 |
| 1 | vertex x edge | 300 |
| 2 | vertex x face, edge x edge | 180 + 315 |
| 3 | edge x face | 360 |
| 4 | face x face | 100 |

Every ordered count is twice the corresponding unordered count. The
certificate checks integral `d^2=0` in both complexes, and that the quotient
projection is an integral chain map. Its first three layers match every
existing native vertex, movement edge and oriented local polygon, up to
the harmless choice of each polygon's orientation. Thus the old native
complex X is exactly the two-skeleton of `D(K)/swap`.

The cells of dimension three and four are an auxiliary mathematical
completion used in the proof. They do not introduce native operations or
new relations on paths: attaching cells of dimension at least three
preserves the fundamental group. No integer or floating matrix rank is
being used to infer a fundamental group from homology.

### Equivariant deformation proof

There is a direct deformation from the ordered distinct-point space
`F_2(|K|)` onto D(K), without subdividing K. Write the two points in their
global barycentric coordinates `l_i,r_i`, and put

```
m_i = min(l_i,r_i),  M = sum_i m_i,
H_t(l,r) = ((l-t*m)/(1-t*M), (r-t*m)/(1-t*M)),  0 <= t <= 1.
```

Distinctness implies `M<1`, so the denominator stays positive. Each new
coordinate is nonnegative; each vector sums to one; and its support stays
inside the original simplex. The difference of the two coordinate vectors
is `(l-r)/(1-t*M)`, which never vanishes. At `t=1` the supports are
disjoint. If the original supports were disjoint then `m=0` and the pair
is fixed throughout. The formula is continuous, agrees across simplex
boundaries, and commutes with particle swap. It is therefore an equivariant
strong deformation retraction. The exact rational evaluator checks this
formula even for pairs arbitrarily close to the excluded diagonal; sampled
checks support the implementation, while the displayed inequalities prove
the statement for all real coordinates.

Transport the pair by any homeomorphism `|K| -> S^2`. For unit vectors x,y
with `x != y`, let `m=(x+y)/2` and `d=(x-y)/2`. Then `m.d=0` and `d != 0`.
The second equivariant strong deformation retraction is

```
G_t(x,y) = (((1-t)*m+d)/N_t, ((1-t)*m-d)/N_t),
N_t^2 = (1-t)^2 * ||m||^2 + ||d||^2 > 0.
```

Both outputs remain unit vectors and distinct. At `t=1` they are antipodal;
antipodal inputs stay fixed. Swap sends d to -d, hence exchanges the two
outputs. The antipodal-pair space is `S^2`, with swap its antipodal map.
Its quotient is `RP^2`. Equivariant deformation retractions descend to
quotients, so the full unordered deleted product has the homotopy type of
`RP^2`, and the ordered deleted product has the homotopy type of `S^2`.
Their fundamental groups are consequently `Z/2` and the trivial group.
The same groups hold for their two-skeleta.

This does **not** identify X itself with `RP^2` up to homotopy. X has Euler
characteristic `66-300+495=261`; the full unordered deleted product has
Euler characteristic one. The theorem upgrades the path classification,
without discarding the higher topology of the existing two-complex.

The connected ordered cover is the actual ordered-particle lift. The
certificate checks all 600 ordered movement edges against the original
sorted-pair chart parity. The native triangular exchange loop lifts to a
path ending at the swapped ordered pair, so it represents the nontrivial
generator. In particular, a closed native path is null-homotopic precisely
when this parity is even. The topology proof establishes existence of a
homotopy through the admitted local cells; it does not supply a short
rewrite algorithm or a bound on its number of elementary deformations.

### Exhaustive flat unitary representations

A flat rank-r unitary transport assigns a unitary between the fibers at
the endpoints of every edge, assigns its inverse on the reversed edge,
and has identity product around every admitted two-cell. Choosing frames
by transport along a spanning tree reduces its classification to unitary
representations of the fundamental group, up to conjugacy at a basepoint.

Since `pi_1(X)=Z/2`, the only remaining datum is a unitary U satisfying
`U^2=I`. Unitarity then gives `U*=U^{-1}=U`; its two orthogonal spectral
projections are `(I+U)/2` and `(I-U)/2`. Thus

```
U is unitarily equivalent to I_(r+) direct-sum (-I_(r-)),  r+ + r- = r.
```

There are exactly two one-dimensional characters, +1 and -1, and no
irreducible higher-dimensional unitary representation for this domain.
This excludes a free continuous exchange angle under these flatness
relations. It does not exclude curved connections, a changed deformation
domain, more particles, or additional quantum control factors. Nor does
it turn a classical path sign into a phase of a specified quantum channel.
Those operational choices require their own complete instrument and readout
account.

The exact implementation is
[statistics_topology.py](statistics_topology.py); its tests use independent
tensor boundaries, oriented deck transformations, all local cell lifts,
exact rational deformation identities, and malformed-domain controls.

```
python -m pytest -q model/exchange_sector/test_statistics_topology.py
```
